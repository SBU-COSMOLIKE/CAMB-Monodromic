    ! Equations module allowing for fairly general k-essence models
    ! Following notation and conventions from https://arxiv.org/pdf/1709.01544.pdf
    ! by João Rebouças, April 2024

    module KEssence
    use DarkEnergyInterface
    use results
    use constants
    use classes
    use Interpolation
    implicit none
    
    private

    real(dl), parameter :: Tpl = sqrt(kappa*hbar/c**5)  ! sqrt(8 pi G hbar/c^5), reduced planck time

    ! General base class. Specific implemenetations should inherit, defining Vofphi and setting up
    ! initial conditions and interpolation tables
    type, extends(TDarkEnergyModel) :: TKEssence
        integer :: DebugLevel = 0 !higher then zero for some debug output to console
        real(dl) :: astart = 1e-7_dl
        real(dl) :: integrate_tol = 1e-6_dl
        real(dl), dimension(:), allocatable :: sampled_a, phi_a, X_a
        ! Steps for log a and linear spacing, switching at max_a_log (set by Init)
        integer, private :: npoints_linear, npoints_log
        real(dl), private :: dloga, da, log_astart, max_a_log
        real(dl), private, dimension(:), allocatable :: ddphi_a, ddX_a
        class(CAMBdata), pointer, private :: State
    contains
        procedure :: Vofphi ! V(phi) potential [+ any cosmological constant]
        procedure :: ValsAta !get phi and phi' at scale factor a, e.g. by interpolation in precomputed table
        procedure :: Init => TKEssence_Init
        procedure :: PerturbedStressEnergy => TKEssence_PerturbedStressEnergy
        procedure :: PerturbationEvolve => TKEssence_PerturbationEvolve
        procedure :: BackgroundDensityAndPressure => TKEssence_BackgroundDensityAndPressure
        procedure :: EvolveBackground
        procedure :: EvolveBackgroundLog
        procedure, private :: X_start => TKEssence_X_start
        procedure :: GetOmegaFromInitial
    end type TKEssence

    ! JVR Modification: adding a new type, TMonodromicKEssence
    ! Check Equation (1) from https://arxiv.org/pdf/1709.01544.pdf
    type, extends(TKEssence) :: TMonodromicKEssence
        real(dl) :: alpha = 0.2_dl
        real(dl) :: C = 1.0_dl
        real(dl) :: A = 0.05_dl
        real(dl) :: nu = 50.0_dl 
        integer :: npoints = 200
    contains
        procedure :: Vofphi => TMonodromicKEssence_VofPhi
        procedure :: Init => TMonodromicKEssence_Init
        procedure :: ReadParams =>  TMonodromicKEssence_ReadParams
        procedure, nopass :: PythonClass => TMonodromicKEssence_PythonClass
        procedure, nopass :: SelfPointer => TMonodromicKEssence_SelfPointer
        procedure, private :: check_error => TMonodromicKEssence_check_error
        procedure :: get_initial_phi
        procedure :: get_initial_X
    end type TMonodromicKEssence

    procedure(TClassDverk) :: dverk

    public TKEssence,  TMonodromicKEssence
    
    contains

    ! -----------------------------------------------------------------------
    ! Base procedures
    ! -----------------------------------------------------------------------

    real(dl) function VofPhi(this, phi, deriv) result (Vofphii)
        !Get the KEssence potential as function of phi
        !The input variable phi is sqrt(8*Pi*G)*psi, where psi is the field
        !Returns (8*Pi*G)^(1-deriv/2)*d^{deriv}V(psi)/d^{deriv}psi evaluated at psi
        !return result is in 1/Mpc^2 units [so times (Mpc/c)^2 to get units in 1/Mpc^2]
        class(TKEssence), intent(in) :: this
        real(dl), intent(in) :: phi
        integer, intent(in) :: deriv
        real(dl), parameter :: units = MPC_in_sec**2 /Tpl**2  !convert to units of 1/Mpc^2

        call MpiStop('KEssence classes must override to provide VofPhi')
        VofPhii = 0
        !if (deriv==0) then
        !    Vofphi= norm*this%m*exp(-this%sigma_model*phi)
        !else if (deriv ==1) then
        !    Vofphi=-norm*this%m*sigma_model*exp(-this%sigma_model*phi)
        !else if (deriv ==2) then
        !    Vofphi=norm*this%m*sigma_model**2*exp(-this%sigma_model*phi)
        !else
        !    stop 'Invalid deriv in Vofphi'
        !end if
        !VofPhi = VOfPhi* MPC_in_sec**2 /Tpl**2  !convert to units of 1/Mpc^2
    end function VofPhi

    subroutine TKEssence_Init(this, State)
        !Make interpolation table, etc,
        !At this point massive neutrinos have been initialized
        !so grho_no_de can be used to get density and pressure of other components at scale factor a
        class(TKEssence), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State

        select type(State)
        class is (CAMBdata)
            this%State => State
        end select

        this%is_cosmological_constant = .false.
        this%num_perturb_equations = 2
        this%log_astart = log(this%astart)
    end subroutine  TKEssence_Init

    subroutine TKEssence_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
        !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
        class(TKEssence), intent(inout) :: this
        real(dl), intent(in) :: grhov, a
        real(dl), intent(out) :: grhov_t
        real(dl), optional, intent(out) :: w
        real(dl) V, a2, grhov_lambda, phi, X

        if (this%is_cosmological_constant) then
            grhov_t = grhov * a * a
            if (present(w)) w = -1_dl
        elseif (a >= this%astart) then
            call this%ValsAta(a, phi, X)
            V = this%Vofphi(phi, 0)
            grhov_t = V*X*(-1._dl + 3*X)*a*a
            if (present(w)) then
                w = (-1._dl + X)/(-1._dl + 3*X)
            end if
        else
            phi = this%phi_a(2)
            X = this%X_a(2)
            V = this%Vofphi(phi, 0)
            grhov_t = V*X*(-1._dl + 3*X)*a*a
            if (present(w)) then
                w = (-1._dl + X)/(-1._dl + 3*X)
            end if
        end if
    end subroutine TKEssence_BackgroundDensityAndPressure

    subroutine EvolveBackgroundLog(this,num,loga,y,yprime)
        ! Evolve the background equation in terms of loga.
        ! Variables are phi=y(1), a^2 phi' = y(2)
        ! Assume otherwise standard background components
        class(TKEssence) :: this
        integer num
        real(dl) y(num),yprime(num)
        real(dl) loga, a

        a = exp(loga)
        call this%EvolveBackground(num, a, y, yprime)
        yprime = yprime*a
    end subroutine EvolveBackgroundLog

    subroutine EvolveBackground(this,num,a,y,yprime)
        ! Evolve the background equation in terms of a.
        ! Variables are phi=y(1), a^2 phi' = y(2)
        ! Assume otherwise standard background components
        class(TKEssence) :: this
        integer :: num
        real(dl), intent(inout) :: y(num),yprime(num)
        real(dl) :: a, a2, tot
        real(dl) :: phi, grhode, X, H

        a2 = a**2
        phi = y(1)
        X = y(2)

        grhode = this%Vofphi(phi, 0)*X*(-1._dl + 3*X)
        tot = this%state%grho_no_de(a)/(a2*a2) + grhode ! 8*pi*G*rho

        H = sqrt(tot/3.0d0)
        yprime(1) = sqrt(2*X*this%State%grhocrit/3)/(a*H) ! dphi/da
        
        yprime(2) = -sqrt(this%State%grhocrit/3)*this%Vofphi(phi, 1)*sqrt(2*X)*(-X + 3*X**2)/(this%Vofphi(phi, 0)*(6*X - 1)) - 6*H*X*(2*X - 1)/(6*X - 1) ! dX/dt
        yprime(2) = yprime(2)/(a*H) ! dX/da
    end subroutine EvolveBackground

    real(dl) function TKEssence_X_start(this,phi)
        class(TKEssence) :: this
        real(dl) :: phi
        TKEssence_X_start = 0
    end function TKEssence_X_start

    subroutine ValsAta(this,a,aphi,aX)
        class(TKEssence) :: this
        !Do interpolation for background phi and X at a (precomputed in Init)
        real(dl) a, aphi, aX
        real(dl) a0,b0,ho2o6,delta,da
        integer ix

        if (a >= 0.9999999d0) then
            aphi= this%phi_a(this%npoints_linear+this%npoints_log)
            aX= this%X_a(this%npoints_linear+this%npoints_log)
            return
        elseif (a < this%astart) then
            aphi = this%phi_a(2)
            aX = this%X_a(2) ! JVR NOTE: this makes things slower but it works for A = 0!
            return
        elseif (a > this%max_a_log) then
            delta= a-this%max_a_log
            ix = this%npoints_log + int(delta/this%da)
        else
            delta= log(a)-this%log_astart
            ix = int(delta/this%dloga)+1
        end if
        da = this%sampled_a(ix+1) - this%sampled_a(ix)
        a0 = (this%sampled_a(ix+1) - a)/da
        b0 = 1 - a0
        ho2o6 = da**2/6._dl
        aphi=b0*this%phi_a(ix+1) + a0*(this%phi_a(ix)-b0*((a0+1)*this%ddphi_a(ix)+(2-a0)*this%ddphi_a(ix+1))*ho2o6)
        aX=b0*this%X_a(ix+1) + a0*(this%X_a(ix)-b0*((a0+1)*this%ddX_a(ix)+(2-a0)*this%ddX_a(ix+1))*ho2o6)
    end subroutine ValsAta

    subroutine TKEssence_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
        !Get density perturbation and heat flux
        class(TKEssence), intent(inout) :: this
        real(dl), intent(out) :: dgrhoe, dgqe
        real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
        real(dl), intent(in) :: ay(*)
        real(dl), intent(inout) :: ayprime(*)
        integer, intent(in) :: w_ix
        real(dl) :: phi, X, delta_phi, delta_phi_prime, delta_X, phidot, V, V_prime
        
        ! Following Kunhao's notes
        call this%ValsAta(a, phi, X)
        V = this%Vofphi(phi, 0)
        V_prime = this%Vofphi(phi, 1)
        phidot = sqrt(2*X*this%State%grhocrit/3)
        delta_phi = ay(w_ix)
        delta_phi_prime = ay(w_ix+1)
        delta_X = phidot*delta_phi_prime/a/(this%state%grhocrit/3)
        dgrhoe = V*(-1._dl + 2*X)*delta_X - V_prime*X*(-1._dl + X)*delta_phi + 4*X*V*delta_X + 2*X*V_prime*(-1._dl + 2*X)*delta_phi
        dgqe = V*(-1._dl + 2*X)*k*phidot*delta_phi/a
    end subroutine TKEssence_PerturbedStressEnergy

    subroutine TKEssence_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
        !Get conformal time derivatives of the density perturbation and velocity
        class(TKEssence), intent(in) :: this
        real(dl), intent(inout) :: ayprime(:)
        real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
        integer, intent(in) :: w_ix
        real(dl) :: phi, X, X_prime, X_dot, V, V_prime, V_primeprime, phidot, phi_prime, phi_primeprime, phi_dotdot, H_curly, a2
        real(dl) :: delta_phi, delta_phi_prime
        real(dl) :: grhode, tot
        real(dl) :: P_X, P_XX, P_Xphi, P_phiphi, P_phiphiX
        real(dl) :: A_tilde, B_tilde, C_tilde, D_tilde

        ! Following Kunhao's thesis, Section 5.3.3, Equation 5.4.1
        a2 = a*a
        call this%ValsAta(a, phi, X)
        phidot = sqrt(2*X*this%State%grhocrit/3)
        phi_prime = a*phidot
        delta_phi = y(w_ix)
        delta_phi_prime = y(w_ix+1)
        V = this%Vofphi(phi, 0)
        V_prime = this%Vofphi(phi, 1)
        V_primeprime = this%Vofphi(phi, 2)
        
        grhode = V*X*(-1._dl + 3*X)*a2
        tot = this%state%grho_no_de(a)/a2 + grhode ! 8*pi*G*rho*a2
        H_curly = sqrt(tot/3._dl)

        ! P = V(-X + X^2)
        P_X = V*(-1._dl + 2*X)
        P_XX = 2*V
        P_Xphi = V_prime*(-1._dl + 2*X)
        P_phiphi = V_primeprime*X*(-1._dl + X)
        P_phiphiX = V_primeprime*(-1._dl + 2*X)

        X_dot = -sqrt(this%State%grhocrit/3)*V_prime*sqrt(2*X)*(-X + 3*X**2)/(V*(6*X - 1)) - 6*(H_curly/a)*X*(2*X - 1)/(6*X - 1)
        X_prime = a*X_dot
        phi_dotdot = (this%State%grhocrit/3) * X_dot/phidot
        phi_primeprime = a2*(phi_dotdot + H_curly*phidot/a)
        A_tilde = P_XX*phidot**2 + P_X
        B_tilde = 2*H_curly*P_X + 2*V_prime*a*phidot**3 + P_XX*(3*phi_prime*phi_primeprime/a2 - H_curly*phidot**2) + phi_prime*P_Xphi
        C_tilde = -a2*P_phiphi + k*k*P_X + 2*V_prime*X_prime*phi_prime + P_Xphi*phi_primeprime + P_phiphiX*phi_prime**2 + 2*H_curly*P_Xphi*phi_prime
        D_tilde = k*z*P_X*phi_prime
        ayprime(w_ix) = delta_phi_prime ! delta_phi'
        ayprime(w_ix+1) = -(D_tilde + C_tilde*delta_phi + B_tilde*delta_phi_prime)/A_tilde ! delta_phi''
    end subroutine TKEssence_PerturbationEvolve

    real(dl) function GetOmegaFromInitial(this, astart, phi, X, atol)
        ! Get omega_de today given particular conditions phi and X at a = astart
        class(TKEssence) :: this
        real(dl), intent(IN) :: astart, phi, X, atol
        integer, parameter ::  NumEqs = 2
        real(dl) :: c(24), w(NumEqs, 9), y(NumEqs), a_switch, y_prime(2)
        integer :: ind, i
        integer, parameter :: nsteps_log = 1000, nsteps_linear = 2000
        
        real(dl) :: da, dloga, loga, a
        
        ind = 1
        a_switch = 1e-3
        dloga = (log(a_switch) - log(astart))/nsteps_log
        y(1) = phi
        y(2) = X
        
        do i = 1, nsteps_log
            loga = log(astart) + i*dloga
            call this%EvolveBackgroundLog(NumEqs, loga, y, y_prime)
            y(1) = y(1) + y_prime(1)*dloga
            y(2) = y(2) + y_prime(2)*dloga
        end do

        da = (1._dl - a_switch)/nsteps_linear
        do i = 1, nsteps_linear
            a = a_switch + i*da
            call this%EvolveBackground(NumEqs, a, y, y_prime)
            y(1) = y(1) + y_prime(1)*da
            y(2) = y(2) + y_prime(2)*da
        end do

        GetOmegaFromInitial = this%Vofphi(y(1), 0)*y(2)*(-1._dl + 3*y(2))/this%State%grhocrit
    end function GetOmegaFromInitial

    ! -----------------------------------------------------------------------
    ! Implementations for Monodromic KEssence model
    ! -----------------------------------------------------------------------

    real(dl) function TMonodromicKEssence_VofPhi(this, phi, deriv) result(VofPhi)
        !The input variable phi is sqrt(8*Pi*G)*psi
        !Returns (8*Pi*G)^(1-deriv/2)*d^{deriv}V(psi)/d^{deriv}psi evaluated at psi
        !return result is in 1/Mpc^2 units [so times (Mpc/c)^2 to get units in 1/Mpc^2]
        class(TMonodromicKEssence), intent(in) :: this
        real(dl), intent(in) :: phi
        integer, intent(in) :: deriv
        real(dl), parameter :: units = MPC_in_sec**2 /Tpl**2  !convert to units of 1/Mpc^2

        ! Assume f = sqrt(kappa)*f_theory = f_theory/M_pl
        ! m = m_theory/M_Pl
        ! JVR Note: removing the `units` factor that appears in both other implementations
        ! Assuming C has the units of `grhocrit`
        ! Because w is close to -1 and rho_de ~ grho_crit ~ V
        if (deriv == 0) then
            if (phi < 0) then
                ! print *, "WARNING: negative field value: this shouldn't happen for this model. Using the second phi value in the interpolation table instead"
                Vofphi = this%C*this%phi_a(2)**(-this%alpha)*(1.0_dl - this%A*sin(this%nu*this%phi_a(2)))
                return
            end if
            Vofphi = this%C*phi**(-this%alpha)*(1.0_dl - this%A*sin(this%nu*phi))
            
            if (isnan(Vofphi)) then
                print*, "ERROR: for phi =", phi, "V is NaN"
                stop
            end if
        else if (deriv == 1) then
            if (phi < 0) then
                ! print *, "WARNING: negative field value: this shouldn't happen for this model. Using the second phi value in the interpolation table instead"
                Vofphi = this%C*(-this%alpha*this%phi_a(2)**(-this%alpha-1)*(1 - this%A*sin(this%nu*this%phi_a(2))) - this%phi_a(2)**(-this%alpha)*this%A*this%nu*cos(this%nu*this%phi_a(2)))
                return
            end if

            Vofphi = this%C*(-this%alpha*phi**(-this%alpha-1)*(1 - this%A*sin(this%nu*phi)) - phi**(-this%alpha)*this%A*this%nu*cos(this%nu*phi))
            if (isnan(Vofphi)) then
                print*, "ERROR: for phi =", phi, "V' is NaN"
                stop
            end if
        else if (deriv == 2) then
            if (phi < 0) then
                ! print *, "WARNING: negative field value: this shouldn't happen for this model. Using the second phi value in the interpolation table instead"
                Vofphi = this%C*this%phi_a(2)**(-this%alpha-2)*(this%A*sin(this%nu*this%phi_a(2))*(this%phi_a(2)**2*this%nu**2 - this%alpha**2 - this%alpha) + 2._dl*this%A*this%nu*this%alpha*this%phi_a(2)*cos(this%nu*this%phi_a(2)) + this%alpha*(1._dl + this%alpha))
                return
            end if
            Vofphi = this%C*phi**(-this%alpha-2)*(this%A*sin(this%nu*phi)*(phi**2*this%nu**2 - this%alpha**2 - this%alpha) + 2._dl*this%A*this%nu*this%alpha*phi*cos(this%nu*phi) + this%alpha*(1._dl + this%alpha))

            if (isnan(Vofphi)) then
                print*, "ERROR: for phi =", phi, "V'' is NaN"
                stop
            end if
        end if
    end function TMonodromicKEssence_VofPhi

    real(dl) function get_initial_X(this, a) result(initial_X)
        ! Initial conditions from tracking solution
        ! See equation (18) from https://arxiv.org/pdf/1709.01544.pdf
        class(TMonodromicKEssence), intent(in) :: this
        real(dl), intent(in) :: a
        real(dl) :: H_ini, t_ini, p, grho_no_de_at_initial_a, phi_tilde_0

        initial_X = (4._dl - this%alpha)/(8._dl - 3*this%alpha)
    end function get_initial_X

    real(dl) function get_initial_phi(this, a) result(initial_phi)
        ! Initial conditions from tracking solution
        ! See equation (19) from https://arxiv.org/pdf/1709.01544.pdf
        class(TMonodromicKEssence), intent(in) :: this
        real(dl), intent(in) :: a
        real(dl) :: H_ini, t_ini, grho_no_de_at_initial_a, X_initial

        grho_no_de_at_initial_a = this%State%grho_no_de(a)/a**4
        X_initial = this%get_initial_X(a)
        H_ini = sqrt(grho_no_de_at_initial_a/3._dl)
        t_ini = 1._dl/(2._dl*H_ini)
        
        initial_phi = sqrt(2*X_initial*this%state%grhocrit/3._dl) * t_ini
    end function get_initial_phi

    subroutine TMonodromicKEssence_Init(this, State)
        class(TMonodromicKEssence), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State
        integer,  parameter :: NumEqs = 2, max_iters = 20
        integer,  parameter :: nsteps_linear = 1500, nsteps_log = 1500, nsteps = nsteps_log + nsteps_linear
        real(dl), parameter :: omega_de_tol = 1e-4
        real(dl), parameter :: splZero = 0._dl
        real(dl), parameter :: a_start = 1e-5, a_switch = 1e-3
        real(dl), parameter :: dloga = (log(a_switch) - log(a_start))/nsteps_log, da = (1._dl - a_switch)/nsteps_linear
        real(dl)            :: y(NumEqs), y_prime(NumEqs)
        real(dl)            :: omega_de_target, om, om1, om2
        real(dl)            :: C_1, C_2, new_C, a, loga, atol, initial_phi, initial_X, a_line, b_line, error
        real(dl)            :: phi, X
        integer             :: i
        Type(TTimer)        :: Timer
        real(dl)            :: grho_no_de, grho_de, fde

        !Make interpolation table, etc,
        !At this point massive neutrinos have been initialized
        !so grho_no_de can be used to get density and pressure of other components at scale factor a
        
        call this%TKEssence%Init(State)
        select type(State)
        class is (CAMBdata)
            omega_de_target = State%Omega_de
        end select

        if (allocated(this%phi_a)) then
            print*, "WARNING: the interpolation table is already allocated. This shouldn't be happening but we are deallocating anyway"
            deallocate(&
                this%phi_a,      &
                this%X_a,   &
                this%ddphi_a,    &
                this%ddX_a, &
                this%sampled_a   &
            )
        end if

        allocate(&
            this%phi_a(nsteps),      &
            this%X_a(nsteps),   &
            this%ddphi_a(nsteps),    &
            this%ddX_a(nsteps), &
            this%sampled_a(nsteps)   &
        )
        
        ! Binary search for C
        C_1 = this%State%grhov * 0.5_dl
        C_2 = this%State%grhov * 1.3_dl
        print*, "Shooting for C with tentative values: ", C_1, C_2
        
        ! See if current C is giving correct omega_de now
        atol = 1d-8
        initial_X = this%get_initial_X(a_start)
        this%C = C_1
        initial_phi = this%get_initial_phi(a_start)
        om1 = this%GetOmegaFromInitial(a_start, initial_phi, initial_X, atol)
        this%C = C_2
        om2 = this%GetOmegaFromInitial(a_start, initial_phi, initial_X, atol)
        print*, "Target Omega_de:", omega_de_target
        print*, "C = ", C_1, "=> omega_de = ", om1
        print*, "C = ", C_2, "=> omega_de = ", om2
        
        do i = 1, max_iters
            if (om1 > omega_de_target .or. om2 < omega_de_target) then
                write (*,*) 'WARNING: initial guesses for C did not bracket the required value'
            end if
            a_line = (om2 - om1)/(C_2 - C_1)
		    b_line = om2 - a_line*C_2
            new_C = (omega_de_target - b_line)/a_line
            this%C = new_C
            om = this%GetOmegaFromInitial(a_start, initial_phi, initial_X, atol)
            error = (om - omega_de_target)/omega_de_target
            print*, "C = ", new_C, "=> omega_de = ", om, "(error = ", error, ")"
            
            if (abs(error) < omega_de_tol) then 
                print*, "Finished shooting successfully after ", i, "iterations"
                exit
            end if

            if (om < omega_de_target) then
                om1 = om
                C_1 = new_C
            else
                om2 = om
                C_2 = new_C
            end if
        end do

        y(1) = initial_phi
        y(2) = initial_X

        do i = 1, nsteps_log
            loga = log(a_start) + i*dloga
            call this%EvolveBackgroundLog(NumEqs, loga, y, y_prime)
            y(1) = y(1) + y_prime(1)*dloga
            y(2) = y(2) + y_prime(2)*dloga
            this%sampled_a(i) = exp(loga)
            this%phi_a(i) = y(1)
            this%X_a(i) = y(2)
        end do

        do i = 1, nsteps_linear
            a = a_switch + i*da
            call this%EvolveBackground(NumEqs, a, y, y_prime)
            y(1) = y(1) + y_prime(1)*da
            y(2) = y(2) + y_prime(2)*da
            this%sampled_a(nsteps_log + i) = a
            this%phi_a(nsteps_log + i) = y(1)
            this%X_a(nsteps_log + i) = y(2)
        end do

        ! JVR NOTE: we need to deallocate phi_a, X_a, sampled_a
        ! this might be causing memory leaks and subsequent segmentation faults in original CAMB!
        ! deallocate(phi_a, X_a, sampled_a)
        
        ! Must set the fields
        this%astart = a_start
        this%npoints_linear = nsteps_linear
        this%npoints_log = nsteps_log
        this%da = da
        this%dloga = dloga
        this%log_astart = log(a_start)
        this%max_a_log = a_switch

        call spline(this%sampled_a, this%phi_a, nsteps, splZero, splZero, this%ddphi_a)
        call spline(this%sampled_a, this%X_a, nsteps, splZero, splZero, this%ddX_a)

        call this%ValsAta(a, phi, X)

    end subroutine TMonodromicKEssence_Init

    subroutine TMonodromicKEssence_ReadParams(this, Ini)
        use IniObjects
        class(TMonodromicKEssence) :: this
        class(TIniFile), intent(in) :: Ini
        call this%TDarkEnergyModel%ReadParams(Ini)
    end subroutine TMonodromicKEssence_ReadParams

    function TMonodromicKEssence_PythonClass()
        character(LEN=:), allocatable :: TMonodromicKEssence_PythonClass
        TMonodromicKEssence_PythonClass = 'MonodromicKEssence'
    end function TMonodromicKEssence_PythonClass

    subroutine TMonodromicKEssence_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TMonodromicKEssence), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TMonodromicKEssence_SelfPointer

    logical function TMonodromicKEssence_check_error(this, afrom, aend) result (check_error)
        class(TMonodromicKEssence) :: this
        real(dl) afrom, aend

        if (global_error_flag/=0) then
            write(*,*) 'MonodromicKEssence error integrating', afrom, aend
            write(*,*) this%A, this%C, this%nu, this%alpha
            stop
            check_error = .false.
            return
        end if
        check_error = .true.
    end function TMonodromicKEssence_check_error

    end module KEssence
