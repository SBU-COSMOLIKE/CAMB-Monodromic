    ! Equations module allowing for fairly general quintessence models
    !
    ! by Antony Lewis (http://cosmologist.info/)

    !!FIX March 2005: corrected update to next treatment of tight coupling
    !!Fix Oct 2011: a^2 factor in ayprime(EV%w_ix+1) [thanks Raphael Flauger]
    ! Oct 2013: update for latest CAMB, thanks Nelson Lima, Martina Schwind
    ! May 2020: updated for CAMB 1.x+

    ! Notes at http://antonylewis.com/notes/CAMB.pdf

    !This module is not well tested, use at your own risk!

    !Need to specify Vofphi function, and also initial_phi
    !You may also need to change other things to get it to work with different types of quintessence model

    !It works backwards, in that it assumes Omega_de is Omega_Q today, then does a binary search on the
    !initial conditions to find what is required to give that Omega_Q today after evolution.

    module Quintessence
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
    type, extends(TDarkEnergyModel) :: TQuintessence
        integer :: DebugLevel = 0 !higher then zero for some debug output to console
        real(dl) :: astart = 1e-7_dl
        real(dl) :: integrate_tol = 1e-6_dl
        real(dl), dimension(:), allocatable :: sampled_a, phi_a, phidot_a
        ! Steps for log a and linear spacing, switching at max_a_log (set by Init)
        integer, private :: npoints_linear, npoints_log
        real(dl), private :: dloga, da, log_astart, max_a_log
        real(dl), private, dimension(:), allocatable :: ddphi_a, ddphidot_a
        class(CAMBdata), pointer, private :: State
    contains
        procedure :: Vofphi !V(phi) potential [+ any cosmological constant]
        procedure :: ValsAta !get phi and phi' at scale factor a, e.g. by interpolation in precomputed table
        procedure :: Init => TQuintessence_Init
        procedure :: PerturbedStressEnergy => TQuintessence_PerturbedStressEnergy
        procedure :: PerturbationEvolve => TQuintessence_PerturbationEvolve
        procedure :: BackgroundDensityAndPressure => TQuintessence_BackgroundDensityAndPressure
        procedure :: EvolveBackground
        procedure :: EvolveBackgroundLog
        procedure, private :: phidot_start => TQuintessence_phidot_start
        procedure :: GetOmegaFromInitial
    end type TQuintessence

    ! Specific implementation for early quintessence + cosmologial constant, assuming the early component
    ! energy density fraction is negligible at z=0.
    ! The specific parameterization of the potential implemented is the axion model of arXiv:1908.06995
    type, extends(TQuintessence) :: TEarlyQuintessence
        real(dl) :: n = 3._dl
        real(dl) :: f =0.05 ! sqrt(8*pi*G)*f
        real(dl) :: m = 5d-54 !m in reduced Planck mass units
        real(dl) :: theta_i = 3.1_dl !initial value of phi/f
        real(dl) :: frac_lambda0 = 1._dl !fraction of dark energy density that is cosmological constant today
        logical :: use_zc = .true. !adjust m to fit zc
        real(dl) :: zc, fde_zc !readshift for peak f_de and f_de at that redshift
        integer :: npoints = 5000 !baseline number of log a steps; will be increased if needed when there are oscillations
        integer :: min_steps_per_osc = 10
        real(dl), dimension(:), allocatable :: fde, ddfde
    contains
        procedure :: Vofphi => TEarlyQuintessence_VofPhi
        procedure :: Init => TEarlyQuintessence_Init
        procedure :: ReadParams =>  TEarlyQuintessence_ReadParams
        procedure, nopass :: PythonClass => TEarlyQuintessence_PythonClass
        procedure, nopass :: SelfPointer => TEarlyQuintessence_SelfPointer
        procedure, private :: fdeAta
        procedure, private :: fde_peak
        procedure, private :: check_error
        procedure :: calc_zc_fde
    end type TEarlyQuintessence

    ! JVR Modification: adding a new type, TMonodromicQuintessence
    ! Check Equation (1) from https://arxiv.org/pdf/1709.01544.pdf
    type, extends(TQuintessence) :: TMonodromicQuintessence
        real(dl) :: alpha = 0.2_dl
        real(dl) :: C = 1.0_dl
        real(dl) :: A = 0.05_dl
        real(dl) :: nu = 50.0_dl 
        integer :: npoints = 200
    contains
        procedure :: Vofphi => TMonodromicQuintessence_VofPhi
        procedure :: Init => TMonodromicQuintessence_Init
        procedure :: ReadParams =>  TMonodromicQuintessence_ReadParams
        procedure, nopass :: PythonClass => TMonodromicQuintessence_PythonClass
        procedure, nopass :: SelfPointer => TMonodromicQuintessence_SelfPointer
        procedure, private :: check_error => TMonodromicQuintessence_check_error
        procedure :: get_initial_phi
    end type TMonodromicQuintessence

    procedure(TClassDverk) :: dverk

    public TQuintessence, TEarlyQuintessence, TMonodromicQuintessence
    
    contains

    ! -----------------------------------------------------------------------
    ! Base procedures
    ! -----------------------------------------------------------------------

    real(dl) function VofPhi(this, phi, deriv) result (Vofphii)
        !Get the quintessence potential as function of phi
        !The input variable phi is sqrt(8*Pi*G)*psi, where psi is the field
        !Returns (8*Pi*G)^(1-deriv/2)*d^{deriv}V(psi)/d^{deriv}psi evaluated at psi
        !return result is in 1/Mpc^2 units [so times (Mpc/c)^2 to get units in 1/Mpc^2]
        class(TQuintessence), intent(in) :: this
        real(dl), intent(in) :: phi
        integer, intent(in) :: deriv
        real(dl), parameter :: units = MPC_in_sec**2 /Tpl**2  !convert to units of 1/Mpc^2

        call MpiStop('Quintessence classes must override to provide VofPhi')
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

    subroutine TQuintessence_Init(this, State)
        !Make interpolation table, etc,
        !At this point massive neutrinos have been initialized
        !so grho_no_de can be used to get density and pressure of other components at scale factor a
        class(TQuintessence), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State

        select type(State)
        class is (CAMBdata)
            this%State => State
        end select

        this%is_cosmological_constant = .false.
        this%num_perturb_equations = 2
        this%log_astart = log(this%astart)
    end subroutine  TQuintessence_Init

    subroutine TQuintessence_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
        !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
        class(TQuintessence), intent(inout) :: this
        real(dl), intent(in) :: grhov, a
        real(dl), intent(out) :: grhov_t
        real(dl), optional, intent(out) :: w
        real(dl) V, a2, grhov_lambda, phi, phidot

        if (this%is_cosmological_constant) then
            grhov_t = grhov * a * a
            if (present(w)) w = -1_dl
        elseif (a >= this%astart) then
            a2 = a**2
            call this%ValsAta(a, phi, phidot)
            V = this%Vofphi(phi, 0)
            grhov_t = phidot**2/2 + a2*V
            if (present(w)) then
                w = (phidot**2/2 - a2*V)/grhov_t
            end if
        else
            grhov_t = 0
            if (present(w)) w = -1
        end if
    end subroutine TQuintessence_BackgroundDensityAndPressure

    subroutine EvolveBackgroundLog(this,num,loga,y,yprime)
        ! Evolve the background equation in terms of loga.
        ! Variables are phi=y(1), a^2 phi' = y(2)
        ! Assume otherwise standard background components
        class(TQuintessence) :: this
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
        class(TQuintessence) :: this
        integer :: num
        real(dl) y(num),yprime(num)
        real(dl) a, a2, tot
        real(dl) phi, grhode, phi_prime, a_prime

        a2 = a**2
        phi = y(1)
        phi_prime = y(2)/a2

        grhode = a2*(0.5d0*phi_prime**2 + a2*this%Vofphi(phi, 0))
        tot = this%state%grho_no_de(a) + grhode

        a_prime = sqrt(tot/3.0d0)
        yprime(1) = phi_prime/a_prime !d phi /d a
        yprime(2) = -a2**2*this%Vofphi(phi,1)/a_prime
    end subroutine EvolveBackground

    real(dl) function TQuintessence_phidot_start(this,phi)
        class(TQuintessence) :: this
        real(dl) :: phi
        TQuintessence_phidot_start = 0
    end function TQuintessence_phidot_start

    subroutine ValsAta(this,a,aphi,aphidot)
        class(TQuintessence) :: this
        !Do interpolation for background phi and phidot at a (precomputed in Init)
        real(dl) a, aphi, aphidot
        real(dl) a0,b0,ho2o6,delta,da
        integer ix

        if (a >= 0.9999999d0) then
            aphi= this%phi_a(this%npoints_linear+this%npoints_log)
            aphidot= this%phidot_a(this%npoints_linear+this%npoints_log)
            return
        elseif (a < this%astart) then
            aphi = this%phi_a(1)
            aphidot = 0
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
        aphidot=b0*this%phidot_a(ix+1) + a0*(this%phidot_a(ix)-b0*((a0+1)*this%ddphidot_a(ix)+(2-a0)*this%ddphidot_a(ix+1))*ho2o6)
    end subroutine ValsAta

    subroutine TQuintessence_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
        !Get density perturbation and heat flux
        class(TQuintessence), intent(inout) :: this
        real(dl), intent(out) :: dgrhoe, dgqe
        real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
        real(dl), intent(in) :: ay(*)
        real(dl), intent(inout) :: ayprime(*)
        integer, intent(in) :: w_ix
        real(dl) phi, phidot, clxq, vq

        call this%ValsAta(a,phi,phidot)
        clxq=ay(w_ix)
        vq=ay(w_ix+1)
        dgrhoe= phidot*vq +clxq*a**2*this%Vofphi(phi,1)
        dgqe= k*phidot*clxq
    end subroutine TQuintessence_PerturbedStressEnergy

    subroutine TQuintessence_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
        !Get conformal time derivatives of the density perturbation and velocity
        class(TQuintessence), intent(in) :: this
        real(dl), intent(inout) :: ayprime(:)
        real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
        integer, intent(in) :: w_ix
        real(dl) clxq, vq, phi, phidot

        call this%ValsAta(a,phi,phidot) !wasting time calling this again..
        clxq=y(w_ix)
        vq=y(w_ix+1)
        ayprime(w_ix)= vq
        ayprime(w_ix+1) = - 2*adotoa*vq - k*z*phidot - k**2*clxq - a**2*clxq*this%Vofphi(phi,2)
    end subroutine TQuintessence_PerturbationEvolve

    real(dl) function GetOmegaFromInitial(this, astart, phi, phidot, atol)
        !Get omega_de today given particular conditions phi and phidot at a = astart
        class(TQuintessence) :: this
        real(dl), intent(IN) :: astart, phi, phidot, atol
        integer, parameter ::  NumEqs = 2
        real(dl) :: c(24), w(NumEqs, 9), y(NumEqs), a_switch, y_prime(2)
        integer :: ind, i
        integer, parameter :: nsteps_log = 2000, nsteps_linear = 2000
        
        real(dl) :: da, dloga, loga, a
        
        ind = 1
        a_switch = 1e-3
        dloga = (log(a_switch) - log(astart))/nsteps_log
        y(1) = phi
        y(2) = phidot*astart**2
        
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
        
        GetOmegaFromInitial = (0.5d0*y(2)**2 + this%Vofphi(y(1),0))/this%State%grhocrit
    end function GetOmegaFromInitial

    ! -----------------------------------------------------------------------
    ! Implementations for Early Quintessence example, axion potential from e.g. arXiv: 1908.06995
    ! -----------------------------------------------------------------------

    real(dl) function TEarlyQuintessence_VofPhi(this, phi, deriv) result(Vofphi)
        !The input variable phi is sqrt(8*Pi*G)*psi
        !Returns (8*Pi*G)^(1-deriv/2)*d^{deriv}V(psi)/d^{deriv}psi evaluated at psi
        !return result is in 1/Mpc^2 units [so times (Mpc/c)^2 to get units in 1/Mpc^2]
        class(TEarlyQuintessence    ), intent(in) :: this
        real(dl), intent(in) :: phi
        integer, intent(in) :: deriv
        real(dl), parameter :: units = MPC_in_sec**2 /Tpl**2  !convert to units of 1/Mpc^2
        real(dl) :: theta, costheta

        ! Assume f = sqrt(kappa)*f_theory = f_theory/M_pl
        ! m = m_theory/M_Pl
        theta = phi/this%f
        if (deriv==0) then
            Vofphi = units*this%m**2*this%f**2*(1 - cos(theta))**this%n + this%frac_lambda0*this%State%grhov
        else if (deriv ==1) then
            Vofphi = units*this%m**2*this%f*this%n*(1 - cos(theta))**(this%n-1)*sin(theta)
        else if (deriv ==2) then
            costheta = cos(theta)
            Vofphi = units*this%m**2*this%n*(1 - costheta)**(this%n-1)*(this%n*(1+costheta) -1)
        end if
    end function TEarlyQuintessence_VofPhi

    subroutine TEarlyQuintessence_Init(this, State)
        use Powell
        class(TEarlyQuintessence), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State
        real(dl) aend, afrom
        integer, parameter ::  NumEqs=2
        real(dl) c(24),w(NumEqs,9), y(NumEqs)
        integer ind, i, ix
        real(dl), parameter :: splZero = 0._dl
        real(dl) lastsign, da_osc, last_a, a_c
        real(dl) initial_phi, initial_phidot, a2
        real(dl), dimension(:), allocatable :: sampled_a, phi_a, phidot_a, fde
        integer npoints, tot_points, max_ix
        logical has_peak
        real(dl) fzero, xzero
        integer iflag, iter
        Type(TTimer) :: Timer
        Type(TNEWUOA) :: Minimize
        real(dl) log_params(2), param_min(2), param_max(2)

        !Make interpolation table, etc,
        !At this point massive neutrinos have been initialized
        !so grho_no_de can be used to get density and pressure of other components at scale factor a

        call this%TQuintessence%Init(State)

        if (this%use_zc) then
            !Find underlying parameters m,f to give specified zc and fde_zc (peak early dark energy fraction)
            !Input m,f are used as starting values for search, which is done by brute force
            !(so should generalize easily, but not optimized for this specific potential)
            log_params(1) = log(this%f)
            log_params(2) = log(this%m)

            if (.false.) then
                ! Can just iterate linear optimizations when nearly orthogonal
                call Timer%Start()
                do iter = 1, 2
                    call brentq(this,match_fde,log(0.01_dl),log(10._dl), 1d-3,xzero,fzero,iflag)
                    if (iflag/=0) print *, 'BRENTQ FAILED f'
                    this%f = exp(xzero)
                    print *, 'match to m, f =', this%m, this%f, fzero
                    call brentq(this,match_zc,log(1d-55),log(1d-52), 1d-3,xzero,fzero,iflag)
                    if (iflag/=0) print *, 'BRENTQ FAILED m'
                    this%m = exp(xzero)
                    print *, 'match to m, f =', this%m, this%f, fzero
                    call this%calc_zc_fde(fzero, xzero)
                    print *, 'matched outputs', fzero, xzero
                end do
                call Timer%WriteTime('Timing for fitting')
            end if
            if (this%DebugLevel>0) call Timer%Start()
            !Minimize in log f, log m
            ! param_min(1) = log(0.001_dl)
            ! param_min(2) = log(1d-58)
            ! param_max(1) = log(1e5_dl)
            ! param_max(2) = log(1d-50)
            ! if (Minimize%BOBYQA(this, match_fde_zc, 2, 5, log_params,param_min, &
            !           param_max, 0.8_dl,1e-4_dl,this%DebugLevel,2000)) then

            if (Minimize%NEWUOA(this, match_fde_zc, 2, 5, log_params,&
                0.8_dl,1e-4_dl,this%DebugLevel,500)) then

                if (Minimize%Last_bestfit > 1e-3) then
                    global_error_flag = error_darkenergy
                    global_error_message= 'TEarlyQuintessence ERROR converging solution for fde, zc'
                    write(*,*) 'last-bestfit= ', Minimize%Last_bestfit
                    return
                end if
                this%f = exp(log_params(1))
                this%m = exp(log_params(2))
                if (this%DebugLevel>0) then
                    call this%calc_zc_fde(fzero, xzero)
                    write(*,*) 'matched outputs Bobyqa zc, fde = ', fzero, xzero
                end if
            else
                global_error_flag = error_darkenergy
                global_error_message= 'TEarlyQuintessence ERROR finding solution for fde, zc'
                return
            end if
            if (this%DebugLevel>0) call Timer%WriteTime('Timing for parameter fitting')
        end if

        this%dloga = (-this%log_astart)/(this%npoints-1)

        !use log spacing in a up to max_a_log, then linear. Switch where step matches
        this%max_a_log = 1.d0/this%npoints/(exp(this%dloga)-1)
        npoints = (log(this%max_a_log)-this%log_astart)/this%dloga + 1

        if (allocated(this%phi_a)) then
            deallocate(this%phi_a,this%phidot_a)
            deallocate(this%ddphi_a,this%ddphidot_a, this%sampled_a)
        end if
        allocate(phi_a(npoints),phidot_a(npoints), sampled_a(npoints), fde(npoints))

        !initial_phi  = 10  !  0.3*grhom/m**3
        !initial_phi2 = 100 !   6*grhom/m**3
        !
        !!           initial_phi  = 65 !  0.3*grhom/m**3
        !!           initial_phi2 = 65 !   6*grhom/m**3
        !
        !astart=1d-9
        !
        !!See if initial conditions are giving correct omega_de now
        !atol=1d-8
        !initial_phidot =  astart*this%phidot_start(this%initial_phi)
        !om1= this%GetOmegaFromInitial(astart,initial_phi,initial_phidot, atol)
        !
        !print*, State%omega_de, 'first trial:', om1
        !if (abs(om1-State%omega_de > this%omega_tol)) then
        !    !if not, do binary search in the interval
        !    OK=.false.
        !    initial_phidot = astart*this%phidot_start(initial_phi2)
        !    om2= this%GetOmegaFromInitial(astart,initial_phi2,initial_phidot, atol)
        !    if (om1 > State%omega_de .or. om2 < State%omega_de) then
        !        write (*,*) 'initial phi values must bracket required value.  '
        !        write (*,*) 'om1, om2 = ', real(om1), real(om2)
        !        stop
        !    end if
        !    do iter=1,100
        !        deltaphi=initial_phi2-initial_phi
        !        phi =initial_phi + deltaphi/2
        !        initial_phidot =  astart*Quint_phidot_start(phi)
        !        om = this%GetOmegaFromInitial(astart,phi,initial_phidot,atol)
        !        if (om < State%omega_de) then
        !            om1=om
        !            initial_phi=phi
        !        else
        !            om2=om
        !            initial_phi2=phi
        !        end if
        !        if (om2-om1 < 1d-3) then
        !            OK=.true.
        !            initial_phi = (initial_phi2+initial_phi)/2
        !            if (FeedbackLevel > 0) write(*,*) 'phi_initial = ',initial_phi
        !            exit
        !        end if
        !
        !    end do !iterations
        !    if (.not. OK) stop 'Search for good intial conditions did not converge' !this shouldn't happen
        !
        !end if !Find initial

        initial_phi = this%theta_i*this%f

        y(1)=initial_phi
        initial_phidot =  this%astart*this%phidot_start(initial_phi)
        y(2)= initial_phidot*this%astart**2

        phi_a(1)=y(1)
        phidot_a(1)=y(2)/this%astart**2
        sampled_a(1)=this%astart
        da_osc = 1
        last_a = this%astart
        max_ix =0

        ind=1
        afrom=this%log_astart
        do i=1, npoints-1
            aend = this%log_astart + this%dloga*i
            ix = i+1
            sampled_a(ix)=exp(aend)
            a2 = sampled_a(ix)**2
            call dverk(this,NumEqs,EvolveBackgroundLog,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)
            if (.not. this%check_error(exp(afrom), exp(aend))) return
            call EvolveBackgroundLog(this,NumEqs,aend,y,w(:,1))
            phi_a(ix)=y(1)
            phidot_a(ix)=y(2)/a2
            if (i==1) then
                lastsign = y(2)
            elseif (y(2)*lastsign < 0) then
                !derivative has changed sign. Use to probe any oscillation scale:
                da_osc = min(da_osc, exp(aend) - last_a)
                last_a = exp(aend)
                lastsign= y(2)
            end if

            !Define fde as ratio of early dark energy density to total
            fde(ix) = 1/((this%state%grho_no_de(sampled_a(ix)) +  this%frac_lambda0*this%State%grhov*a2**2) &
                /(a2*(0.5d0* phidot_a(ix)**2 + a2*this%Vofphi(y(1),0))) + 1)
            if (max_ix==0 .and. ix > 2 .and. fde(ix)< fde(ix-1)) then
                max_ix = ix-1
            end if
            if (sampled_a(ix)*(exp(this%dloga)-1)*this%min_steps_per_osc > da_osc) then
                !Step size getting too big to sample oscillations well
                exit
            end if
        end do

        ! Do remaining steps with linear spacing in a, trying to be small enough
        this%npoints_log = ix
        this%max_a_log = sampled_a(ix)
        this%da = min(this%max_a_log *(exp(this%dloga)-1), &
            da_osc/this%min_steps_per_osc, (1- this%max_a_log)/(this%npoints-this%npoints_log))
        this%npoints_linear = int((1- this%max_a_log)/ this%da)+1
        this%da = (1- this%max_a_log)/this%npoints_linear

        tot_points = this%npoints_log+this%npoints_linear
        allocate(this%phi_a(tot_points),this%phidot_a(tot_points))
        allocate(this%ddphi_a(tot_points),this%ddphidot_a(tot_points))
        allocate(this%sampled_a(tot_points), this%fde(tot_points), this%ddfde(tot_points))
        this%sampled_a(1:ix) = sampled_a(1:ix)
        this%phi_a(1:ix) = phi_a(1:ix)
        this%phidot_a(1:ix) = phidot_a(1:ix)
        this%sampled_a(1:ix) = sampled_a(1:ix)
        this%fde(1:ix) = fde(1:ix)

        ind=1
        afrom = this%max_a_log
        do i=1, this%npoints_linear
            ix = this%npoints_log + i
            aend = this%max_a_log + this%da*i
            a2 =aend**2
            this%sampled_a(ix)=aend
            call dverk(this,NumEqs,EvolveBackground,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)
            if (.not. this%check_error(afrom, aend)) return
            call EvolveBackground(this,NumEqs,aend,y,w(:,1))
            this%phi_a(ix)=y(1)
            this%phidot_a(ix)=y(2)/a2

            this%fde(ix) = 1/((this%state%grho_no_de(aend) +  this%frac_lambda0*this%State%grhov*a2**2) &
                /(a2*(0.5d0* this%phidot_a(ix)**2 + a2*this%Vofphi(y(1),0))) + 1)
            if (max_ix==0 .and. this%fde(ix)< this%fde(ix-1)) then
                max_ix = ix-1
            end if
        end do

        call spline(this%sampled_a,this%phi_a,tot_points,splZero,splZero,this%ddphi_a)
        call spline(this%sampled_a,this%phidot_a,tot_points,splZero,splZero,this%ddphidot_a)
        call spline(this%sampled_a,this%fde,tot_points,splZero,splZero,this%ddfde)
        has_peak = .false.
        if (max_ix >0) then
            ix = max_ix
            has_peak = this%fde_peak(a_c, this%sampled_a(ix), this%sampled_a(ix+1), this%fde(ix), &
                this%fde(ix+1), this%ddfde(ix), this%ddfde(ix+1))
            if (.not. has_peak) then
                has_peak = this%fde_peak(a_c, this%sampled_a(ix-1), this%sampled_a(ix), &
                    this%fde(ix-1), this%fde(ix), this%ddfde(ix-1), this%ddfde(ix))
            end if
        end if
        if (has_peak) then
            this%zc = 1/a_c-1
            this%fde_zc = this%fdeAta(a_c)
        else
            if (this%DebugLevel>0) write(*,*) 'TEarlyQuintessence: NO PEAK '
            this%zc = -1
        end if
        if (this%DebugLevel>0) then
            write(*,*) 'TEarlyQuintessence zc, fde used', this%zc, this%fde_zc
        end if
    end subroutine TEarlyQuintessence_Init

    logical function check_error(this, afrom, aend)
        class(TEarlyQuintessence) :: this
        real(dl) afrom, aend

        if (global_error_flag/=0) then
            write(*,*) 'TEarlyQuintessence error integrating', afrom, aend
            write(*,*) this%n, this%f, this%m, this%theta_i
            stop
            check_error = .false.
            return
        end if
        check_error= .true.
    end function check_error

    logical function fde_peak(this, peak, xlo, xhi, Flo, Fhi, ddFlo, ddFhi)
        class(TEarlyQuintessence) :: this
        real(dl), intent(out) :: peak
        real(dl) Delta
        real(dl), intent(in) :: xlo, xhi, ddFlo, ddFhi,Flo, Fhi
        real(dl) a, b, c, fac

        !See if derivative has zero in spline interval xlo .. xhi

        Delta = xhi - xlo

        a = 0.5_dl*(ddFhi-ddFlo)/Delta
        b = (xhi*ddFlo-xlo*ddFhi)/Delta
        c = (Fhi-Flo)/Delta+ Delta/6._dl*((1-3*xhi**2/Delta**2)*ddFlo+(3*xlo**2/Delta**2-1)*ddFhi)
        fac = b**2-4*a*c
        if (fac>=0) then
            fac = sqrt(fac)
            peak = (-b + fac)/2/a
            if (peak >= xlo .and. peak <= xhi) then
                fde_peak = .true.
                return
            else
                peak = (-b - fac)/2/a
                if (peak >= xlo .and. peak <= xhi) then
                    fde_peak = .true.
                    return
                end if
            end if
        end if
        fde_peak = .false.
    end function fde_peak

    function match_zc(this, logm)
        class(TEarlyQuintessence), intent(inout) :: this
        real(dl), intent(in) :: logm
        real(dl) match_zc, zc, fde_zc

        this%m = exp(logm)
        call this%calc_zc_fde(zc, fde_zc)
        match_zc = zc - this%zc

        end function match_zc

        function match_fde(this, logf)
        class(TEarlyQuintessence), intent(inout) :: this
        real(dl), intent(in) :: logf
        real(dl) match_fde, zc, fde_zc

        this%f = exp(logf)
        call this%calc_zc_fde(zc, fde_zc)
        match_fde = fde_zc - this%fde_zc
    end function match_fde

    function match_fde_zc(this, x)
        class(TEarlyQuintessence) :: this
        real(dl), intent(in) :: x(:)
        real(dl) match_fde_zc, zc, fde_zc

        this%f = exp(x(1))
        this%m = exp(x(2))
        call this%calc_zc_fde(zc, fde_zc)

        match_fde_zc = (log(this%fde_zc)-log(fde_zc))**2 + (log(zc)-log(this%zc))**2
        if (this%DebugLevel>1) then
            write(*,*) 'search f, m, zc, fde_zc, chi2', this%f, this%m, zc, fde_zc, match_fde_zc
        end if
    end function match_fde_zc

    subroutine calc_zc_fde(this, z_c, fde_zc)
        class(TEarlyQuintessence), intent(inout) :: this
        real(dl), intent(out) :: z_c, fde_zc
        real(dl) aend, afrom
        integer, parameter ::  NumEqs=2
        real(dl) c(24),w(NumEqs,9), y(NumEqs)
        integer ind, i, ix
        real(dl), parameter :: splZero = 0._dl
        real(dl) a_c
        real(dl) initial_phi, initial_phidot, a2
        real(dl), dimension(:), allocatable :: sampled_a, fde, ddfde
        integer npoints, max_ix
        logical has_peak
        real(dl) a0, b0, da

        ! Get z_c and f_de(z_c) where z_c is the redshift of (first) peak of f_de (de energy fraction)
        ! Do this by forward propagating until peak, then get peak values by cubic interpolation

        initial_phi = this%theta_i*this%f
        this%log_astart = log(this%astart)
        this%dloga = (-this%log_astart)/(this%npoints-1)

        npoints = (-this%log_astart)/this%dloga + 1
        allocate(sampled_a(npoints), fde(npoints), ddfde(npoints))

        y(1)=initial_phi
        initial_phidot =  this%astart*this%phidot_start(initial_phi)
        y(2)= initial_phidot*this%astart**2
        sampled_a(1)=this%astart
        max_ix =0
        ind=1
        afrom=this%log_astart
        do i=1, npoints-1
            aend = this%log_astart + this%dloga*i
            ix = i+1
            sampled_a(ix)=exp(aend)
            a2 = sampled_a(ix)**2
            call dverk(this,NumEqs,EvolveBackgroundLog,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)
            if (.not. this%check_error(exp(afrom), exp(aend))) return
            call EvolveBackgroundLog(this,NumEqs,aend,y,w(:,1))
            fde(ix) = 1/((this%state%grho_no_de(sampled_a(ix)) +  this%frac_lambda0*this%State%grhov*a2**2) &
                /((0.5d0*y(2)**2/a2 + a2**2*this%Vofphi(y(1),0))) + 1)
            if (max_ix==0 .and. ix > 2 .and. fde(ix)< fde(ix-1)) then
                max_ix = ix-1
            end if
            if (max_ix/=0 .and. ix > max_ix+4) exit
        end do

        call spline(sampled_a,fde,ix,splZero,splZero,ddfde)
        has_peak = .false.
        if (max_ix >0) then
            has_peak = this%fde_peak(a_c, sampled_a(max_ix), sampled_a(max_ix+1), fde(max_ix), &
                fde(max_ix+1), ddfde(max_ix), ddfde(max_ix+1))
            if (.not. has_peak) then
                has_peak = this%fde_peak(a_c, sampled_a(max_ix-1), sampled_a(max_ix), &
                    fde(max_ix-1), fde(max_ix), ddfde(max_ix-1), ddfde(max_ix))
            end if
        end if
        if (has_peak) then
            z_c = 1/a_c-1
            ix = int((log(a_c)-this%log_astart)/this%dloga)+1
            da = sampled_a(ix+1) - sampled_a(ix)
            a0 = (sampled_a(ix+1) - a_c)/da
            b0 = 1 - a0
            fde_zc=b0*fde(ix+1) + a0*(fde(ix)-b0*((a0+1)*ddfde(ix)+(2-a0)*ddfde(ix+1))*da**2/6._dl)
        else
            write(*,*) 'calc_zc_fde: NO PEAK'
            z_c = -1
            fde_zc = 0
        end if
    end subroutine calc_zc_fde

    function fdeAta(this,a)
        class(TEarlyQuintessence) :: this
        real(dl), intent(in) :: a
        real(dl) fdeAta, aphi, aphidot, a2

        call this%ValsAta(a, aphi, aphidot)
        a2 = a**2
        fdeAta = 1/((this%state%grho_no_de(a) +  this%frac_lambda0*this%State%grhov*a2**2) &
            /(a2*(0.5d0* aphidot**2 + a2*this%Vofphi(aphi,0))) + 1)
    end function fdeAta

    subroutine TEarlyQuintessence_ReadParams(this, Ini)
        use IniObjects
        class(TEarlyQuintessence) :: this
        class(TIniFile), intent(in) :: Ini
        call this%TDarkEnergyModel%ReadParams(Ini)
    end subroutine TEarlyQuintessence_ReadParams

    function TEarlyQuintessence_PythonClass()
        character(LEN=:), allocatable :: TEarlyQuintessence_PythonClass
        TEarlyQuintessence_PythonClass = 'EarlyQuintessence'
    end function TEarlyQuintessence_PythonClass

    subroutine TEarlyQuintessence_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TEarlyQuintessence), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TEarlyQuintessence_SelfPointer

    ! -----------------------------------------------------------------------
    ! Implementations for Monodromic Quintessence model
    ! -----------------------------------------------------------------------

    real(dl) function TMonodromicQuintessence_VofPhi(this, phi, deriv) result(VofPhi)
        !The input variable phi is sqrt(8*Pi*G)*psi
        !Returns (8*Pi*G)^(1-deriv/2)*d^{deriv}V(psi)/d^{deriv}psi evaluated at psi
        !return result is in 1/Mpc^2 units [so times (Mpc/c)^2 to get units in 1/Mpc^2]
        class(TMonodromicQuintessence), intent(in) :: this
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
            Vofphi = this%C*(-this%alpha*phi**(-this%alpha-1)*(1 - this%A*sin(this%nu*phi)) - phi**(-this%alpha)*this%A*this%nu*cos(this%nu*phi))
            ! print *, "phi = ", phi, "V' = ", Vofphi
            ! if (isnan(Vofphi)) stop "V' is NaN"
        else if (deriv == 2) then
            if (phi < 0) then
                ! print *, "WARNING: negative field value: this shouldn't happen for this model. Using the second phi value in the interpolation table instead"
                Vofphi = this%C*this%phi_a(2)**(-this%alpha-2)*(this%A*sin(this%nu*this%phi_a(2))*(this%phi_a(2)**2*this%nu**2 - this%alpha**2 - this%alpha) + 2._dl*this%A*this%nu*this%alpha*this%phi_a(2)*cos(this%nu*this%phi_a(2)) + this%alpha*(1._dl + this%alpha))
                return
            end if
            Vofphi = this%C*phi**(-this%alpha-2)*(this%A*sin(this%nu*phi)*(phi**2*this%nu**2 - this%alpha**2 - this%alpha) + 2._dl*this%A*this%nu*this%alpha*phi*cos(this%nu*phi) + this%alpha*(1._dl + this%alpha))

            if (isnan(Vofphi)) stop "V'' is NaN"
        end if
    end function TMonodromicQuintessence_VofPhi

    real(dl) function get_initial_phi(this, a) result(initial_phi)
        ! Initial conditions from tracking solution
        ! See equations (4) from https://arxiv.org/pdf/1709.01544.pdf
        class(TMonodromicQuintessence), intent(in) :: this
        real(dl), intent(in) :: a
        real(dl) :: H_ini, t_ini, p, grho_no_de_at_initial_a, phi_tilde_0

        grho_no_de_at_initial_a = this%State%grho_no_de(a)/a**4
        H_ini = this%State%CP%H0 * sqrt(grho_no_de_at_initial_a/this%State%grhocrit)
        t_ini = 1._dl/(2._dl*H_ini)
        p = 2._dl/(2._dl + this%alpha)
        phi_tilde_0 = ((p**2 + p)/(this%alpha*this%C))**(-1/(2+this%alpha))
        initial_phi = phi_tilde_0 * t_ini**p
    end function get_initial_phi

    subroutine TMonodromicQuintessence_Init(this, State)
        class(TMonodromicQuintessence), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State
        integer,  parameter :: NumEqs = 2, max_iters = 20
        integer,  parameter :: nsteps_linear = 1000, nsteps_log = 1000, nsteps = nsteps_log + nsteps_linear
        real(dl), parameter :: omega_de_tol = 1e-4
        real(dl), parameter :: splZero = 0._dl
        real(dl), parameter :: a_start = 1e-7, a_switch = 1e-3
        real(dl), parameter :: dloga = (log(a_switch) - log(a_start))/nsteps_log, da = (1._dl - a_switch)/nsteps_linear
        real(dl)            :: y(NumEqs), y_prime(NumEqs)
        real(dl)            :: omega_de_target, om, om1, om2
        real(dl)            :: C_1, C_2, new_C, a, loga, atol, initial_phi, initial_phidot, a_line, b_line, error
        real(dl)            :: phi, phidot
        integer             :: i
        Type(TTimer)        :: Timer

        !Make interpolation table, etc,
        !At this point massive neutrinos have been initialized
        !so grho_no_de can be used to get density and pressure of other components at scale factor a
        
        call this%TQuintessence%Init(State)
        select type(State)
        class is (CAMBdata)
            omega_de_target = State%Omega_de
        end select

        if (allocated(this%phi_a)) then
            print*, "WARNING: the interpolation table is already allocated. This shouldn't be happening but we are deallocating anyway"
            deallocate(&
                this%phi_a,      &
                this%phidot_a,   &
                this%ddphi_a,    &
                this%ddphidot_a, &
                this%sampled_a   &
            )
        end if

        allocate(&
            this%phi_a(nsteps),      &
            this%phidot_a(nsteps),   &
            this%ddphi_a(nsteps),    &
            this%ddphidot_a(nsteps), &
            this%sampled_a(nsteps)   &
        )
        
        ! Binary search for C
        C_1 = this%State%grhov * 0.5_dl
        C_2 = this%State%grhov * 1.3_dl
        print*, "Shooting for C with tentative values: ", C_1, C_2
        
        ! See if current C is giving correct omega_de now
        atol = 1d-8
        initial_phidot = 0d0
        this%C = C_1
        initial_phi = this%get_initial_phi(a_start)
        om1 = this%GetOmegaFromInitial(a_start, initial_phi, initial_phidot, atol)
        this%C = C_2
        initial_phi = this%get_initial_phi(a_start)
        om2 = this%GetOmegaFromInitial(a_start, initial_phi, initial_phidot, atol)
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
            initial_phi = this%get_initial_phi(a_start)
            om = this%GetOmegaFromInitial(a_start, initial_phi, initial_phidot, atol)
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
        y(2) = 0d0
        loga = log(a_start)
        do i = 1, nsteps_log
            call this%EvolveBackgroundLog(NumEqs, loga, y, y_prime)
            y(1) = y(1) + y_prime(1)*dloga
            y(2) = y(2) + y_prime(2)*dloga
            this%sampled_a(i) = exp(loga)
            this%phi_a(i) = y(1)
            this%phidot_a(i) = y(2)/this%sampled_a(i)**2
            loga = loga + dloga
        end do
        
        a = a_switch
        do i = 1, nsteps_linear
            call this%EvolveBackground(NumEqs, a, y, y_prime)
            y(1) = y(1) + y_prime(1)*da
            y(2) = y(2) + y_prime(2)*da
            this%sampled_a(nsteps_log + i) = a
            this%phi_a(nsteps_log + i) = y(1)
            this%phidot_a(nsteps_log + i) = y(2)/a**2
            a = a + da
        end do

        ! JVR NOTE: we need to deallocate phi_a, phidot_a, sampled_a
        ! this might be causing memory leaks and subsequent segmentation faults in original CAMB!
        ! deallocate(phi_a, phidot_a, sampled_a)
        
        ! Must set the fields
        this%astart = a_start
        this%npoints_linear = nsteps_linear
        this%npoints_log = nsteps_log
        this%da = da
        this%dloga = dloga
        this%log_astart = log(a_start)
        this%max_a_log = a_switch

        call spline(this%sampled_a, this%phi_a, nsteps, splZero, splZero, this%ddphi_a)
        call spline(this%sampled_a, this%phidot_a, nsteps, splZero, splZero, this%ddphidot_a)

    end subroutine TMonodromicQuintessence_Init

    subroutine TMonodromicQuintessence_ReadParams(this, Ini)
        use IniObjects
        class(TMonodromicQuintessence) :: this
        class(TIniFile), intent(in) :: Ini
        call this%TDarkEnergyModel%ReadParams(Ini)
    end subroutine TMonodromicQuintessence_ReadParams

    function TMonodromicQuintessence_PythonClass()
        character(LEN=:), allocatable :: TMonodromicQuintessence_PythonClass
        TMonodromicQuintessence_PythonClass = 'MonodromicQuintessence'
    end function TMonodromicQuintessence_PythonClass

    subroutine TMonodromicQuintessence_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TMonodromicQuintessence), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TMonodromicQuintessence_SelfPointer

    logical function TMonodromicQuintessence_check_error(this, afrom, aend) result (check_error)
        class(TMonodromicQuintessence) :: this
        real(dl) afrom, aend

        if (global_error_flag/=0) then
            write(*,*) 'MonodromicQuintessence error integrating', afrom, aend
            write(*,*) this%A, this%C, this%nu, this%alpha
            stop
            check_error = .false.
            return
        end if
        check_error = .true.
    end function TMonodromicQuintessence_check_error

    end module Quintessence
