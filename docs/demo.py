import numpy as np
import matplotlib.pyplot as plt
import camb
print(f"Using CAMB {camb.__version__} installed at {camb.__file__}")

base_cosmology = camb.set_params(# Background
	H0=67.56, ombh2=0.02238280, omch2=0.1201075, TCMB=2.7255,
	# Dark Energy
	dark_energy_model = 'MonodromicQuintessence', alpha=0.2, A=0.05, nu=50,
	# Neutrinos
	num_nu_massless=3.044, num_nu_massive=0, nu_mass_numbers=[0],
	# Initial Power Spectrum
	As = 2.100549e-09, ns = 0.9660499, 
	YHe = 0.246, WantTransfer=True
)

base_results = camb.get_results(base_cosmology)