import numpy as np

from .common import time_at_z, print_warning, logConversion, H, nan_clean
from .__init__ import redshift, logEnergies, DarkAgesError
from scipy.integrate import odeint as solve_ode
from scipy.integrate import trapz

#particle_dict has structure 'particle_name':[mass,spin_label,spin_orientations, flavorfactor, antifactor, colorfactor]
particle_dict = {'gamma': [0,'1',2,1,1,1],
			  'neutrino': [0,'1/2 N',1,3,2,1],
			  'electron': [5.110e-6,'1/2 C',2,1,2,1],
				  'muon': [1.057e-1,'1/2 C',2,1,2,1],
				   'tau': [1.777,'1/2 C',2,1,2,1],
					'up': [2.2e-3,'1/2 C',2,1,2,3],
				  'down': [4.7e-3,'1/2 C',2,1,2,3],
				 'charm': [1.28,'1/2 C',2,1,2,3],
			   'strange': [9.6e-2,'1/2 C',2,1,2,3],
				   'top': [173.1,'1/2 C',2,1,2,3],
				'bottom': [4.18,'1/2 C',2,1,2,3],
			   'w-boson': [80.39,'1',3,1,2,1],
			   'z-boson': [91.19,'1',3,1,1,1],
				 'gluon': [6e-1,'1',2,8,1,1],
				 'higgs': [125.09,'0',1,1,1,1],
				   #'pi0': [1.350e-1,'0',1,1,1,1],
				  #'piCh': [1.396e-1,'0',1,1,2,1]
}

def get_spin(particle):
	spin_label = particle_dict.get(particle)[1]
	tmp_dict = {'0':0.0,
			'1/2 N':0.5,
			'1/2 C':0.5,
				'1':1.0}
	return tmp_dict.get(spin_label)

def PBH_spectrum_at_m( mass, *particles):
	ret = np.zeros((len(logEnergies)), dtype=np.float64)
	E = logConversion(logEnergies - 9)
	if ('ALL' in particles):
		particles = particle_dict.keys()
	if particles:
		for particle in particles:
			spin = get_spin(particle)
			fraction = fraction_at_M( mass, particle )
			ret[:] += fraction*np.asarray(np.vectorize(Primary_spectrum, excluded = ['PBH_mass','spin']).__call__( energy=E[:], PBH_mass=mass, spin=spin))
	else:
		raise DarkAgesError('There is no particles given')
	return ret

def PBH_spectrum( PBH_mass_ini, *particles):
	ret = np.zeros((len(logEnergies), len(redshift)), dtype=np.float64)
	mass_at_z = PBH_mass_at_z(PBH_mass_ini, z_start =1e4)[-1,:]
	for idx, mass in enumerate(mass_at_z):
		ret[:,idx] = PBH_spectrum_at_m( mass, *particles )
	return ret

def full_spectrum( PBH_mass_ini ):
	ret = np.zeros((3, len(logEnergies), len(redshift)), dtype=np.float64)
	E = logConversion(logEnergies - 9)
	ret[0,:,:] = PBH_spectrum( PBH_mass_ini, 'electron')
	ret[1,:,:] = PBH_spectrum( PBH_mass_ini, 'gamma')
	ret[2,:,:] = PBH_spectrum( PBH_mass_ini, 'ALL')
	return E, redshift, ret


def get_temperature_from_mass(mass):
	return 1.06e13 / mass

def get_mass_from_temperature(temperature):
	return 1.06e13 / temperature

def Primary_spectrum( energy, PBH_mass, spin ):
	if abs(PBH_mass) < np.inf and PBH_mass > 0.:
		PBH_temperature = get_temperature_from_mass(PBH_mass)
		Gamma = 27 * ( 6.70861e-39)**2 * (energy * (PBH_mass * 5.61e23))**2
		return (1/(2*np.pi)) * Gamma / (np.exp( energy / PBH_temperature ) - (-1)**int(np.round(2*spin)))
	else:
		return 0.

def fraction_at_M( PBH_mass, channel, *all_channels):
	if all_channels:
		return  PBH_F_of_M( PBH_mass, channel ) / PBH_F_of_M( PBH_mass, *all_channels )
	else:
		return  PBH_F_of_M( PBH_mass, channel ) / PBH_F_of_M( PBH_mass )

def PBH_F_of_M( PBH_mass, *particles ):
	# factor_dict has the structure 'spin_label' : [f, beta]
	factor_dict = {'0':[0.267, 2.66],
			   '1/2 N':[0.142, 4.53],
			   '1/2 C':[0.147, 4.53],
			       '1':[0.060, 6.04]}

	def single_contribution(PBH_mass, mass_of_particle, multiplicity, beta):
		if mass_of_particle > 0:
			equivalent_BH_mass = get_mass_from_temperature(mass_of_particle)
			return multiplicity * np.exp(-PBH_mass / (beta*equivalent_BH_mass) )
		else:
			return multiplicity

	ret = 0.
	if particles and ('ALL' not in particles):
		for particle in particles:
			if particle in particle_dict:
				values = particle_dict.get(particle)
				mass_of_particle = values[0]
				multiplicity = factor_dict[values[1]][0] * values[2] * values[3] * values[4] * values[5]
				beta = factor_dict[values[1]][1]
				ret += single_contribution(PBH_mass, mass_of_particle, multiplicity, beta)
			else:
				print_warning('The particle "{:s}" is not recognized'.format(particle))
	else:
		for particle, values in particle_dict.iteritems():
			mass_of_particle = values[0]
			multiplicity = factor_dict[values[1]][0] * values[2] * values[3] * values[4] * values[5]
			beta = factor_dict[values[1]][1]
			ret += single_contribution(PBH_mass, mass_of_particle, multiplicity, beta)

	return ret

def PBH_dMdt(PBH_mass, time, scale):
	if PBH_mass > 0:
		ret = -5.34e25*(scale**(-3))*PBH_F_of_M( scale*PBH_mass ) * (PBH_mass)**(-2)
		return ret
	else:
		return 0

def PBH_mass_at_t(initial_PBH_mass, t_start = 1e10):
	def jac(mass, time, scale):
		out = np.zeros((1,1))
		out[0,0] = 0 # partial_(dMdt) / partial_t
		#out[0] = -2*PBH_dMdt(mass, time, scale) / (mass*scale)
		return out

	times = np.zeros((redshift.shape[0]+2), dtype=np.float64)
	times[0] = t_start
	times[1:-1] = time_at_z(redshift)
	times[-1] = time_at_z(1.)
	times.sort()

	log10_ini_mass = np.log10(initial_PBH_mass)
	scale = 10**(np.floor(log10_ini_mass)+5)
	initial_PBH_mass *= 1/scale

	temp_t = 10**np.linspace(np.log10(times[0]), np.log10(times[-1]), 1e4)

	temp_mass, full_info = solve_ode(PBH_dMdt, initial_PBH_mass, temp_t, Dfun=jac, args=(scale,), full_output=1,mxstep=10000)
	PBH_mass_at_t = np.interp(times, temp_t, temp_mass[:,0])

	#PBH_mass_at_t, full_info = solve_ode(PBH_dMdt, initial_PBH_mass, times, args=(scale,), full_output=1,mxstep=10000)

	#print('\n\n')
	#for key, val in full_info.iteritems():
	#	print('{0}\t>>>>\t{1}'.format(key,val))
	#print('\n\n')

	out = np.array([times, scale*PBH_mass_at_t])
	for idx in xrange(out.shape[-1]):
		if out[1,idx] <= 0:
			out[1,idx:] = np.nan
			break

	return out[:,:]

def PBH_mass_at_z(initial_PBH_mass, redshift=redshift, z_start = 1e4):
	times =  time_at_z(redshift)
	temp_grid = PBH_mass_at_t(initial_PBH_mass, t_start = time_at_z(z_start))[:,1:-1]
	out = np.zeros((2,len(times)), dtype=np.float64)
	out[1,:] = nan_clean(np.interp(times, temp_grid[0,:], temp_grid[1,:]))
	out[0,:] = redshift
	return out
