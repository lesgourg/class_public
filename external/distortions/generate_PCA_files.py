#!/usr/bin/env python

import numpy as np
import sys
import scipy.interpolate as sciint
from numpy.linalg import norm as vector_norm
from numpy.linalg import eigh as eigen_vals_vecs
import os
import matplotlib.pyplot as plt

# Read inputs
if(len(sys.argv)==14):
  sd_detector_name = sys.argv[1]
  sd_detector_nu_min = eval(sys.argv[2])
  sd_detector_nu_max = eval(sys.argv[3])
  sd_detector_nu_delta = eval(sys.argv[4])
  sd_detector_bin_number = eval(sys.argv[5])
  sd_z_min = eval(sys.argv[6])
  sd_z_max = eval(sys.argv[7])
  sd_z_size = eval(sys.argv[8])
  sd_detector_delta_Ic = eval(sys.argv[9])
  sd_PCA_size = eval(sys.argv[10])
  z_th = eval(sys.argv[11])
  DI_units = eval(sys.argv[12])   # = 2.70062634e-18
  x_to_nu = eval(sys.argv[13])    # = 56.7798
  has_noisefile = False
elif(len(sys.argv)==11):
  sd_detector_name = sys.argv[1]
  sd_external_path = sys.argv[2]
  sd_noisefile_name = sys.argv[3]
  sd_z_min = eval(sys.argv[4])
  sd_z_max = eval(sys.argv[5])
  sd_z_size = eval(sys.argv[6])
  sd_PCA_size = eval(sys.argv[7])
  z_th = eval(sys.argv[8])
  DI_units = eval(sys.argv[9])   # = 2.70062634e-18
  x_to_nu = eval(sys.argv[10])    # = 56.7798
  has_noisefile = True
else:
  raise Exception("generate_PCA_files.py received invalid input arguments")

def PCA_string_to_array(line,delimiter=" "):
  line = line.replace("\n","")
  if delimiter is not "\t":
    line = line.replace("\t","")
  return np.array([float(x) for x in line.split(delimiter) if (x is not "" and x is not " ")])

def read_noisefile(filename):
  with open(filename) as det_noise:
    header = True
    while(header):
      line = det_noise.readline()
      if(line.startswith("#")):
        continue
      header=False
    Nrows,Ncols = PCA_string_to_array(line)
    #Skip first line containing Nx,Ncols
    line = det_noise.readline()
    cols = []
    while(line):
      cols.append(PCA_string_to_array(line))
      line = det_noise.readline()
    cols = np.array(cols).T
    assert(int(Ncols)==len(cols))
    assert(int(Nrows)==len(cols[0]))
    return len(cols[0]),cols[0]/x_to_nu,cols[1]*1e-26

dir_path = os.path.dirname(os.path.realpath(__file__))

# Read external file Greens_data.dat
readfile = "Greens_data.dat"
with open(os.path.join(dir_path,readfile)) as f:
  # Read the header first
  header = True
  while(header):
    line = f.readline()
    if(line.startswith("#")):
      continue
    # The first line of the header without the "#" is still part of the header
    header=False

  # Read the first line specifying z
  Greens_z = PCA_string_to_array(f.readline())
  Greens_Nz = len(Greens_z)
  Greens_lnz = np.log(Greens_z+1.)

  # Read T_ini,T_last and rho
  Greens_T_ini = PCA_string_to_array(f.readline())
  Greens_T_last = PCA_string_to_array(f.readline())
  Greens_drho = PCA_string_to_array(f.readline())

  # Calculate the difference in Temperature
  Greens_dT = (Greens_T_last-Greens_T_ini)/Greens_T_ini

  # Read the rest of the file
  done = False
  Greens_data_full = []
  while(not done):
    line = f.readline()
    if(not line):
      done = True
    else:
      Greens_data_full.append(PCA_string_to_array(line))
  Greens_data_full = np.array(Greens_data_full).T

  # Seperate the rest of the data into x, Green(z,x) and the blackbody
  Greens_x = Greens_data_full[0]
  Greens_Nx = len(Greens_x)
  Greens_G_th = Greens_data_full[1:Greens_Nz+1]
  Greens_blackbody = Greens_data_full[Greens_Nz+1]

  # Spline Greens function for interpolation
  Greens_G_th_Spline = [None for index_x_old in range(Greens_Nx)]
  for index_x_old in range(Greens_Nx):
    Greens_G_th_Spline[index_x_old] = sciint.CubicSpline(Greens_lnz,Greens_G_th[:,index_x_old])

  # Spline Greens dT for interpolation
  Greens_T_ini_Spline = sciint.CubicSpline(Greens_lnz,Greens_T_ini)
  Greens_T_last_Spline = sciint.CubicSpline(Greens_lnz,Greens_T_last)
  Greens_dT_Spline = sciint.CubicSpline(Greens_lnz,Greens_dT)
  Greens_drho_Spline = sciint.CubicSpline(Greens_lnz,Greens_drho)

  # Define new z and x arrays
  Nz_arr = sd_z_size
  z_arr = np.logspace(np.log10(sd_z_min),np.log10(sd_z_max),Nz_arr)
  lnz_arr = np.log(z_arr+1.)

  if has_noisefile:
    Nx_arr,x_arr,deltaIc_arr = read_noisefile(os.path.join(sd_external_path,sd_noisefile_name))
  else:
    Nx_arr = sd_detector_bin_number+1
    x_arr = np.linspace(sd_detector_nu_min/x_to_nu,sd_detector_nu_max/x_to_nu,Nx_arr)

  # Define visibility function
  #bb_vis = np.exp(-(z_arr/2.021e6)**2.5)
  bb_vis = np.exp(-(z_arr/z_th)**2.5)

  # The Gth file of Chluba subtracts away some part of the G_T distortion into a shift from T_ini to T_last
  # Here we calculate backwards, and obtain the shift of f_g due to the internal dT
  df_g = Greens_dT_Spline(lnz_arr)/Greens_drho_Spline(lnz_arr)

  # Initialize spectral shapes
  G_th = np.zeros((Nx_arr,Nz_arr))
  DI_T_shift = np.zeros((Nx_arr,Nz_arr))
  Gdist = np.zeros(Nx_arr)
  Ydist = np.zeros(Nx_arr)
  Mdist = np.zeros(Nx_arr)

  # Interpolate Green's function
  index_x_old = 0
  for index_x_new,x in enumerate(x_arr):
    # Define spectral shapes
    Gdist[index_x_new] = (x**4*np.exp(x)/(np.exp(x)-1)**2)*DI_units*1.0e18
    Ydist[index_x_new] = Gdist[index_x_new]*(x/np.tanh(x/2.)-4.)
    Mdist[index_x_new] = Gdist[index_x_new]*(1./2.19229-1./x)

    x_s = Greens_T_ini_Spline(lnz_arr)/Greens_T_last_Spline(lnz_arr)*x
    x_z = x*lnz_arr/lnz_arr
    DI_T_shift[index_x_new,:] = DI_units*1.0e26*x_z**3.*(np.exp(-x_s)/(1.-np.exp(-x_s))-np.exp(-x_z)/(1.-np.exp(-x_z)))/Greens_drho_Spline(lnz_arr)

    try:
      # Find position in xarray
      while(x>Greens_x[index_x_old]):
        index_x_old += 1
      # Linear interpolation in x
      frac = (x-Greens_x[index_x_old])/(Greens_x[index_x_old+1]-Greens_x[index_x_old])

      # Cubic interpolation for all values of z
      lowx_vals = Greens_G_th_Spline[index_x_old](lnz_arr)
      highx_vals = Greens_G_th_Spline[index_x_old+1](lnz_arr)

      G_th[index_x_new,:] = (lowx_vals*(1.-frac)+highx_vals*frac)
      G_th[index_x_new,:] *= bb_vis*1.e-8
      G_th[index_x_new,:] += DI_T_shift[index_x_new,:]*1.e-8
      #G_th[index_x_new,:] += Gdist[index_x_new]*df_g
    except:
      raise ValueError("{} is not in the file range [{},{}] for file '{}'".format(x,Greens_x[0],Greens_x[-1],readfile))

  # Begin orthonormlization
  # Y distortion
  e_Y = Ydist/vector_norm(Ydist)
  M_Y = np.dot(e_Y,Mdist)
  G_Y = np.dot(e_Y,Gdist)
  # Mu distortion
  Mperp = Mdist-M_Y*e_Y
  e_M = Mperp/vector_norm(Mperp)
  G_M = np.dot(e_M,Gdist)
  # G distortion
  Gperp = Gdist-G_Y*e_Y-G_M*e_M
  e_G = Gperp/vector_norm(Gperp)

  f_g = np.zeros(Nz_arr)
  f_mu = np.zeros(Nz_arr)
  f_y = np.zeros(Nz_arr)
  # Now, factorize G into orthonormal subspace
  for index_z in range(Nz_arr):
    # Compute non-normalized components
    f_g[index_z]  = (np.dot(G_th[:,index_z],e_G))/vector_norm(Gperp)
    f_mu[index_z] = (np.dot(G_th[:,index_z],e_M)-G_M*f_g[index_z])/vector_norm(Mperp)
    f_y[index_z]  = (np.dot(G_th[:,index_z],e_Y)-M_Y*f_mu[index_z]-G_Y*f_g[index_z])/vector_norm(Ydist)

  # Now we can re-normalize our functions and add the shift
  J_g = 4.*f_g
  J_mu = f_mu/1.401
  J_y  = 4.*f_y

  # Calculate non-normalized residual
  Residual = np.zeros((Nx_arr,Nz_arr))
  for index_x in range(Nx_arr):
    for index_z in range(Nz_arr):
      Residual[index_x,index_z] = G_th[index_x,index_z]-Gdist[index_x]*f_g[index_z]-Ydist[index_x]*f_y[index_z]-Mdist[index_x]*f_mu[index_z]

  # Calculate Fisher matrix
  Fisher = np.zeros((Nz_arr,Nz_arr))
  delta_ln_z = np.log(z_arr[1])-np.log(z_arr[0])
  for index_za in range(Nz_arr):
    for index_zb in range(Nz_arr):
      if has_noisefile:
        Fisher[index_za,index_zb] = np.sum(Residual[:,index_za]*Residual[:,index_zb]*pow(delta_ln_z/deltaIc_arr[:]*1.e8,2.))
      else:
        Fisher[index_za,index_zb] = np.sum(Residual[:,index_za]*Residual[:,index_zb]*pow(delta_ln_z/sd_detector_delta_Ic*1.e8,2.))

  # Solve eigenvalue problem
  eigvals,eigvecs = eigen_vals_vecs(Fisher)
  eigvals = eigvals[::-1]
  eigvecs = eigvecs[:,::-1]

  E_vecs = np.real(eigvecs[:,:sd_PCA_size]).T
  S_vecs = np.zeros((sd_PCA_size,Nx_arr))
  for index_pca in range(sd_PCA_size):
    for index_x in range(Nx_arr):
      S_vecs[index_pca][index_x] = np.dot(E_vecs[index_pca],Residual[index_x,:]*delta_ln_z)

  # Create output files
  form = "%.6e" #Output formatting

  # Write file for branching ratio (Evec)
  with open(os.path.join(dir_path,sd_detector_name+"_branching_ratios.dat"),"w") as brfile:
    brfile.write("# In the file there is: z, J_T, J_y, J_mu, E_i (i=1-{})\n".format(sd_PCA_size))
    brfile.write("# The first line contains the number of lines and the number of columns.\n".format(sd_PCA_size))
    brfile.write("{} {}\n".format(Nz_arr,sd_PCA_size))
    for index_z in range(Nz_arr):
      brfile.write((form+" ") % z_arr[index_z])
      brfile.write((form+" ") % f_g[index_z])
      brfile.write((form+" ") % f_y[index_z])
      brfile.write((form    ) % f_mu[index_z])
      for index_pca in range(sd_PCA_size):
        brfile.write((" "+form) % E_vecs[index_pca][index_z])
      brfile.write("\n")

  # Write file for distortion shapes (Svec)
  with open(os.path.join(dir_path,sd_detector_name+"_distortions_shapes.dat"),"w") as dsfile:
    dsfile.write("# In the file there is: nu, G_T, Y_SZ, M_mu, S_i (i=1-{})\n".format(sd_PCA_size))
    dsfile.write("# The first line contains the number of lines and the number of columns.\n".format(sd_PCA_size))
    dsfile.write("{} {}\n".format(Nx_arr,sd_PCA_size))
    for index_x in range(Nx_arr):
      dsfile.write((form+" ") % (x_arr[index_x]*x_to_nu))
      dsfile.write((form+" ") % Gdist[index_x])
      dsfile.write((form+" ") % Ydist[index_x])
      dsfile.write((form    ) % Mdist[index_x])
      for index_pca in range(sd_PCA_size):
        dsfile.write((" "+form) % S_vecs[index_pca][index_x])
      dsfile.write("\n")

  # Update list of detectors
  # Open and read already present list
  with open(os.path.join(dir_path,"detectors_list.dat"),"a") as detector_file:
    if has_noisefile:
      detector_file.write('%s  %s\n' % (sd_detector_name, sd_noisefile_name))
    else:
      detector_file.write('%s  %.6e  %.6e  %.6e  %i  %.6e\n' % (sd_detector_name, sd_detector_nu_min, sd_detector_nu_max, sd_detector_nu_delta, sd_detector_bin_number, sd_detector_delta_Ic))


