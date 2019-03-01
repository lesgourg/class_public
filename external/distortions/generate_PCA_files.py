#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sciint
from numpy.linalg import norm as vector_norm
from numpy.linalg import eig as eigen_vals_vecs # eigh is different from eig

def PCA_string_to_array(line,delimiter=" "):
  line = line.replace("\n","")
  if delimiter is not "\t":
    line = line.replace("\t","")
  return np.array([float(x) for x in line.split(delimiter) if (x is not "" and x is not " ")])

def gen_PCA_files(detector_name, nu_min, nu_max, delta_nu, delta_I_c, N_PCA):

  readfile = "Greens_data.dat"

  DI_units = 2.70062634e-18
  x_to_nu = 56.7798

  # Read external file Greens_data.dat
  with open(readfile) as f:
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

    # Read (and ignore) T_ini,T_last and rho
    Greens_T_ini = PCA_string_to_array(f.readline())
    Greens_T_last = PCA_string_to_array(f.readline())
    Greens_rho = PCA_string_to_array(f.readline())

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
      Greens_G_th_Spline[index_x_old] = sciint.CubicSpline(Greens_z,Greens_G_th[:,index_x_old])

    # Define new z and x arrays
    z_min = 1.02e3
    z_max = 5.0e6
    z_size = 1000
    z_arr = np.logspace(np.log10(z_min),np.log10(z_max),z_size)
    Nz_arr = len(z_arr)

    x_size = (nu_max-nu_min)/delta_nu
    x_arr = np.linspace(nu_min/x_to_nu,nu_max/x_to_nu,x_size)
    Nx_arr = len(x_arr)

    # Define visibility function
    z_th = 1.98e6
    bb_vis = np.exp(-(z_arr/z_th)**2.5)

    # Initialize spectral shapes
    G_th = np.zeros((Nx_arr,Nz_arr))
    Gdist = np.zeros(Nx_arr)
    Ydist = np.zeros(Nx_arr)
    Mdist = np.zeros(Nx_arr)

    # Interpolate Green's function
    index_x_old = 0
    for index_x_new,x in enumerate(x_arr):
      try:
        # Find position in xarray
        while(x>Greens_x[index_x_old]):
          index_x_old += 1
        # Linear interpolation in x
        frac = (x-Greens_x[index_x_old])/(Greens_x[index_x_old+1]-Greens_x[index_x_old])

        # Cubic interpolation for all values of z
        lowx_vals = Greens_G_th_Spline[index_x_old](z_arr)
        highx_vals = Greens_G_th_Spline[index_x_old+1](z_arr)

        G_th[index_x_new,:] = (lowx_vals*(1.-frac)+highx_vals*frac)
        G_th[index_x_new,:] *= 1.0e-8*bb_vis # Units, and energy deposition
	#print x, G_th[index_x_new,400], Greens_x[index_x_old]
      except:
        raise ValueError("{} is not in the file range [{},{}] for file '{}'".format(x,Greens_x[0],Greens_x[-1],readfile))

      # Define spectral shapes
      Gdist[index_x_new] = (x**4*np.exp(x)/(np.exp(x)-1)**2)*DI_units*1.0e18
      Ydist[index_x_new] = Gdist[index_x_new]*(x/np.tanh(x/2.)-4.)
      Mdist[index_x_new] = Gdist[index_x_new]*(1./2.19229-1./x)

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
      # Normalize, and re-instate energy conservation at early times
      f_g[index_z]  = 4.*f_g[index_z] + (1.-bb_vis[index_z])
      f_mu[index_z] = f_mu[index_z]/1.401
      f_y[index_z]  = 4.*f_y[index_z]

    # Calculate non-normalized residual
    Residual = np.zeros((Nz_arr,Nz_arr))
    for index_x in range(Nx_arr):
      for index_z in range(Nz_arr):
        Residual[index_x,index_z] = G_th[index_x,index_z]-Gdist[index_x]*f_g[index_z]/4-Ydist[index_x]*f_y[index_z]/4-Mdist[index_x]*f_mu[index_z]*1.401

    # Calculate non-normalized fisher matrix
    Fisher = np.zeros((Nz_arr,Nz_arr))
    for index_za in range(Nz_arr):
      for index_zb in range(Nz_arr):
        Fisher[index_za,index_zb] = np.sum(Residual[:,index_za]*Residual[:,index_zb])

    # Normalize fisher matrix
    delta_ln_z = np.log(z_arr[1])-np.log(z_arr[0])
    normalization = (delta_ln_z/(delta_I_c*1.e-18))**2
    Fisher/=normalization

    # Solve eigenvalue problem
    eigvals,eigvecs = eigen_vals_vecs(Fisher)

    E_vecs = np.real(eigvecs[:,:N_PCA]).T
    S_vecs = np.zeros((N_PCA,Nx_arr))
    for index_pca in range(N_PCA):
      for index_x in range(Nx_arr):
        S_vecs[index_pca][index_x] = np.dot(E_vecs[index_pca],Residual[index_x,:])

    # Create output files
    form = "%.10e" #Output formatting

    # Write file for branching ratio (Evec)
    with open(detector_name+"_branching_ratios.dat","w") as brfile:
      brfile.write("# In the file there is: z, J_T, J_y, J_mu, E_i (i=1-{})\n".format(N_PCA))
      brfile.write("{} {}\n".format(Nz_arr,N_PCA))
      for index_z in range(Nz_arr):
        brfile.write((form+" ") % z_arr[index_z])
        brfile.write((form+" ") % f_g[index_z])
        brfile.write((form+" ") % f_y[index_z])
        brfile.write((form    ) % f_mu[index_z])
        for index_pca in range(N_PCA):
          brfile.write((" "+form) % E_vecs[index_pca][index_z])
        brfile.write("\n")

    # Write file for distortion shapes (Svec)
    with open(detector_name+"_distortions_shapes.dat","w") as dsfile:
      dsfile.write("# In the file there is: nu, G_T, Y_SZ, M_mu, S_i (i=1-{})\n".format(N_PCA))
      dsfile.write("{} {}\n".format(Nx_arr,N_PCA))
      for index_x in range(Nx_arr):
        dsfile.write((form+" ") % x_arr[index_x])
        dsfile.write((form+" ") % Gdist[index_x])
        dsfile.write((form+" ") % Ydist[index_x])
        dsfile.write((form    ) % Mdist[index_x])
        for index_pca in range(N_PCA):
          dsfile.write((" "+form) % S_vecs[index_pca][index_x])
        dsfile.write("\n")

    # Update list of detectors
    # Open and read already present list
#    with open("detector_list.dat") as detector_list:
      # Read the header first
#      header = True
#      while(header):
#        line = detector_list.readline()
#        if(line.startswith("#")):
#          continue
        # The first line of the header without the "#" is still part of the header
#        header=False

#     done = False
#     old_list = []
#     while(not done):
#      line = detector_list.readline()
#      if(not line):
#        done = True
#      else:
#        old_list.append(PCA_string_to_array(line))
#     print old_list


gen_PCA_files("PIXIE_TEST", 30., 1000., 15., 5.0e-26, 8)

