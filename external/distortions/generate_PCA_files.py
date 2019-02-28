import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sciint
from numpy.linalg import norm as vector_norm
from numpy.linalg import eig as eigen_vals_vecs #eigh is different from eig
def PCA_string_to_array(line,delimiter=" "):
  line = line.replace("\n","")
  if delimiter is not "\t":
    line = line.replace("\t","")
  return np.array([float(x) for x in line.split(delimiter) if (x is not "" and x is not " ")])

def gen_PCA_files(readfile,xarray,zarray,zth,DI_units,delta_I_c,N_PCA):
  with open(readfile) as f:
    #Read the header first
    header = True
    while(header):
      line = f.readline()
      if(line.startswith("#")):
        continue
      #The first line of the header without the "#" is still part of the header
      header=False

    #Read the first line specifying z
    greens_z = PCA_string_to_array(f.readline())
    Nz = len(greens_z)

    #Read (and ignore) T_ini,T_last and rho
    T_ini = PCA_string_to_array(f.readline())
    T_last = PCA_string_to_array(f.readline())
    rho = PCA_string_to_array(f.readline())

    #Read the rest of the file
    done = False
    greens_data_full = []
    while(not done):
      line = f.readline()
      if(not line):
        done = True
      else:
        greens_data_full.append(PCA_string_to_array(line))
    greens_data_full = np.array(greens_data_full).T

    #Seperate the rest of the data into x, Green(z,x) and the blackbody
    greens_x = greens_data_full[0]
    Nx = len(greens_x)
    greens_Green = greens_data_full[1:Nz+1]
    blackbody = greens_data_full[Nz+1]

    #Define visibility function
    bb_vis = np.exp(-(zarray/zth)**2.5)

    #Spline Greens function for interpolation
    greens_GreenSpline = [None for ioldx in range(Nx)]
    for ioldx in range(Nx):
      greens_GreenSpline[ioldx] = sciint.CubicSpline(greens_z,greens_Green[:,ioldx])

    #Define new arrays
    Nxnew = len(xarray)
    Nznew = len(zarray)
    newGreen = np.zeros((Nxnew,Nznew))
    Gdist = np.zeros(Nxnew)
    Ydist = np.zeros(Nxnew)
    Mdist = np.zeros(Nxnew)

    #Interpolate Greens function
    ioldx = 0
    for inewx,x in enumerate(xarray):
      try:
        #Find position in xarray
        while(x<greens_x[ioldx]):
          ioldx+=1
        #Linear interpolation in x
        frac = (x-greens_x[ioldx])/(greens_x[ioldx+1]-greens_x[ioldx])
        #Cubic interpolation for all values of z
        lowx_vals = greens_GreenSpline[ioldx](zarray)
        highx_vals = greens_GreenSpline[ioldx+1](zarray)
        newGreen[inewx,:] = (lowx_vals *(1.-frac) + highx_vals * frac)
        newGreen[inewx,:] *= 1.0e-8*bb_vis # Units, and energy deposition
      except:
        raise ValueError("{} is not in the file range [{},{}] for file '{}'".format(x,greens_x[0],greens_x[-1],readfile))
      #Define additional distortions
      Gdist[inewx] = (x**4*np.exp(x)/(np.exp(x)-1)**2)*DI_units*1.0e18
      Ydist[inewx] = Gdist[inewx]*(x/np.tanh(x/2.)-4.)
      Mdist[inewx] = Gdist[inewx]*(1./2.19229-1./x)

    #Begin orthonormlization
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


    f_g = np.zeros(Nznew)
    f_mu = np.zeros(Nznew)
    f_y = np.zeros(Nznew)
    #Now, factorize G into orthonormal subspace
    for iz in range(Nznew):
      #Compute non-normalized components
      f_g[iz]  = (np.dot(newGreen[:,iz],e_G))/vector_norm(Gperp)
      f_mu[iz] = (np.dot(newGreen[:,iz],e_M) - G_M*f_g[iz])/vector_norm(Mperp)
      f_y[iz]  = (np.dot(newGreen[:,iz],e_Y) - M_Y*f_mu[iz] - G_Y*f_g[iz])/vector_norm(Ydist)
      #Normalize, and re-instate energy conservation at early times
      f_g[iz]  = 4.*f_g[iz] + (1.-bb_vis[iz])
      f_mu[iz] = f_mu[iz]/1.401
      f_y[iz]  = 4.*f_y[iz]

    #Calculate non-normalized residual
    Residual = np.zeros((Nxnew,Nznew))
    for ix in range(Nxnew):
      for iz in range(Nznew):
        Residual[ix,iz] = newGreen[ix,iz]-Gdist[ix]*f_g[iz]/4-Ydist[ix]*f_y[iz]/4-Mdist[ix]*f_mu[iz]*1.401

    #Calculate non-normalized fisher matrix
    Fisher = np.zeros((Nznew,Nznew))
    for iza in range(Nznew):
      for izb in range(Nznew):
        Fisher[iza,izb] = np.sum(Residual[:,iza]*Residual[:,izb])

    #Normalize fisher matrix
    delta_ln_z = np.log(zarray[1])-np.log(zarray[0])
    normalization = (delta_ln_z/delta_I_c)**2
    Fisher/=normalization

    #Solve eigenvalue problem
    eigvals,eigvecs = eigen_vals_vecs(Fisher)

    Evecs = np.real(eigvecs[:,:N_PCA]).T
    Svecs = np.zeros((N_PCA,Nxnew))
    for npca in range(N_PCA):
      for ix in range(Nxnew):
        Svecs[npca][ix] = np.dot(Evecs[npca],Residual[ix,:])

    form = "%.10e" #Output formatting

    #Write file for branching ratio (Evec)
    with open("DET_branching_ratios.dat","w") as brfile:
      brfile.write("# In the file there is: z, J_T, J_y, J_mu, E_i (i=1-{})\n".format(N_PCA))
      brfile.write("{} {}\n".format(Nznew,N_PCA))
      for iz in range(Nznew):
        brfile.write((form+" ") % zarray[iz])
        brfile.write((form+" ") % f_g[iz])
        brfile.write((form+" ") % f_y[iz])
        brfile.write((form    ) % f_mu[iz])
        for npca in range(N_PCA):
          brfile.write((" "+form) % Evecs[npca][iz])
        brfile.write("\n")

    #Write file for distortion shapes (Svec)
    with open("DET_distortions_shapes.dat","w") as dsfile:
      dsfile.write("# In the file there is: nu, G_T, Y_SZ, M_mu, S_i (i=1-{})\n".format(N_PCA))
      dsfile.write("{} {}\n".format(Nxnew,N_PCA))
      for ix in range(Nxnew):
        dsfile.write((form+" ") % xarray[ix])
        dsfile.write((form+" ") % Gdist[ix])
        dsfile.write((form+" ") % Ydist[ix])
        dsfile.write((form    ) % Mdist[ix])
        for npca in range(N_PCA):
          dsfile.write((" "+form) % Svecs[npca][ix])
        dsfile.write("\n")

    print(Evecs)
    print(Svecs)
xnew = np.linspace(0.01,5.0,500)
znew = np.logspace(np.log10(1.02e3),np.log10(5.0e6),1000)
N_PCA = 4
gen_PCA_files("Greens_data.dat",xnew,znew,1.98e6,2.70062634e-18,5.0e-8,N_PCA)

"""

  /** Calculate unity vectors with relative dot products and moduli */
  M_y = 0.;
  G_y = 0.;
  for(index_x=0; index_x<psd->x_size; ++index_x){
    e_y[index_x] = Y[index_x]/mod_Y;
    M_y += e_y[index_x]*M[index_x];
    G_y += e_y[index_x]*G[index_x];
  }

  mod_M_per_squared = 0.;
  for(index_x=0; index_x<psd->x_size; ++index_x){
    M_per[index_x] = M[index_x]-M_y*e_y[index_x];
    mod_M_per_squared += pow(M_per[index_x],2.);
  }
  mod_M_per = sqrt(mod_M_per_squared);

  G_mu = 0.;
  for(index_x=0; index_x<psd->x_size; ++index_x){
    e_mu[index_x] = M_per[index_x]/mod_M_per;
    G_mu += e_mu[index_x]*G[index_x];
  }
  mod_G_per_squared = 0.;
  for(index_x=0; index_x<psd->x_size; ++index_x){
    G_per[index_x] = G[index_x]-G_y*e_y[index_x]-G_mu*e_mu[index_x]; 
    mod_G_per_squared += pow(G_per[index_x],2.);
  }
  mod_G_per = sqrt(mod_G_per_squared);
  for(index_x=0; index_x<psd->x_size; ++index_x){
    e_g[index_x] = G_per[index_x]/mod_G_per;
  }
  /** Calculate branching ratios */
  for(index_z=0; index_z<psd->z_size; ++index_z){
    bb_vis = exp(-pow(psd->z[index_z]/psd->z_th,2.5));

    e_g_dot_G_th = 0.;
    e_mu_dot_G_th = 0.;
    e_y_dot_G_th = 0.;
    for(index_x=0; index_x<psd->x_size; ++index_x){
      e_g_dot_G_th += e_g[index_x]*G_th[index_x][index_z];
      e_mu_dot_G_th += e_mu[index_x]*G_th[index_x][index_z];
      e_y_dot_G_th += e_y[index_x]*G_th[index_x][index_z]; 
    }

    f_g[index_z] = 4.*e_g_dot_G_th/mod_G_per;
    f_mu[index_z] = (e_mu_dot_G_th-G_mu*f_g[index_z]/4.)/mod_M_per/1.401;
    f_y[index_z] = (e_y_dot_G_th-1.401*M_y*f_mu[index_z]-G_y*f_g[index_z]/4.)*4./mod_Y;
    f_g[index_z] += 1.-bb_vis;
  }
"""
