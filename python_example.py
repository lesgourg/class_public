# import classy module
from classy import Class
# create instance of the class "Class"
LambdaCDM = Class()
# pass input parameters
LambdaCDM.set({’omega_b’:0.022032,’omega_cdm’:0.12038,’h’:0.67556,’A_s’:2.215e -9,’n_s’:0.9619,’tau_reio’:0.0925})
LambdaCDM.set({’output’:’tCl,pCl,lCl,mPk’,’lensing’:’yes’,’P_k_max_1/Mpc’:3.0}) # run class
LambdaCDM.compute()
# get all C_l output
cls = LambdaCDM.lensed_cl(2500)
# To check the format of cls
cls.viewkeys()

ll = cls[’ell’][2:] clTT = cls[’tt’][2:] clEE = cls[’ee’][2:] clPP = cls[’pp’][2:

                                                                               
