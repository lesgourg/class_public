//KLASS
#include"ClassEngine.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

using namespace std;


// example run: one can specify as a second argument a preicison file
int main(int argc,char** argv){

  const int l_max_scalars=1200;

  //CLASS config
  ClassParams pars;

  pars.add("100*theta_s",1.04);
  pars.add("omega_b",0.0220);
  pars.add("omega_cdm",0.1116);
  pars.add("A_s",2.42e-9);
  pars.add("n_s",.96);
  pars.add("tau_reio",0.09);

  pars.add("k_pivot",0.05);
  pars.add("YHe",0.25);
  pars.add("output","tCl,pCl,lCl"); //pol +clphi

  pars.add("l_max_scalars",l_max_scalars);
  pars.add("lensing",true); //note boolean

  pars.add("background_verbose",1);
  pars.add("thermodynamics_verbose",1);
  pars.add("perturbations_verbose",1);
  pars.add("transfer_verbose",1);
  pars.add("primordial_verbose",1);
  pars.add("spectra_verbose",1);
  pars.add("nonlinear_verbose",1);
  pars.add("lensing_verbose",1);

  ClassEngine* tKlass(0);

  try{
    //le calculateur de spectres
    if (argc==2){
      string pre=string(argv[1]);
      tKlass=new ClassEngine(pars,pre);
    }
    else{
      tKlass=new ClassEngine(pars);
    }

    //cout.precision( 16 );
    //tKlass->writeCls(cout);

    ofstream outfile;
    const char* outfile_name = "testKlass_Cl_lensed.dat";
    outfile.open(outfile_name, ios::out | ios::trunc );
    tKlass->writeCls(outfile);
    cout << "Cl's written in file " << outfile_name << endl;
  }
  catch (std::exception &e){
    cout << "GOSH" << e.what() << endl;
  }

  delete tKlass;

}
