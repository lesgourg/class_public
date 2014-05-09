
//KLASS
#include"ClassEngine.hh"

#include <iostream>
#include<string>
#include <stdexcept>

using namespace std;


// example run: one can specify as a second argument a preicison file
int main(int argc,char** argv){
  
  //jusqu'a ou en l
  const int l_max_scalars=1200;

  //CLASS config
  ClassParams pars;
  
  //pars.add("H0",70.3);
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


  ClassEngine* KKK(0);

  try{
    //le calculateur de spectres
    if (argc==2){
      string pre=string(argv[1]);
      KKK=new ClassEngine(pars,pre);
    }
    else{
      KKK=new ClassEngine(pars);
    }
    
    cout.precision( 16 );
    KKK->writeCls(cout);
  }
  catch (std::exception &e){
    cout << "GIOSH" << e.what() << endl;
  }

  delete KKK;

}
