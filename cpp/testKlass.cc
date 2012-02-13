
//KLASS
#include"ClassEngine.hh"

#include <iostream>
#include<string>

using namespace std;

int main(){
  
  //jusqu'a ou en l
  const int lmax_scalars=2000;

  //CLASS config
  ClassParams pars;
  
  pars.add("H0",72);
  pars.add("omega_b",2e-2);
  pars.add("omega_cdm",0.10);
  pars.add("A_s",2.8e-9);
  pars.add("n_s",1.);
  pars.add("tau_reio",0.08);

  pars.add("k_pivot",0.002);
  pars.add("YHe",0.25);
  pars.add("output","tCl,pCl,lCl"); //pol +clphi

  pars.add("lmax_scalars",lmax_scalars);
  pars.add("lensing",false); //note boolean

  pars.add("modes","s"); //scalars
  pars.add("ic","ad"); //adiabatic

  try{
    //le calculateur de spectres
    ClassEngine* KKK=new ClassEngine(pars);
    
    cout.precision( 16 );
    for (int l=2;l<=lmax_scalars;l++) {
      cout << l << "\t" 
	   << KKK->getCl(ClassEngine::TT,l) << "\t" 
	   << KKK->getCl(ClassEngine::TE,l) << "\t" 
	   << KKK->getCl(ClassEngine::EE,l) << "\t" 
	   << KKK->getCl(ClassEngine::BB,l)<< "\t" 
	   << KKK->getCl(ClassEngine::PP,l)<< "\t" 
	   << KKK->getCl(ClassEngine::TP,l)<< "\t"
	   << KKK->getCl(ClassEngine::EP,l)  
	   << endl;
    }
    delete KKK;
  }
  //class engine throws std:exceptions
  catch(exception& e){
    cerr << e.what() << endl;
  }

}
