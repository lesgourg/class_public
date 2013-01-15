
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
  
  pars.add("H0",70.3);
  pars.add("omega_b",0.0220);
  //pars.add("omega_b",0.0215);
  pars.add("omega_cdm",0.1116);
  pars.add("A_s",2.42e-9);
  pars.add("n_s",.96);
  pars.add("z_reio",10.4);

  pars.add("k_pivot",0.002);
  pars.add("YHe",0.25);
  pars.add("output","tCl,pCl,lCl"); //pol +clphi

  pars.add("l_max_scalars",l_max_scalars);
  pars.add("lensing",false); //note boolean

  pars.add("modes","s"); //scalars
  pars.add("ic","ad"); //adiabatic


  try{
    //le calculateur de spectres
    ClassEngine* KKK(0);
    if (argc==2){
      string pre=string(argv[1]);
      KKK=new ClassEngine(pars,pre);
    }
    else{
      KKK=new ClassEngine(pars);
    }
    
    cout.precision( 16 );

    for (int l=2;l<=l_max_scalars;l++) {
      try{
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
      catch (std::exception &e){
	cout << "GIOSH" << e.what() << endl;
      }

    }
    delete KKK;
  }
  //class engine throws std:exceptions
  catch(exception& e){
    cerr << e.what() << endl;
  }

}
