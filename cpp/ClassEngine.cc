//--------------------------------------------------------------------------
//
// Description:
// 	class ClassEngine : see header file (ClassEngine.hh) for description.
//
//------------------------------------------------------------------------
//-----------------------
// This Class's Header --
//-----------------------
#include "ClassEngine.hh"
// ------------------------
// Collaborating classes --
//-------------------------
// C++
//--------------------
#include<iostream>
#include <iomanip>
#include<string>
#include<cmath>
#include <stdexcept>
#include<sstream>

using namespace std;

template<typename T> std::string str(const T &x){
  std::ostringstream os;
  os << x;
  return os.str();
};
//specilization
template<> std::string str (const float &x){
  std::ostringstream os;
  os << setprecision(8) << x;
  return os.str();
}
template<> std::string str (const double &x){
  std::ostringstream os;
  os << setprecision(16) << x;
  return os.str();
}
template<> std::string str (const bool &x){
  { return x ? "yes" : "no"; }
}

template<> std::string str (const std::string &x) {return x;}

std::string str (const char* s){return string(s);}

//instanciations
template string str(const int &x);
template string str(const signed char &x);
template string str(const unsigned char &x);
template string str(const short &x);
template string str(const unsigned short &x);
template string str(const unsigned int &x);
template string str(const long &x);
template string str(const unsigned long &x);
template string str(const long long &x);
template string str(const unsigned long long &x);

//---------------
// Constructors --
//----------------
ClassEngine::ClassEngine(const ClassParams& pars): cl(0),dofree(true){

  //prepare fp structure
  size_t n=pars.size();
  //
  parser_init(&fc,n,_errmsg);
  
  //config
  for (size_t i=0;i<pars.size();i++){
    strcpy(fc.name[i],pars.key(i).c_str());
    strcpy(fc.value[i],pars.value(i).c_str());
  }

    //input
  if (input_init(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,_errmsg) == _FAILURE_) 
    throw invalid_argument(_errmsg);

  //init bessel
  if (bessel_init(&pr,&bs) == _FAILURE_) throw  invalid_argument(_errmsg);

  //calcul class
  computeCls();
  
  //cout <<"creating " << sp.ct_size << " arrays" <<endl;
  cl=new double[sp.ct_size];

  printFC();

}

//--------------
// Destructor --
//--------------
ClassEngine::~ClassEngine()
{

  printFC();
  dofree && freeStructs();

  //clean
  if (bessel_free(&bs) == _FAILURE_) throw domain_error(bs.error_message);

  delete [] cl;

}

//-----------------
// Member functions --
//-----------------
int ClassEngine::updateParValues(const std::vector<double>& par){
  dofree && freeStructs();
  for (size_t i=0;i<par.size();i++) {
    double val=par[i];
    strcpy(fc.value[i],str(val).c_str());
  }
  return computeCls();
}

//print content of file_content
void ClassEngine::printFC() {
  printf("FILE_CONTENT SIZE=%d\n",fc.size);
  for (int i=0;i<fc.size;i++) printf("%d : %s = %s\n",i,fc.name[i],fc.value[i]);


}
int ClassEngine::class_assuming_bessels_computed(
				    struct file_content *pfc,
				    struct precision * ppr,
				    struct background * pba,
				    struct thermo * pth,
				    struct perturbs * ppt,
				    struct bessels * pbs,
				    struct transfers * ptr,
				    struct primordial * ppm,
				    struct spectra * psp,
				    struct nonlinear * pnl,
				    struct lensing * ple,
				    struct output * pop,
				    ErrorMsg errmsg) {
  


  if (input_init(pfc,ppr,pba,pth,ppt,pbs,ptr,ppm,psp,pnl,ple,pop,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    dofree=false;
    return _FAILURE_;
  }

  if (background_init(ppr,pba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",pba->error_message);
    dofree=false;
    return _FAILURE_;
  }
    
  if (thermodynamics_init(ppr,pba,pth) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",pth->error_message);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (perturb_init(ppr,pba,pth,ppt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",ppt->error_message);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  /*
  if (bessel_init(ppr,pbs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n=>%s\n",pbs->error_message);
    dofree=false;
    return _FAILURE_;
  } 
  */

  if (transfer_init(ppr,pba,pth,ppt,pbs,ptr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",ptr->error_message);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (primordial_init(ppr,ppt,ppm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",ppm->error_message);
    transfer_free(&tr);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (spectra_init(ppr,pba,ppt,ptr,ppm,psp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",psp->error_message);
    primordial_free(&pm);
    transfer_free(&tr);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (nonlinear_init(ppr,pba,pth,ppt,pbs,ptr,ppm,psp,pnl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",pnl->error_message);
    spectra_free(&sp);
    primordial_free(&pm);
    transfer_free(&tr);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (lensing_init(ppr,ppt,psp,pnl,ple) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",ple->error_message);
    nonlinear_free(&nl);
    spectra_free(&sp);
    primordial_free(&pm);
    transfer_free(&tr);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }
  
  //fprintf(stderr,"%d %e %e %e\n",l,cl[l][0],cl[l][1],cl[l][2]);

  dofree=true;
  return _SUCCESS_;
}


int ClassEngine::computeCls(){

  //printFC();
  //new call
  return this->class_assuming_bessels_computed(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,_errmsg);
}

int
ClassEngine::freeStructs(){
  
  
  if (lensing_free(&le) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",le.error_message);
    return _FAILURE_;
  }
  
  if (nonlinear_free(&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }
  
  if (spectra_free(&sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }
    
  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }
  
  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (perturb_free(&pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;
}


double
ClassEngine::getCl(ClassEngine::cltype t,const long &l){

  if (!dofree) throw out_of_range("no Cl available because CLASS failed");

  if (output_total_cl_at_l(&sp,&le,&op,static_cast<double>(l),cl) == _FAILURE_)
    throw out_of_range(sp.error_message);

  double zecl=-1;

  
  switch(t)
    {
    case TT:
      (sp.has_tt==_TRUE_) ? zecl=cl[sp.index_ct_tt] : throw invalid_argument("no ClTT available");
      break;
    case TE:
      (sp.has_te==_TRUE_) ? zecl=cl[sp.index_ct_te] : throw invalid_argument("no ClTE available");
      break; 
    case EE:
      (sp.has_ee==_TRUE_) ? zecl=cl[sp.index_ct_ee] : throw invalid_argument("no ClEE available");
      break;
    case BB:
      (sp.has_bb==_TRUE_) ? zecl=cl[sp.index_ct_bb] : throw invalid_argument("no ClBB available");
      break;
    case PP:
      (sp.has_pp==_TRUE_) ? zecl=cl[sp.index_ct_pp] : throw invalid_argument("no ClPhi-Phi available");
      break;
    case TP:
      (sp.has_tp==_TRUE_) ? zecl=cl[sp.index_ct_tp] : throw invalid_argument("no ClT-Phi available");
      break;
    case EP:
      (sp.has_ep==_TRUE_) ? zecl=cl[sp.index_ct_ep] : throw invalid_argument("no ClE-Phi available");
      break;
    }
  
  //delete [] cl;
  return zecl*l*(l+1)/(2*M_PI)*1e12*T_cmb()*T_cmb();

}
