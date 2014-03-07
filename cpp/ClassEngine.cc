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
#include<numeric>

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
  if (input_init(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,_errmsg) == _FAILURE_)
    throw invalid_argument(_errmsg);

  //proetction parametres mal defini
  for (size_t i=0;i<pars.size();i++){
    if (fc.read[i] !=_TRUE_) throw invalid_argument(string("invalid CLASS parameter: ")+fc.name[i]);
  }

  //calcul class
  computeCls();

  //cout <<"creating " << sp.ct_size << " arrays" <<endl;
  cl=new double[sp.ct_size];

  //printFC();

}


ClassEngine::ClassEngine(const ClassParams& pars,const string & precision_file): cl(0),dofree(true){

  struct file_content fc_precision;
  fc_precision.size = 0;
  //decode pre structure
  if (parser_read_file(const_cast<char*>(precision_file.c_str()),&fc_precision,_errmsg) == _FAILURE_){
    throw invalid_argument(_errmsg);
  }

  //pars
  struct file_content fc_input;
  fc_input.size = 0;

 //prepare fc par structure
  size_t n=pars.size();
  parser_init(&fc_input,n,_errmsg);
  //config
  for (size_t i=0;i<pars.size();i++){
    strcpy(fc_input.name[i],pars.key(i).c_str());
    strcpy(fc_input.value[i],pars.value(i).c_str());
  }

  //concatenate both
  if (parser_cat(&fc_input,&fc_precision,&fc,_errmsg) == _FAILURE_) throw invalid_argument(_errmsg);

  //parser_free(&fc_input);
  parser_free(&fc_precision);

  //input
  if (input_init(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,_errmsg) == _FAILURE_)
    throw invalid_argument(_errmsg);

  //proetction parametres mal defini
  for (size_t i=0;i<pars.size();i++){
    if (fc.read[i] !=_TRUE_) throw invalid_argument(string("invalid CLASS parameter: ")+fc.name[i]);
  }

  //calcul class
  computeCls();

  //cout <<"creating " << sp.ct_size << " arrays" <<endl;
  cl=new double[sp.ct_size];

  //printFC();


}



//--------------
// Destructor --
//--------------
ClassEngine::~ClassEngine()
{

  //printFC();
  dofree && freeStructs();

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
int ClassEngine::class(
                       struct file_content *pfc,
                       struct precision * ppr,
                       struct background * pba,
                       struct thermo * pth,
                       struct perturbs * ppt,
                       struct primordial * ppm,
                       struct nonlinear * pnl,
                       struct transfers * ptr,
                       struct spectra * psp,
                       struct lensing * ple,
                       struct output * pop,
                       ErrorMsg errmsg) {



  if (input_init(pfc,ppr,pba,pth,ppt,ptr,ppm,psp,pnl,ple,pop,errmsg) == _FAILURE_) {
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

  if (primordial_init(ppr,ppt,ppm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",ppm->error_message);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (nonlinear_init(ppr,pba,pth,ppt,ppm,pnl) == _FAILURE_)  {
    printf("\n\nError in nonlinear_init \n=>%s\n",pnl->error_message);
    primordial_free(&pm);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (transfer_init(ppr,pba,pth,ppt,pnl,ptr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",ptr->error_message);
    nonlinear_free(&nl);
    primordial_free(&pm);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (spectra_init(ppr,pba,ppt,ppm,pnl,ptr,psp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",psp->error_message);
    transfer_free(&tr);
    nonlinear_free(&nl);
    primordial_free(&pm);
    perturb_free(&pt);
    thermodynamics_free(&th);
    background_free(&ba);
    dofree=false;
    return _FAILURE_;
  }

  if (lensing_init(ppr,ppt,psp,pnl,ple) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",ple->error_message);
    spectra_free(&sp);
    transfer_free(&tr);
    nonlinear_free(&nl);
    primordial_free(&pm);
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
  return this->class_main(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,_errmsg);
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

// int
// ClassEngine::l_size(Engine::cltype t){
//   int lmax(-1);

//   switch(t)
//     {
//     case TT:
//       if (sp.has_tt==_TRUE_) lmax=sp.l_size[sp.index_ct_tt];
//       break;
//     case TE:
//       if (sp.has_te==_TRUE_) lmax=sp.l_size[sp.index_ct_te] ;
//       break;
//     case EE:
//       if (sp.has_ee==_TRUE_) lmax=sp.l_size[sp.index_ct_ee] ;
//       break;
//     case BB:
//       if (sp.has_bb==_TRUE_) lmax=sp.l_size[sp.index_ct_bb] ;
//       break;
//     case PP:
//       if (sp.has_pp==_TRUE_) lmax=sp.l_size[sp.index_ct_pp] ;
//       break;
//     case TP:
//       if (sp.has_tp==_TRUE_) lmax=sp.l_size[sp.index_ct_tp] ;
//       break;
//     case EP:
//       if (sp.has_ep==_TRUE_) lmax=sp.l_size[sp.index_ct_ep] ;
//       break;
//     }
//   return lmax;
// }



double
ClassEngine::getCl(Engine::cltype t,const long &l){

  if (!dofree) throw out_of_range("no Cl available because CLASS failed");

  if (output_total_cl_at_l(&sp,&le,&op,static_cast<double>(l),cl) == _FAILURE_){
    cerr << ">>>fail getting Cl type=" << (int)t << " @l=" << l <<endl;
    throw out_of_range(sp.error_message);
  }

  double zecl=-1;

  double tomuk=1e6*Tcmb();
  double tomuk2=tomuk*tomuk;

  switch(t)
    {
    case TT:
      (sp.has_tt==_TRUE_) ? zecl=tomuk2*cl[sp.index_ct_tt] : throw invalid_argument("no ClTT available");
      break;
    case TE:
      (sp.has_te==_TRUE_) ? zecl=tomuk2*cl[sp.index_ct_te] : throw invalid_argument("no ClTE available");
      break;
    case EE:
      (sp.has_ee==_TRUE_) ? zecl=tomuk2*cl[sp.index_ct_ee] : throw invalid_argument("no ClEE available");
      break;
    case BB:
      (sp.has_bb==_TRUE_) ? zecl=tomuk2*cl[sp.index_ct_bb] : throw invalid_argument("no ClBB available");
      break;
    case PP:
      (sp.has_pp==_TRUE_) ? zecl=cl[sp.index_ct_pp] : throw invalid_argument("no ClPhi-Phi available");
      break;
    case TP:
      (sp.has_tp==_TRUE_) ? zecl=tomuk*cl[sp.index_ct_tp] : throw invalid_argument("no ClT-Phi available");
      break;
    case EP:
      (sp.has_ep==_TRUE_) ? zecl=tomuk*cl[sp.index_ct_ep] : throw invalid_argument("no ClE-Phi available");
      break;
    }

  return zecl;

}
void
ClassEngine::getCls(const std::vector<unsigned>& lvec, //input
		      std::vector<double>& cltt,
		      std::vector<double>& clte,
		      std::vector<double>& clee,
		      std::vector<double>& clbb)
{
  cltt.resize(lvec.size());
  clte.resize(lvec.size());
  clee.resize(lvec.size());
  clbb.resize(lvec.size());

  for (size_t i=0;i<lvec.size();i++){
    try{
      cltt[i]=getCl(ClassEngine::TT,lvec[i]);
      clte[i]=getCl(ClassEngine::TE,lvec[i]);
      clee[i]=getCl(ClassEngine::EE,lvec[i]);
      clbb[i]=getCl(ClassEngine::BB,lvec[i]);
    }
    catch(exception &e){
      throw e;
    }
  }

}


bool
ClassEngine::getLensing(const std::vector<unsigned>& lvec, //input
		std::vector<double>& clpp    ,
		std::vector<double>& cltp  ,
		std::vector<double>& clep  ){


  clpp.resize(lvec.size());
  cltp.resize(lvec.size());
  clep.resize(lvec.size());

  for (size_t i=0;i<lvec.size();i++){
    try{
      clpp[i]=getCl(ClassEngine::PP,lvec[i]);
      cltp[i]=getCl(ClassEngine::TP,lvec[i]);
      clep[i]=getCl(ClassEngine::EP,lvec[i]);
    }
    catch(exception &e){
      cout << "plantage!" << endl;
      cout << __FILE__ << e.what() << endl;
      return false;
    }
  }
  return true;
}


void
ClassEngine::writeCls(std::ostream &of,int ttmax){

  vector<unsigned> lvec(ttmax-1,1);
  lvec[0]=2;
  partial_sum(lvec.begin(),lvec.end(),lvec.begin());

  vector<double> cltt,clte,clee,clbb,clpp,cltp,clep;
  bool hasLensing=false;
  try{
    getCls(lvec,cltt,clte,clee,clbb);
    hasLensing=getLensing(lvec,clpp,cltp,clep);
  }
  catch (std::exception &e){
    cout << "GIOSH" << e.what() << endl;
  }

  //cout.precision( 16 );
  for (size_t i=0;i<lvec.size();i++) {
    of << lvec[i] << "\t"
       << cltt[i] << "\t"
       << clte[i] << "\t"
       << clee[i] << "\t"
       << clbb[i];
    if (hasLensing){
      of << "\t" << clpp[i] << "\t" << cltp[i] << "\t" << clep[i];
    }
    of << "\n";
  }



}
