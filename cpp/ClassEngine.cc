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
#include<cassert>

//#define DBUG

using namespace std;

template<typename T> std::string str(const T &x){
  std::ostringstream os;
  os << x;
  return os.str();
}
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
ClassEngine::ClassEngine(const ClassParams& pars, bool verbose): cl(0),dofree(true){

  //prepare fp structure
  size_t n=pars.size();
  //
  parser_init(&fc,n,(char*)"pipo",_errmsg);

  //config
  for (size_t i=0;i<pars.size();i++){
    strcpy(fc.name[i],pars.key(i).c_str());
    strcpy(fc.value[i],pars.value(i).c_str());
    //store
    parNames.push_back(pars.key(i));
    //identify lmax
    if(verbose) cout << pars.key(i) << "\t" << pars.value(i) <<endl;
    if (pars.key(i)=="l_max_scalars") {
      istringstream strstrm(pars.value(i));
      strstrm >> _lmax;
    }
  }
  if( verbose ) cout << __FILE__ << " : using lmax=" << _lmax <<endl;
  // assert(_lmax>0); // this collides with transfer function calculations

    //input
  if (input_read_from_file(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&sd,&op,_errmsg) == _FAILURE_)
    throw invalid_argument(_errmsg);

  //proetction parametres mal defini
  for (size_t i=0;i<pars.size();i++){
    if (fc.read[i] !=_TRUE_) throw invalid_argument(string("invalid CLASS parameter: ")+fc.name[i]);
  }

  //calcul class
  computeCls();

  //cout <<"creating " << sp.ct_size << " arrays" <<endl;
  if( pt.has_cl_cmb_temperature || pt.has_cl_cmb_polarization || pt.has_cl_lensing_potential ){
    cl=new double[sp.ct_size];
  }

  //printFC();

}


ClassEngine::ClassEngine(const ClassParams& pars,const string & precision_file, bool verbose): cl(0),dofree(true){

  struct file_content fc_precision;
  fc_precision.size = 0;
  //decode pre structure
  if (parser_read_file(const_cast<char*>(precision_file.c_str()),&fc_precision,_errmsg) == _FAILURE_){
    throw invalid_argument(_errmsg);
  }

  //pars
  struct file_content fc_input;
  fc_input.size = 0;
  fc_input.filename=new char[1];
 //prepare fc par structure
  size_t n=pars.size();
  parser_init(&fc_input,n,(char*)"pipo",_errmsg);
  //config
  for (size_t i=0;i<pars.size();i++){
    strcpy(fc_input.name[i],pars.key(i).c_str());
    strcpy(fc_input.value[i],pars.value(i).c_str());
    if (pars.key(i)=="l_max_scalars") {
      istringstream strstrm(pars.value(i));
      strstrm >> _lmax;
    }
  }
  if( verbose ) cout << __FILE__ << " : using lmax=" << _lmax <<endl;
  assert(_lmax>0);



  //concatenate both
  if (parser_cat(&fc_input,&fc_precision,&fc,_errmsg) == _FAILURE_) throw invalid_argument(_errmsg);

  //parser_free(&fc_input);
  parser_free(&fc_precision);

  //input
  if (input_read_from_file(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&sd,&op,_errmsg) == _FAILURE_)
    throw invalid_argument(_errmsg);

  //proetction parametres mal defini
  for (size_t i=0;i<pars.size();i++){
    if (fc.read[i] !=_TRUE_) throw invalid_argument(string("invalid CLASS parameter: ")+fc.name[i]);
  }

  //calcul class
  computeCls();

  //cout <<"creating " << sp.ct_size << " arrays" <<endl;
  if( pt.has_cl_cmb_temperature || pt.has_cl_cmb_polarization || pt.has_cl_lensing_potential ){
    cl=new double[sp.ct_size];
  }
  //printFC();


}



//--------------
// Destructor --
//--------------
ClassEngine::~ClassEngine()
{

  //printFC();
  dofree && freeStructs();

  if( pt.has_cl_cmb_temperature || pt.has_cl_cmb_polarization || pt.has_cl_lensing_potential ){
    delete [] cl;
  }

}

//-----------------
// Member functions --
//-----------------
bool ClassEngine::updateParValues(const std::vector<double>& par){
  dofree && freeStructs();
  for (size_t i=0;i<par.size();i++) {
    double val=par[i];
    strcpy(fc.value[i],str(val).c_str());
    strcpy(fc.name[i],parNames[i].c_str());
#ifdef DBUG
    cout << "update par values #" << i << "\t" <<  val << "\t" << str(val).c_str() << endl;
#endif
  }
  int status=computeCls();
#ifdef DBUG
  cout << "update par status=" << status << " succes=" << _SUCCESS_ << endl;
#endif

  return (status==_SUCCESS_);
}

//print content of file_content
void ClassEngine::printFC() {
  printf("FILE_CONTENT SIZE=%d\n",fc.size);
  for (int i=0;i<fc.size;i++) printf("%d : %s = %s\n",i,fc.name[i],fc.value[i]);


}
int ClassEngine::class_main(
			    struct file_content *pfc,
			    struct precision * ppr,
			    struct background * pba,
			    struct thermo * pth,
			    struct perturbs * ppt,
			    struct transfers * ptr,
			    struct primordial * ppm,
			    struct spectra * psp,
			    struct nonlinear * pnl,
			    struct lensing * ple,
			    struct distortions * psd,
			    struct output * pop,
			    ErrorMsg errmsg) {


  if (input_read_from_file(pfc,ppr,pba,pth,ppt,ptr,ppm,psp,pnl,ple,psd,pop,errmsg) == _FAILURE_) {
    printf("\n\nError running input_read_from_file \n=>%s\n",errmsg);
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

  if (distortions_init(ppr,pba,pth,ppt,ppm,psd) == _FAILURE_) {
    printf("\n\nError in distortions_init \n=>%s\n",psd->error_message);
    lensing_free(&le);
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

  dofree=true;
  return _SUCCESS_;
}


int ClassEngine::computeCls(){

#ifdef DBUG
  cout <<"call computecls" << endl;
  //printFC();
#endif

  int status=this->class_main(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&sd,&op,_errmsg);
#ifdef DBUG
  cout <<"status=" << status << endl;
#endif
  return status;

}

int
ClassEngine::freeStructs(){


  if (distortions_free(&sd) == _FAILURE_) {
    printf("\n\nError in distortions_free \n=>%s\n",sd.error_message);
    return _FAILURE_;
  }

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

void ClassEngine::call_perturb_sources_at_tau(
                           int index_md,
                           int index_ic,
                           int index_tp,
                           double tau,
                           double * psource
                           ) {
  if( perturb_sources_at_tau( &pt, index_md, index_ic, index_tp, tau, psource ) == _FAILURE_){
    cerr << ">>>fail getting Tk type=" << (int)index_tp <<endl;
    throw out_of_range(pt.error_message);
  }
}

void
ClassEngine::getTk( double z,
   std::vector<double>& k,
   std::vector<double>& d_cdm,
   std::vector<double>& d_b,
   std::vector<double>& d_ncdm,
   std::vector<double>& d_tot,
   std::vector<double>& t_cdm,
   std::vector<double>& t_b,
  std::vector<double>& t_ncdm,
   std::vector<double>& t_tot )
{

  if (!dofree) throw out_of_range("no Tk available because CLASS failed");

  double tau;
  int index;
  //transform redshift in conformal time
  background_tau_of_z(&ba,z,&tau);

  if(log(tau) < pt.ln_tau[0]){
    cerr << "Asking sources at a z bigger than z_max_pk, something probably went wrong\n";
    throw out_of_range(pt.error_message);
  }

  double *pvecback=new double[ba.bg_size];
  background_at_tau(&ba,tau,long_info,inter_normal, &index, pvecback);
  double fHa = pvecback[ba.index_bg_f] * (pvecback[ba.index_bg_a]*pvecback[ba.index_bg_H]);
  delete[] pvecback;

  // copy transfer func data to temporary
  const size_t index_md = pt.index_md_scalars;
  d_cdm.assign( pt.k_size[index_md], 0.0 );
  d_b.assign( pt.k_size[index_md], 0.0 );
  d_ncdm.assign( pt.k_size[index_md], 0.0 );
  d_tot.assign( pt.k_size[index_md], 0.0 );
  t_cdm.assign( pt.k_size[index_md], 0.0 );
  t_b.assign( pt.k_size[index_md], 0.0 );
  t_ncdm.assign( pt.k_size[index_md], 0.0 );
  t_tot.assign( pt.k_size[index_md], 0.0 );

  if( pt.ic_size[index_md] > 1 ){
    cerr << ">>>have more than 1 ICs, will use first and ignore others" << endl;
  }

  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_delta_cdm, tau, &d_cdm[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_delta_b, tau, &d_b[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_delta_ncdm1, tau, &d_ncdm[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_delta_tot, tau, &d_tot[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_theta_b, tau, &t_b[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_theta_ncdm1, tau, &t_ncdm[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_theta_tot, tau, &d_tot[0]);

  //
  std::vector<double> h_prime(pt.k_size[index_md],0.0), eta_prime(pt.k_size[index_md],0.0);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_eta_prime, tau, &eta_prime[0]);
  call_perturb_sources_at_tau(index_md, 0, pt.index_tp_h_prime, tau, &h_prime[0]);

  // gauge trafo velocities, store k-vector
  for (int index_k=0; index_k<pt.k_size[index_md]; index_k++)
  {
    auto ak = pt.k[index_md][index_k];

    // write data to vectors
    k.push_back( ak );

    // use the conformal Newtonian gauge for velocities
    // not correct, but N-body gauge currently not implemented
    double alphak2 = (h_prime[index_k]+6*eta_prime[index_k])/2;

    t_cdm[index_k]  = (-alphak2) / fHa;
    t_b[index_k]    = (-alphak2 + t_b[index_k]) / fHa;
    t_ncdm[index_k] = (-alphak2 + t_ncdm[index_k]) / fHa;
    t_tot[index_k]  = (-alphak2 + t_tot[index_k]) / fHa;
  }
}

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


double ClassEngine::get_f(double z)
{
  double tau;
  int index;
  double *pvecback;
  //transform redshift in conformal time
  background_tau_of_z(&ba,z,&tau);

  //pvecback must be allocated
  pvecback=(double *)malloc(ba.bg_size*sizeof(double));

  //call to fill pvecback
  background_at_tau(&ba,tau,long_info,inter_normal, &index, pvecback);



  double f_z=pvecback[ba.index_bg_f];
#ifdef DBUG
  cout << "f_of_z= "<< f_z <<endl;
#endif
  return f_z;
}


double ClassEngine::get_sigma8(double z)
{
  double tau;
  int index;
  double *pvecback;
  double sigma8 = 0.;
  //transform redshift in conformal time
  background_tau_of_z(&ba,z,&tau);

  //pvecback must be allocated
  pvecback=(double *)malloc(ba.bg_size*sizeof(double));

  //call to fill pvecback
  background_at_tau(&ba,tau,long_info,inter_normal, &index, pvecback);
  //background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback);
  spectra_sigma(&ba,&pm,&sp,8./ba.h,z,&sigma8);

#ifdef DBUG
  cout << "sigma_8= "<< sigma8 <<endl;
#endif
  return sigma8;
}

// ATTENTION FONCTION BIDON - GET omegam ! -------------------
double ClassEngine::get_Az(double z)
{
  double Dv = get_Dv(z);
  // A(z)=100DV(z)sqrt(~mh2)/cz
  double omega_bidon = 0.12 ;
  double Az = 100.*Dv*sqrt(omega_bidon)/(3.e8*z); // is there speed of light somewhere ?
  return Az;
}
//      --------------------------

double ClassEngine::get_Dv(double z)
{
  double tau;
  int index;
  double *pvecback;
  //transform redshift in conformal time
  background_tau_of_z(&ba,z,&tau);

  //pvecback must be allocated
  pvecback=(double *)malloc(ba.bg_size*sizeof(double));

  //call to fill pvecback
  background_at_tau(&ba,tau,long_info,inter_normal, &index, pvecback);

  double H_z=pvecback[ba.index_bg_H];
  double D_ang=pvecback[ba.index_bg_ang_distance];
#ifdef DBUG
  cout << "H_z= "<< H_z <<endl;
  cout << "D_ang= "<< D_ang <<endl;
#endif
  double D_v;

  D_v=pow(D_ang*(1+z),2)*z/H_z;
  D_v=pow(D_v,1./3.);
#ifdef DBUG
  cout << D_v << endl;
#endif
  return D_v;
}

double ClassEngine::get_Fz(double z)
{
  double tau;
  int index;
  double *pvecback;
  //transform redshift in conformal time
  background_tau_of_z(&ba,z,&tau);

  //pvecback must be allocated
  pvecback=(double *)malloc(ba.bg_size*sizeof(double));

  //call to fill pvecback
  background_at_tau(&ba,tau,long_info,inter_normal, &index, pvecback);


  double H_z=pvecback[ba.index_bg_H];
  double D_ang=pvecback[ba.index_bg_ang_distance];
#ifdef DBUG
  cout << "H_z= "<< H_z <<endl;
  cout << "D_ang= "<< D_ang <<endl;
#endif
  double F_z = (1.+z) * D_ang * H_z /(3.e8) ; // is there speed of light somewhere ?
  return F_z;
}

double ClassEngine::get_Hz(double z)
{
  double tau;
  int index;
  double *pvecback;
  //transform redshift in conformal time
  background_tau_of_z(&ba,z,&tau);

  //pvecback must be allocated
  pvecback=(double *)malloc(ba.bg_size*sizeof(double));

  //call to fill pvecback
  background_at_tau(&ba,tau,long_info,inter_normal, &index, pvecback);


  double H_z=pvecback[ba.index_bg_H];

  return(H_z);

}


double ClassEngine::get_Da(double z)
{
  double tau;
  int index;
  double *pvecback;
  //transform redshift in conformal time
  background_tau_of_z(&ba,z,&tau);

  //pvecback must be allocated
  pvecback=(double *)malloc(ba.bg_size*sizeof(double));

  //call to fill pvecback
  background_at_tau(&ba,tau,long_info,inter_normal, &index, pvecback);


  double H_z=pvecback[ba.index_bg_H];
  double D_ang=pvecback[ba.index_bg_ang_distance];
#ifdef DBUG
  cout << "H_z= "<< H_z <<endl;
  cout << "D_ang= "<< D_ang <<endl;
#endif
  return D_ang;
}
