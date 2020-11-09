//--------------------------------------------------------------------------
//
// Description:
// 	class ClassEngine :
// encapsulation of class calls
//
//
// Author List:
//	Stephane Plaszczynski (plaszczy@lal.in2p3.fr)
//
// History (add to end):
//	creation:   ven. nov. 4 11:02:20 CET 2011
//
//-----------------------------------------------------------------------

#ifndef ClassEngine_hh
#define ClassEngine_hh

//CLASS
#include"class.h"

#include"Engine.hh"
//STD
#include<string>
#include<vector>
#include<utility>
#include<ostream>

using std::string;

//general utility to convert safely numerical types to string
template<typename T> std::string str(const T &x);
//specialisations
template<> std::string str (const float &x);
template<> std::string str (const double &x);
template<> std::string str (const bool &x); //"yes" or "no"
template<> std::string str (const std::string &x);

std::string str(const char* x);
//////////////////////////////////////////////////////////////////////////
//class to encapsulate CLASS parameters from any type (numerical or string)
class ClassParams{
public:

  ClassParams(){};
  ClassParams( const ClassParams& o):pars(o.pars){};

  //use this to add a CLASS variable
  template<typename T> unsigned add(const string& key,const T& val){
  pars.push_back(make_pair(key,str(val)));
  return pars.size();
  }

  //accesors
  inline unsigned size() const {return pars.size();}
  inline string key(const unsigned& i) const {return pars[i].first;}
  inline string value(const unsigned& i) const {return pars[i].second;}


private:
  std::vector<std::pair<string,string> > pars;
};

///////////////////////////////////////////////////////////////////////////
class ClassEngine : public Engine
{

  friend class ClassParams;

public:
  //constructors
  ClassEngine(const ClassParams& pars, bool verbose=true );
  //with a class .pre file
  ClassEngine(const ClassParams& pars, const string & precision_file, bool verbose=true);


  // destructor
  ~ClassEngine();

  //modfiers: _FAILURE_ returned if CLASS pb:
  bool updateParValues(const std::vector<double>& par);


  //get value at l ( 2<l<lmax): in units = (micro-K)^2
  //don't call if FAILURE returned previously
  //throws std::execption if pb

  double getCl(Engine::cltype t,const long &l);
  void getCls(const std::vector<unsigned>& lVec, //input
	      std::vector<double>& cltt,
	      std::vector<double>& clte,
	      std::vector<double>& clee,
	      std::vector<double>& clbb);
  bool getLensing(const std::vector<unsigned>& lVec, //input
	      std::vector<double>& clphiphi,
	      std::vector<double>& cltphi,
	      std::vector<double>& clephi);

  void call_perturb_sources_at_tau(
                           int index_md,
                           int index_ic,
                           int index_tp,
                           double tau,
                           double * psource
                           );

  void getTk( double z,
        std::vector<double>& k,
        std::vector<double>& d_cdm,
        std::vector<double>& d_b,
        std::vector<double>& d_ncdm,
        std::vector<double>& d_tot,
        std::vector<double>& t_cdm,
        std::vector<double>& t_b,
        std::vector<double>& t_ncdm,
        std::vector<double>& t_tot );

 //for BAO
  inline double z_drag() const {return th.z_d;}
  inline double rs_drag() const {return th.rs_d;}
  double get_Dv(double z);

  double get_Da(double z);
  double get_sigma8(double z);
  double get_f(double z);

  double get_Fz(double z);
  double get_Hz(double z);
  double get_Az(double z);

  double getTauReio() const {return th.tau_reio;}

  //may need that
  inline int numCls() const {return sp.ct_size;};
  inline double Tcmb() const {return ba.T_cmb;}

  inline int l_max_scalars() const {return _lmax;}

  //print content of file_content
  void printFC();

private:
  //structures class en commun
  struct file_content fc;
  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct distortions sd;      /* for spectral distortions */
  struct output op;           /* for output files */

  ErrorMsg _errmsg;            /* for error messages */
  double * cl;

  //helpers
  bool dofree;
  int freeStructs();

  //call once /model
  int computeCls();

  int class_main(
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
		 ErrorMsg errmsg);
  //parnames
  std::vector<std::string> parNames;

protected:


};

#endif

