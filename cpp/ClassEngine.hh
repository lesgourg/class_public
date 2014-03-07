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
//  updated to class v2.1 by Julien Lesgourgues (7 March 2014)
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
  ClassEngine(const ClassParams& pars);
  //with a class .pre file
  ClassEngine(const ClassParams& pars,const string & precision_file);


  // destructor
  ~ClassEngine();

  //modfiers: _FAILURE_ returned if CLASS pb:
  int updateParValues(const std::vector<double>& par);


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


  //may need that
  inline int numCls() const {return sp.ct_size;};
  inline double Tcmb() const {return ba.T_cmb;}

  //print content of file_content
  void printFC();
  void writeCls(std::ostream &o,int lmax);


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
                 struct primordial * ppm,
                 struct nonlinear * pnl,
                 struct transfers * ptr,
                 struct spectra * psp,
                 struct lensing * ple,
                 struct output * pop,
                 ErrorMsg errmsg);

protected:


};


;
#endif
