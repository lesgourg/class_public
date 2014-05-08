//--------------------------------------------------------------------------
//
// Description:
// 	class Engine :
//base class for Boltzmann code
//
//
// Author List:
//	Stephane Plaszczynski (plaszczy@lal.in2p3.fr)
//
// History (add to end):
//	creation:   Tue Mar 13 15:28:50 CET 2012 
//
//------------------------------------------------------------------------

#ifndef Engine_hh
#define Engine_hh

#include<vector>
#include<ostream>

class Engine
{

public:

  enum cltype {TT=0,EE,TE,BB,PP,TP,EP}; //P stands for phi (lensing potential)

  //constructors
  Engine();

  //pure virtual:
  virtual bool updateParValues(const std::vector<double>& cosmopars)=0;

  // units = (micro-K)^2
  virtual void getCls(const std::vector<unsigned>& lVec, //input 
		      std::vector<double>& cltt, 
		      std::vector<double>& clte, 
		      std::vector<double>& clee, 
		      std::vector<double>& clbb)=0;
  

  virtual bool getLensing(const std::vector<unsigned>& lVec, //input 
			  std::vector<double>& clpp, 
			  std::vector<double>& cltp, 
			  std::vector<double>& clep)=0;


  virtual double z_drag() const=0;
  virtual double rs_drag() const =0;

  virtual double get_Dv(double z)=0;
  
  virtual double get_Da(double z)=0;
  virtual double get_sigma8(double z)=0;
  virtual double get_f(double z)=0;
  virtual double get_Fz(double z)=0;
  virtual double get_Az(double z)=0;
  virtual double get_Hz(double z)=0;

  virtual double getTauReio() const=0;

  // destructor
  virtual ~Engine(){};

  //write Cl model+lensing in ostream
  virtual void writeCls(std::ostream &o);
  inline int lmax() {return _lmax;}

protected:
  int _lmax;

};

#endif

