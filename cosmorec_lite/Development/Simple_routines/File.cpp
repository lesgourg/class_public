//-----------------------------
// Author Jens Chluba October 2003
//-----------------------------
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "File.h"

using namespace std;

void inputFile::open(const char* name)
{
  nm=name;
  file.open(name, ios::in );
  // go to beginning of file
  file.seekg(0, ios::beg);
  
  if (!file)
    {
      cerr << " Fehler beim oeffnen des Files: " << name 
       << ". ueberpruefen Sie bitte den Pfad.\n\n "; 
      exit(0);
    }
  return;
}

void inputFile::getNchar(char* cstr, int N)
{
  int j=0;
  
  for (int i=0; i<N && file.eof()==0; ++i )
    {
      char ch;
      file.get(ch);
      
      if (!(ch==' '))
    {
      cstr[j] = ch;
      j++;
    }
      
      cstr[j] = '\0';
    }
}

void inputFile::read_xy(double x[], double y[], const long N)
{
  if(N>0)
    {         
      file.seekg(0, ios::beg);
      
      for(l=0; l<N && file.eof()==0; l++)
    {
      file >> x[l];
      file >> y[l];       
    }
      
      if(l!= N && file.eof()==1) 
    cout << " Class inputFile :: read_xy(): " 
         << "End of file was reached before arrays were completelly filled!" 
         << "\n Only " << l-1 << " points have been read." << endl; 
    }
  else 
    cout << " Class inputFile :: read_xy(): You requested 0 values to be read." << endl;
  
  return;
}

void inputFile::read_cols(double x[], double y[], int colx, int coly, int ncols, const long N)
{
  int col=0;
  double dummy;
  
  if(N>0)
    {         
      file.seekg(0, ios::beg);
      
      for(l=0; l<N && file.eof()==0; l++)
    {
      while(col < colx){ file >> x[l]; col++;}
      while(col < coly){ file >> y[l]; col++;}
      while(col < ncols){ file >> dummy; col++;}
      col=0;
    }
      
      if(l!= N && file.eof()==1) 
    cout << " Class inputFile :: read_cols(): " 
         << "End of file was reached before arrays were completelly filled!" 
         << "\n Only " << l-1 << " points have been read." << endl; 
    }
  else 
    cout << " Class inputFile :: read_cols(): You requested 0 values to be read." << endl;
  
  return;
}

void inputFile::read_col(double x[], int colx, int ncols, const long N)
{
  int col=0;
  double dummy;
  
  if(N>0)
    {         
      file.seekg(0, ios::beg);
      
      for(l=0; l<N && file.eof()==0; l++)
    {
      while(col < colx){ file >> x[l]; col++;}
      while(col < ncols){ file >> dummy; col++;}
      col=0;
    }
      
      if(l!= N && file.eof()==1) 
    cout << " Class inputFile :: read_col(): " 
         << "End of file was reached before arrays were completelly filled!" 
         << "\n Only " << l-1 << " points have been read." << endl; 
    }
  else 
    cout << " Class inputFile :: read_col(): You requested 0 values to be read." << endl;
  
  return;
}

// Public Funktionen....
void inputFile::read_x(double x[], const long N)
{
  double dummy;
  if(N>0)
    {         
      file.seekg(0, ios::beg);
      
      for(l=0; l<N && file.eof()==0; l++)
    {
      file >> x[l];
      file >> dummy;          
    }
      
      if(l!= N && file.eof()==1) 
    cout << " Class inputFile :: read_x(): " 
         << "End of file was reached before arrays were completelly filled!" 
         << "\n Only " << l-1 << " points have been read." << endl; 
    }
  else 
    cout << " Class inputFile :: read_x(): You requested 0 values to be read." << endl;
  
  return;
}

void inputFile::read_y(double y[], const long N)
{
  double dummy;
  if(N>0)
    {         
      file.seekg(0, ios::beg);
      
      for(l=0; l<N && file.eof()==0; l++)
    {
      file >> dummy;
      file >> y[l];       
    }
      
      if(l!= N && file.eof()==1) 
    cout << " Class inputFile :: read_y(): " 
         << "End of file was reached before arrays were completelly filled!" 
         << "\n Only " << l-1 << " points have been read." << endl; 
    }
  else 
    cout << " Class inputFile :: read_y(): You requested 0 values to be read." << endl;
  
  return;
}

double inputFile::get_valN(const long NR)
{
  r=0.0;
  
  if(NR>0)
    {         
      file.seekg(0, ios::beg);
      
      for(l=0; l<NR && file.eof()==0; l++) file >> r;
      
      if(l!=NR && file.eof()==1) 
    {
      cout << " Class inputFile :: get_valN(): " 
           << " End of file was reached before arrays were completelly filled! "
           << "\n Only " << l-1 << " points have been read." << endl; 
      
      r=0.0;
    }
    }
  else 
    cout << " Class inputFile :: get_valN(): Numbering starts at NR==1!"
     << " You requested 0. value to be read." << endl;
  
  return r; 
}

int inputFile::get_intN(const long NR)
{
  return (int)get_valN(NR);
  /*
    ir=0;
    
    if(NR>0)
    {         
    file.seekg(0, ios::beg);
    
    for(l=0; l<NR && file.eof()==0; l++) file >> ir;
    
    if(l!=NR && file.eof()==1) 
    {
    cout << " Class inputFile :: get_intN(): " 
    << " End of file was reached before arrays were completelly filled! "
    << "\n Only " << l-1 << " points have been read." << endl; 
    
    ir=0;
    }
    }
    else 
    cout << " Class inputFile :: get_intN(): Numbering starts at NR==1!"
    << " You requested 0. value to be read." << endl;
    
    return ir; 
  */
}

string inputFile::get_stringN(const long NR)
{
  str="";
  
  if(NR>0)
    {         
      file.seekg(0, ios::beg);
      
      for(l=0; l<NR && file.eof()==0; l++) file >> str;
      
      if(l!=NR && file.eof()==1) 
    {
      cout << " Class inputFile :: get_stringN(): " 
           << " End of file was reached before arrays were completelly filled! "
           << "\n Only " << l-1 << " points have been read." << endl; 
      
      str="";
    }
    }
  else 
    cout << " Class inputFile :: get_stringN(): Numbering starts at NR==1!"
     << " You requested 0. value to be read." << endl;
  
  return str; 
}

double inputFile::get_next()
{
  if(file.eof()==0) file >> r;
  else
    {
      cout << " Class inputFile :: get_next(): " 
       << " End of file was already reached! " << endl; 
      
      r=0.0;
    }
  
  return r; 
}

int inputFile::get_next_int()
{
  if(file.eof()==0) file >> ir;
  else
    {
      cout << " Class inputFile :: get_next_int(): " 
       << " End of file was already reached! " << endl; 
      
      ir=0;
    }
  
  return ir; 
}

char inputFile::get_next_char()
{
  if(file.eof()==0) file >> ch;
  else
    {
      cout << " Class inputFile :: get_next_char(): " 
       << " End of file was already reached! " << endl; 
      
      //ch;
    }
  
  return ch; 
}

string inputFile::get_next_string()
{
  if(file.eof()==0) file >> str;
  else
    {
      cout << " Class inputFile :: get_next_string(): " 
       << " End of file was already reached! " << endl; 
      
      str="";
    }
  
  return str; 
}

string inputFile::get_next_line()
{
  string dum;
  str="";

  if(file.eof()==0) getline(file, str);
  else
    {
      cout << " Class inputFile :: get_next_line(): " 
       << " End of file was already reached! " << nm<< endl; 
    }
  
  return str; 
}
