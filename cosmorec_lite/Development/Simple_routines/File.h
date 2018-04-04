#ifndef _File_
#define _File_

#include <stdlib.h>
#include <fstream>
#include <string>

using namespace std;

class inputFile
{
private:
    ifstream file;
    long l, ir;
    double r;
    string str, nm;
    char ch;
    
    void getNchar(char* cstr, int N);
    
public:
    // Konstruktoren der Klasse File....
    inputFile (const string& name){ nm=name; open(name);}
    inputFile (const char* name){ nm=name; open(name);}
    inputFile (){};
    ~inputFile (){ file.close(); }
    
    void open(const char* name);
    void open(const string& name){nm=name; open(name.c_str());}
    void close(){ file.close(); return;}
    
    // Informationen ueber status....
    bool iseof(){ return file.eof(); }
    
    // Public functions....
    void read_xy(double x[], double y[], const long N);
    void read_cols(double x[], double y[], int colx, int coly, int ncols, const long N);
    void read_col(double x[], int colx, int ncols, const long N);
    void read_x(double x[], const long N);
    void read_y(double y[], const long N);
    
    double get_valN(const long NR);
    int get_intN(const long NR);
    string get_stringN(const long NR);
    double get_next();
    int get_next_int();
    char get_next_char();
    string get_next_string();
    string get_next_line();
};

#endif
