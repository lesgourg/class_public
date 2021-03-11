#include "common.h"

void class_protect_sprintf(char* dest, char* tpl,...) {
  va_list args;
  va_start(args,tpl);
  vsnprintf(dest, 2048,tpl,args);
  va_end(args);
}

void class_protect_fprintf(FILE* stream, char* tpl,...) {
  va_list args;
  char dest[6000];
  va_start(args,tpl);
  vsnprintf(dest, 2048,tpl,args);
  va_end(args);
  fprintf(stream,"%s",dest);
}

void* class_protect_memcpy(void* dest, void* from, size_t sz) {
  return memcpy(dest, from,sz);
}

int get_number_of_titles(char * titlestring){
  int i;
  int number_of_titles=0;

  for (i=0; i<strlen(titlestring); i++){
    if (titlestring[i] == '\t')
      number_of_titles++;
  }
  return number_of_titles;
}

/**
 * Finds wether or not a file exists.
 *
 * @param fname  Input: File name
 * @return boolean
 */

int file_exists(const char *fname){
  FILE *file = fopen(fname, "r");
  if (file != NULL){
    fclose(file);
    return _TRUE_;
  }

  return _FALSE_;

}

/**
 * Finds whether two doubles are equal or which one is bigger
 *
 * @param a Input: first number
 * @param b Input: second number
 * @return -1, 1 or 0
 */

int compare_doubles(const void *a,
                    const void *b){
  double *x = (double *) a;
  double *y = (double *) b;
  if (*x < *y)
    return -1;
  else if
    (*x > *y) return 1;
  return 0;
}

/**
 * This function detects if a string begins with a character,
 * ignoring whitespaces during its search
 *
 * returns the result, NOT the _SUCCESS_ or _FAILURE_ codes.
 * (This is done such that it can be used inside of an if statement)
 *
 * @param thestring  Input: string to test
 * @param beginchar  Input: the character by which the string begins or not
 * @return boolean
 */

int string_begins_with(char* thestring, char beginchar){

  /** Define temporary variables */
  int int_temp=0;
  int strlength = strlen((thestring));
  int result = _FALSE_;

  /** Check through the beginning of the string to see if the beginchar is met */
  for(int_temp=0;int_temp<strlength;++int_temp){
    /* Skip over whitespaces (very important) */
    if(thestring[int_temp]==' ' || thestring[int_temp]=='\t'){continue;}
    /* If the beginchar is met, everything is good */
    else if(thestring[int_temp]==beginchar){result=_TRUE_;}
    /* If something else is met, cancel */
    else{break;}
  }

  return result;
}

/**
 * Get version number
 *
 * @param version  Output: The string to write the version number into
 * @return the error status
 */

int class_version(
                  char * version
                  ) {

  sprintf(version,"%s",_VERSION_);
  return _SUCCESS_;
}
