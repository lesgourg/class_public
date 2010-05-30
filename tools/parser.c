#include "parser.h"

int parser_read_file(
		     char * filename,
		     struct file_content * pfc,
		     ErrorMsg errmsg
		     ){
  FILE * inputfile;
  char line[_LINE_LENGTH_MAX_];
  int counter;
  int is_data;
  FileArg name;
  FileArg value;

  class_open(inputfile,filename,"r",errmsg);

  counter = 0;
  while (fgets(line,_LINE_LENGTH_MAX_,inputfile) != NULL) {
    class_call(parser_read_line(line,&is_data,name,value,errmsg),errmsg,errmsg);
    if (is_data == _TRUE_) counter++;
  }

  class_test(counter == 0,
	     errmsg,
	     "No readable input in file %s",filename);

  class_alloc(pfc->filename,(strlen(filename)+1)*sizeof(char),errmsg);
  strcpy(pfc->filename,filename);
  pfc->size = counter++;
  class_alloc(pfc->name,pfc->size*sizeof(FileArg),errmsg);
  class_alloc(pfc->value,pfc->size*sizeof(FileArg),errmsg);

  rewind(inputfile);

  counter = 0;
  while (fgets(line,_LINE_LENGTH_MAX_,inputfile) != NULL) {
    class_call(parser_read_line(line,&is_data,name,value,errmsg),errmsg,errmsg);
    if (is_data == _TRUE_) {
      strcpy(pfc->name[counter],name);
      strcpy(pfc->value[counter],value);
      counter++;
    }
  }

  fclose(inputfile);

  return _SUCCESS_;

}

int parser_free(
		struct file_content * pfc
		) {
  free(pfc->filename);
  free(pfc->name);
  free(pfc->value);
}

int parser_read_line(
		char * line,
		int * is_data,
		char * name,
		char * value,
		ErrorMsg errmsg
		) {

  char * phash; 
  char * pequal;
  char * left;
  char * right;

  /* check that there is an '=' */

  pequal=strchr(line,'=');
  if (pequal == NULL) {*is_data = _FALSE_; return _SUCCESS_;}

  /* if yes, check that there is not an '#' before the '=' */

  phash=strchr(line,'#');
  if ((phash != NULL) && (phash-pequal<2)) {*is_data = _FALSE_; return _SUCCESS_;}
  
  /* get the name, i.e. the block before the '=' */

  left=line;
  while (left[0]==' ') {
    left++;
  }

  right=pequal-1;
  while (right[0]==' ') {
    right--;
  }

  if (right-left < 0) {*is_data = _FALSE_; return _SUCCESS_;}

  class_test(right-left+1 >= _ARGUMENT_LENGTH_MAX_,
	     errmsg,
	     "name starting by '%s' too long; shorten it or increase _ARGUMENT_LENGTH_MAX_",strncpy(name,left,(_ARGUMENT_LENGTH_MAX_-1)));
  
  strncpy(name,left,right-left+1);
  name[right-left+1]='\0';

  /* get the value, i.e. the block after the '=' */

  left = pequal+1;
  while (left[0]==' ') {
    left++;
  }

  if (phash == NULL)
    right = line+strlen(line)-1;
  else 
    right = phash-1;

  while (right[0]<=' ') {
    right--;
  }
 
  if (right-left < 0) {*is_data = _FALSE_; return _SUCCESS_;}

  class_test(right-left+1 >= _ARGUMENT_LENGTH_MAX_,
	     errmsg,
	     "value starting by '%s' too long; shorten it or increase _ARGUMENT_LENGTH_MAX_",strncpy(value,left,(_ARGUMENT_LENGTH_MAX_-1)));

  strncpy(value,left,right-left+1);
  value[right-left+1]='\0';

  *is_data = _TRUE_;

  return _SUCCESS_;

}

int parser_read_int(
		    struct file_content * pfc,
		    char * name,
		    int * value,
		    ErrorMsg errmsg
		    ) {
  int index;
  int i;

  index=0;

  while ((strcmp(pfc->name[index],name) != 0) && (index < pfc->size))
    index++;

  class_test(index == pfc->size,
	     errmsg,
	     "did not find parameter %s in file %s\n",name,pfc->filename);

  class_test(sscanf(pfc->value[index],"%d",value) != 1,
	     errmsg,
	     "could not read value of parameter %s in file %s\n",name,pfc->filename);

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
	       errmsg,
	       "multiple entry of parameter %s in file %s\n",name,pfc->filename);
  }

  return _SUCCESS_;

}

int parser_read_double(
		    struct file_content * pfc,
		    char * name,
		    double * value,
		    ErrorMsg errmsg
		    ) {
  int index;
  int i;

  index=0;

  while ((strcmp(pfc->name[index],name) != 0) && (index < pfc->size))
    index++;

  class_test(index == pfc->size,
	     errmsg,
	     "did not find parameter %s in file %s\n",name,pfc->filename);

  class_test(sscanf(pfc->value[index],"%lf",value) != 1,
	     errmsg,
	     "could not read value of parameter %s in file %s\n",name,pfc->filename);

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
	       errmsg,
	       "multiple entry of parameter %s in file %s\n",name,pfc->filename);
  }

  return _SUCCESS_;

}

int parser_read_string(
		       struct file_content * pfc,
		       char * name,
		       FileArg * value,
		       ErrorMsg errmsg
		       ) {
  int index;
  int i;

  index=0;

  while ((strcmp(pfc->name[index],name) != 0) && (index < pfc->size))
    index++;

  class_test(index == pfc->size,
	     errmsg,
	     "did not find parameter %s in file %s\n",name,pfc->filename);

  strcpy(*value,pfc->value[index]);

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
	       errmsg,
	       "multiple entry of parameter %s in file %s\n",name,pfc->filename);
  }

  return _SUCCESS_;

}
