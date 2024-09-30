#include "parser.h"

int parser_read_file(char * filename,
                     struct file_content * pfc,
                     ErrorMsg errmsg){

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

  class_call(parser_init(pfc,counter,filename,errmsg),
             errmsg,
             errmsg);

  rewind(inputfile);

  counter = 0;
  while (fgets(line,_LINE_LENGTH_MAX_,inputfile) != NULL) {
    class_call(parser_read_line(line,&is_data,name,value,errmsg),errmsg,errmsg);
    if (is_data == _TRUE_) {
      strcpy(pfc->name[counter],name);
      strcpy(pfc->value[counter],value);
      pfc->read[counter]=_FALSE_;
      counter++;
    }
  }

  fclose(inputfile);

  return _SUCCESS_;

}

int parser_init(struct file_content * pfc,
                int size,
                char * filename,
                ErrorMsg errmsg) {

  if (size > 0) {
    pfc->size=size;
    class_alloc(pfc->filename,(strlen(filename)+1)*sizeof(char),errmsg);
    strcpy(pfc->filename,filename);
    class_alloc(pfc->name,size*sizeof(FileArg),errmsg);
    class_alloc(pfc->value,size*sizeof(FileArg),errmsg);
    class_alloc(pfc->read,size*sizeof(short),errmsg);
  }

  return _SUCCESS_;
}

int parser_free(struct file_content * pfc) {

  if (pfc->size > 0) {
    free(pfc->name);
    free(pfc->value);
    free(pfc->read);
    free(pfc->filename);
  }

  return _SUCCESS_;
}

int parser_read_line(char * line,
                     int * is_data,
                     char * name,
                     char * value,
                     ErrorMsg errmsg) {

  char * phash;
  char * pequal;
  char * left;
  char * right;

  /* check that there is an '=' (if you want the role of '=' to be
     played by ':' you only need to substitute it in the next line and
     recompile) */

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
  /* Ignore any of " ' */
  if(left[0]=='\'' || left[0]=='\"'){
    left++;
  }

  right=pequal-1;
  while (right[0]==' ') {
    right--;
  }
  if(right[0]=='\'' || right[0]=='\"'){
    right--;
  }

  /* deal with missing variable names */

  class_test(right-left < 0,
             errmsg,
             "Syntax error in the input line '%s': no name passed on the left of the '=' sign",line);

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

int parser_read_int(struct file_content * pfc,
                    char * name,
                    int * value,
                    int * found,
                    ErrorMsg errmsg) {
  int index;
  int i;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* read parameter value. If this fails, return an error */

  class_test(sscanf(pfc->value[index],"%d",value) != 1,
             errmsg,
             "could not read value of parameter '%s' in file '%s'\n",name,pfc->filename);

  /* if parameter read correctly, set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */

  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is found,
     return an error. */

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
               errmsg,
               "multiple entry of parameter '%s' in file '%s'\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */

  return _SUCCESS_;

}

int parser_read_double(struct file_content * pfc,
                       char * name,
                       double * value,
                       int * found,
                       ErrorMsg errmsg) {
  int index;
  int i;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* read parameter value. If this fails, return an error */

  class_test(sscanf(pfc->value[index],"%lg",value) != 1,
             errmsg,
             "could not read value of parameter '%s' in file '%s'\n",name,pfc->filename);

  /* if parameter read correctly, set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */

  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is found,
     return an error. */

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
               errmsg,
               "multiple entry of parameter '%s' in file '%s'\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */

  return _SUCCESS_;

}

int parser_read_double_and_position(struct file_content * pfc,
                                    char * name,
                                    double * value,
                                    int * position,
                                    int * found,
                                    ErrorMsg errmsg) {
  int index;
  int i;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* read parameter value. If this fails, return an error */

  class_test(sscanf(pfc->value[index],"%lg",value) != 1,
             errmsg,
             "could not read value of parameter '%s' in file '%s'\n",name,pfc->filename);

  /* if parameter read correctly, set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */

  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is found,
     return an error. */

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
               errmsg,
               "multiple entry of parameter '%s' in file '%s'\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */

  * position = index;

  return _SUCCESS_;

}

int parser_read_string(struct file_content * pfc,
                       char * name,
                       FileArg * value,
                       int * found,
                       ErrorMsg errmsg) {
  int index;
  int i;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* read parameter value. */

  strcpy(*value,pfc->value[index]);

  /* Set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */

  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is found,
     return an error. */

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
               errmsg,
               "multiple entry of parameter '%s' in file '%s'\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */

  return _SUCCESS_;

}

int parser_read_list_of_doubles(struct file_content * pfc,
                                char * name,
                                int * size,
                                double ** pointer_to_list,
                                int * found,
                                ErrorMsg errmsg) {
  int index;
  int i;

  char * string;
  char * substring;
  FileArg string_with_one_value;

  double * list;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* count number of comas and compute size = number of comas + 1 */
  i = 0;
  string = pfc->value[index];
  do {
    i ++;
    substring = strchr(string,',');
    string = substring+1;
  } while(substring != NULL);

  *size = i;

  /* free and re-allocate array of values */
  class_alloc(list,*size*sizeof(double),errmsg);
  *pointer_to_list = list;

  /* read one double between each comas */
  i = 0;
  string = pfc->value[index];
  do {
    i ++;
    substring = strchr(string,',');
    if (substring == NULL) {
      strcpy(string_with_one_value,string);
    }
    else {
      strncpy(string_with_one_value,string,(substring-string));
      string_with_one_value[substring-string]='\0';
    }
    class_test(sscanf(string_with_one_value,"%lg",&(list[i-1])) != 1,
               errmsg,
               "could not read %dth value of list of parameters '%s' in file '%s'\n",
               i,name,pfc->filename);
    string = substring+1;
  } while(substring != NULL);

  /* if parameter read correctly, set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */

  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is found,
     return an error. */

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
               errmsg,
               "multiple entry of parameter '%s' in file '%s'\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */

  return _SUCCESS_;

}

int parser_read_list_of_integers(struct file_content * pfc,
                                 char * name,
                                 int * size,
                                 int ** pointer_to_list,
                                 int * found,
                                 ErrorMsg errmsg) {
  int index;
  int i;

  char * string;
  char * substring;
  FileArg string_with_one_value;

  int * list;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* count number of comas and compute size = number of comas + 1 */
  i = 0;
  string = pfc->value[index];
  do {
    i ++;
    substring = strchr(string,',');
    string = substring+1;
  } while(substring != NULL);

  *size = i;

  /* free and re-allocate array of values */
  class_alloc(list,*size*sizeof(int),errmsg);
  *pointer_to_list = list;

  /* read one integer between each comas */
  i = 0;
  string = pfc->value[index];
  do {
    i ++;
    substring = strchr(string,',');
    if (substring == NULL) {
      strcpy(string_with_one_value,string);
    }
    else {
      strncpy(string_with_one_value,string,(substring-string));
      string_with_one_value[substring-string]='\0';
    }
    class_test(sscanf(string_with_one_value,"%d",&(list[i-1])) != 1,
               errmsg,
               "could not read %dth value of list of parameters '%s' in file '%s'\n",
               i,name,pfc->filename);
    string = substring+1;
  } while(substring != NULL);

  /* if parameter read correctly, set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */

  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is found,
     return an error. */

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
               errmsg,
               "multiple entry of parameter '%s' in file '%s'\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */
  return _SUCCESS_;

}

int parser_read_list_of_strings(struct file_content * pfc,
                                char * name,
                                int * size,
                                char ** pointer_to_list,
                                int * found,
                                ErrorMsg errmsg) {
  int index;
  int i;

  char * string;
  char * substring;
  FileArg string_with_one_value;

  char * list;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* count number of comas and compute size = number of comas + 1 */
  i = 0;
  string = pfc->value[index];
  do {
    i ++;
    substring = strchr(string,',');
    string = substring+1;
  } while(substring != NULL);

  *size = i;

  /* free and re-allocate array of values */
  class_alloc(list,*size*sizeof(FileArg),errmsg);
  *pointer_to_list = list;

  /* read one string between each comas */
  i = 0;
  string = pfc->value[index];
  do {
    i ++;
    substring = strchr(string,',');
    if (substring == NULL) {
      strcpy(string_with_one_value,string);
    }
    else {
      strncpy(string_with_one_value,string,(substring-string));
      string_with_one_value[substring-string]='\0';
    }
    strcpy(list+(i-1)*_ARGUMENT_LENGTH_MAX_,string_with_one_value);
    //Insert EOL character:
    *(list+i*_ARGUMENT_LENGTH_MAX_-1) = '\n';
    string = substring+1;
  } while(substring != NULL);

  /* if parameter read correctly, set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */
  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is
     found,
     return an error. */
  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
               errmsg,
               "multiple entry of parameter '%s' in file '%s'\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */
  return _SUCCESS_;
}

int parser_cat(struct file_content * pfc1,
               struct file_content * pfc2,
               struct file_content * pfc3,
               ErrorMsg errmsg) {

  int i;

  class_test(pfc1->size < 0.,
             errmsg,
             "size of file_content structure probably not initialized properly\n");

  class_test(pfc2->size < 0.,
             errmsg,
             "size of file_content structure probably not initialized properly\n");

  if (pfc1->size == 0) {
    class_alloc(pfc3->filename,(strlen(pfc2->filename)+1)*sizeof(char),errmsg);
    class_sprintf(pfc3->filename,"%s",pfc2->filename);
  }
  if (pfc2->size == 0) {
    class_alloc(pfc3->filename,(strlen(pfc1->filename)+1)*sizeof(char),errmsg);
    class_sprintf(pfc3->filename,"%s",pfc1->filename);
  }
  if ((pfc1->size !=0) && (pfc2->size != 0)) {
    class_alloc(pfc3->filename,(strlen(pfc1->filename)+strlen(pfc2->filename)+5)*sizeof(char),errmsg);
    class_sprintf(pfc3->filename,"%s or %s",pfc1->filename,pfc2->filename);
  }

  pfc3->size = pfc1->size + pfc2->size;
  class_alloc(pfc3->value,pfc3->size*sizeof(FileArg),errmsg);
  class_alloc(pfc3->name,pfc3->size*sizeof(FileArg),errmsg);
  class_alloc(pfc3->read,pfc3->size*sizeof(short),errmsg);

  for (i=0; i < pfc1->size; i++) {
    strcpy(pfc3->value[i],pfc1->value[i]);
    strcpy(pfc3->name[i],pfc1->name[i]);
    pfc3->read[i]=pfc1->read[i];
  }

  for (i=0; i < pfc2->size; i++) {
    strcpy(pfc3->value[i+pfc1->size],pfc2->value[i]);
    strcpy(pfc3->name[i+pfc1->size],pfc2->name[i]);
    pfc3->read[i+pfc1->size]=pfc2->read[i];
  }

  return _SUCCESS_;

}



int parser_check_options(char * strinput, char ** options, int N_options, int* valid){

  int i, j, n_option, string_length, option_length;
  int found;
  char str[_ARGUMENT_LENGTH_MAX_];
  strcpy(str,strinput);

  *valid = _TRUE_;

  for(n_option = 0; n_option < N_options; ++n_option){

    string_length = strlen(str);                // Length of string
    option_length = strlen(options[n_option]);  // Length of option to remove from string

    for(i=0; i <= string_length - option_length; i++){

      /* Match first character */
      if(str[i] == options[n_option][0]){
        found = _TRUE_;
      }
      else{
        found = _FALSE_;
        continue;
      }

      /* Only if first character fits, try all the others */
      for(j=0; j<option_length; j++){
        if(str[i + j] != options[n_option][j]){
          found = _FALSE_;
          break;
        }
      }

      /* If it is not seperated by anything afterwards, it is not its own word */
      if(str[i+j] != ',' && str[i+j] != ' ' && str[i+j] != '\t' && str[i+j] != '\n' && str[i+j] != '\0' && str[i+j] != '.'  && str[i+j] != '&'){
        found = _FALSE_;
      }

      /* If it is not seperated by anything before, it is not its own word */
      if(i>0 && str[i-1] != ',' && str[i-1] != ' ' && str[i-1] != '\t' && str[i-1] != '\n' && str[i-1] != '\0' && str[i-1] != '.'  && str[i-1] != '&'){
        found = _FALSE_;
      }

      /* *
       * If word is found then shift all characters to left
       * and decrement the string length
       * */
      if(found == _TRUE_){

        for(j=i; j<=string_length - option_length; j++){
            str[j] = str[j+option_length];
        }

        string_length = string_length - option_length;

        // We will match the next occurrence of the word from the current index.
        i--;

      }
    }
  }

  string_length = strlen(str);
  /* Check that the remaining characters are exclusively seperators */
  for(i=0;i<string_length;i++){
    if( str[i] != ',' && str[i] != ' ' && str[i] != '\t' && str[i] != '\n' && str[i] != '\0' && str[i] != '.'  && str[i] != '&'){
      *valid = _FALSE_;
    }
  }

  return _SUCCESS_;
}

int parser_extend(struct file_content * pfc, int N_extend, ErrorMsg errmsg) {
  // Append N_extend empty entries in the vectors of pfc
  pfc->size += N_extend;
  class_realloc(pfc->name,  pfc->size*sizeof(FileArg), errmsg);
  class_realloc(pfc->value, pfc->size*sizeof(FileArg), errmsg);
  class_realloc(pfc->read,  pfc->size*sizeof(short),   errmsg);
  for (int index = pfc->size - N_extend; index < pfc->size; index++) {
    pfc->name[index][0] = '\0';
    pfc->value[index][0] = '\0';
    pfc->read[index] = _FALSE_;
  }
  return _SUCCESS_;
}

int parser_copy(struct file_content * pfc_source, struct file_content * pfc_destination, int index_start, int index_end) {
  // Copy the entries from index_start to index_end from pfc_source to pfc_destination
  memcpy(pfc_destination->name  + index_start, pfc_source->name  + index_start, (index_end - index_start)*sizeof(FileArg));
  memcpy(pfc_destination->value + index_start, pfc_source->value + index_start, (index_end - index_start)*sizeof(FileArg));
  memcpy(pfc_destination->read  + index_start, pfc_source->read  + index_start, (index_end - index_start)*sizeof(short));
  return _SUCCESS_;
}

int parser_init_from_pfc(struct file_content * pfc_source, struct file_content * pfc_destination, ErrorMsg errmsg) {
  class_call(parser_init(pfc_destination,
                         pfc_source->size,
                         pfc_source->filename,
                         errmsg),
             errmsg,errmsg);
  parser_copy(pfc_source, pfc_destination, 0, pfc_source->size);
  return _SUCCESS_;
}
