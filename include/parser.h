#ifndef __PARSER__
#define __PARSER__

#include "common.h"

#define _LINE_LENGTH_MAX_ 1024 /**< size of the string read in each line of the file (extra characters not taken into account) */
#define _ARGUMENT_LENGTH_MAX_ 1024 /**< maximum size of each argument (name or value), including the final null character */

typedef char FileArg[_ARGUMENT_LENGTH_MAX_];

/* after reading a given file, all relevant information stored in this structure, in view of being processed later*/
struct FileContent {
  FileContent() {}
  ~FileContent() {
    if (filename) free(filename);
    if (name) free(name);
    if (value) free(value);
    if (read) free(read);
  }
  char* filename = nullptr;
  int size = 0;
  FileArg* name = nullptr;
  FileArg* value = nullptr;
  short* read = nullptr;
  bool is_shooting = false;
};

/**************************************************************/

/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int parser_read_file(
		     char * filename,
		     FileContent* pfc,
		     ErrorMsg errmsg
		     );

int parser_init(
		FileContent* pfc,
		int size,
        char * filename,
		ErrorMsg errmsg
		);

int parser_read_line(
		char * line,
		int * is_data,
		char * name,
		char * value,
		ErrorMsg errmsg
		);

int parser_read_int(
		    FileContent* pfc,
		    char * name,
		    int * value,
		    int * found,
		    ErrorMsg errmsg
		    );

int parser_read_double(
		    FileContent* pfc,
		    const char* name,
		    double * value,
		    int * found,
		    ErrorMsg errmsg
		    );

  int parser_read_double_and_position(
                                      FileContent* pfc,
                                      char * name,
                                      double * value,
                                      int * position,
                                      int * found,
                                      ErrorMsg errmsg
                                      );

int parser_read_string(
		       FileContent* pfc,
		       char * name,
		       FileArg * value,
		       int * found,
		       ErrorMsg errmsg
		       );

int parser_read_list_of_doubles(
				FileContent* pfc,
				char * name,
				int * size,
				double ** pointer_to_list,
				int * found,
				ErrorMsg errmsg
				);

int parser_read_list_of_integers(
				FileContent* pfc,
				char * name,
				int * size,
				int ** pointer_to_list,
				int * found,
				ErrorMsg errmsg
				);

int parser_read_list_of_strings(
				FileContent* pfc,
				char * name,
				int * size,
				char ** pointer_to_list,
				int * found,
				ErrorMsg errmsg
				);

int parser_cat(
	       const FileContent* pfc1,
	       const FileContent* pfc2,
	       FileContent* pfc3,
	       ErrorMsg errmsg
	       );

#ifdef __cplusplus
}
#endif

#endif
