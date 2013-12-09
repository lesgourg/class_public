#include "common.h"

void class_protect_sprintf(char* dest, char* tpl,...) {
  va_list args;
  int sz;
  va_start(args,tpl);
  sz=vsnprintf(dest, 2048,tpl,args);
  va_end(args);
}

void class_protect_fprintf(FILE* stream, char* tpl,...) {
  va_list args;
  int sz;
  char dest[6000];
  va_start(args,tpl);
  sz=vsnprintf(dest, 2048,tpl,args);
  va_end(args);
  fprintf(stream,dest);
}

void* class_protect_memcpy(void* dest, void* from, size_t sz) {
  return memcpy(dest, from,sz);
}
