/**
 * This file should not be modified
 * */

#ifdef __ASSIGN_DEFAULT_PRECISION__      
#define class_precision_parameter(NAME,TYPE,DEF_VALUE)          \
ppr->NAME = DEF_VALUE;                       
#endif                                        
#ifdef __ALLOCATE_PRECISION_PARAMETER__      
#define class_precision_parameter(NAME,TYPE,DEF_VALUE)          \
TYPE NAME;                                    
#endif
#ifdef __PARSE_PRECISION_PARAMETER__
#define class_precision_parameter(NAME,TYPE,DEF_VALUE)          \
class_read_ ## TYPE(#NAME,ppr->NAME);
#endif


#ifdef __ASSIGN_DEFAULT_PRECISION__      
#define class_string_parameter(NAME,DIR,STRING)   \
sprintf(ppr->NAME,__CLASSDIR__);                  \
strcat(ppr->NAME,DIR);
#endif                                        
#ifdef __ALLOCATE_PRECISION_PARAMETER__      
#define class_string_parameter(NAME,DIR,STRING)    \
FileName NAME;                                    
#endif
#ifdef __PARSE_PRECISION_PARAMETER__
#define class_string_parameter(NAME,DIR,STRING)     \
class_read_string(STRING,ppr->NAME);
#endif


#ifdef __ASSIGN_DEFAULT_PRECISION__      
#define class_type_parameter(NAME,READ_TP,REAL_TP,DEF_VAL) \
ppr->NAME = DEF_VAL;                       
#endif                                        
#ifdef __ALLOCATE_PRECISION_PARAMETER__      
#define class_type_parameter(NAME,READ_TP,REAL_TP,DEF_VAL) \
REAL_TP NAME;                                    
#endif
#ifdef __PARSE_PRECISION_PARAMETER__
#define class_type_parameter(NAME,READ_TP,REAL_TP,DEF_VAL) \
class_read_ ## READ_TP(#NAME,ppr->NAME);
#endif
