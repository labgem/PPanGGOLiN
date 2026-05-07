#ifndef GENMEMO_H
#define GENMEMO_H

/*\
    genmemo.h

    Prototypes of generic routines for memory allocation

    Jan 1998

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

    Vers-mod  Date         Description

    1.05-a    26-JAN-1998  Creation
    1.06-a    20-SEP-1998  Add macro freenull
\*/


#include <stdlib.h>  /* size_t */


/* ------------------------------------------------------------------- */
/* #define freenull( p ) do { \ 
   GenFree( p ); \
   p = NULL; \
   } while(0) 
*/


/* ------------------------------------------------------------------- */
void* GenAlloc
(
 size_t       nelem,        /* I : number of elements to allocate */ 
 size_t       elsize,       /* I : size in bytes of each element */
 int          doexit,       /* I : 1 if failure exits, 0 if only return NULL */
 const char*  where,        /* I : name of calling function */
 const char*  what          /* I : name of allocated object */
) ;

/* ------------------------------------------------------------------- */
void GenFree( void* ptr ) ;




#endif
