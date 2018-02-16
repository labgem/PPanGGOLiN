/*\

    main.c

    A main() interface to compile a standalone program which is
    executable from the operating system.

    Jan 1998

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

Vers-mod  Date         Who  Description

1.05-a    25-JAN-1998  MD   Create to make an interface to operating system
1.05-b    30-JAN-1998  MD   Prototype of called mainfunc() in mainfunc.h
\*/


#include "mainfunc.h"   /* mainfunc() */

/* ------------------------------------------------------------------- */
main( int argc, const char *argv[] )
/*\
    Interface to a main routine.
\*/
/* ------------------------------------------------------------------- */
{
  return mainfunc( argc, argv ) ;
}
/* ~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
