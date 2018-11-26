#ifndef NEM_EXE_H
#define NEM_EXE_H
/*
 * nem_exe.h
 *
 * Prototype of main function called either by mexFunction() or main()
 *
 */

#include "nem_typ.h"    /* DataT, ... */
#include "nem_alg.h"    /* ClassifyByNem, ... */
#include "nem_rnd.h"    /* RandomInteger, ... */
#include "lib_io.h"     /* ReadOpeningComments, ... */
#include "genmemo.h"    /* GenAlloc, ... */

#include <stdio.h>      /* printf, ... */
#ifdef __TURBOC__
#include <alloc.h>      /* coreleft, ... */ 
#endif
#include <string.h>     /* strncpy, ... */
#include <math.h>       /* sqrt, ... */

extern int nem(const char* Fname,
        const int nk,
        const char* algo,
        const float beta,
        const char* convergence,
        const float convergence_th,
        const char* format,
        const int it_max,
        const int dolog,
        const char* model_family,
        const char* proportion,
        const char* dispersion,
        const int init_mode,
        const char* init_file,
        const char* out_file_prefix,
        const int seed);
#endif
