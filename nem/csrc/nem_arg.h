#include "nem_typ.h"    /* DataT, ... */

/*
    V1.06-a   17-JUN-1998  NoiseModel -> StatModel
 */

int NemArgs
    (
          int           Argc,
          const char    *Argv[],
          char          *Fname,         /* O */
          StatModelT    *StatModelP,    /* O */
          NemParaT      *NemParaP,      /* O */
	  char*         DatadescS,      /* O */
	  DataT*        DataP,          /* O */
	  ImageNeighT*  ImageP,         /* O */
	  TypeET*       SpatialTypeP    /* O */	  
    ) ;


void PrintUsage( const char* CmdName ) ;

extern const char *CritStrVC[ CRIT_NB ] ;
