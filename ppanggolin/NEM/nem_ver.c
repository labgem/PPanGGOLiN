/*\
    NEM_VER.C

    Informations about software successive versions.
    To be updated at each modification.
\*/

#include "nem_ver.h"  /* Exported prototypes */

const char *NemVersionStrC = "1.08" ; /* Current version of NEM */

void PrintVersions( FILE* F )         /* Describes successive versions */
{
    fprintf( F , "\n" ) ;
    fprintf( F , " Vers  Date      Description\n" ) ;
    fprintf( F , " 0.00  01.02.96  First version released on WEB\n" ) ;
    fprintf( F , " 1.00  31.05.96  Random start, image defaults, long help\n" ) ;
    fprintf( F , " 1.01  27.06.96  NCEM, sequential E-step\n" ) ;
    fprintf( F , " 1.02  17.10.96  Partially known labels\n" ) ;
    fprintf( F , " 1.03  02.10.97  Init centers from known labels, log is optional, augmented .mf\n" ) ;
    fprintf( F , " 1.04  11.01.98  Estimation of beta, Gibbsian EM, longer help\n" ) ;
    fprintf( F , " 1.05  09.04.98  Missing data, initialization modified\n" ) ;
    fprintf( F , " 1.06  26.02.99  Pseudo-likelihood beta, Laplace distributions\n" ) ;
    fprintf( F , " 1.07  08.04.99  Bernoulli distributions\n" ) ;
    fprintf( F , " 1.08  21.07.17  Add param input by file rather than by arguments \n" ) ;
    fprintf( F , "\n" ) ;
}

