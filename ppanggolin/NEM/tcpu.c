/*
 *  tcpu.c
 *
 *  This utility computes the cpu time used by a given command
 *  and writes this time into a specified file.
 */

#include <stdio.h>      /* FILE */
#include <stdlib.h>     /* malloc() */
#include <sys/times.h>  /* times() */
#include <time.h>  /* CLK_TCK() CLOCKS_PER_SEC()*/

main(int argc, char *argv[] ) 
{
   clock_t begintime, endtime;
   struct tms *a_tms;

   FILE* fid ;
   int   sts ;


   if ( argc < 3 )
     {
       fprintf( stderr , "Usage : %s  cmd  timefile\n" , argv[ 0 ] ) ;
       return 1 ;
     }

   a_tms = ( struct tms *) malloc( sizeof (struct tms));

   times(a_tms); begintime = a_tms->tms_cutime;

   sts = system( argv[1] );

   times(a_tms); endtime = a_tms->tms_cutime;

   fid = fopen( argv[2], "w" ) ;
   fprintf( fid , " %5.3f \n", ((double)(endtime-begintime)/CLOCKS_PER_SEC));
   fclose( fid ) ;

   return sts ;
}
