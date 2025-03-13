/*\

    NEM_RND.C

    Programme NEM (Neighborhood EM) : routines de tirage aleatoire

    Van Mo DANG       Avril 96

Vers-mod  Date         Who  Description

1.04-a    10-OCT-1997  MD   add RandomPermutationAlgo()
1.04-b    05-NOV-1997  MD   use random/srandom instead of lrand48/srand48
1.04-c    12-JAN-1997  MD   add RandomReal()
\*/


#include <sys/types.h>   /* time_t */
#include <time.h>        /* time() */
#define _GNU_SOURCE      // FIX: Required for musl build because srandom is not POSIX
#include <stdlib.h>      /* srand48(), lrand48() */

#include "nem_rnd.h"

#define MAXRAND   0x7fffffff   /* 2**31 - 1 = maximum value of random() */


void  RandomSeedByTime( void ) 
{
    time_t   t ;

    t = time( 0 ) ;

#ifdef __TURBOC__
    srand( (unsigned) t ) ;
#else
    srandom( t ) ; /*V1.04-b*/
#endif
}


int   RandomInteger( int Mini, int Maxi ) 
{
    int       span ;
    long int  nrandom ;
    int       result ;


    if ( Mini >= Maxi )
      {
        return Maxi ;
      }

    span = Maxi - Mini + 1 ;

#ifdef __TURBOC__
    nrandom = rand( ) ;
#else
    nrandom = random( ) ; /*V1.04-b*/
#endif

    result = ( (int) ( nrandom % span ) ) + Mini  ;

    return  result ;
}


/*V1.05-a*/
float   RandomFloat( float Mini, float Maxi )
{
    float     span ;
    long int  nrandom ;
    float     result ;


    if ( Mini >= Maxi )
      {
        return Maxi ;
      }

    span = Maxi - Mini ;

#ifdef __TURBOC__
    nrandom = rand( ) ;
#else
    nrandom = random( ) ; /*V1.04-b*/
#endif

    result = ( (float) nrandom / MAXRAND ) * span + Mini  ;

    return  result ;
}





/* =========================== */

void RandomPermutationAlgo( int* TabV , int Nb )      /*V1.04-a*/

/* =========================== */
{
  int icou ;
  int iech ;
  int valech ;

  for ( icou = 0 ; icou < Nb ; icou ++ )
    {
      iech = RandomInteger( 0 , Nb - 1 ) ;

      valech       = TabV[ iech ] ;

      TabV[ iech ] = TabV[ icou ] ;

      TabV[ icou ] = valech ;
    }
}



/* Tests 

#include <stdio.h>

main( int argc, char *argv[] )
{
  int       n ;
  int       mini, maxi ;
  int       i ;

  if ( argc < 4 )
    {
      printf( "Au moins 3 args\n" );
      return 2 ;
    }
  n = atoi( argv[ 1 ] ) ;

  mini = atoi( argv[ 2 ] ) ;
  maxi = atoi( argv[ 3 ] ) ;

  RandomSeedByTime( ) ;

  for ( i = 0 ; i < n ; i ++ )
    {
      int r = RandomInteger( mini , maxi  ) ;

      fprintf( stdout, "%8d  ", r ) ;
    }

  fprintf( stdout, "\n" ) ;
  return 0 ;
}
*/
