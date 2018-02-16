/*\

    randord.c

    This utility outputs a random permutation of set {1,...,N},
    where N is given on the command line. The result is printed
    on standard output.

    05-NOV-1997

    Mo Van Dang
    UMR CNRS 6599
    Universite de Technologie de Compiegne, France

\*/


#include <stdio.h>     /* printf() */
#include <stdlib.h>    /* srandom() */

#include "nem_rnd.h"   /* RandomSeedByTime(), RandomInteger() */



/* =========================== */

main( int argc, char *argv[] )

/* =========================== */
{
  int  nbElts ;
  int* tabEltsV ;
  int  i ;

  void PrintUsage( const char* CmdS ) ;
  void Error( const char *MsgS , int ExitCode ) ;
  void GiveRandomSeed( int  MySeed ) ;


  /* Check arguments */
  if ( argc < 2 )
    {
      PrintUsage( argv[0] ) ;
      return 1 ;
    }

  nbElts = atoi( argv[ 1 ] ) ;
  if ( nbElts <= 0 )
    Error( "Number of elements must be greater than zero" , 1 ) ;


  if ( argc < 3 )
    RandomSeedByTime( ) ;
  else
    GiveRandomSeed( atoi( argv[ 2 ] ) ) ;


  /* Allocate and intialize the vector of integers */
  if ( ( tabEltsV = calloc( nbElts , sizeof( int ) ) ) == NULL )
    Error( "Not enough memory for that many elements" , 2 ) ;

  for ( i = 0 ; i < nbElts ; i ++ )
    tabEltsV[ i ] = i + 1 ;

  /* Run random permutation algorithm */
  RandomPermutationAlgo( tabEltsV , nbElts ) ;


  /* Print result to standard output */
  for ( i = 0 ; i < nbElts ; i ++ )
    fprintf( stdout , "%3d " , tabEltsV[ i ] ) ;
  fprintf( stdout , "\n" ) ;

  free( tabEltsV ) ;
  return 0 ;
}


/* =========================== */

void PrintUsage( const char* CmdS )

/* =========================== */
{
  fprintf( stderr , "\nSyntax :\n\n" ) ;
  fprintf( stderr , "  %s  N  [ seed ]\n\n" , CmdS ) ;
  fprintf( stderr , "  This utility prints to standard output a random permutation\n" ) ;
  fprintf( stderr , "  of set {1,...,N}.  A seed of the random generator may be specified\n" ) ;
  fprintf( stderr , "  (by default, seed = system time).\n\n" ) ;
}


/* =========================== */

void Error( const char *MsgS , int ExitCode )

/* =========================== */
{

  switch( ExitCode )
    {
    case 1 : 
      fprintf( stderr , "\n*** Input Error : %s\n\n" , MsgS ) ;
      break ;

    case 2 :
      fprintf( stderr , "\n*** Runtime Error : %s\n\n" , MsgS ) ;
      break ;

    default :
      fprintf( stderr , "\n*** Fatal Error : %s\n\n" , MsgS ) ;
    }

  exit( ExitCode ) ;

}



/* =========================== */

void GiveRandomSeed( int  MySeed )

/* =========================== */
{
#ifdef __TURBOC__
    srand( (unsigned) MySeed ) ;
#else
    srandom( MySeed ) ;
#endif
}



