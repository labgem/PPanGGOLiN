/*
 * program to compute permuted classification error in 2 cluster case.
 *
 * Mo Dang, oct 1997.
 */

#include <stdio.h>
#include <stdlib.h>


main(int argc, char *argv[])
{
  int    nbmax ;
  int*   clas1V ;
  int*   clas2V ;
  int    nbelts ;
  int    nbelts2 ;
  int    ielt ;
  int    nbdif ;
  float  percent ;

  int ReadIntVector( const char *fname, int* VP[], int nbmax ) ;

  if ( argc < 4 )
    {
      fprintf( stderr, "Usage : %s file1 file2 nbmax\n" , argv[0] ) ;
      fprintf( stderr, "  computes the percentage of misclassification\n" ) ;
      return 1 ;
    }

  nbmax = atoi( argv[ 3 ] ) ;
  nbelts = ReadIntVector( argv[1] , &clas1V , nbmax ) ;
  if ( nbelts <= 0 )
    return 2 ;

  nbelts2 = ReadIntVector( argv[2], &clas2V, nbmax ) ;
  if ( nbelts2 != nbelts )
    {
      fprintf( stderr, "Error : %d elts in %s and %d elts in %s\n", 
	       nbelts, argv[1], nbelts2, argv[2] ) ;
      return 3 ;
    }

  for ( ielt = 1, nbdif = 0 ; ielt <= nbelts ; ielt ++ )
    {
      if ( clas1V[ ielt ] != clas2V[ ielt ] )
	nbdif ++ ;
    }

  if ( nbdif > nbelts / 2 ) 
    nbdif = nbelts - nbdif ;

  percent = ( 100.0 * nbdif ) / nbelts ;

  fprintf( stdout , "Error = %3.1f %%  (%d on %d)\n", percent ,
	   nbdif , nbelts ) ;

  return 0 ;
}





int ReadIntVector( const char *fname, int* VP[], int nbmax )
{
  FILE* fid ;
  int   nbelts ;

  if ( ( fid = fopen( fname , "r" ) ) == NULL )
    {
      fprintf( stderr , "Cannot open file %s\n" ,  fname ) ;
      *VP = NULL ;
      return 2 ;
    }

  if( ( *VP = calloc( nbmax + 1 , sizeof( int ) ) ) == NULL )
    {
      fprintf( stderr , "Not enough memory\n" ) ;
      fclose( fid ) ;
      return 2 ;
    }

  for ( nbelts = 0 ; ( nbelts <= nbmax ) && (! feof( fid ) ) ; )
    {
      int n ;

      if ( fscanf( fid , "%d", &n ) >= 1 )
	{
	  nbelts ++ ;
	  (*VP)[ nbelts ] = n ;
	}
    }

  if ( nbelts == 0 )
    {
      free( *VP ) ;
      *VP = NULL ;
    }

  fclose( fid ) ;

  return nbelts ;
}
