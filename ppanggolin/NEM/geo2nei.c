/*\

    geo2nei.c

    This utility takes as input a set of geographic coordinates
    and computes the neighborhood system induced by the distances
    between the objects.

    November 1997

    Mo Van Dang
    Universite de Technologie de Compiegne
    URA CNRS 817

\*/

#include <stdio.h>      /* FILE , printf ... */
#include <stdlib.h>     /* malloc, ... */
#include "values.h"     /* MAXFLOAT */
#include <math.h>       /* sqrt */

#include "lib_io.h"     /* AskFileToRead ... */

#define VERSION   "1.00"


#define MAT( A , nl , nc , l , c ) ( ( A[ ( l * nc ) + c ] ) )


#define TRUE        1
#define FALSE       0

#define MAX_DIM     4
#define LEN_DESC    500

#define EXPFAC      3.0

#define SEPAR_STR   " \t"

#ifndef MAXFLOAT
 #ifdef FLT_MAX
   #define MAXFLOAT FLT_MAX
 #else
   #define MAXFLOAT 3.40282347e+38F
 #endif
#endif


typedef enum {
  WEI_NONE,
  WEI_CONST,
  WEI_EXP,
  WEI_NB
}
WeiModT ;


typedef struct
{
    int     Nelts;                /* number of elements */
    float   AveSumWei ;           /* ave sum of weights */
    float   MinSumWei ;           /* min sum of weights */
    float   MaxSumWei ;           /* max sum of weights */
    float   AveDisNearest ;       /* ave distance to 4th nearest neighbor */
    float   MinDisNearest ;       /* min distance to 4th nearest neighbor */
    float   MaxDisNearest ;       /* max distance to 4th nearest neighbor */
    float   AveNbNei ;            /* ave number of neighbors */
    int     MinNbNei ;            /* minimum number of neighbors */
    int     MaxNbNei ;            /* maximum number of neighbors */
}
StatsT ;



main( int argc , char *argv[] )
{
    int     sts = 0 ;                   /* status returned by the program */
    char    namedat[ LEN_FILE + 1 ] ;   /* input data file name */
    char    namenei[ LEN_FILE + 1 ] ;   /* output neighbors file name */
    char    descnei[ LEN_DESC + 1 ] ;   /* comment of neigbors file */
    int     nd ;                        /* number of spatial variables */
    int     d ;                         /* current spatial variable */
    int     spacol[ MAX_DIM ] ;         /* columns of spatial variables */
    int     mincol , ncol ;             /* input file : nb of columns */
    int     npt ;                       /* input file : nb of lines */
    float*  locM ;                      /* spatial locations : npt * nd */

    int     ranknei ;                   /* rank of nearest neighbor */
    float   weifac ;                    /* factor -> average sum weights = 4 */
    float   disthres ;                  /* pts neighbors if dis<=threshold */
    WeiModT weimod ;                    /* weight computation mode */

    StatsT  stats ;                     /* neighbor stats */

    void PrintHelp( const char* CmdName ) ;
    void PrintVersion( const char* CmdName ) ;

    int ComputeNeighbors
        (
         const float*   PtsM , 
         int            Npt , 
         int            Nc , 
	 int            RankNei ,
         float          Dth, 
         WeiModT        WeiMod , 
	 float          WeiFac ,
         const char*    NameNei , 
         const char*    Desc ,
         int            DoSave , 
	 StatsT*        NeiStatsP
        ) ;

    void DisplayStats
      (
       int              RankNei, /* I : rank of nearest neighbor */
       const StatsT*    StatsP   /* I : distances/neighbors statistics */
      ) ;


    if ( argc > 1 )
      {
	if ( strcmp( argv[ 1 ] , "-h" ) == 0 )
	  {
	    PrintHelp( argv[ 0 ] ) ;
	    return 1 ;
	  }
	else if ( strcmp( argv[ 1 ] , "-v" ) == 0 )
	  {
	    PrintVersion( argv[ 0 ] ) ;
	    return 1 ;
	  }
	else
	  {
	    printf( "Type %s -h to get help or -v to know current version\n" ,
		    argv[ 0 ] ) ;
	    return 2 ;
	  }
      }

    printf( "* * *  Welcome to GEO2NEI program V%s  * * *\n\n" ,
	    VERSION ) ;

    /* Ask for all parameters : 
       input file name, number and columns of spatial variables ;
    */
    if ( AskFileToRead( "input objects-variables", namedat ) != 0 )
       return 2 ;
    if ( AskFileToWrite( "output neighborhood" , FALSE , namenei ) != 0 )
       return 2 ;
    printf( "Enter comment in neighbors file : " ) ;
    gets( descnei ) ;

    if ( AskInteger( "number of spatial coordinates" , 2 , 1 , MAX_DIM , 
                     & nd ) != 0 )
       return 2 ;

    for ( d = 0 ; d < nd ; d ++ )
    {
        char    msg[ 120 ] ;

        sprintf( msg , "column of spatial coordinate %d" , d + 1 ) ;
        if ( AskInteger( msg , d + 1 , 1 , 100 , & spacol[ d ] ) != 0 )
           return 2 ;
        spacol[ d ] -- ; /* C-like indices start from 0 */
    }

    if ( AskInteger( "Average number of neighbors" , 4 , 1 , 20 , 
                     & ranknei ) != 0 )
      return 2 ;

    if ( AskInteger( "Weight computation {0=>none, 1=>const 2=>exp[(-d/dmax)^2]}" , 
		     WEI_NONE , WEI_NONE , WEI_NB-1 , (int*) & weimod ) != 0 )
      return 2 ;
    
/*    if ( AskFloat( "Distance threshold" , 1.00 , 0.01 , 100.00 , 
/*		     & disthres ) != 0 )
/*	 return 2 ;
 */

    printf( "Analyzing file %s ...\n" , namedat ) ;
    switch( CountLinesColumns( namedat , SEPAR_STR , 
                               &mincol , &ncol , &npt ) )
    {
        case 0 : printf( "File %s : %d points, %d columns\n", 
                         namedat, npt, ncol ) ;
                 break ;
        case 1 : printf( "Error : '%s' unequal columns (%d to %d)\n" ,
                         namedat , mincol , ncol ) ;
                 return 2 ;
        default :
                printf( "Error while analyzing file %s\n" , namedat ) ;
                return 2 ;
    }

    for ( d = 0 ; d < nd ; d ++ )
    {
        if ( spacol[ d ] >= ncol )
        {
            printf( "Error : spatial var. %d's column = %d > nb columns\n" ,
                    d + 1 , spacol[ d ] + 1 ) ;
            return 2 ;
        }
    }


    /* Allocate and read spatial locations from file to memory 
    */
    if ( ( locM = malloc( npt * nd * sizeof( float ) ) ) ==  NULL )
    {
        printf( "Error : out of memory for %d site locations\n" , npt ) ; 
        return 2 ;
    }

    if ( ReadSelectedColumns( namedat , npt , ncol , nd , spacol , 
                              locM ) != 0 )
    {
        printf( "Error while reading spatial locations from file\n" ) ;
        free( locM ) ;
        return 2 ;
    }

    /* Compute and save neighborhood system
    */
    printf( "Computing neighborhood graph ...\n" ) ;

    /* Step 1 : compute distance stats to set distance threshold */
    disthres = 1.0 ;
    weifac = 1.0 ;
    if ( ComputeNeighbors( locM , npt , nd , ranknei , disthres , 
			   WEI_NONE , weifac , namenei , descnei , FALSE , 
			   &stats ) != -1 )
      {
	/* Step 2 : compute neighborhood graph with non-normalized weights */
	disthres = stats.AveDisNearest * 1.1 ;
	if ( ComputeNeighbors( locM , npt , nd , ranknei , disthres , 
			       weimod , weifac , namenei , descnei , FALSE , 
			       &stats ) != -1 )	
	  {
	    /* Step 3 : save neighborhood graph with normalized weights */
	    weifac = 4.0 / stats.MaxSumWei ;
	    printf( "weifac = %g\n", weifac ) ;
	    if ( ComputeNeighbors( locM , npt , nd , ranknei , disthres , 
			       weimod , weifac , namenei , descnei , TRUE , 
			       &stats ) != -1 )	
	      DisplayStats( ranknei, &stats ) ;
	    else
	      printf( "Error while computing neighborhood graph 2\n" ) ;
	  }
	else
	  printf( "Error while computing neighborhood graph 1\n" ) ;
	    
      }
    else
      printf( "Error while computing distance statistics\n" ) ;


    /* Free allocated memory */
    free( locM ) ;

    return sts ;

} /* main() */




/* ------------------------------------------------------------------- */
int ComputeNeighbors
        (
         const float*   PtsM , 
         int            Npt , 
         int            Nd , 
	 int            RankNei ,
         float          Dth , 
         WeiModT        WeiMod , 
	 float          WeiFac ,
         const char*    NameNei , 
         const char*    Desc ,
         int            DoSave , 
	 StatsT*        StatsP
        ) 
/* ------------------------------------------------------------------- */
{
    FILE*  fnei ;       /* File to write into */
    float  sqrDth ;     /* Squared distance threshold */

    float* sqrDisV ;    /* [Npt] Current point's distance to other points */
    float* weiV ;       /* [Npt] Current point's neighbors weights > 0.0 */
    int*   neiV ;       /* [Npt] Current point's neighbors indices 0..Npt-1 */

    int    pt ;         /* Index of current point 0..Npt-1 */
    int    ok ;         /* FALSE if problem in writing file */

    void StartStats
      (
       int      Npt ,          /* I : number of points */
       StatsT*  StatsP         /* O : distances/neighbors statistics */
      ) ;

    void CompuSqrDisV
      (
       int           Pt ,      /* I : index of current point 0..Npt-1 */
       int           RankNei,  /* I : neighbor rank to compute distance */
       const float*  PtsM ,    /* I : [Npt*Nd] matrix of point coordinates */
       int           Npt ,     /* I : number of points */
       int           Nd ,      /* I : number of coordinates */
       float*        SqrDisV , /* O : [Npt] distances to other points */
       float*        SqrDnearP /* O : distance to nearest point */
      ) ;

    void SqrDisVToNeiV
      ( 
       int           Pt ,      /* I : index of current point 0..Npt-1 */
       float         SqrDth ,  /* I : squared distance threshold */ 
       const float*  SqrDisV , /* I : [Npt] distances to other points */
       int           Npt ,     /* I : number of points */
       WeiModT       WeiMod ,  /* I : weight computation mode */
       float         WeiFac ,  /* I : weight normalizing factor */
       int*          NbNeiP ,  /* O : number of neighbors 0..Npt-1 */
       int*          NeiV ,    /* O : [Npt] neighbors indices 0..Npt-1 */
       float*        WeiV ,    /* O : [Npt] neighbors weights > 0.0 */
       float*        SumWeiP   /* O : sum of weights > 0.0 */
      ) ;

    void UpdateStats
      (
       int           Pt ,      /* I : index of current point 0..Npt-1 */
       float         SqrDnear, /* I : distance to nearest point */
       int           NbNei ,   /* I : number of neighbors 0..Npt-1 */
       float         SumWei ,  /* I : sum of neighbor weights */
       StatsT*       StatsP    /* I/O : distances/neighbors statistics */
      ) ;

    int SaveNeiV
      ( 
       FILE*         Fnei ,    /* I/O : file to write into */
       int           Pt ,      /* I : index of current point 0..Npt-1 */
       int           NbNei ,   /* I : number of neighbors 0..Npt-1 */
       const int*    NeiV      /* I : [Npt] neighbors indices 0..Npt-1 */
      ) ;

    int SaveWeiV
      ( 
       FILE*         Fnei ,    /* I/O : file to write into */
       int           NbNei ,   /* I : number of neighbors 0..Npt-1 */
       const float*  WeiV      /* I : [Npt] neighbors weights > 0.0 */
      ) ;

    void CloseStats
      (
       StatsT*  StatsP         /* I/O : distances/neighbors statistics */
      ) ;



    /* Convert parameters to squared distances */
    sqrDth = Dth * Dth ;

    /* If save to file requested, create file and write header */
    if ( DoSave ) 
      {
	if ( ( fnei = fopen( NameNei , "w" ) ) == NULL )
	  {
	    printf( "Error : cannot create file %s\n", NameNei ) ;
	    return -1 ;
	  }

	fprintf( fnei , "# %s\n" , Desc ) ;
	if ( WeiMod == WEI_NONE )
	  fprintf( fnei , "0\n" ) ;
	else
	  fprintf( fnei , "1\n" ) ;
      }

    /* Allocate memory for local computations */
    sqrDisV = malloc( Npt * sizeof( float ) ) ;
    weiV    = malloc( Npt * sizeof( float ) ) ;
    neiV    = malloc( Npt * sizeof( int ) ) ;
    if ( ( sqrDisV == NULL ) || ( weiV == NULL ) || ( neiV == NULL ) )
      {
	printf( "Could not allocate distance/weight/neigbor vectors [%d]\n" ,
		 Npt ) ;
	if ( sqrDisV != NULL ) free( sqrDisV ) ;
	if ( weiV    != NULL ) free( weiV ) ;

	return -1 ;
      }


    /* For each point */
    for ( pt = 0 , 
	    ok = TRUE , 
	    StartStats( Npt , StatsP ) ; 
	  ( pt < Npt ) && ok ; 
	  pt ++ )
      {
	float   sqrdisNearest ;   /* Distance to nearest point */
	int     nbnei ;           /* Number of neighbors */
	float   sumwei ;          /* Sum of neighbor weights */

	/* Compute its distance to other points */
	CompuSqrDisV( pt , RankNei, PtsM , Npt , Nd , 
		      sqrDisV , &sqrdisNearest ) ;

	/* Threshold the computed distances to get neighbors and weights */
	SqrDisVToNeiV( pt , sqrDth , sqrDisV , Npt , WeiMod , WeiFac ,
		       &nbnei , neiV , weiV , &sumwei ) ;

	/* printf( "sumwei %d = %g\n" , pt , sumwei ) ; */

	/* Update distance/neighbor statistics */
	UpdateStats( pt , sqrdisNearest , nbnei , sumwei , StatsP ) ;

	/* Eventually save neighbors and weights to file */
	if ( DoSave )
	  {
	    if ( SaveNeiV( fnei , pt , nbnei , neiV ) != 0 )
	      ok = FALSE ;
	    else
	      if ( WeiMod != WEI_NONE )
		if ( SaveWeiV( fnei , nbnei , weiV ) != 0 )
		  ok = FALSE ;
	    fprintf( fnei , "\n" ) ;
	  }
      }
    CloseStats( StatsP ) ;


    /* Free memory used for local computations */
    free( sqrDisV ) ;
    free( neiV ) ;
    free( weiV ) ;

    /* If save to file requested, close file */
    if ( DoSave ) 
      {
	fclose( fnei ) ;
      }

    if ( ok ) 
      return 0 ;
    else
      return -1 ;

} /* end of ComputeNeighbors() */



/* ------------------------------------------------------------------- */
void StartStats
 (
  int      Npt ,          /* I : number of points */
  StatsT*  StatsP         /* O : distances/neighbors statistics */
 ) 
/* ------------------------------------------------------------------- */
{
    StatsP->Nelts = 0 ;

    StatsP->AveSumWei = 0.0 ;
    StatsP->MinSumWei = MAXFLOAT ;
    StatsP->MaxSumWei = 0.0 ;

    StatsP->AveDisNearest = 0.0 ;
    StatsP->MinDisNearest = MAXFLOAT ;
    StatsP->MaxDisNearest = MINFLOAT ;

    StatsP->AveNbNei     = 0.0 ;
    StatsP->MinNbNei      = Npt ;
    StatsP->MaxNbNei      = 0 ;
}


/* ------------------------------------------------------------------- */
void UpdateStats
 (
  int           Pt ,      /* I : index of current point 0..Npt-1 */
  float         SqrDnear, /* I : distance to nearest point */
  int           NbNei ,   /* I : number of neighbors 0..Npt-1 */
  float         SumWei ,  /* I : sum of neighbor weights */
  StatsT*       StatsP    /* I/O : distances/neighbors statistics */
 ) 
/* ------------------------------------------------------------------- */
{

    StatsP->Nelts ++ ;


    StatsP->AveSumWei += SumWei ;

    if ( SumWei < StatsP->MinSumWei )
      StatsP->MinSumWei = SumWei ;

    if ( SumWei > StatsP->MaxSumWei )
      StatsP->MaxSumWei = SumWei ;


    StatsP->AveDisNearest += sqrt( SqrDnear ) ;

    if ( SqrDnear < StatsP->MinDisNearest )
      StatsP->MinDisNearest = SqrDnear ;

    if ( SqrDnear > StatsP->MaxDisNearest )
      StatsP->MaxDisNearest = SqrDnear ;


    StatsP->AveNbNei += NbNei ;

    if ( NbNei < StatsP->MinNbNei ) 
      StatsP->MinNbNei = NbNei ;

    if ( NbNei > StatsP->MaxNbNei ) 
      StatsP->MaxNbNei = NbNei ;

}


/* ------------------------------------------------------------------- */
void CloseStats
 (
  StatsT*  StatsP         /* I/O : distances/neighbors statistics */
 ) 
/* ------------------------------------------------------------------- */
{
  if ( StatsP->MinDisNearest >= 0 )
    StatsP->MinDisNearest = sqrt( StatsP->MinDisNearest ) ;

  if ( StatsP->MaxDisNearest >= 0 )
    StatsP->MaxDisNearest = sqrt( StatsP->MaxDisNearest ) ;

  if ( StatsP->Nelts > 0 )
    {
      StatsP->AveSumWei     /= StatsP->Nelts ;
      StatsP->AveDisNearest /= StatsP->Nelts ;
      StatsP->AveNbNei      /= StatsP->Nelts ;
    }
}


/* ------------------------------------------------------------------- */
void DisplayStats
 (
  int              RankNei, /* I : rank of nearest neighbor */
  const StatsT*    StatsP   /* I : distances/neighbors statistics */
 ) 
/* ------------------------------------------------------------------- */
{
    printf( "Statistics (ave, min, max) : \n" ) ;
    printf( "  Distance to %dth nearest site :  %8.2f  (%8.2f to %8.2f)\n" ,
	    RankNei,
	    StatsP->AveDisNearest, 
	    StatsP->MinDisNearest , StatsP->MaxDisNearest );
    printf( "  Number of neighbors          :  %8.2f  (%8d to %8d)\n" , 
	    StatsP->AveNbNei ,
	    StatsP->MinNbNei , StatsP->MaxNbNei );
    printf( "  Sum of neighbor weights      :  %8.2f  (%8.2f to %8.2f)\n" ,
	    StatsP->AveSumWei, 
	    StatsP->MinSumWei , StatsP->MaxSumWei );
}

/* ------------------------------------------------------------------- */
static int floatcompare(const void *x, const void *y)
/* ------------------------------------------------------------------- */
{
  const float* xx = x;
  const float* yy = y;

  if (*xx > *yy)
    return (1);
  if (*xx < *yy)
    return (-1);
  return (0);
}

/* ------------------------------------------------------------------- */
/* Computes the distances from a given point to other points 
 */
void CompuSqrDisV
 (
  int           Pt ,      /* I : index of current point 0..Npt-1 */
  int           RankNei,  /* I : neighbor rank to compute distance */
  const float*  PtsM ,    /* I : [Npt*Nd] matrix of point coordinates */
  int           Npt ,     /* I : number of points */
  int           Nd ,      /* I : number of coordinates */
  float*        SqrDisV , /* O : [Npt] distances to other points */
  float*        SqrDnearP /* O : distance to RankNei nearest point */
 )
/* ------------------------------------------------------------------- */
{
    int   j ;
    float dnear ;

    for ( j = 0 , dnear = MAXFLOAT ; j < Npt ; j ++ )
      {
	if ( j != Pt )
	  {
	    int   d ;
	    float sqrdis ;

	    for ( d = 0 , sqrdis = 0.0 ; d < Nd ; d ++ )
	      {
	        float cpt = MAT( PtsM , Npt , Nd , Pt , d ) ;
		float cj  = MAT( PtsM , Npt , Nd , j , d ) ;
		float dif = cpt - cj ;

		sqrdis = sqrdis + dif * dif ;
	      }

	    /* +++ */ if ( sqrdis == 0.0 ) 
	      printf( "*** Warning : Points %d and %d have same location\n" ,
		      Pt+1 , j+1 ) ;

	    SqrDisV[ j ] = sqrdis ;
	    if ( sqrdis < dnear ) 
	      dnear = sqrdis ;
	  }
	else
	  SqrDisV[ j ] = 0.0 ;
      }

    /* Compute distance to RankNei'th nearest location */
    {
      float* sortdis_1n = calloc( Npt , sizeof( float ) ) ;

      if ( sortdis_1n == NULL )
	{
	  (*SqrDnearP) = dnear ;
	  return ;
	}

      /* Sort the distances to other locations 
       */
      memcpy( sortdis_1n , SqrDisV , Npt * sizeof( float ) ) ;

      qsort( sortdis_1n , Npt, sizeof( float ), floatcompare ) ;

      (*SqrDnearP) = sortdis_1n[ RankNei - 1 ] ;

      free( sortdis_1n ) ;
    }
}


/* ------------------------------------------------------------------- */
/* Thresholds the computed distances to get neighbors and weights */
void SqrDisVToNeiV
 ( 
  int           Pt ,      /* I : index of current point 0..Npt-1 */
  float         SqrDth ,  /* I : squared distance threshold */ 
  const float*  SqrDisV , /* I : [Npt] distances to other points */
  int           Npt ,     /* I : number of points */
  WeiModT       WeiMod ,  /* I : weight computation mode */
  float         WeiFac ,  /* I : weight normalizing factor */
  int*          NbNeiP ,  /* O : number of neighbors 0..Npt-1 */
  int*          NeiV ,    /* O : [Npt] neighbors indices 0..Npt-1 */
  float*        WeiV ,    /* O : [Npt] neighbors weights > 0.0 */
  float*        SumWeiP   /* O : sum of weights > 0.0 */
 ) 
/* ------------------------------------------------------------------- */
{
    int   j ;

    for ( j = 0 , 
	    (*NbNeiP) = 0,
	    (*SumWeiP) = 0.0 ; 
	  j < Npt ; 
	  j ++ )
      {
	if ( j != Pt )
	  {
	    if ( SqrDisV[ j ] <= SqrDth )
	      {
		NeiV[ (*NbNeiP) ] = j ;

		switch( WeiMod )
		  {
		  case WEI_CONST:
		    WeiV[ (*NbNeiP) ] = 1.0 ;
		    break ;

		  case WEI_EXP:
		    WeiV[ (*NbNeiP) ] = exp( - EXPFAC * 
					     SqrDisV[ j ] / SqrDth ) ;
		    break ;

		  default:
		    WeiV[ (*NbNeiP) ] = 1.0 ;
		  }
		WeiV[ (*NbNeiP) ] *= WeiFac ;

		(*SumWeiP) += WeiV[ (*NbNeiP) ] ;

		(*NbNeiP) ++ ;
	      }
	  }
      }
}



/* ------------------------------------------------------------------- */
int SaveNeiV
 ( 
  FILE*         Fnei ,    /* I/O : file to write into */
  int           Pt ,      /* I : index of current point 0..Npt-1 */
  int           NbNei ,   /* I : number of neighbors 0..Npt-1 */
  const int*    NeiV      /* I : [Npt] neighbors indices 0..Npt-1 */
 ) 
/* ------------------------------------------------------------------- */
{
    int nei ;
    int ok ;

    fprintf( Fnei , "%4d  %3d  " , Pt + 1 , NbNei ) ;
    for ( nei = 0 , ok = TRUE ; ( nei < NbNei ) && ok ; nei ++ )
      {
	ok = ( fprintf( Fnei , "%4d " , NeiV[ nei ] + 1 ) > 0 ) ;
      }

    return ( ok ? 0 : -1 ) ;
}



int SaveWeiV
 ( 
  FILE*         Fnei ,    /* I/O : file to write into */
  int           NbNei ,   /* I : number of neighbors 0..Npt-1 */
  const float*  WeiV      /* I : [Npt] neighbors weights > 0.0 */
 ) 
{
    int nei ;
    int ok ;

    for ( nei = 0 , ok = TRUE ; ( nei < NbNei ) && ok ; nei ++ )
      {
	ok = ( fprintf( Fnei , " %4.2f" , WeiV[ nei ] ) > 0 ) ;
      }

    return ( ok ? 0 : -1 ) ;
}



void PrintHelp( const char* CmdName )
{
  printf( "\n" ) ;
  printf( "This program computes a neighborhood graph given a set of \n" ) ;
  printf( "spatial coordinates. It saves the resulting graph in a \n" ) ;
  printf( "neighborhood file which may be used as input to the \n" ) ;
  printf( "program nem_exe. It uses a simple thresholding of the \n" ) ;
  printf( "euclidean distances between the objects.\n" ) ;
  printf( "\n" ) ;
  printf( "Before running %s, the spatial coordinates should be given in \n" ,
	  CmdName ) ;
  printf( "an ASCII file, 1 line/object, 1 column/spatial coordinate.\n" ) ;
  printf( "Example (spatial coordinates in columns 2 and 3 of file) : \n" ) ;
  printf( "   x 15 50 x x   => 1st object is at position x = 15  y = 50\n" ) ;
  printf( "   x 30 20 x x   => 2nd object is at position x = 30  y = 20\n" ) ;
  printf( "The elements on a line should be separated by spaces or tabs\n" ) ;
  printf( "\n" ) ;
  printf( "--- Press ENTER for more ---\n" ) ;
  getchar( ) ;
  printf( "The program will prompt you for :\n" ) ;
  printf( "1 - The name of the input spatial coordinates file\n" ) ;
  printf( "2 - The name of the output neighborhood file\n" ) ;
  printf( "3 - How many spatial coordinates and their column numbers\n" ) ;
  printf( "4 - The desired average number of neighbors of an object ;\n" ) ;
  printf( "    this parameter is used to compute the distance threshold\n" ) ;
  printf( "5 - How to calculate the weights : no weigths or exponential\n" ) ;
  printf( "    exponential means  w_ij = A * exp - 3 (d_ij / d_th)^2\n" ) ;
  printf( "    i.e. weights decreasing as  exp (- squared distance).\n" ) ;
  printf( "    Constant A is computed to have max_i( sum_j w_ij ) = 4,\n" ) ;
  printf( "    as in 4 nearest neighbor unweighted graphs\n" ) ;
  printf( "\n" ) ;
}


void PrintVersion( const char* CmdName )
{
  printf( "\n" ) ;
  printf( "Version 1.00 (14-NOV-1997)\n" ) ;
  printf( "==========================\n" ) ;
  printf( "14-NOV-1997\n" ) ;
  printf( "First complete version. \n" ) ;
  printf( "Added help, computation of distance threshold and weights.\n" ) ;
  printf( "\n" ) ;
}
