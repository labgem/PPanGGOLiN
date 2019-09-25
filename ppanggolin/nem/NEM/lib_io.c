/*\
    lib_io.c

    Routines for input/output and memory allocation

    June 96

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

    Vers-mod  Date         Description

    1.05-a    26-JAN-1998  Remove temporarily added alloc_exit to lib_io.c
\*/

#include <stdio.h>   /* FILE */
#include <string.h>  /* strcpy */
#include <stdlib.h>  /* FILE */

#include "lib_io.h"   /* prototypes of exported functions */

#define TRUE       1
#define FALSE      0

#define LEN_LINE   500
#define LEN_FIELD  30



/* ------------------------------------------------------------------- */
int   ReadOpeningComments  /* 0 : OK, -1 : can't open, 1 : LenComment small */
      (
       const char*   FileName,    /* I : name of file to open */
       const char*   MarkerS,     /* I : comment marker of beginning */
       int           LenComment,  /* I : length of allocated CommentS */
       FILE**        FP ,         /* O : opened file, NULL if impossible */
       char*         CommentS     /* O : read comment */
      )
/* ------------------------------------------------------------------- */
{
    int     mlen = strlen( MarkerS ) ;
    int     iscomment ;
    int     nblines ;
    int     sts ;
    char    line[ LEN_LINE + 1 ] ;
    int     iline ;

    if ( ( (*FP) = fopen( FileName , "r" ) ) == NULL )
      return -1 ;

    /* Read comment lines */
    strcpy( CommentS, "" ) ;
    for ( iscomment = TRUE , nblines = 0 , sts = 0 ; 
          iscomment && (! feof( (*FP) )) ; 
          nblines ++ )
    {
        if ( fgets( line, LEN_LINE, (*FP) ) != NULL )
        {
            iscomment = ( strstr( line , MarkerS ) == line ) ;
            if ( iscomment )
            {
              if ( sts == 0 )
                {
                  if ( (int)( strlen( CommentS ) + strlen( &line[ mlen ] ) ) > 
                      LenComment )
                    sts = 1 ;
                
                  strncat( CommentS, &(line[ mlen ]) , LenComment ) ;
                }
            }
        }
    }
    nblines -- ;

    fclose( (*FP) ) ;

    /* Skip comment lines */
    (*FP) = fopen( FileName, "r" ) ;
    for ( iline = 0 ; iline < nblines ; iline ++ )
      {
        fgets( line , LEN_LINE , (*FP) ) ;
      }

    return sts ;

} /* end of ReadOpeningComments() */

/* ------------------------------------------------------------------- */
int AskFileToRead           /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of file to read */
        char*       NameF   /* O : name of (readable) file */
    )
/* ------------------------------------------------------------------- */
{
    int     exists ;    /* TRUE if last entered file name exists */
    int     nbask ;     /* nb of repeated asks */

    for ( exists = FALSE , nbask = 1 ; 
          ( ! exists ) && ( nbask <= MAX_ASK ) ; 
          nbask ++ )
    {
        if ( nbask == 1 )
            printf( "Name of  %s  file  (RETURN to quit) : ", Desc ) ;

        gets( NameF ) ;

        if ( strlen(NameF) != 0 )
        {
            FILE*   f ;

            if ( ( f = fopen( NameF , "r" ) ) != NULL )
            {
                fclose( f ) ;
                exists = TRUE ;
            }
            else
            {
                printf( " '%s' does not exist. " , NameF ) ;
                if ( nbask < MAX_ASK )
                    printf( "Please type again : " ) ;
                else
                    printf( "\n" ) ;
                exists = FALSE ;
            }
        }
        else nbask = MAX_ASK ;
    }

    if ( exists )
       return 0 ;   /* last file typed in exists */
    else
        return -1 ; /* MAX_ASK or more unsuccessful tries */

} /* end of AskFileToRead() */

/* ------------------------------------------------------------------- */
int AskFileToWrite          /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of file to write */
        int         Conf ,  /* I : TRUE if ask confirmation for overwrite */
        char*       NameF   /* O : name of (readable) file */
    )
/* ------------------------------------------------------------------- */
{
    int     writeok ; /* TRUE if entered file name creation ok */
    int     nbask ;   /* nb of repeated asks */

    for ( writeok = FALSE , nbask = 1 ; 
          ( ! writeok ) && ( nbask <= MAX_ASK ) ; 
          nbask ++ )
    {
        FILE*   f ;
        int     accept ;

        printf( "Name of  %s  file to create : " , Desc ) ;
        gets( NameF ) ;

        if ( strlen(NameF) != 0 )
        {
            /* By default, accept file overwriting ; this acceptation
               will be unvalidated if confirmation asked, and
               file already exists, and user refuses to overwrite file
            */
            accept = TRUE ;
            if ( ( Conf ) && ( ( f = fopen( NameF , "r" ) ) != NULL ) )
            {
                char c ;

                fclose( f ) ;   /* after successful open "r" , close file */

                printf( "File %s already exists. Overwrite it ? (y/n/q) " ,
                        NameF ) ;

                c = getchar() ;
                getchar() ; /* to empty buffer */
                switch( c )
                {
                    case 'y' : accept = TRUE ; break ;
                    case 'q' : accept = FALSE ; nbask = MAX_ASK ; break ;
                    default  : accept = FALSE ;
                }
            }

            if ( accept )
            {
                /* Try to create file, if fails : directory may be invalid */
                if ( ( f = fopen( NameF , "w" ) ) != NULL )
                {
                    fclose( f ) ;   /* after open "w" , close file */
                    remove( NameF ) ; /* remove file */
                    writeok = TRUE ;
                }
                else
                {
                    printf( " Cannot create '%s' (check name/permission)\n" ,
                            NameF ) ;
                }
            }
        }
        else    nbask = MAX_ASK ;
    }

    if ( writeok )
       return 0 ;   /* last file typed in could be created */
    else
        return -1 ; /* MAX_ASK or more unsuccessful tries */

} /* end of AskFileToWrite() */


/* ------------------------------------------------------------------- */
int AskInteger              /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of number to type in */
        int         Def ,   /* I : default value */
        int         Min ,   /* I : minimum value */
        int         Max ,   /* I : maximum value */
        int*        NbReadP /* O : number read */
    )
/* ------------------------------------------------------------------- */
{
    int     numberok ; /* TRUE if entered file name creation ok */
    int     nbask ;    /* nb of repeated asks */

    for ( numberok = FALSE , nbask = 1 ; 
          ( ! numberok ) && ( nbask <= MAX_ASK ) ; 
          nbask ++ )
    {
        char    stringread[ 132 + 1 ] ;

        printf( "Enter  %s  ( %d <= n <= %d )  [%d]  : " , 
                Desc , Min , Max , Def ) ;
        gets( stringread ) ;

        if ( strlen( stringread ) != 0 )
        {
            if ( ( sscanf( stringread , "%d" , NbReadP ) == 1 ) &&
                 ( Min <= (*NbReadP) ) && ( (*NbReadP) <= Max )  )
            {
                numberok = TRUE ;
            }
            else printf( " Invalid number\n" ) ;
        }
        else
        {
            (*NbReadP) = Def ;
            numberok = TRUE ;
        }
    }

    if ( numberok )
       return 0 ;   /* last file typed in could be created */
    else
        return -1 ; /* MAX_ASK or more unsuccessful tries */

} /* end of AskInteger() */


/* ------------------------------------------------------------------- */
int AskFloat                /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of number to type in */
        float       Def ,   /* I : default value */
        float       Min ,   /* I : minimum value */
        float       Max ,   /* I : maximum value */
        float*      NbReadP /* O : number read */
    )
/* ------------------------------------------------------------------- */
{
    int     numberok ; /* TRUE if entered file name creation ok */
    int     nbask ;    /* nb of repeated asks */

    for ( numberok = FALSE , nbask = 1 ; 
          ( ! numberok ) && ( nbask <= MAX_ASK ) ; 
          nbask ++ )
    {
        char    stringread[ 132 + 1 ] ;

        printf( "Enter  %s  ( %g <= n <= %g )  [%g]  : " , 
                Desc , Min , Max , Def ) ;
        gets( stringread ) ;

        if ( strlen( stringread ) != 0 )
        {
            if ( ( sscanf( stringread , "%f" , NbReadP ) == 1 ) &&
                 ( Min <= (*NbReadP) ) && ( (*NbReadP) <= Max )  )
            {
                numberok = TRUE ;
            }
            else printf( " Invalid number\n" ) ;
        }
        else
        {
            (*NbReadP) = Def ;
            numberok = TRUE ;
        }
    }

    if ( numberok )
       return 0 ;   /* last file typed in could be created */
    else
        return -1 ; /* MAX_ASK or more unsuccessful tries */

} /* end of AskFloat() */


/* ------------------------------------------------------------------- */
int CountTokens               /* ret : nb of tokens in Line */
    (
        const char* Line ,    /* I : line to analyze */
        const char* SeparS    /* I : separator between tokens */
    ) 
/* ------------------------------------------------------------------- */
{
    int            NbTokens ;                /* to be returned */
    static char    myline[ LEN_LINE + 1 ] ;  /* static to avoid stacking */
    int            len ;
    char*          p;

    /* Copy to local string, and eventually strip off newline char */
    strncpy( myline , Line , LEN_LINE ) ;
    len = strlen( myline ) ;
    if ( myline[ len - 1 ] == '\n' )
       myline[ len - 1 ] = '\0' ;

    for ( NbTokens = 0 , p = strtok( myline , SeparS ) ;
          p != NULL ;
          p = strtok( NULL , SeparS ) )
    {
        NbTokens ++ ;
    }

    return NbTokens ;

} /* end of CountTokens() */


/* ------------------------------------------------------------------- */
int CountLinesColumns          /* 0/1 = same/dif. nb of col, -1 = problem */
    (
        const char* NameF ,    /* I : file name to analyze */
        const char* SeparS ,   /* I : separator between columns */
        int*        MinColP ,  /* O : minimum number of columns*/
        int*        MaxColP ,  /* O : maximum number of columns*/
        int*        NbLinesP   /* O : number of lines */
    )
/* ------------------------------------------------------------------- */
{
    FILE*   finp ;                  /* input file handler */
    char    line[ LEN_LINE + 1 ] ;  /* last line read */
    int     nblines ;
    int     mincols=0 , maxcols=0 ;
    int     noteq ;


    if ( ( finp = fopen( NameF , "r" ) ) == NULL )
    {
        printf( "Error : can't open file %s\n" , NameF ) ;
        return -1 ;
    }

    /* Read first line */
    nblines = 0 ;
    if ( fgets( line , LEN_LINE , finp ) != NULL )
    {
        maxcols = CountTokens( line , SeparS ) ;
        mincols = maxcols ;
        if ( maxcols > 0 )
           nblines ++ ;
    }

    /* Read following lines */
    noteq = FALSE ;
    while ( ! feof( finp ) )
    {

        if ( fgets( line , LEN_LINE , finp ) != NULL )
        {
            int nbcols = CountTokens( line , SeparS ) ;

            if ( nbcols > 0 )
            {
                nblines ++ ;
                if ( nbcols != maxcols )
                {
                    noteq = TRUE ;
                    if ( nbcols > maxcols )
                       maxcols = nbcols ;
                    else
                        mincols = nbcols ;
                }
            }
        }
    }

    (*NbLinesP) = nblines ;
    (*MinColP) = mincols ;
    (*MaxColP) = maxcols ;

    /* Close file */
    fclose( finp ) ;

    if ( noteq )
       return 1 ;
    else
        return 0 ;

} /* end of CountLinesColumns() */


/* ------------------------------------------------------------------- */
int ReadSelectedColumns        /* 0 = OK , -1 = problem */
    (
        const char* NameF ,    /* I : name of file to read */
        int         Npt ,      /* I : number of lines to read */
        int         Ntot ,     /* I : total number of columns per line */
        int         Nsel ,     /* I : number of selected columns */
        const int*  SelCol ,   /* I : selected columns [ at least Nsel ] */
        float*      PtsM       /* O : points read [ Npt * Nsel ] */
    )
/* ------------------------------------------------------------------- */
{
    FILE*   finp ;  /* input file handler */
    int     ok ;    /* TRUE if file format is correct */
    int     i ;     /* current line :      0..Npt-1 */
    int     c ;     /* current column :    0..Ntot-1 */
    int     sel ;   /* current selection : 0..Nsel-1 */
    char    field[ LEN_FIELD + 1 ] ;    /* last token read */


    /* Open file */
    if ( ( finp = fopen( NameF , "r" ) ) == NULL )
    {
        printf( " Error : can't open file %s\n" , NameF ) ;
        return -1 ;
    }

    /* Read field by field */
    for ( i = 0 , ok = TRUE ; ( i < Npt ) && ok ; i ++ )
    {
        for ( c = 0 ; ( c < Ntot ) && ok ; c ++ )
        {
            if ( fscanf( finp , "%s" , field ) == 1 )
            {
                float x ;
                int   isfloat = ( sscanf( field , "%f" , &x ) == 1 ) ;

                for ( sel = 0 ; ( sel < Nsel ) && ok ; sel ++ )
                {
                    if ( SelCol[ sel ] == c )
                    {
                        if ( isfloat )
                        {
                            PtsM[ ( i * Nsel ) + sel ] = x ;
                        }
                        else
                        {
                          printf( " In '%s', [%d,%d] = '%s' not a number\n" ,
                                    NameF , i + 1 , c + 1 , field ) ;
                          ok = FALSE ;
                        }
                    }
                }
            } /* end - if read successful for element [i,c] */
            else 
            {
                printf( " File '%s', cannot read line %d, column %d\n" ,
                        NameF , i + 1 , c + 1 ) ;
                ok = FALSE ;
            }
        } /* end - for each column c */
    } /* end - for each line i */

    /* Close file */
    fclose( finp ) ;

    if ( ok )
        return 0 ;
    else
        return -1 ;

} /* end of ReadSelectedColumns() */


/* ------------------------------------------------------------------- */
/* Test of routine ReadOpeningComments 

#define FNAME         "t.dat"
#define LEN_COMMENT   1000

main()
{
  int   sts ;
  FILE* f ;
  char  comS[ LEN_COMMENT + 1 ] ;
  char  line[ LEN_LINE + 1 ] ;

  sts = ReadOpeningComments( FNAME , "//" , LEN_COMMENT , &f , comS ) ;

  if ( sts != -1 )
    {
      printf( "Comments %s of file %s :\n%s\n" , 
             (sts == 1) ? "(shortened)" : "",  FNAME , comS ) ;
      printf( "Remaining lines : \n" ) ;

      while( !feof( f ) )
        {
          if ( fgets( line , LEN_LINE , f ) != NULL ) 
            printf( line ) ;
        }
      fclose( f ) ;
    }
  return 0 ;
}

*/

/* ------------------------------------------------------------------- */
/* Test of routine AskFileToRead 

main()
{
    int     sts ;
    char    fname[ LEN_FILE + 1 ] ;

    sts = AskFileToRead( "test" , fname ) ;

    printf( "*** AskFileToRead returned file name '%s' (status %d)\n" ,
            fname , sts ) ;
    return sts ;
}

*/


/* ------------------------------------------------------------------- */
/* Test of routine AskFileToWrite

main()
{
    int     sts ;
    char    fname[ LEN_FILE + 1 ] ;

    sts = AskFileToWrite( "test" , TRUE, fname ) ;

    printf( "*** AskFileToWrite returned file name '%s' (status %d)\n" ,
            fname , sts ) ;
    return sts ;
}

*/

/* ------------------------------------------------------------------- */
/* Test of routine AskInteger 
#include <stdlib.h>
main( int argc , char *argv[] )
{
    int     sts ;
    int     n ;

    if ( argc < 4 ) return 1 ;

    sts = AskInteger( "test number" , atoi( argv[1] ) , atoi( argv[2] ) ,
                      atoi( argv[3] ) , & n ) ;

    printf( "*** AskInteger returned '%d' (status %d)\n" ,
            n , sts ) ;
    return sts ;
}

*/

/* ------------------------------------------------------------------- */
/* Test of routine AskFloat 
#include <stdlib.h>
main( int argc , char *argv[] )
{
    int     sts ;
    float   x ;

    if ( argc < 4 ) return 1 ;

    sts = AskFloat( "test number" , atof( argv[1] ) , atof( argv[2] ) ,
                      atoi( argv[3] ) , & x ) ;

    printf( "*** AskFloat returned '%f' (status %d)\n" ,
            x , sts ) ;
    return sts ;
}

*/

/* ------------------------------------------------------------------- */
/* Test of routine CountTokens : 2 command line args = args of function 
main( int argc , char *argv[] )
{
    int     sts ;

    if ( argc < 3 ) return 1 ;

    sts = CountTokens( argv[1] , argv[2] ) ;

    printf( "*** CountTokens( '%s' , '%s' )  returned '%d'\n" , 
            argv[1] , argv[2] , sts ) ;
    return sts ;
}
*/

/* ------------------------------------------------------------------- */
/* Test of routine CountLinesColumns : 
   2 command line args = args of function 

main( int argc , char *argv[] )
{
    int     sts ;
    int     minc , maxc , nbl ;

    if ( argc < 3 ) return 1 ;

    sts = CountLinesColumns( argv[1] , argv[2] , &minc , &maxc , &nbl ) ;

    printf( "*** CountLinesColumns( '%s' , '%s' ) returned '%d'\n" , 
            argv[1] , argv[2] , sts ) ;
    printf( "*** minc = %d   maxc = %d    nbl = %d\n" , minc,maxc,nbl ) ;

    return sts ;
}

*/


/* ------------------------------------------------------------------- */
/* Test of routine ReadSelectedColumns : 
   command line args : 
   1 = filename , 2 = nb sel. col , 3, 4, ... = sel. col

#include <stdlib.h>
main( int argc , char *argv[] )
{
    int     sts ;
    int     minc , nbl , nbc , nbsel , s ;
    int     selcV[ 5 ] ;
    float*  ptM ;

    if ( argc < 4 ) return 1 ;
    nbsel = atoi( argv[2] ) ;
    if ( nbsel > 5 ) return 1 ;
    for ( s = 0 ; s < nbsel ; s ++ )
        selcV[ s ] = atoi( argv[ 3 + s ] ) - 1 ;

    sts = CountLinesColumns( argv[1] , " " , &minc , &nbc , &nbl ) ;
    if ( sts != 0 )
    { printf( "Error : CountLinesColumns returned %d\n" , sts ) ;
      return 2 ; }

    if ( ( ptM = malloc( nbl * nbsel * sizeof( float ) ) ) == NULL )
       return 3 ;

    sts = ReadSelectedColumns( argv[1] , nbl , nbc , nbsel , selcV , ptM ) ;

    printf( "*** ReadSelectedColumns( '%s' , %d , %d , %d , [%d,%d,%d] ) returned '%d'\n" , 
            argv[1] , nbl , nbc , nbsel , selcV[ 0 ] , selcV[ 1 ] , 
            selcV[ 2 ] , sts ) ;

    {
        int i , s ;

        for ( i = 0 ; i < nbl ; i ++ )
        {
            printf( "*** Pt %d : " , i + 1 ) ;
            for ( s = 0 ; s < nbsel ; s ++ )
                printf( "  %4.2f" , ptM[ i * nbsel + s ] ) ;
            printf( "\n" ) ;
        }
    }

    free( ptM ) ;
    return sts ;
}

*/


/* ------------------------------------------------------------------- */
