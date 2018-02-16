/*\
    txt2hlp.c

    This utility converts a text file into a C file consisting of a
    'help' function. "Batch" version (parameters given by command line
    arguments, no keyboard input requested).

    New version 2.0 11-JAN-1998 : allows several help functions.
\*/

#include <stdio.h>      /* printf, FILE */
#include <string.h>     /* strcpy */
#include <time.h>       /* time */


#define TRUE        1
#define FALSE       0

#define LEN_FILE    100
#define LEN_FUNC    32 
#define LEN_LINE    100
#define LEN_COMMENT 500

#define NB_PARA     4

#define MARKER_CHAR '%'

main( int Argc , char* Argv[] )
{
    char    ftextName[ LEN_FILE + 1 ] ;
    char    fhelpName[ LEN_FILE + 1 ] ;
    char    funcName[ LEN_FUNC + 1 ] ;
    char    CommentS[ LEN_COMMENT + 1 ] ;
    char    rep ;
    FILE*   fh ;

    int     hFileYes ;

    FILE*   ftext ;

    void PrintUsage( const char* Cmd ) ;


    printf( "\nWelcome to program TXT2HLP\n\n" ) ;


    /* Fetch command line arguments */
    if ( ( Argc - 1 ) == 0 )
      {
	/* No args -> print help */
	PrintUsage( Argv[ 0 ] ) ;
	return 1 ;
      }

    if ( ( Argc - 1 ) < NB_PARA )
      {
	printf( "Error : only %d parameters (%d requested)\n" ,
	         ( Argc - 1 ) ,  NB_PARA ) ;
	PrintUsage( Argv[ 0 ] ) ;
	return 2 ;
      }


    strncpy( ftextName , Argv[ 1 ] , LEN_FILE ) ;
    strncpy( fhelpName , Argv[ 2 ] , LEN_FILE ) ;
    rep                = Argv[ 3 ][ 0 ] ;
    strncpy( CommentS  , Argv[ 4 ] , LEN_COMMENT ) ;

    printf( "Input text file name  : %s\n" , ftextName ) ;
    printf( "Output help file name : %s.c\n" , fhelpName ) ;
    printf( "___.h file ? (y/n)    : %c\n" , rep ) ;
    hFileYes = ( rep == 'y' ) ;
    printf( "Comment               : %s\n" , CommentS ) ;
    printf( "\n" );


    /* Open input file */
    if ( ( ftext = fopen( ftextName , "r" ) ) == NULL )
    {
        printf( "Error : cannot read file %s\n" , ftextName ) ;
        return 2 ;
    }


    /* If ___.h file requested, create it */
    if ( hFileYes )
    {
        char    fhName[ LEN_FILE + 1 ] ;

        strcpy( fhName , fhelpName ) ;
        strncat( fhName ,  ".h" , LEN_FILE ) ;

        /* Open ___.h file */
        if ( ( fh = fopen( fhName , "w" ) ) == NULL )
        {
            printf( "Error : cannot write file %s\n" , fhName ) ;
            fclose( ftext ) ;
            return 2 ;
        }

        printf( "Writing file %s ...\n" , fhName ) ;

	/* Start conditional inclusion of .h file */
        fprintf( fh , "#ifndef %s_H\n", fhelpName ) ;
	fprintf( fh , "#define %s_H\n\n", fhelpName ) ;

        /* Include necessary ___.h files */
        fprintf( fh , "#include <stdio.h>  /* FILE */\n\n" ) ;

    }


    /* Treat C file */
    {
        char    fcName[ LEN_FILE + 1 ] ;
        FILE*   fc ;
	time_t  timer = time( NULL ) ;
	int     in_function ;  /* 1 if in a function, 0 else */

        strcpy( fcName , fhelpName ) ;
        strncat( fcName ,  ".c" , LEN_FILE ) ;

        /* Open ___.c file */
        if ( ( fc = fopen( fcName , "w" ) ) == NULL )
        {
            printf( "Error : cannot write file %s\n" , fcName ) ;
            fclose( ftext ) ;
            return 2 ;
        }

        printf( "Writing file %s ...\n" , fcName ) ;

        /* Start ___.c file */
        fprintf( fc , "/*\\\n\n" ) ;

        fprintf( fc , "    %s.c\n\n" , fhelpName ) ;

        fprintf( fc , "    %s\n\n" , CommentS ) ;

	fprintf( fc , "    %s\n", asctime( localtime( &timer ) ) ) ;

        fprintf( fc , "\\*/\n\n" ) ; 

        fprintf( fc , "#include \"%s.h\"  /* Exported prototypes */\n\n" , 
                 fhelpName ) ;

        /* Scan lines from text file to ___.c file */
	in_function = 0 ;
	funcName[ 0 ] = '\0' ;
        while ( ! feof( ftext ) )
        {
            char    line[ LEN_LINE + 1 ] ;

            if ( fgets( line , LEN_LINE , ftext ) != NULL )
	      {
		/* Suppress newline character */
                int     len     = strlen( line ) ;
                line[ len - 1 ] = '\0' ;

		/* If this line is a structure marker */
		if ( line[ 0 ] == MARKER_CHAR )
		  {
		    /* If not currently in function, start a new function 
		       and include it in header file */
		    if ( ! in_function )
		      {
			in_function = 1 ;
			strncpy( funcName, &line[ 1 ], LEN_FUNC ) ;
			fprintf( fc , "void %s( FILE* F )\n{\n" , funcName ) ;
			fprintf( fh , "extern void %s( FILE* F ) ;\n\n" , 
				 funcName ) ;
		      }
		    else /* currently in a function, end it */
		      {
			in_function = 0 ;
			fprintf( fc , "} /* end of %s() */\n\n\n", funcName ) ;
		      }
		  }
		else /* not a structure marker, if in function, 
			make print out instruction */
		  {
		    if ( in_function )
		      fprintf( fc , "    fprintf( F , \"%s\\n\" ) ;\n" , 
			       line ) ;
		  }

	      } /* end if fgets() != NULL */
        }

	/* If current function not ended, end it */
	if ( in_function )
	  {
	    printf( "Warning : function %s not ended, I will do it\n", 
		    funcName ) ;
	    fprintf( fc , "} /* end of %s() */\n\n\n", funcName ) ;
	  }

        /* Close ___.c file */
        fclose( fc ) ;
    }


    /* If ___.h file requested, close it */
    if ( hFileYes )
      {
	/* Close conditional inclusion of .h file */
        fprintf( fh , "#endif\n" ) ;

        /* Close ___.h file */
        fclose( fh ) ;
      }

    /* Close input file */
    fclose( ftext ) ;

    printf( "\nBye\n" ) ;

    return 0 ;
} /* end of main() */



/* -------------------------------------------------------------- */
void PrintUsage( const char* Cmd )
{
    printf( "\nSyntax :\n" ) ;
    printf( "   %s   in_file  out_file  hfile_yn  comment\n\n" , 
	    Cmd ) ;
}



