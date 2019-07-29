/*\
    lib_io.h

    Prototypes of routines for input/output

    June 96

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

    Vers-mod  Date         Description

    1.05-a    26-JAN-1998  Add ExitET
\*/

#include <stdio.h>   /* FILE */

#define LEN_FILE     132     /* maximum length of file name */
#define MAX_ASK      5       /* maximum number of repeated asked inputs */

typedef enum
{

  EXIT_OK,           /* Good:    program achieved processing normally */
  EXIT_W_RESULT,     /* Warning: program achieved with unusable result */
  EXIT_E_ARGS,       /* Error:   user gave wrong argument syntax */ 
  EXIT_E_FILE,       /* Error:   files not found or wrong format */
  EXIT_E_MEMORY,     /* Error:   program ran out of memory */
  EXIT_E_SYSTEM,     /* Error:   a system call failed */
  EXIT_E_BUG,        /* Error:   internal program inconsistency */
  EXIT_NB

} ExitET ;

/* ------------------------------------------------------------------- */
int   ReadOpeningComments  /* 0 : OK, -1 : can't open, 1 : LenComment small */
      (
       const char*   FileName,    /* I : name of file to open */
       const char*   MarkerS,     /* I : comment marker of beginning */
       int           LenComment,  /* I : length of allocated CommentS */
       FILE**        FP ,         /* O : opened file, NULL if impossible */
       char*         CommentS     /* O : read comment */
      ) ;

/* ------------------------------------------------------------------- */
int AskFileToRead           /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of file to read */
        char*       NameF   /* O : name of (readable) file */
    ) ;

/* ------------------------------------------------------------------- */
int AskFileToWrite          /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of file to write */
        int         Conf ,  /* I : TRUE if ask confirmation for overwrite */
        char*       NameF   /* O : name of (readable) file */
    ) ;

/* ------------------------------------------------------------------- */
int AskInteger              /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of number to type in */
        int         Def ,   /* I : default value */
        int         Min ,   /* I : minimum value */
        int         Max ,   /* I : maximum value */
        int*        NbReadP /* O : number read */
    ) ;

/* ------------------------------------------------------------------- */
int AskFloat                /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of number to type in */
        float       Def ,   /* I : default value */
        float       Min ,   /* I : minimum value */
        float       Max ,   /* I : maximum value */
        float*      NbReadP /* O : number read */
    ) ;

/* ------------------------------------------------------------------- */
int CountLinesColumns          /* 0/1 = same/dif. nb of col, -1 = problem */
    (
        const char* NameF ,    /* I : file name to analyze */
        const char* SeparS ,   /* I : separator between columns */
        int*        MinColP ,  /* O : minimum number of columns*/
        int*        MaxColP ,  /* O : maximum number of columns*/
        int*        NbLinesP   /* O : number of lines */
    ) ;

/* ------------------------------------------------------------------- */
int CountTokens               /* ret : nb of tokens in Line */
    (
        const char* Line ,    /* I : line to analyze */
        const char* SeparS    /* I : separator between tokens */
    ) ;

/* ------------------------------------------------------------------- */
int ReadSelectedColumns        /* 0 = OK , -1 = problem */
    (
        const char* NameF ,    /* I : name of file to read */
        int         Npt ,      /* I : number of lines to read */
        int         Ntot ,     /* I : total number of columns per line */
        int         Nsel ,     /* I : number of selected columns */
        const int*  SelCol ,   /* I : selected columns [ at least Nsel ] */
        float*      PtsM       /* O : points read [ Npt * Nsel ] */
    ) ;
