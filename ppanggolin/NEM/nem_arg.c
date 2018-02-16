/*\

    NEM_ARG.C

    Programme NEM (Neighborhood EM) : Traitement des arguments

    Van Mo DANG       Janvier 96

    Steps to add an option :
    1 - check unused option characters by running 'nem_exe' without args
    2 - add option fetching in NemArgs()
    3 - add field to NemParaT in "nem_typ.h"  : OptMeaning
    4 - add eventual enum in "nem_typ.h"      : OptET
    5 - add eventual string array             : OptStrVC
    6 - add option default value              : DEFAULT_OPT
    7 - add comment to this header            : 1.06-k blabla
    8 - add short help in PrintUsage()        : -o [%s ] blabla {%s|%s}
    9 - add long help/vers in "nem_user.txt"

Vers-mod  Date         Who Description

1.03-a    22-AUG-1997  MD  '-l' to output logfile (default now is NO log)
1.04-a    04-OCT-1997  MD  '-s f myfile.zz' instead of just '-s f'
1.04-b    05-OCT-1997  MD  '-B heu <bstep>' implemented
1.04-c    09-OCT-1997  MD  if "-o -" log name is input base name
1.04-d    10-OCT-1997  MD  '-c cvthres' implemented
1.04-e    10-OCT-1997  MD  '-S seed' implemented
1.04-f    10-OCT-1997  MD  '-O order' implemented
1.04-g    13-OCT-1997  MD  '-D' implemented
1.04-h    13-OCT-1997  MD  '-R' implemented (reference class)
1.04-i    13-OCT-1997  MD  '-C crit'
1.04-j    04-NOV-1997  MD  '-I eiter'
1.04-k    11-JAN-1997  MD  now '-h general/filein/fileout/examples/vers'
1.05-a    12-JAN-1998  MD  '-M miss' implemented
1.05-b    05-FEB-1998  MD  add 'M' in '-C crit'
1.06-a    17-JUN-1998  MD  ReadStrFile and alloc parameters <- from nem_exe.c
1.06-b    23-JUN-1998  MD  use new StatModel structure (old NoiseModel)
1.06-b    01-AUG-1998  MD  model switch now '-m norm/lapl p* s**'
1.06-c    03-AUG-1998  MD  add '-s mi' and '-s mf'
1.06-d    10-SEP-1998  MD  add '-U para|seq'
1.06-e    10-SEP-1998  MD  add '-t random|first'
1.06-f    15-SEP-1998  MD  add '-G nbit conv step'
1.06-g    21-SEP-1998  MD  change '-c 0.01' to '-c {none|[clas]|crit} 0.01'
1.06-h    05-OCT-1998  MD  default nb of random inits = D * 5 
1.06-i    05-OCT-1998  MD  default crit = CRIT_M
1.07-a    08-APR-1999  MD  Add Bernoulli family
1.08-a    20-JUI-2017  GG   Add param input by file rather than by arguments
\*/

#include "nem_typ.h"    /* DataT, ... */
#include "nem_ver.h"    /* PrintVersions */
#include "nem_hlp.h"    /* PrintHelpGeneral ... */
#include "lib_io.h"     /* ReadOpeningComments */ /*V1.06-a*/
#include "genmemo.h"    /* GenAlloc, ... */ /*V1.06-a*/
#include <stdio.h>      /* printf, ... */
#include <stdlib.h>     /* atof, ... */
#include <string.h>     /* strncpy, ... */
#include <time.h>       /* time() */          /*V1.04-e*/

#include "nem_arg.h"    /* prototypes */

#define DEFAULT_ALGO         ALGO_NEM
#define DEFAULT_BETA         1.0
#define DEFAULT_BTAMODE      BETA_FIX           /*V1.04-b*/
#define DEFAULT_BTAHEUSTEP   0.1                /*V1.04-b*/
#define DEFAULT_BTAHEUMAX    2.0                /*V1.04-b*/
#define DEFAULT_BTAHEUDDROP  0.8             	/*V1.04-b*/
#define DEFAULT_BTAHEUDLOSS  0.5             	/*V1.04-b*/
#define DEFAULT_BTAHEULLOSS  0.02            	/*V1.04-b*/
#define DEFAULT_BTAGRADNIT   1               	/*V1.06-f*/
#define DEFAULT_BTAGRADCVTH  0.001           	/*V1.06-f*/
#define DEFAULT_BTAGRADSTEP  0.0            	/*V1.06-f*/
#define DEFAULT_BTAGRADRAND  0            	/*V1.06-f*/
#define DEFAULT_CRIT         CRIT_M             /*V1.04-h*/
#define DEFAULT_CVTHRES      0.01               /*V1.04-d*/
#define DEFAULT_CVTEST       CVTEST_CLAS        /*V1.06-g*/
#define DEFAULT_FAMILY       FAMILY_NORMAL      /*V1.06-b*/
#define DEFAULT_DISPER       DISPER___          /*V1.06-b*/
#define DEFAULT_PROPOR       PROPOR__           /*V1.06-b*/
#define DEFAULT_FORMAT       FORMAT_HARD
#define DEFAULT_INIT         INIT_SORT
#define DEFAULT_MISSING      MISSING_REPLACE    /*V1.05-a*/
#define DEFAULT_NO_PARAM_FILE NO_PARAM_FILE     /*V1.08-a*/
#define DEFAULT_SORTEDVAR    0
#define DEFAULT_NEIGHSPEC    NEIGH_FOUR
#define DEFAULT_NBITERS      100
#define DEFAULT_NBEITERS     1
#define DEFAULT_NBRANDINITS  5                  /*V1.06-h*/
#define DEFAULT_ORDER        ORDER_DIRECT       /*V1.04-f*/
#define DEFAULT_UPDATE       UPDATE_SEQ         /*V1.06-d*/
#define DEFAULT_TIE          TIE_RANDOM         /*V1.06-e*/


typedef enum
{
  HELP_GENERAL,
  HELP_OPTIONS,
  HELP_EXAMPLES,
  HELP_FILEIN,
  HELP_FILEOUT,
  HELP_VERSIONS,
  HELP_NB

} HelpET ;


static const char *TypeStrC[ TYPE_NB ] = { "S", "I", "N" } ; /*V1.06-a*/

static const char *AlgoStrVC[ ALGO_NB ] = { "nem" , "ncem" , "gem" } ;
static const char *BtaStrVC[ BETA_NB ] = { "fix" , "psgrad", 
					   "heu_d" , "heu_l" } ;
       const char *CritStrVC[ CRIT_NB ] = { "U" , "M", "D" , "L" } ;
static const char *CvTestStrVC[ CVTEST_NB ] = { "none", "clas" , "crit" } ;
static const char *FormatStrVC[ FORMAT_NB ] = { "hard", "fuzzy" } ;
static const char *HelpStrVC[ HELP_NB ] = {
  "general",
  "options",
  "examples",
  "filein",
  "fileout",
  "versions"
} ; /* V1.04-k*/

static const char *InitStrVC[ INIT_NB ] = { "s", "r", "m", "f", "l" } ;
/*V1.05-a*/
static const char *MissStrVC[ MISSING_NB ] = { "replace" , "ignore" } ;
/*V1.06-b*/
static const char *FamilyStrVC[ FAMILY_NB ] = { "norm", "lapl", "bern" } ;
static const char *DisperStrVC[ DISPER_NB ] = { "s__", "sk_", "s_d", "skd" } ;
static const char *ProporStrVC[ PROPOR_NB ] = { "p_", "pk" } ;
static const char *NeighStrVC[ NEIGH_NB ] = { "4", "f" } ;
static const char *OrderStrVC[ ORDER_NB ] = { "direct", "random" } ;/*V1.04-f*/
static const char *UpdateStrVC[ UPDATE_NB ] = { "seq", "para" } ;/*V1.06-d*/
static const char *TieStrVC[ TIE_NB ] = { "random", "first" } ;/*V1.06-e*/



/* ==================== LOCAL FUNCTION PROTOTYPING =================== */


/* ==================== GLOBAL FUNCTION DEFINITION =================== */


/* ------------------------------------------------------------------- */
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
        ) 
/* ------------------------------------------------------------------- */
{
  StatusET     err = STS_OK ;
  const char*  func = "NemArgs" ;
  int          nk ;
  int          nbopt ;
  int          iopt ;
  const char** opts;        

  int  ReadStrFile  /*V1.06-a*/
         (
             const char*  BaseName,     /* I */
             char*        Comment,      /* O */
             DataT*       DataP,        /* O */
             ImageNeighT* ImageP,       /* O */
             TypeET*      TypeP         /* O */
         );                             

  char GetOptionSwitch( const char *Arg ) ;
  int  GetEnum( const char* S, const char* SV[], int SizeV );/*V1.04-b*/
  StatusET PrintHelp( const char* HelpType ) ;


  switch ( Argc ) /*V1.04-k*/
  {
  case 1:  /* No args : print usage */
    PrintUsage( Argv[ 0 ] ) ;
    return STS_I_DONE ;

  case 2:  /* 1 arg : only possibility is -v to ask version information */
    if ( GetOptionSwitch( Argv[ 1 ] ) == 'v' )
      {
	PrintVersions( stdout ) ;
	return STS_I_DONE ;
      }
    else
      return STS_E_ARG ;

  case 3:  /* 2 args : -h help or file nk */
    switch( GetOptionSwitch( Argv[ 1 ] ) )
      {
      case '?':
      case 'h':
	return PrintHelp( Argv[2] ) ;

      default: /* continue processing args */
             break;
      }

  default:  /* 3 or more args : continue processing args */
    break ;
  }

  strncpy( Fname, Argv[1], LEN_FILENAME ) ;
  sscanf( Argv[2], "%d", &nk ) ;
  StatModelP->Spec.K = nk ;
  if ( nk <= 0 )
  {
      fprintf( stderr, "Nb of classes must be > 0 (here %s)\n",
               Argv[2] ) ;
      return STS_E_ARG ;
  }


  /* !!! Read structure file */ /*V1.06-a*/
  if ( ( err = ReadStrFile( Fname, 
			    DatadescS,
			    DataP, 
			    ImageP,
			    SpatialTypeP ) ) != STS_OK )
    return err ;

  /* !!! Allocate model parameters */ /*V1.06-a*/
  StatModelP->Para.Prop_K    = GenAlloc( nk, sizeof(float), 
					 1, func, "Prop_K" ) ;
  StatModelP->Para.Disp_KD   = GenAlloc( nk * DataP->NbVars, sizeof(float), 
				       1, func, "Disp_KD" ) ;
  StatModelP->Para.Center_KD = GenAlloc( nk * DataP->NbVars, sizeof(float), 
					 1, func, "Center_KD" ) ;

  StatModelP->Para.NbObs_K   = GenAlloc( nk, sizeof(float), 
					 1, func, "NbObs_K" ) ;
  StatModelP->Para.NbObs_KD  = GenAlloc( nk * DataP->NbVars, sizeof(float), 
					 1, func, "NbObs_KD" ) ;
  StatModelP->Para.Iner_KD   = GenAlloc( nk * DataP->NbVars, sizeof(float), 
					 1, func, "NbObs_KD" ) ;

  StatModelP->Desc.DispSam_D = GenAlloc( DataP->NbVars, sizeof(float), 
					  1, func, "DispSam_D" );
  StatModelP->Desc.MiniSam_D = GenAlloc( DataP->NbVars, sizeof(float), 
					  1, func, "MiniSam_D" );
  StatModelP->Desc.MaxiSam_D = GenAlloc( DataP->NbVars, sizeof(float), 
					  1, func, "MaxiSam_D" );

  /* Set default value of optional parameters */
  StatModelP->Spec.ClassFamily = DEFAULT_FAMILY ;
  StatModelP->Spec.ClassDisper = DEFAULT_DISPER ;
  StatModelP->Spec.ClassPropor = DEFAULT_PROPOR ;


  NemParaP->Algo          = DEFAULT_ALGO ;
  StatModelP->Para.Beta   = DEFAULT_BETA ;          /*V1.06-b*/
  StatModelP->Spec.BetaModel = DEFAULT_BTAMODE ;       /*V1.04-b*/
  NemParaP->BtaHeuStep    = DEFAULT_BTAHEUSTEP ;    /*V1.04-b*/
  NemParaP->BtaHeuMax     = DEFAULT_BTAHEUMAX ;
  NemParaP->BtaHeuDDrop   = DEFAULT_BTAHEUDDROP ;
  NemParaP->BtaHeuDLoss   = DEFAULT_BTAHEUDLOSS ;
  NemParaP->BtaHeuLLoss   = DEFAULT_BTAHEULLOSS ;
  NemParaP->BtaPsGrad.NbIter    = DEFAULT_BTAGRADNIT  ;/*V1.06-g*/
  NemParaP->BtaPsGrad.ConvThres = DEFAULT_BTAGRADCVTH ;
  NemParaP->BtaPsGrad.Step      = DEFAULT_BTAGRADSTEP ;
  NemParaP->BtaPsGrad.RandInit  = DEFAULT_BTAGRADRAND ;
  NemParaP->Crit          = DEFAULT_CRIT ;          /*V1.04-h*/
  NemParaP->CvThres       = DEFAULT_CVTHRES ;       /*V1.04-d*/
  NemParaP->CvTest        = CVTEST_CLAS ;           /*V1.06-g*/
  NemParaP->DoLog         = FALSE ;                 /*V1.03-a previously TRUE*/
  NemParaP->NbIters       = DEFAULT_NBITERS ;
  NemParaP->NbEIters      = DEFAULT_NBEITERS ;
  NemParaP->NbRandomInits = DEFAULT_NBRANDINITS*DataP->NbVars ;  /*V1.06-h*/
  NemParaP->Seed          = time( NULL ) ;          /*V1.04-e*/
  NemParaP->Format        = DEFAULT_FORMAT ;
  NemParaP->InitMode      = DEFAULT_INIT ;
  NemParaP->ParamFileMode = DEFAULT_NO_PARAM_FILE ;
  NemParaP->SortedVar     = DEFAULT_SORTEDVAR ;
  NemParaP->NeighSpec     = DEFAULT_NEIGHSPEC ;
  NemParaP->VisitOrder    = DEFAULT_ORDER ;         /*V1.04-f*/
  NemParaP->SiteUpdate    = DEFAULT_UPDATE ;        /*V1.06-d*/
  NemParaP->TieRule       = DEFAULT_TIE ;           /*V1.06-e*/
  NemParaP->Debug         = FALSE ;                 /*V1.04-g*/
  strncpy( NemParaP->OutBaseName, Fname, LEN_FILENAME ) ;
  strncpy( NemParaP->NeighName, Fname, LEN_FILENAME ) ;
  strncat( NemParaP->NeighName, ".nei", LEN_FILENAME ) ;
  strncpy( NemParaP->RefName, "", LEN_FILENAME ) ;
  
  /* Treat options */
  nbopt = Argc - 3 ;
  opts  = & ( Argv[ 3 ] ) ;

  for ( iopt = 0, err = STS_OK ; 
        ( iopt < nbopt ) && ( err == STS_OK ) ;
        iopt ++ )
    {
      char copt = GetOptionSwitch( opts[ iopt ] ) ;

      if ( copt != '\0' ) 
        switch( copt )
        {
        case 'a' :    /* next arg is type of algorithm */
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->Algo = GetEnum( opts[ iopt ], AlgoStrVC, ALGO_NB ) ;
              if ( NemParaP->Algo == -1 )
                {
                  fprintf( stderr, " Unknown type of algorithm %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'b' :    /* next arg is beta */
          iopt ++ ;
          if ( iopt < nbopt )
            {
              StatModelP->Para.Beta = atof( opts[ iopt ] ) ;
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'c' :    /* next arg : conv test */  /*V1.04-d*/ /*V1.06-g*/
	  iopt ++ ;
          if ( iopt < nbopt ) /* conv test given */ {
	    NemParaP->CvTest=GetEnum( opts[ iopt ], CvTestStrVC, CVTEST_NB );
	    if ( NemParaP->CvTest == -1 ) {
	      fprintf( stderr, " Unknown convergence test %s\n", opts[iopt] ) ;
	      err = STS_E_ARG ;
	    }
	    else if ( NemParaP->CvTest != CVTEST_NONE ) /* get threshold */ {
	      iopt ++ ;
	      if ( iopt < nbopt ) /* threshold given */ {
		NemParaP->CvThres = atof( opts[ iopt ] ) ;
		if ( NemParaP->CvThres <= 0 ) {
		  fprintf( stderr, " Conv threshold must be > 0 (here %s)\n",
			   opts[ iopt ] ) ;
		  err = STS_E_ARG ;
		} /* else threshold > 0 : OK */
	      } 
	      else /* threshold not given */ {
		fprintf( stderr, " Expected threshold for conv test %s\n", 
			 opts[ iopt - 1 ] ) ;
		err = STS_E_ARG ;
	      }
            } /* threshold checked */
	  } /* conv test given */
          else {
	    fprintf( stderr, " Expected value after switch %s\n", 
		     opts[ iopt - 1 ] ) ;
	    err = STS_E_ARG ;
	  }
          break ;


        case 'f' :    /* next arg is partition input/output format */
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->Format=GetEnum( opts[ iopt ], FormatStrVC, FORMAT_NB );
              if ( NemParaP->Format == -1 )
                {
                  fprintf( stderr, " Unknown format %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'i' :    /* next arg is number of iterations */
          iopt ++ ;
          if ( iopt < nbopt )
            {
              NemParaP->NbIters = atoi( opts[ iopt ] ) ;
              if ( NemParaP->NbIters < 0 )
                {
                  fprintf( stderr, "Nb iterations must be >= 0 (here %s)\n",
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'l' :    /* next arg is 'y' if log requested */  /*V1.03-a*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
              if ( opts[ iopt ][ 0 ] == 'y' )
      	  NemParaP->DoLog = TRUE ;
      	else
      	  NemParaP->DoLog = FALSE ;
            }
          else
            {
              fprintf( stderr, " Expected y or n after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'm' :    /* next 3 args is stat model */ /*V1.06-b*/
          if ( iopt + 3 < nbopt )
            {
	      iopt ++ ;
	      StatModelP->Spec.ClassFamily = 
		GetEnum( opts[ iopt ], FamilyStrVC, FAMILY_NB );
              if ( StatModelP->Spec.ClassFamily == -1 )
                {
                  fprintf( stderr, " Unknown family %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
	      iopt ++ ;
	      StatModelP->Spec.ClassPropor = 
		GetEnum( opts[ iopt ], ProporStrVC, PROPOR_NB );
              if ( StatModelP->Spec.ClassPropor == -1 )
                {
                  fprintf( stderr, " Unknown proportion %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
	      iopt ++ ;
	      StatModelP->Spec.ClassDisper = 
		GetEnum( opts[ iopt ], DisperStrVC, DISPER_NB );
              if ( StatModelP->Spec.ClassDisper == -1 )
                {
                  fprintf( stderr, " Unknown dispersion %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected 3 values after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'n' :    /* next arg is neighborhood specification */
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->NeighSpec=GetEnum( opts[ iopt ], NeighStrVC, NEIGH_NB);
              if ( NemParaP->NeighSpec == -1 )
                {
                  fprintf( stderr, " Unknown neighborhood specification %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'o' :    /* next arg is output file name */
          iopt ++ ;
          if ( iopt < nbopt )
            {
              strncpy( NemParaP->OutBaseName, opts[ iopt ], LEN_FILENAME ) ;
            }
          else
            {
              fprintf( stderr, " Expected file name after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 's' :    /* next arg is initialization mode */ 
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->InitMode = GetEnum( opts[ iopt ], InitStrVC, INIT_NB );

              switch( NemParaP->InitMode )  {

                case INIT_SORT:  /* next arg is sorted var */
                  iopt ++ ;
                  if ( iopt < nbopt )
                    {
                      NemParaP->SortedVar = atoi( opts[ iopt ] ) - 1 ;
                    } 
                  else
                    {
                      fprintf( stderr, 
                               " Expected variable index after switch %s\n", 
                               opts[ iopt - 1 ] ) ;
                      err = STS_E_ARG ;
                    }
                  break ;

                case INIT_RANDOM:  /* next arg is number of initializations */
                  iopt ++ ;
                  if ( iopt < nbopt )
                    {
                      NemParaP->NbRandomInits = atoi( opts[ iopt ] ) ;
                    } 
                  else
                    {
                      fprintf( stderr, 
                               " Expected nb of trials after switch %s\n", 
                               opts[ iopt - 1 ] ) ;
                      err = STS_E_ARG ;
                    }
                  break ;

                case INIT_FILE:  /* next arg is name of initial partition */
		  iopt ++ ;   /*V1.04-a*/
		  if ( iopt < nbopt )
		    {
		      strncpy( NemParaP->StartName, opts[ iopt ], 
			       LEN_FILENAME ) ;
		    }
		  else
		    {
		      fprintf( stderr, 
			       " Expected file name after switch %s\n", 
			       opts[ iopt - 1 ] ) ;
		      err = STS_E_ARG ;
		    }
                  break ;

		case INIT_LABEL:/* next arg is name of partial labels */
		  iopt ++ ;
		  if ( iopt < nbopt )
		    {
		      strncpy( NemParaP->LabelName, opts[ iopt ], 
			       LEN_FILENAME ) ;
		    }
		  else
		    {
		      fprintf( stderr, 
			       " Expected file name after switch %s\n", 
			       opts[ iopt - 1 ] ) ;
		      err = STS_E_ARG ;
		    }
		  break ;

		case INIT_PARAM_FILE: /* next args are parameter values */
		    
		  iopt ++ ;   /*V1.08-a*/
		  if ( iopt < nbopt )
		    {
		      strncpy( NemParaP->ParamName, opts[ iopt ], 
			       LEN_FILENAME ) ;
		    }
		  else
		    {
		      fprintf( stderr, 
			       " Expected file name after switch %s\n", 
			       opts[ iopt - 1 ] ) ;
		      err = STS_E_ARG ;
		    }
                  break ;
		  break ;

                default:
                  fprintf( stderr, " Unknown init mode %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected init mode after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 't' :    /* next arg is MAP tie rule */  /*V1.06-e*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->TieRule=GetEnum(opts[ iopt ],TieStrVC,TIE_NB);
              if ( NemParaP->TieRule == -1 )
                {
                  fprintf( stderr, " Unknown tie rule %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;



        case 'B' :    /* next arg is beta estimation method */ /*V1.04-b*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
              StatModelP->Spec.BetaModel = 
		GetEnum( opts[ iopt ], BtaStrVC, BETA_NB );
	      if ( StatModelP->Spec.BetaModel == -1 )
		{
                  fprintf( stderr, " Unknown type of beta estimation %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
		}
	      else {
		if ( ( StatModelP->Spec.BetaModel == BETA_HEUL ) |
		     ( StatModelP->Spec.BetaModel == BETA_HEUD ) )
		  NemParaP->Crit = CRIT_U ;  /* force to choose U criterion */
	      }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;



        case 'C' :    /* next arg is crit to choose local max */ /*V1.04-i*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
              NemParaP->Crit = GetEnum( opts[ iopt ], CritStrVC, CRIT_NB );
	      if ( NemParaP->Crit == -1 )
		{
                  fprintf( stderr, " Unknown criterion name %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
		}
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;



	case 'G' :   /* next 4 args are parameters of beta gradient */
	  if ( iopt + 4 < nbopt )
	    {
	      NemParaP->BtaPsGrad.NbIter    = atof( opts[ iopt + 1 ] ) ;
	      NemParaP->BtaPsGrad.ConvThres = atof( opts[ iopt + 2 ] ) ;
	      NemParaP->BtaPsGrad.Step      = atof( opts[ iopt + 3 ] ) ;
	      NemParaP->BtaPsGrad.RandInit  = atoi( opts[ iopt + 4 ] ) ;
	      iopt += 4 ;
	    }
	  else
	    {
	      fprintf( stderr, " Expected 4 values after switch %s\n", 
		       opts[ iopt - 1 ] ) ;
	      err = STS_E_ARG ;
	    }
	  break ;




	case 'H' :   /* next 5 args are parameters of beta heuristic */
	  if ( iopt + 5 < nbopt )
	    {
	      NemParaP->BtaHeuStep = atof( opts[ iopt + 1 ] ) ;
	      NemParaP->BtaHeuMax  = atof( opts[ iopt + 2 ] ) ;
	      NemParaP->BtaHeuDDrop = atof( opts[ iopt + 3 ] ) ;
	      NemParaP->BtaHeuDLoss = atof( opts[ iopt + 4 ] ) ;
	      NemParaP->BtaHeuLLoss = atof( opts[ iopt + 5 ] ) ;
	      iopt += 5 ;
	    }
	  else
	    {
	      fprintf( stderr, " Expected 4 values after switch %s\n", 
		       opts[ iopt - 1 ] ) ;
	      err = STS_E_ARG ;
	    }
	  break ;


        case 'I' :    /* next arg is number of E-step iterations */ /*V1.04-j*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
              NemParaP->NbEIters = atoi( opts[ iopt ] ) ;
              if ( NemParaP->NbEIters < 0 )
                {
                  fprintf( stderr, "Nb iterations must be >= 0 (here %s)\n",
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'M' :    /* next arg is missing data mode */     /*V1.05-a*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->MissMode=GetEnum(opts[ iopt ], MissStrVC, MISSING_NB);
              if ( NemParaP->MissMode == -1 )
                {
                  fprintf( stderr, " Unknown missing mode %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;



        case 'O' :    /* next arg is visit order type */     /*V1.04-f*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->VisitOrder=GetEnum(opts[ iopt ], OrderStrVC, ORDER_NB);
              if ( NemParaP->VisitOrder == -1 )
                {
                  fprintf( stderr, " Unknown visit order type %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;



        case 'R' :    /* next arg is reference class file */  /*V1.04-h*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
              strncpy( NemParaP->RefName, opts[ iopt ], LEN_FILENAME ) ;
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;



        case 'S' :    /* next arg is random generator seed */  /*V1.04-e*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
              NemParaP->Seed = atol( opts[ iopt ] ) ;
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'T' :    /* set debug mode to be true */  /*V1.04-g*/
	  NemParaP->Debug = TRUE ;
          break ;


        case 'U' :    /* next arg is site update mode */  /*V1.06-d*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->SiteUpdate=GetEnum(opts[ iopt ],UpdateStrVC,UPDATE_NB);
              if ( NemParaP->SiteUpdate == -1 )
                {
                  fprintf( stderr, " Unknown site update mode %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'v' :    /* print out version information */
          PrintVersions( stdout ) ;
          err = STS_I_DONE ;
          break ;


        case '?' :    /* long help requested */
        case 'h' :    /*V1.04-k*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      err = PrintHelp( opts[ iopt ] ) ;
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;

        default :
          fprintf( stderr, "Unknown option switch '%s'\n", opts[ iopt ] ) ;
          err = STS_E_ARG ;
        } /* end if ( copt != '\0' ) switch ( copt ) */
      else
        {
          fprintf( stderr, " Expected an option switch at arg %d (%s)\n",
                   iopt + 3 , opts[ iopt ] ) ;
          err = STS_E_ARG ;
        }
    } /* end for ( iopt ... ) */

  strncpy( NemParaP->OutName, NemParaP->OutBaseName, LEN_FILENAME ) ;
  strncat( NemParaP->OutName, 
           NemParaP->Format == FORMAT_HARD ? EXT_OUTNAMEHARD : EXT_OUTNAMEFUZZY,
           LEN_FILENAME ) ;

  if ( NemParaP->DoLog )
  {
    if ( strcmp( NemParaP->OutBaseName , "-" ) != 0 ) /*V1.04-c*/
      {
	strncpy( NemParaP->LogName, NemParaP->OutBaseName, LEN_FILENAME ) ;
	strncat( NemParaP->LogName, EXT_LOGNAME, LEN_FILENAME ) ;
      }
    else
      {
	strncpy( NemParaP->LogName, Fname , LEN_FILENAME ) ;
	strncat( NemParaP->LogName, EXT_LOGNAME, LEN_FILENAME ) ;	
      }
  }
  else
  {
      strcpy( NemParaP->LogName, "" ) ;
  }

  return err ;

} /* end of NemArgs() */



/* ==================== LOCAL FUNCTION DEFINITION =================== */



/* ------------------------------------------------------------------- */
int  ReadStrFile
            (
                const char*  BaseName,          /* I */
                char*        CommentS,          /* O */
                DataT*       DataP,             /* O */
                ImageNeighT* ImageP,            /* O */
                TypeET*      TypeP              /* O */
            )
/*\
  Read structure file
       expected format : 
         type size nbvars
           type : N | I | S
           size : (if type == I) lx ly (if type == N or S) nbpts
           nbvars : integer
\*/
/* ------------------------------------------------------------------- */
{
    char    infname[ LEN_FILENAME + 1 ] ;
    FILE    *fstr ;
    int     nbelts ;
    char    type[100] ;  /* N | I | S */
    int     n1, n2, n3 ;

    void my_strupr( char *s ) ;

    strncpy( infname, BaseName, LEN_FILENAME ) ;
    strncat( infname, ".str" , LEN_FILENAME ) ;

    if ( ReadOpeningComments( infname , "#" , LEN_LINE , 
                              & fstr , CommentS ) == -1 )
    {
        fprintf( stderr, "File %s does not exist\n", infname ) ;
        return STS_E_FILEIN ;
    }


    /* Now we expect a line having structure : "I nl nc d" or "S n d" or
       "N n d"
    */
    nbelts = fscanf( fstr, "%s %d %d %d", type, &n1, &n2, &n3 ) ;
    fclose( fstr ) ;  /*V1.05-h*/

    if ( nbelts < 3 )
    {
        fprintf( stderr, "Structure file (%s) not enough fields\n",
                 infname ) ;
        return STS_E_FILE ;        
    }

    my_strupr( type ) ;

    if ( ! strcmp( type, TypeStrC[ TYPE_NONSPATIAL ] ) )
    {
        /* No spatial information */
        DataP->NbPts  = n1 ;
        DataP->NbVars = n2 ;
        *TypeP        = TYPE_NONSPATIAL ;
    }
    else if ( ! strcmp( type, TypeStrC[ TYPE_SPATIAL ] ) )
    {
        /* Spatial irregular */
        DataP->NbPts  = n1 ;
        DataP->NbVars = n2 ;
        *TypeP        = TYPE_SPATIAL ;
    }
    else if ( ! strcmp( type, TypeStrC[ TYPE_IMAGE ] ) )
    {
        /* Image */
        if ( nbelts < 4 )
        {
            fprintf( stderr, "Structure file (%s) not enough fields\n",
                     infname ) ;
            return STS_E_FILE ;        
        }             
        ImageP->Nl    = n1 ;
        ImageP->Nc    = n2 ;
        DataP->NbPts  = ImageP->Nl * ImageP->Nc ;
        DataP->NbVars = n3 ;
        *TypeP = TYPE_IMAGE ;
    }
    else
    {
        /* default : Error in structure file */
        fprintf( stderr, "Data type %s unknown in file %s\n", 
                 type, infname ) ;
        return STS_E_FILE ;
    }

    return STS_OK ;

}  /* end of ReadStrFile() */




/* ------------------------------------------------------------------- */
char GetOptionSwitch( const char *Arg ) 
/* ------------------------------------------------------------------- */
{
    if ( ( Arg[ 0 ] == '-' ) || ( Arg[ 0 ] == '/' ) )
      {
        return Arg[ 1 ] ;
      }
    else
      {
        return '\0' ;
      }
} /* end of GetOptionSwitch() */


/* ------------------------------------------------------------------- */
int GetEnum( const char* S, const char* SV[], int SizeV ) /*V1.04-b*/
/*\

    Returns the position of string S in the array SV consisting of
    SizeV strings (indiced from 0 to SizeV - 1).  Returns -1 if S is
    not found in SV.

\*/
/* ------------------------------------------------------------------- */
{
  int pos ;
  int found ;

  for ( pos = 0 , found = FALSE ;
	( pos < SizeV ) && ( ! found ) ;
	pos ++ )
    {
      if ( ! strcmp( S , SV[ pos ] ) )
	found = TRUE ;
    }

  if ( found )
    return ( pos - 1 ) ;  /* '- 1' because for() adds 1 after finding */
  else
    return -1 ;
}

/* ------------------------------------------------------------------- */
void PrintUsage( const char* CmdName )
/* ------------------------------------------------------------------- */
{
    fprintf( stderr, "Usage :   " ) ;
    fprintf( stderr, "  %s   file  nk  [ option1 option2 ... ]\n\n",
            CmdName ) ;
    fprintf( stderr, "  file       base name of input and output files\n" ) ;
    fprintf( stderr, "  nk         number of classes\n\n" ) ;

    fprintf( stderr, "Options :   [ default  ]\n\n" ) ;

    fprintf( stderr, "  -a algo   [ %-8s ]   type of algorithm { ",
            AlgoStrVC[ DEFAULT_ALGO ] ) ;
    {
      AlgoET ialgo ;
      
      for ( ialgo = 0 ; ialgo < ALGO_NB ; ialgo ++ )
        fprintf( stderr, "%s ", AlgoStrVC[ ialgo ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -b beta   [ %-8g ]   spatial coefficient, >= 0.0\n",
            DEFAULT_BETA ) ;

    fprintf( stderr, "  -c wh thr [%4s %5.3g]   convergence test { ",
            CvTestStrVC[ DEFAULT_CVTEST ], DEFAULT_CVTHRES ) ;
    {
      CvemET icvtest ;
      
      for ( icvtest = 0 ; icvtest < CVTEST_NB ; icvtest ++ )
        fprintf( stderr, "%s ", CvTestStrVC[ icvtest ] ) ;
    }
    fprintf( stderr, "} + threshold (>0)\n" ) ; /*V1.04-d*//*V1.06-g*/

    fprintf( stderr, "  -f format [ %-8s ]   partition format, { ",
            FormatStrVC[ DEFAULT_FORMAT ] ) ;
    {
      FormET ipart ;
      
      for ( ipart = 0 ; ipart < FORMAT_NB ; ipart ++ )
        fprintf( stderr, "%s ", FormatStrVC[ ipart ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -i itmax  [ %-8d ]   number of NEM iterations (>= 0)\n",
             DEFAULT_NBITERS ) ;

    fprintf( stderr, "  -l dolog  [ n        ]   log file or not { y n }\n" ) ;  /*V1.03-a*/

    fprintf( stderr, "  -m f p d  [%-4s %-2s %-3s]  model  {",
	     FamilyStrVC[ DEFAULT_FAMILY ],
	     ProporStrVC[ DEFAULT_PROPOR ],
	     DisperStrVC[ DEFAULT_DISPER ] ) ;
    {
      int imodel ;
      
      for ( imodel = 0 ; imodel < FAMILY_NB - 1 ; imodel ++ )
        fprintf( stderr, "%s ", FamilyStrVC[ imodel ] ) ;
      fprintf( stderr, "%s", FamilyStrVC[ FAMILY_NB - 1 ] ) ;
      fprintf( stderr, "} {" ) ;
      for ( imodel = 0 ; imodel < PROPOR_NB - 1 ; imodel ++ )
        fprintf( stderr, "%s ", ProporStrVC[ imodel ] ) ;
      fprintf( stderr, "%s", ProporStrVC[ PROPOR_NB - 1 ] ) ;
      fprintf( stderr, "} {" ) ;
      for ( imodel = 0 ; imodel < DISPER_NB - 1 ; imodel ++ )
        fprintf( stderr, "%s ", DisperStrVC[ imodel ] ) ;
      fprintf( stderr, "%s", DisperStrVC[ DISPER_NB - 1 ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -n neigh  [ %-8s ]   neighbour specification (image), { ",
            NeighStrVC[ DEFAULT_NEIGHSPEC ] ) ;
    {
      NeighET inei ;
      
      for ( inei = 0 ; inei < NEIGH_NB ; inei ++ )
        fprintf( stderr, "%s ", NeighStrVC[ inei ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -o fout   [ file     ]   output files basename\n" ) ;

    fprintf( stderr, "  -s init   [ s 1      ]   init : " ) ;
    fprintf( stderr, "s <v> = sort var <v>\n" ) ;
    fprintf( stderr, "%33s f <ini.uf> = from <ini.uf>\n", " " ) ; /*V1.04-a*/
    fprintf( stderr, "%33s r <n> = <n> random inits\n", " " ) ;
    fprintf( stderr, "%33s l <file> = use known labels from <file>\n", " " ) ;
    fprintf( stderr, "%33s m <file> = initial/fixed parameters\n", " " ) ; /*V1.08-a*/

    fprintf( stderr, "--- Press return for more options ---" ) ;
    getchar( ) ;

    fprintf( stderr, "  -t tie    [ %-8s ]   MAP tie rule  { ",
            TieStrVC[ DEFAULT_TIE ] ) ;
    {
      TieET itie ;
      
      for ( itie = 0 ; itie < TIE_NB ; itie ++ )
        fprintf( stderr, "%s ", TieStrVC[ itie ] ) ;
    }
    fprintf( stderr, "}\n" ) ;


    fprintf( stderr, "  -B bmod   [ %-8s ]   b estimation mode { ",
	     BtaStrVC[ DEFAULT_BTAMODE ] ) ;
    {
      BetaET bmod ;
      for (bmod = 0; bmod < BETA_NB; bmod ++)
	fprintf( stderr, "%s ", BtaStrVC[ bmod ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -C crit   [ %-8s ]   local maximum criterion { ",
	     CritStrVC[ DEFAULT_CRIT ] ) ;
    {
      CritET  crit ;
      for (crit = 0; crit < CRIT_NB; crit ++)
	fprintf( stderr, "%s ", CritStrVC[ crit ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -G nit conv step rand [ %2d %5.3f %3.1f %1d ]  %s\n",
	     DEFAULT_BTAGRADNIT, DEFAULT_BTAGRADCVTH, DEFAULT_BTAGRADSTEP,
	     DEFAULT_BTAGRADRAND,
	     "parameters of beta gradient estimation") ;

    fprintf( stderr, "  -H bstep bmax ddrop dloss lloss [ %3.2f %3.1f %3.1f %3.1f %3.2f ]\n %33s %s\n",
	     DEFAULT_BTAHEUSTEP, DEFAULT_BTAHEUMAX, DEFAULT_BTAHEUDDROP,
	     DEFAULT_BTAHEUDLOSS, DEFAULT_BTAHEULLOSS, " ",
	     "parameters of beta heuristic") ;

    fprintf( stderr, "  -I eiter  [ %-8d ]   number of E-step iterations (>= 1)\n",
             DEFAULT_NBEITERS ) ;/*V1.04-j*/


    fprintf( stderr, "  -M miss   [ %-8s ]   missing data processing { ",
            MissStrVC[ DEFAULT_MISSING ] ) ;
    {
      MissET imissing ;
      
      for ( imissing = 0 ; imissing < MISSING_NB ; imissing ++ )
        fprintf( stderr, "%s ", MissStrVC[ imissing ] ) ;
    }
    fprintf( stderr, "}\n" ) ;


    fprintf( stderr, "  -O order  [ %-8s ]   order of site visit { ",
            OrderStrVC[ DEFAULT_ORDER ] ) ;
    {
      OrderET iorder ;
      
      for ( iorder = 0 ; iorder < ORDER_NB ; iorder ++ )
        fprintf( stderr, "%s ", OrderStrVC[ iorder ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -S seed   [ <time>   ]   random generator seed \n" ) ;     /*V1.04-e*/

    fprintf( stderr, "  -T                       print debugging information\n" ) ;     /*V1.04-f*/

    fprintf( stderr, "  -U update [ %-8s ]   site update scheme  { ",
            UpdateStrVC[ DEFAULT_UPDATE ] ) ;
    {
      UpdET iupdate ;
      
      for ( iupdate = 0 ; iupdate < UPDATE_NB ; iupdate ++ )
        fprintf( stderr, "%s ", UpdateStrVC[ iupdate ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "\n\nYou may also just type arguments : \n" ) ;
    fprintf( stderr, "  -v                      versions information\n" ) ;
    fprintf( stderr, "  -h help_topic           longer help - help topics are\n {" ) ;
    {
      HelpET ihelp ;

      for ( ihelp = 0 ; ihelp < HELP_NB ; ihelp ++ )
	fprintf( stderr, " %s", HelpStrVC[ ihelp ] ) ;
    }
    fprintf( stderr, " } \n" ) ;

} /* end of PrintUsage() */


/* ------------------------------------------------------------------- */

StatusET PrintHelp( const char* HelpType ) 

/* ------------------------------------------------------------------- */
{
  StatusET err = STS_I_DONE ;

  switch( GetEnum( HelpType, HelpStrVC, HELP_NB) )
    {
    case HELP_GENERAL:
      PrintHelpGeneral( stdout ) ;
      break ;

    case HELP_OPTIONS:
      PrintHelpOptions( stdout ) ;
      break ;

    case HELP_EXAMPLES:
      PrintHelpExamples( stdout ) ;
      break ;

    case HELP_FILEIN:
      PrintHelpFileIn( stdout ) ;
      break ;

    case HELP_FILEOUT:
      PrintHelpFileOut( stdout ) ;
      break ;

    case HELP_VERSIONS:
      PrintHelpVersions( stdout ) ;
      break ;

    default:
      fprintf( stderr, " Unknown help type %s\n", 
	       HelpType ) ;
      err = STS_E_ARG ;
    }  

  return err ;
}


/*V1.06-a*/
/* ------------------------------------------------------------------- */
void my_strupr( char *s )
/* ------------------------------------------------------------------- */
{
    if ( s == NULL ) return ;
    for ( ; (*s) != '\0' ; s ++ )
    {
        if ( ( (*s) > 'a' ) && ( (*s) < 'z' ) )
        {
            (*s) = (*s) + ( 'A' - 'a' ) ;
        }
    }
}


/* ======================================================================= */
/*V1.06-c*/
/* start test of NemArgs 
int main( int argc, const char* argv[] )
{
  char        fname[ LEN_FILENAME + 1 ] ;
  char        datadescS[ LEN_LINE + 1 ] ;
  NemParaT    NemPara ;
  SpatialT    Spatial ;
  StatModelT  StatModel ;
  DataT       Data ;

  StatusET  err = 
    NemArgs( argc, argv, fname, & StatModel, & NemPara,
	     datadescS, & Data, & Spatial.NeighData.Image,
	     & Spatial.Type ) ;

  return err ;
}
 end test of NemArgs */
