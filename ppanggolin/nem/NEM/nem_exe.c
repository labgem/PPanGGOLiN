/*\

    NEM_EXE.C

    Neighborhood EM project : 
    main program and input/output functions.

    June 96

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

    This program clusters a given data set, using the 
    Neighborhood EM algorithm proposed by Ambroise 1996.
    All parameters are to be specified on the command line.
    Command syntax can be obtained by typing : 'nem_exe'
    with no argument. Longer help is obtained by typing :
    'nem_exe -h'.

    The program takes as input  a set of objects described, 
    on one hand, by quantitative variables (the 'observed data'), 
    and on the other hand (optionnally), by their neighborhood 
    relationships in the geographic space. It outputs  
    a partition of the objects ( fuzzy or hard ) and class 
    parameters. The outputs are 'locally optimal' according to 
    a 'penalized likelihood' criterion (see below). Only 
    two parameters are required: the number of classes, 
    and the 'beta' coefficient that measures the desired 
    spatial smoothing. The processus is otherwise entirely 
    unsupervised. Other optional parameters may be specified 
    to fit better to the problem at hand.

    The 'penalized likelihood' criterion is designed in order to
    produce a partition and class parameters that fit well to 
    the data (likelihood term), and a partition that is spatially
    homogeneous (penalization term). It is possible to give
    more or less importance to the spatial information by
    increasing or decreasing the 'beta' coefficient.

    The algorithm finds a 'local' maximum of the criterion
    through an iterative processus, akin to gradient ascent : 
    starting from an initial arbitrary solution, it computes 
    successive solutions that gradually improve the criterion ; 
    thence, the final partition depends on the initial partition.
    This program can either read a given initial partition,
    compute it by thresholding one variable's histogram,
    or generate random initial partitions.
\*/

/*\
    NEM_EXE.C

Vers-mod  Date         Who  Description

1.03-a    03-NOV-1996  MD   Add prototype of ReadLabelFile() in main()
1.03-b    03-NOV-1996  MD   Save Pk, Vk, and Ck in don.mf
1.03-c    22-AUG-1997  MD   Unknown labels : memberships set to 0
1.03-d    22-AUG-1997  MD   Save criteria to don.mf file
1.03-e    30-SEP-1997  MD   Save criterion L to don.mf file
1.03-f    30-SEP-1997  MD   Save algorithm parameters into don.mf file
1.04-a    04-OCT-1997  MD   Any name for initial partition file
1.04-b    09-OCT-1997  MD   If "-o -" classification saved to standard output
1.04-c    09-OCT-1997  MD   If "-s f -" initial classification from stdin
1.04-d    10-OCT-1997  MD   Random seed set by option
1.04-e    10-OCT-1997  MD   Visit order ("-O random")
1.04-f    13-OCT-1997  MD   Read reference classification
1.04-g    05-NOV-1997  MD   Call to srand48() replaced by srandom()
1.04-h    02-DEC-1997  MD   Add "GEM" to AlgoDesC
1.05-a    12-JAN-1998  MD   Count missing data
1.05-b    17-JAN-1998  MD   Program exit values changed
1.05-c    25-JAN-1998  MD   main() now mainfunc()
1.05-d    26-JAN-1998  MD   DataP->LabelV init to NULL by default
1.05-e    26-JAN-1998  MD   Free allocated data after processing
1.05-f    26-JAN-1998  MD   GenAlloc/GenFree instead of malloc/calloc/free
1.05-g    30-JAN-1998  MD   mainfunc() prototyped in mainfunc.h
1.05-h    30-JAN-1998  MD   Fixed bug : str file was opened and not closed
1.05-i    05-FEB-1998  MD   my_strupr returns void
1.05-j    05-FEB-1998  MD   save M (Markov pseudo-likelihood) into don.mf
1.06-a    17-JUN-1998  MD   use new StatModel structure (old NoisModel)
1.06-b    23-JUN-1998  MD   LEN_LINE -> nem_typ.h
1.06-c    20-SEP-1998  MD   change format of ref / label file (!!!)
1.06-d    20-SEP-1998  MD   struct ErrInfo added to CriterT
1.06-e    21-SEP-1998  MD   TieRule copied in Errinfo
1.07-a    26-FEB-1999  MD   Add Bernoulli family
1.07-b    26-FEB-1999  MD   Add "\n" at end of final classification file
1.08-a    20-JUI-2017  GG   Add param input by file rather than by arguments
\*/

#define _GNU_SOURCE // FIX: Required for musl build because srandom is not POSIX
#include <stdlib.h> 
#include "nem_exe.h"   /* Prototype of exported mainfunc() */

/* ==================== LOCAL FUNCTION PROTOTYPING =================== */


/* Called by ClassifyByNem */

    static int SaveResults
        ( const int          Npt,                   /* I */
          const int          Nd,                    /* I */
          const float*       ClassifM,              /* I */
          const SpatialT*    SpatialP,              /* I */
          const NemParaT*    NemParaP,              /* I */
          const StatModelT*  StatModelP,            /* I */
	  const CriterT*     CriterP                /* I */ /*V1.03-d*/
        ) ;

static void FreeAllocatedData
(
 DataT*       DataP,       /* O and deallocated */
 SpatialT*    SpatialP,    /* O and deallocated */
 ModelParaT*  ModelParaP,  /* O and deallocated */  /*V1.06-a*/
 CriterT*     CriterP,     /* O and deallocated */  /*V1.06-d*/
 float*       ClassifM     /* O and deallocated */
) ;


    static int SetVisitOrder  /*V1.04-e*/
        ( 
	     int         Npt,          /* I */
	     OrderET     VisitOrder,   /* I */
	     int**       SiteVisitVP   /* O and allocated (Npt) */
	) ;
         
    static int  ReadStrFile
            (
                const char*  BaseName,          /* I */
                char*        CommentS,          /* O */
                DataT*       DataP,             /* O */
                ImageNeighT* ImageP,            /* O */
                TypeET*      TypeP              /* O */
            ) ;

    static int  ReadMatrixFile                     
         (                             
             const char  *FileName,     /* I */      /*V1.04-a*/
             int         Nl,            /* I */
             int         Nc,            /* I */
             float       **MatP         /* O and allocated (size : Nl * Nc) */
         ) ;                             
                                     
    static int  ReadLabelFile  /*V1.03-a*/
         ( 
             const char  *LabelName,    /* I */
             int         Npt,           /* I */
             int         *KfileP,       /* O : file # classes */ /*V1.06-c*/
             int         **LabelVP,     /* O and allocated (Npt) */
             float       **ClassifMP    /* O and allocated (Npt*Nk) */
         );
    static int  ReadParamFile  /*V1.08-a*/
         (
             const char  *ParamName,    /* I */
  	     const FamilyET Family,
             const int   K,          /* I : number of classes */
             const int   D,          /* I : number of variables */
             ParamFileET* ParamFileMode, /* O : specified if parameters will be used at begin or throughout the clustering process */
             ModelParaT* ParaP       /* O : read parameters */
         );
    static int  MakeErrinfo
         ( 
             const char* RefName,       /* I : filename of reference class */
             int         N,             /* I : number of objects */
             int         Kc,            /* I : user number of classes */
	     TieET       TieRule,       /* I : specified MAP tie rule */
	     ErrinfoT*   ErrinfoP,      /* O and allocated */
	     ErrcurT*    ErrcurP        /* O and allocated */
         ) ;

    static int  ReadNeiFile                     
         (                              
             const char  *BaseName,     /* I */
             int         NbPts,         /* I */
             NeighET     NeighSpec,     /* I */
             char*       NeiDescS,      /* O [LEN_LINE+1] */
             SpatialT    *SpatialP      /* O and allocated */
         ) ;



/* Called by ReadNeiFile */

    static int  ReadPtsNeighs
                (
                    FILE        *Fnei,          /* I/O */
                    int         NbPts,          /* I */
                    int         *MaxNeiP,       /* O */
                    NeighDataT  *NeighDataP     /* O and allocated */
                ) ;

    static int  ReadImageNeigh
                (
                    FILE        *Fnei,          /* I/O */
                    NeighDataT  *NeighDataP     /* O and allocated */
                ) ;

    static int  SetImageNeigh
                (   
                    NeighET     NeighSpec,          /* I */
                    char*       NeiDescS,           /* O [LEN_LINE+1] */
                    NeighDataT* NeighDataP          /* O and allocated */
                ) ;


/* Called by MakeErrinfo */

static int factorial(int n);

static int compute_permutations    /* ret 0 if OK, -1 if memory error */
(
 const int Startval,         /* I : start value of integer suite */
 const int K,                /* I : size of integer suite to permute > 0 */
 int*      perms_Kfact_K_p[] /* O : matrix to store permuted values */
) ;


/* Called by compute_permutations */

static int rec_permutations        /* ret 0 if OK, -1 if memory error */
(
 const int array_A[],       /* I : remaining array to permute */
 const int A,               /* I : length of the remaining array : 0..K */
 const int K,               /* I : length of original array */
 int       offset,          /* I : first line of storage to use */
 int       perms_Kfact_K[]  /* O : matrix to store permuted values, use
			       lines :   offset -> offset + A! - 1
			       columns : K - A  -> K - 1 */
) ;

int GetEnum( const char* S, const char* SV[], int SizeV ) ;
void my_strupr( char *s ) ;


//VERSION
const char *NemVersionStrC = "1.08-a";

/* ==================== GLOBAL FUNCTION DEFINITION =================== */


/* ------------------------------------------------------------------- */
int nem(const char* Fname,
        const int   nk,
        const char* algo,
        const float beta,
        const char* convergence,
        const float convergence_th,
        const char* format,
        const int   it_max,
        const int   dolog,
        const char* model_family,
        const char* proportion,
        const char* dispersion,
        const int   init_mode,
        const char* init_file,
        const char* out_file_prefix,
        const int   seed)
/*\
    NEM function.
\*/
/* ------------------------------------------------------------------- */
{
    const char*             func = "nem" ;
    StatusET                err ;
    static  DataT           Data = {0} ;
    static  NemParaT        NemPara = {0} ;
    static  SpatialT        Spatial = {{{0}}} ;
    static  StatModelT      StatModel = {{0}} ;
    static  float           *ClassifM = 0;
    CriterT                 Criteria = {0} ; /*V1.03-d*/

    /* main program algorithm :
       - read all necessary data into memory
       - compute partition by NEM method
       - save results to file
    */

    if(!dolog)
    {
        out_stderr = stderr;
    }
    else
    {
        char name_out_stderr[LEN_FILENAME];
        strncpy( name_out_stderr , out_file_prefix , LEN_FILENAME ) ;  /*V1.04-a*/
        strncat( name_out_stderr , ".stderr", LEN_FILENAME ) ;
        out_stderr = fopen(name_out_stderr, "w");
    }

    fprintf( out_stderr , " * * * NEM (spatial data clustering) v%s * * *\n" ,
             NemVersionStrC ) ;
#ifdef __TURBOC__
    fprintf( out_stderr, "\n Initial free memory : %lu bytes\n\n", 
             (unsigned long) coreleft() );
#endif

    /* RandomSeedByTime( ) ;  V1.04-d*/

    char        datadescS[ LEN_LINE + 1 ] ;

      StatModel.Spec.K = nk ;
    if ( nk <= 0 )
    {
        fprintf( out_stderr, "Nb of classes must be > 0 (here %d)\n",
               nk ) ;
        return STS_E_ARG ;
    }

    if ( ( err = ReadStrFile( Fname, 
                datadescS,
                &Data, 
                &Spatial.NeighData.Image,
                &Spatial.Type) ) != STS_OK )
    return err ;

      /* !!! Allocate model parameters */ /*V1.06-a*/
    StatModel.Para.Prop_K    = GenAlloc( nk, sizeof(float), 
                       1, func, "Prop_K" ) ;
    StatModel.Para.Disp_KD   = GenAlloc( nk * Data.NbVars, sizeof(float), 
                         1, func, "Disp_KD" ) ;
    StatModel.Para.Center_KD = GenAlloc( nk * Data.NbVars, sizeof(float), 
                       1, func, "Center_KD" ) ;
    StatModel.Para.NbObs_K   = GenAlloc( nk, sizeof(float), 
                       1, func, "NbObs_K" ) ;
    StatModel.Para.NbObs_KD  = GenAlloc( nk * Data.NbVars, sizeof(float), 
                       1, func, "NbObs_KD" ) ;
    StatModel.Para.Iner_KD   = GenAlloc( nk * Data.NbVars, sizeof(float), 
                       1, func, "NbObs_KD" ) ;
    StatModel.Desc.DispSam_D = GenAlloc( Data.NbVars, sizeof(float), 
                        1, func, "DispSam_D" );
    StatModel.Desc.MiniSam_D = GenAlloc( Data.NbVars, sizeof(float), 
                        1, func, "MiniSam_D" );
    StatModel.Desc.MaxiSam_D = GenAlloc( Data.NbVars, sizeof(float), 
                        1, func, "MaxiSam_D" );
    /* Set default value of optional parameters */
    StatModel.Spec.ClassFamily = DEFAULT_FAMILY ;
    StatModel.Spec.ClassDisper = DEFAULT_DISPER ;
    StatModel.Spec.ClassPropor = DEFAULT_PROPOR ;
    NemPara.Algo          = DEFAULT_ALGO ;
    StatModel.Para.Beta   = DEFAULT_BETA ;          /*V1.06-b*/
    StatModel.Spec.BetaModel = DEFAULT_BTAMODE ;       /*V1.04-b*/
    NemPara.BtaHeuStep    = DEFAULT_BTAHEUSTEP ;    /*V1.04-b*/
    NemPara.BtaHeuMax     = DEFAULT_BTAHEUMAX ;
    NemPara.BtaHeuDDrop   = DEFAULT_BTAHEUDDROP ;
    NemPara.BtaHeuDLoss   = DEFAULT_BTAHEUDLOSS ;
    NemPara.BtaHeuLLoss   = DEFAULT_BTAHEULLOSS ;
    NemPara.BtaPsGrad.NbIter    = DEFAULT_BTAGRADNIT  ;/*V1.06-g*/
    NemPara.BtaPsGrad.ConvThres = DEFAULT_BTAGRADCVTH ;
    NemPara.BtaPsGrad.Step      = DEFAULT_BTAGRADSTEP ;
    NemPara.BtaPsGrad.RandInit  = DEFAULT_BTAGRADRAND ;
    NemPara.Crit          = DEFAULT_CRIT ;          /*V1.04-h*/
    NemPara.CvThres       = DEFAULT_CVTHRES ;       /*V1.04-d*/
    NemPara.CvTest        = CVTEST_CLAS ;           /*V1.06-g*/
    NemPara.DoLog         = FALSE ;                 /*V1.03-a previously TRUE*/
    NemPara.NbIters       = DEFAULT_NBITERS ;
    NemPara.NbEIters      = DEFAULT_NBEITERS ;
    NemPara.NbRandomInits = DEFAULT_NBRANDINITS ;  /*V1.06-h*/
    NemPara.Seed          = seed ;//time( NULL )          /*V1.04-e*/
    NemPara.Format        = DEFAULT_FORMAT ;
    NemPara.InitMode      = DEFAULT_INIT ;
    NemPara.ParamFileMode = DEFAULT_NO_PARAM_FILE ;
    NemPara.SortedVar     = DEFAULT_SORTEDVAR ;
    NemPara.NeighSpec     = DEFAULT_NEIGHSPEC ;
    NemPara.VisitOrder    = DEFAULT_ORDER ;         /*V1.04-f*/
    NemPara.SiteUpdate    = DEFAULT_UPDATE ;        /*V1.06-d*/
    NemPara.TieRule       = DEFAULT_TIE ;           /*V1.06-e*/
    NemPara.Debug         = FALSE ;                 /*V1.04-g*/
    strncpy( NemPara.OutBaseName, out_file_prefix, LEN_FILENAME ) ;
    strncpy( NemPara.NeighName, Fname, LEN_FILENAME ) ;
    strncpy( NemPara.ParamName, init_file, LEN_FILENAME ) ;
    strncat( NemPara.NeighName, ".nei", LEN_FILENAME ) ;
    strncpy( NemPara.RefName, "", LEN_FILENAME ) ;

    //-----
    NemPara.Algo = GetEnum( algo , AlgoStrVC, ALGO_NB ) ;
    if ( NemPara.Algo == -1 )
    {
        fprintf( out_stderr, " Unknown type of algorithm %s\n", algo ) ;
        err = STS_E_ARG ;
    }
    //-----
    
    if (beta < 0)
    {
        StatModel.Spec.BetaModel = BETA_PSGRAD ;
    }
    else{
        StatModel.Para.Beta = beta ;
    }
    //-----
    NemPara.CvTest=GetEnum( convergence, CvTestStrVC, CVTEST_NB );
    if ( NemPara.CvTest == -1 ) {
      fprintf( out_stderr, " Unknown convergence test %s\n", convergence ) ;
      err = STS_E_ARG ;
    }
    else if ( NemPara.CvTest != CVTEST_NONE ) /* get threshold */ {
        NemPara.CvThres = convergence_th ;
        if ( NemPara.CvThres <= 0 ) {
            fprintf( out_stderr, " Conv threshold must be > 0 (here %f)\n", convergence_th ) ;
            err = STS_E_ARG ;
        } /* else threshold > 0 : OK */
    } 
    //-----
    NemPara.Format=GetEnum( format , FormatStrVC, FORMAT_NB );
    if ( NemPara.Format == -1 )
    {
        fprintf( out_stderr, " Unknown format %s\n", format) ;
        err = STS_E_ARG ;
    }
    //-----
    NemPara.NbIters = it_max ;
    if ( NemPara.NbIters < 0 )
    {
        fprintf( out_stderr, "Nb iterations must be >= 0 (here %d)\n",  it_max ) ;
        err = STS_E_ARG ;
    }
    //-----
    if ( dolog )
        NemPara.DoLog = TRUE ;
    else
        NemPara.DoLog = FALSE ;
    //-----
    StatModel.Spec.ClassFamily = GetEnum( model_family, FamilyStrVC, FAMILY_NB );
    if ( StatModel.Spec.ClassFamily == -1 )
    {
        fprintf( out_stderr, " Unknown family %s\n", model_family ) ;
        err = STS_E_ARG ;
    }
    //-----
    StatModel.Spec.ClassPropor = GetEnum( proportion, ProporStrVC, PROPOR_NB );
    if ( StatModel.Spec.ClassPropor == -1 )
    {
        fprintf( out_stderr, " Unknown proportion %s\n", proportion ) ;
        err = STS_E_ARG ;
    }
    //-----
    StatModel.Spec.ClassDisper = GetEnum( dispersion, DisperStrVC, DISPER_NB );
    if ( StatModel.Spec.ClassDisper == -1 )
    {
        fprintf( out_stderr, " Unknown dispersion %s\n", dispersion) ;
        err = STS_E_ARG ;
    }
    //-----
    NemPara.NeighSpec = NEIGH_FILE;
    //-----
    NemPara.InitMode = init_mode;

    strncpy( NemPara.OutName, NemPara.OutBaseName, LEN_FILENAME ) ;
    strncat( NemPara.OutName, 
           NemPara.Format == FORMAT_HARD ? EXT_OUTNAMEHARD : EXT_OUTNAMEFUZZY,
           LEN_FILENAME ) ;

    if ( NemPara.DoLog )
    {
      if ( strcmp( NemPara.OutBaseName , "-" ) != 0 ) /*V1.04-c*/
        {
      strncpy( NemPara.LogName, NemPara.OutBaseName, LEN_FILENAME ) ;
      strncat( NemPara.LogName, EXT_LOGNAME, LEN_FILENAME ) ;
        }
      else
        {
      strncpy( NemPara.LogName, NemPara.OutBaseName , LEN_FILENAME ) ;
      strncat( NemPara.LogName, EXT_LOGNAME, LEN_FILENAME ) ;   
        }
    }
    else
    {
        strcpy( NemPara.LogName, "" ) ;
    }

    strncpy( NemPara.ParamName, init_file, LEN_FILENAME ) ;

    /* Read points */

    char namedat[ LEN_FILENAME + 1 ] ;   /*V1.04-a*/
    char neidescS[ LEN_LINE + 1 ] ;
    int  klabelfile ;

    strncpy( namedat , Fname , LEN_FILENAME ) ;  /*V1.04-a*/
    strncat( namedat , ".dat", LEN_FILENAME ) ;
    fprintf( out_stderr, "Reading points ...\n" ) ;
    if ( ( err = ReadMatrixFile( namedat,        /*V1.04-a*/
                                 Data.NbPts, 
                                 Data.NbVars, 
                                 &Data.PointsM ) ) != STS_OK )
       return err ;
    /* Count missing data */ /*V1.05-a*/
    {
          int i, j ;

          Data.NbMiss = 0 ;
          for ( i = 0 ; i < Data.NbPts ; i ++ )
        {
          for ( j = 0 ; j < Data.NbVars ; j ++ )
            {
              if ( isnan( Data.PointsM[ i * Data.NbVars + j ] ) )
            Data.NbMiss ++ ;
            }
        }
    }

    /* Allocate and set sites visit order */
    if ( ( err = SetVisitOrder( Data.NbPts,   /*V1.04-e*/
                NemPara.VisitOrder,
                & Data.SiteVisitV ) ) != STS_OK )
      return err ;


    /* Allocate and eventually read initial classification matrix */
    Data.LabelV = NULL ; /*V1.05-d*/

    switch( NemPara.InitMode )
    {
    case INIT_FILE:
        fprintf( out_stderr, "Reading initial partition ...\n" ) ;
        if ( ( err = ReadMatrixFile( NemPara.StartName,     /*V1.04-a*/
				                     Data.NbPts,
                                     StatModel.Spec.K, 
                                     &ClassifM ) ) != STS_OK ) // FIX: ReadMatrixFile needs a float**
            return err ;
        break ;

    case INIT_PARAM_FILE:
        fprintf( out_stderr, "Reading parameter file ...\n" ) ;
        if ( ( err = ReadParamFile( NemPara.ParamName,
                                    StatModel.Spec.ClassFamily,
                                    StatModel.Spec.K,
                                    Data.NbVars,
                                    & NemPara.ParamFileMode,
                                    & StatModel.Para ) ) != STS_OK ) return err ;
    case INIT_SORT:    /* allocate partition (will be initialized later) */
    case INIT_RANDOM:  /* allocate partition (will be initialized later) */
        /* Allocate classification matrix */
        if ( ( (ClassifM) = 
	       GenAlloc( Data.NbPts * StatModel.Spec.K, sizeof( float ),
			 0, func, "(ClassifM)" ) ) == NULL )
	  return STS_E_MEMORY ;
        break ;

    case INIT_LABEL:
      fprintf( out_stderr, "Reading known labels file ...\n" ) ;
      if ( ( err = ReadLabelFile( NemPara.LabelName, Data.NbPts, 
				                  & klabelfile,
                                  & Data.LabelV,
                                  &ClassifM ) ) != STS_OK ) // FIX: ReadMatrixFile needs a float**
        return err ;

      if ( klabelfile != StatModel.Spec.K ) {
	fprintf( out_stderr, 
		 "Error : label file %d classes, command line %d classes\n",
		 klabelfile, StatModel.Spec.K ) ;
	return STS_E_FILE ;
      }
      break;

    default: /* error */
        fprintf( out_stderr, "Unknown initialization mode (%d)\n", 
                 NemPara.InitMode );
        return STS_E_FUNCARG ;
    }


    /* Eventually read reference class file */  /*V1.04-f*/
    if ( ( err = MakeErrinfo( NemPara.RefName, Data.NbPts, 
                  StatModel.Spec.K, NemPara.TieRule, 
                  &Criteria.Errinfo, &Criteria.Errcur ) ) != STS_OK )
      return err ;

    /* Read neighborhood file */
    if ( Spatial.Type != TYPE_NONSPATIAL )
    {
        fprintf( out_stderr, "Reading neighborhood information ...\n" ) ;
        if ( ( err = ReadNeiFile( Fname, 
                                  Data.NbPts, 
                                  NemPara.NeighSpec ,
                                  neidescS,
                                  &Spatial ) ) != STS_OK )
                return err ;
    }
    else
    {
        StatModel.Para.Beta = 0.0 ; /*V1.06-a*/
        Spatial.MaxNeighs   = 0 ;
    }

    fprintf( out_stderr, "\nData : " ) ;
    if ( strcmp( datadescS, "" ) ) fprintf( out_stderr, "%s\n", datadescS ) ;
    else                          fprintf( out_stderr, "\n" ) ;

    fprintf( out_stderr, "  file names =  %10s   |   nb points   = %10d\n", 
             Fname, Data.NbPts ) ;
    fprintf( out_stderr, "  type       =  %10s   |   dim         = %10d\n",
             TypeDesC[ Spatial.Type ], Data.NbVars ) ;
    if ( Spatial.Type == TYPE_IMAGE )
    {
    fprintf( out_stderr, "  image size =  (%4d,%4d)\n", 
             Spatial.NeighData.Image.Nl, Spatial.NeighData.Image.Nc ) ;
    }
    if ( Data.NbMiss > 0 ) /*V1.05-a*/
    {
    fprintf( out_stderr, "  %d missing values / %d\n",
         Data.NbMiss, Data.NbPts * Data.NbVars ) ;
    }

    if ( Spatial.Type != TYPE_NONSPATIAL )
    {
    fprintf( out_stderr, "Neighborhood system :\n  max neighb =  %10d\n",
             Spatial.MaxNeighs ) ;
    fprintf( out_stderr, "%s\n", neidescS ) ;
    }

    fprintf( out_stderr, "\n" ) ;
    fprintf( out_stderr, "NEM parameters :\n" ) ;
    fprintf( out_stderr, "Type of algorithm : '%s'\n" , AlgoDesC[ NemPara.Algo ] ) ;
    fprintf( out_stderr, 
         "  beta       =  %10.2f   |   nk                    = %3d\n", 
             StatModel.Para.Beta, StatModel.Spec.K) ;
    fprintf( out_stderr, 
         "                %10s   |   model                 = %s, %s %s\n",
             " ", 
         FamilyDesVC[ StatModel.Spec.ClassFamily ],
         ProporDesVC[ StatModel.Spec.ClassPropor ],
         DisperDesVC[ StatModel.Spec.ClassDisper ]) ;
    fprintf( out_stderr, "\n" ) ;

    if ( err == STS_OK )
    {
#ifdef __TURBOC__
      srand( (unsigned) NemPara.Seed ) ;
#else
      srandom( NemPara.Seed ) ;  /*V1.04-g*/
#endif

        if ( ( err = ClassifyByNem( &NemPara, &Spatial, &Data, 
                                    &StatModel, ClassifM, 
				    &Criteria ) ) == STS_OK )
        {
            fprintf( out_stderr, "Saving results ...\n" ) ;
            SaveResults( Data.NbPts , Data.NbVars , ClassifM, 
                         &Spatial, &NemPara, &StatModel, &Criteria ) ;
        }

	FreeAllocatedData( &Data, &Spatial, &StatModel.Para, 
			   &Criteria, ClassifM ) ;
    }

    switch( err )
    {
        case STS_OK :
	  if ( strcmp( NemPara.OutBaseName , "-" ) != 0 )  /*V1.04-b*/
	    {
	      fprintf( out_stderr, "NEM completed, classification in %s\n",
		       NemPara.OutName ) ;
	      fprintf( out_stderr, " criteria and parameters in %s%s\n",
		       NemPara.OutBaseName, EXT_MFNAME ) ;

	      if ( NemPara.DoLog )
		{
		  fprintf( out_stderr, "Log of detailed running in %s\n",
			   NemPara.LogName );
		}
	    }
	  else
	    {
	      fprintf( out_stderr, "NEM completed, classification to screen\n") ;
	    }
      fclose(out_stderr);
	  return EXIT_OK ;

        case STS_W_EMPTYCLASS :
             fprintf( out_stderr, "*** NEM warning status : empty class\n" );
             fclose(out_stderr);
             return EXIT_W_RESULT ;
             
        case STS_I_DONE :
             fclose(out_stderr);
             return EXIT_OK ;

        case STS_E_ARG :
             fprintf( out_stderr, "*** NEM error status : bad arguments\n" );
             //PrintUsage( argv[ 0 ] ) ;
             fclose(out_stderr);
             return EXIT_E_ARGS ;

        case STS_E_MEMORY :
             fprintf( out_stderr, "*** NEM error status : not enough memory\n" );
#ifdef __TURBOC__
             fprintf( out_stderr, "\n Memory left : %lu bytes\n", 
                      (unsigned long) coreleft() );
#endif
             fclose(out_stderr);
             return EXIT_E_MEMORY ;

        case STS_E_FILEIN :
             fprintf( out_stderr, "*** NEM error status : can't open a file\n" );
             fclose(out_stderr);
             return EXIT_E_FILE ;

        case STS_E_FILE :
             fprintf( out_stderr, "*** NEM error status : wrong file format\n" );
             fclose(out_stderr);
             return EXIT_E_FILE ;

        case STS_E_FUNCARG :
             fprintf( out_stderr, "*** NEM internal error : bad arguments\n" );
             fclose(out_stderr);
             return EXIT_E_BUG ;

        default :
             fprintf( out_stderr, "*** NEM unknown error status (%d)\n", err );
             fclose(out_stderr);
             return EXIT_E_BUG ;
    } 
} /* end of mainfunc() */



/* ==================== LOCAL FUNCTION DEFINITION =================== */

/* ------------------------------------------------------------------- */
static int SetVisitOrder   /*V1.04-e*/
        ( 
	     int         Npt,          /* I */
	     OrderET     VisitOrder,   /* I */
	     int**       SiteVisitVP   /* O and allocated (Npt) */
	) 
/* ------------------------------------------------------------------- */
{
  int ivis ;


  if (( *SiteVisitVP = GenAlloc( Npt, sizeof( int ), 0, 
				 "SetVisitOrder", "SiteVisitVP" ) ) == NULL )
      return STS_E_MEMORY ;

  for ( ivis = 0 ; ivis < Npt ; ivis ++ )
    (*SiteVisitVP)[ ivis ] = ivis ;

  if ( VisitOrder == ORDER_RANDOM )
    {
      RandomPermutationAlgo( (*SiteVisitVP) , Npt ) ;
    }

  return STS_OK ;
}


/* ------------------------------------------------------------------- */
static int  ReadStrFile
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
        fprintf( stderr, "File str %s does not exist\n", infname ) ;
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
static int ReadMatrixFile
         (
             const char  *FileName,    /* I */      /*V1.04-a*/
             int         Nl,           /* I */
             int         Nc,           /* I */
             float       **MatP        /* O and allocated (size : Nl * Nc) */
         ) 
/* ------------------------------------------------------------------- */
{
    float       *Mat ;
    FILE        *fp ;
    int         il ;
    int         ic ;
    float       x ;

    if ( ( *MatP = GenAlloc( Nl *  Nc, sizeof( float ),
			     0, "ReadMatrixFile", FileName ) ) == NULL )
      return STS_E_MEMORY ;
    Mat = *MatP ;

    if ( strcmp( FileName , "-" ) != 0 )   /*V1.04-c*/
      {
	if ( ( fp = fopen( FileName, "r" ) ) == NULL )
	  {
	    fprintf( stderr, "File matrix %s does not exist\n", FileName ) ;
	    GenFree( Mat ) ;
	    return STS_E_FILEIN ;
	  }
      }
    else
      fp = stdin ;

    for ( il = 0, ic = 0 ; ( il < Nl ) && (! feof( fp )) ; il ++ )
    {
        /* printf( "\n pt %d : ", il ) ; */
        for ( ic = 0 ; ( ic < Nc ) && (! feof( fp )) ; ic ++ )
        {
            if ( fscanf( fp, "%f", &x) == 1 )
            {
                Mat[ ( il * Nc ) + ic ] = x ;
                /* printf( " %g", x ) ; */
            }
            else    ic -- ;
        }
    }

    if ( strcmp( FileName , "-" ) != 0 )   /*V1.04-c*/
      fclose( fp ) ;

    if ( ( il < Nl ) || ( ic < Nc ) )
    {
        if ( ic == 0 )
        {
            ic = Nc ;
            il -- ;
        }
        fprintf( stderr, "%s : short file (%d/%d lines and %d/%d columns)\n", 
                 FileName, il, Nl, ic, Nc ) ;
        GenFree( Mat ) ;
        return STS_E_FILE ;
    }
    else
        return STS_OK ;

}  /* end of ReadMatrixFile() */




/* ------------------------------------------------------------------- */
static int  ReadLabelFile
         ( 
             const char  *LabelName,    /* I */
             int         Npt,           /* I */
	     int         *KfileP,       /* O : file # classes */ /*V1.06-c*/
             int         **LabelVP,     /* O and allocated (Npt) */
             float       **ClassifMP    /* O and allocated (Npt*Nk) */
         )
/* ------------------------------------------------------------------- */
{
    FILE  *fp;
    int   pt ;

    if ( ( fp = fopen( LabelName, "r" ) ) == NULL )
    {
        fprintf( stderr, "File label %s does not exist\n", LabelName ) ;
        return STS_E_FILEIN ;
    }

    /* Read number of classes */ /*V1.06-c*/
    fscanf( fp , "%d" , KfileP ) ;

    /* Allocate classification matrix and vector */
    if ( ( *ClassifMP = 
	   GenAlloc( Npt *  (*KfileP), sizeof( float ),
		     0, "ReadLabelFile", LabelName ) ) == NULL )
        return STS_E_MEMORY ;

    if ( ( *LabelVP = 
	   GenAlloc( Npt, sizeof( int ),
		     0, "ReadLabelFile", LabelName ) ) == NULL )
        return STS_E_MEMORY ;

    /* Read labels and store in matrix and vector */
    for ( pt = 0 ; ( pt < Npt ) && (! feof( fp )) ; pt ++ ) {
        /* Read label indication */
        fscanf( fp , "%d" , & (*LabelVP)[ pt ] ) ;

        /* If label is known and valid */
        if ( ( 0 < (*LabelVP)[ pt ] ) && ( (*LabelVP)[ pt ] <= (*KfileP) ) ) 
	  LabelToClassVector( (*KfileP), (*LabelVP)[ pt ] - 1, 
			      & (*ClassifMP)[ pt * (*KfileP) ] ) ;
	/* Else (label unknown or invalid) */
        else {
            /*V1.03-c*/
	    (*LabelVP)[ pt ] = 0 ; /* enforce label to 0 if it is invalid */
	    LabelToClassVector( (*KfileP), (*LabelVP)[ pt ] - 1, 
				& (*ClassifMP)[ pt * (*KfileP) ] ) ;
	}
      }

    fclose( fp ) ;


    if ( pt < Npt )
    {
        fprintf( stderr, "%s : short file (%d/%d labels)\n", 
                 LabelName, pt , Npt ) ;
        GenFree( (*ClassifMP) ) ; (*ClassifMP) = NULL ;
        GenFree( (*LabelVP) ) ;   (*LabelVP)   = NULL ;

        return STS_E_FILE ;
    }
    else
        return STS_OK ;

}  /* end of ReadLabelFile() */

/* ------------------------------------------------------------------- */
static int  ReadParamFile /*V1.08-a*/
         ( 
             const char  *ParamName,    /* I */
  	         const FamilyET Family,
             const int   K,          /* I : number of classes */
             const int   D,          /* I : number of variables */
             ParamFileET* ParamFileMode, /* O : specified if parameters will be used at begin or throughout the clustering process */
             ModelParaT* ParaP       /* O : read parameters */
         )
/* ------------------------------------------------------------------- */
{
   FILE  *fp;

   if ( ( fp = fopen( ParamName, "r" ) ) == NULL )
   {
       fprintf( stderr, "File param %s does not exist\n", ParamName ) ;
       return STS_E_FILEIN ;
   }

  /* Read if file correspond to parameters at beginning (1) or throughout the clustering process (2) */ /*V1.08-a*/

  int nbvalues = 1 + (K - 1) + K * D + K * D ;
  int miOrMf = 0; 
  if (fscanf( fp , "%d" , &miOrMf) == EOF)
  {
     fprintf( stderr, "The file %s needs at least %d values (%d missing)\n",ParamName,1 + (K - 1) + K * D + K * D,nbvalues) ;
     return STS_E_FILEIN;
  }
  nbvalues-- ;
  if (miOrMf==1)
  { 
    *ParamFileMode = PARAM_FILE_INIT ;
  }
  else if (miOrMf==2)
  {
    *ParamFileMode = PARAM_FILE_FIX ;
  }
  else
  {
    fprintf( stderr, "First line of file %s must be 1 (parameters at beginning) or 2 (fixed parameters throughout the clustering process) \n", ParamName ) ;
    return STS_E_FILEIN ;
  }
  StatusET  sts = STS_OK ;
  int       k ;   /* class counter : 0..K-1 */
  int       d ;   /* variable counter : 0..D-1 */
  float     pK ;  /* remaining proportion for class K */
    
  /* Read proportions */
  float value;
  for ( k = 0, pK = 1 ; 
        k < K - 1 ; 
        k ++ ) {
   if (fscanf( fp , "%f" , &value) == EOF)
   {
        fprintf( stderr, "The file %s needs at least %d values (%d missing)\n",ParamName,1 + (K - 1) + K * D + K * D,nbvalues) ;
        return STS_E_FILEIN;
   }
   nbvalues-- ;
   ParaP->Prop_K[ k ] = value ;
   pK = pK - ParaP->Prop_K[ k ] ;
  }
  ParaP->Prop_K[ K - 1 ] = pK ;
  if ( pK <= 0.0 ) {
    //fprintf( stderr, "Last class has pK = %5.2f <= 0\n", pK ) ;
    sts = STS_E_FILE ;
  }

   /* Read centers */
  for ( k = 0 ; k < K ; k ++ ) {
 
    for ( d = 0 ; d < D ; d ++ ) {

      if (fscanf( fp , "%f" , &value) == EOF)
      {
         fprintf( stderr, "The file %s needs at least %d values (%d missing)\n",ParamName,1 + (K - 1) + K * D + K * D,nbvalues) ;
         return STS_E_FILEIN;
      }
      nbvalues-- ;
      ParaP->Center_KD[ k * D + d ] = value ;
    }
  }

   /* Read dispersions */
  for ( k = 0 ; k < K ; k ++ ) {

    for ( d = 0 ; d < D ; d ++ ) {

       if (fscanf( fp , "%f" , &value) == EOF)
       {
         fprintf( stderr, "The file %s needs at least %d values (%d missing)\n",ParamName,1 + (K - 1) + K * D + K * D,nbvalues) ;
         return STS_E_FILEIN;
       }
       
       nbvalues-- ;
       if ( Family == FAMILY_NORMAL ) {
         ParaP->Disp_KD[ k * D + d ] = value * value ;
       }
       else
         ParaP->Disp_KD[ k * D + d ] = value ;
       if ( ParaP->Disp_KD[ k * D + d ] <= 0 ) {

          fprintf( stderr, "Dispersion(k=%d, d=%d) = %5.3f <= 0\n", 
	         k+1, d+1, ParaP->Disp_KD[ k * D + d ] ) ;
         sts = STS_E_FILE ;
      }
    }
  }

  if (fscanf( fp , "%f" , &value) != EOF)
  {
      fprintf( stderr, "The file %s needs do not need more than %d values\n",ParamName,1 + (K - 1) + K * D + K * D) ;
      sts = STS_E_FILEIN;
  }

  fclose( fp ) ;
          
  return sts ;

}  /* end of ReadParamFile() */

/* ------------------------------------------------------------------- */
static int  MakeErrinfo
         ( 
             const char* RefName,       /* I : filename of reference class */
             int         N,             /* I : number of objects */
             int         Kc,            /* I : user number of classes */
	     TieET       TieRule,       /* I : specified MAP tie rule */
	     ErrinfoT*   ErrinfoP,      /* O and allocated */
	     ErrcurT*    ErrcurP        /* O and allocated */
         )
/* ------------------------------------------------------------------- */
{
  StatusET err ;
  int      *tmpV ;
  int      ipt ;

  if ( strcmp( RefName, "" ) != 0 ) {
    ErrinfoP->Kc = Kc ;
    if ( ( err = ReadLabelFile( RefName, N, 
				& ErrinfoP->Kr,
				& tmpV,
				& ErrinfoP->Refclas_N_Kr ) ) != STS_OK )
      return err ;

    /* Check all reference labels ok */
    for ( ipt = 0, err = STS_OK ; 
	  ( ipt < N ) && ( err == STS_OK ) ; ipt ++ ) {
      if ( ( tmpV[ ipt ] <= 0 ) || ( tmpV[ ipt ] > ErrinfoP->Kr ) ) {
	fprintf( stderr, 
		 "Reference class for point %d not in 1..%d \n", 
		 ipt + 1, ErrinfoP->Kr ) ;
	err = STS_E_FILE ;
      }
    }
    GenFree( tmpV ) ; tmpV = NULL ;
    if ( err != STS_OK )
      return err ;

    /* Compute greatest number of classes and permutations of classes */
    ErrinfoP->Km = ( ErrinfoP->Kc > ErrinfoP->Kr ) ? 
      ErrinfoP->Kc : ErrinfoP->Kr ;

    ErrinfoP->Kmfac = factorial( ErrinfoP->Km ) ;

    ErrinfoP->TieRule = TieRule ;

    compute_permutations( 0, ErrinfoP->Km, & ErrinfoP->Perm_Kmfac_Km ) ;

    /* Allocate and initialize later computed stuff */
    if ( ( ErrcurP->Agree_Km_Km = 
	   GenAlloc( ErrinfoP->Km * ErrinfoP->Km, sizeof( float ),
		     0, "MakeErrinfo", "Agree_Km_Km" ) ) == NULL )
        return STS_E_MEMORY ;

    if ( ( ErrcurP->Loclas_N_Kc = 
	   GenAlloc( N * ErrinfoP->Kc, sizeof( float ),
		     0, "MakeErrinfo", "Loclas_N_Kc" ) ) == NULL )
        return STS_E_MEMORY ;

    ErrcurP->Ibestpermut = -1 ;
    ErrcurP->Errorrate   = -2.0 ;

    return STS_OK ;
  }
  else {
    ErrinfoP->Kr = NULL ;
    ErrinfoP->Refclas_N_Kr = NULL ;
    ErrcurP->Errorrate     = -1.0 ;
    return STS_OK ;
  }
}  /* end of MakeErrinfo() */





/* ------------------------------------------------------------------- */
static int factorial(int n)
/* ------------------------------------------------------------------- */
{
  int result = 1;

  for (; n>0; n--)
    result *= n;

  return result;
}


/* ------------------------------------------------------------------- */
static int compute_permutations    /* ret 0 if OK, -1 if memory error */
(
 const int Startval,         /* I : start value of integer suite */
 const int K,                /* I : size of integer suite to permute > 0 */
 int*      perms_Kfact_K_p[] /* O : matrix to store permuted values */
)
/* ------------------------------------------------------------------- */
{
  int   ik ;      /* current index and value of integer suite : 0..K-1 */
  int*  array_K ; /* integer suite to permute : [ 0 ... (K-1) ] */
  int   Kfact ;   /* K! */
  int   err ;     /* 0 if OK, -1 if memory error */

  if ( K <= 0 )
    return 1 ;

  Kfact = factorial( K ) ;

  if( ( (*perms_Kfact_K_p) = malloc( Kfact * K * sizeof( int ) ) ) == NULL )
    return -1 ;

  if( ( array_K = malloc( K * sizeof( int ) ) ) == NULL ) {
    free( (*perms_Kfact_K_p) ) ;
    (*perms_Kfact_K_p) = NULL ;
    return -1 ;
  }

  for ( ik = 0 ; ik < K ; ik ++ )
    array_K[ ik ] = Startval + ik ;

  err = rec_permutations( array_K, K, K, 0, (*perms_Kfact_K_p) ) ;

  free( array_K ) ;

  return err ;
}


/* ------------------------------------------------------------------- */
static int rec_permutations        /* ret 0 if OK, -1 if memory error */
(
 const int array_A[],       /* I : remaining array to permute */
 const int A,               /* I : length of the remaining array : 0..K */
 const int K,               /* I : length of original array */
 int       offset,          /* I : first line of storage to use */
 int       perms_Kfact_K[]  /* O : matrix to store permuted values, use
			       lines :   offset -> offset + A! - 1
			       columns : K - A  -> K - 1 */
)
/* ------------------------------------------------------------------- */
{
  int    err ;         /* 0 if currently OK */
  int    ia ;          /* array element currently removed : 0..A-1 */
  int    ja ;          /* array element currently copied : 0..A-1 (neq ia) */
  int    am1fact ;     /* (A-1)! */
  int    iam1fact ;    /* current perms line (without offset) : 0..am1fact-1 */
  int*   redarr_am1 ;  /* array with one element removed (A-1) */

  am1fact = factorial( A - 1 ) ;

  /* Check against out of bounds offset */
  if ( ( offset < 0 ) || ( factorial( K ) < ( offset + A * am1fact ) ) )
    return 1 ;
  
  if ( ( redarr_am1 = malloc( ( A - 1 ) * sizeof( int ) ) ) == NULL ) 
    return -1;

  /* For each element of given array */
  for ( ia = 0, err = 0 ; ( ia < A ) && ( err == 0 ) ; ia ++ ) {

    /* Copy (A-1)! times this element into the column ( K - A ) of perms, 
       starting from line ( offset + ia * (A-1)! ) (skip previous ia's) */
    for ( iam1fact = 0 ; iam1fact < am1fact ; iam1fact ++ ) 
      perms_Kfact_K[ ( offset + ia * am1fact + iam1fact ) * K + ( K - A ) ] =
	array_A[ ia ] ;

    /* Copy array into reduced array without this element */
    for ( ja = 0 ; ja < ia ; ja ++ )
      redarr_am1[ ja ] = array_A[ ja ] ;
    for ( ja = ia + 1 ; ja < A ; ja ++ )
      redarr_am1[ ja - 1 ] = array_A[ ja ] ;

    /* Recursive call with reduced array */
    err = rec_permutations( redarr_am1, A-1, K, offset + ia * am1fact,
			    perms_Kfact_K ) ;
  }

  free( redarr_am1 ) ;
  return err ;
}




/* ------------------------------------------------------------------- */
static int  ReadNeiFile
         ( 
             const char  *BaseName,          /* I */
             int         NbPts,              /* I */
             NeighET     NeighSpec,          /* I */
	     char*       NeiDescS,           /* O [LEN_LINE+1] */
             SpatialT    *SpatialP           /* I/O and allocated */
         ) 
/* ------------------------------------------------------------------- */
{
    char        infname[ LEN_FILENAME + 1 ] ;
    FILE        *fnei ;
    StatusET    err = STS_OK ;

    /* Read neighbourhood data depending on spatial type */
    if ( SpatialP->Type == TYPE_SPATIAL )
    {
        strncpy( infname, BaseName, LEN_FILENAME ) ;
        strncat( infname, ".nei" , LEN_FILENAME ) ;

        if ( ReadOpeningComments( infname , "#" , LEN_LINE , 
                                  & fnei , NeiDescS ) == -1 )
        {
            fprintf( stderr, "File Neigh %s File does not exist\n", infname ) ;
            return STS_E_FILEIN ;
        }

        err = ReadPtsNeighs( fnei, NbPts, 
                             &SpatialP->MaxNeighs, 
                             &SpatialP->NeighData ) ;
        fclose( fnei ) ;
    }
    else    /* Type is image */
    {
        if ( NeighSpec == NEIGH_FILE )
        {
            strncpy( infname, BaseName, LEN_FILENAME ) ;
            strncat( infname, ".nei" , LEN_FILENAME ) ;
            if ( ReadOpeningComments( infname , "#" , LEN_LINE , 
		                      & fnei , NeiDescS ) == -1 )
            {
                fprintf( stderr, "File Neigh %s does not exist\n", infname ) ;
                return STS_E_FILEIN ;
            }

            err = ReadImageNeigh( fnei, &SpatialP->NeighData ) ;
            fclose( fnei ) ;
        }
        else
        {
            err = SetImageNeigh( NeighSpec, NeiDescS, &SpatialP->NeighData ) ;
        }

        SpatialP->MaxNeighs = SpatialP->NeighData.Image.NbNeigh ;

    }


    return err ;

}  /* end of ReadNeiFile() */


/* ------------------------------------------------------------------- */
static int      ReadPtsNeighs
                (
                    FILE        *Fnei,          /* I/O */
                    int         NbPts,          /* I */
                    int         *MaxNeiP,       /* O */
                    NeighDataT  *NeighDataP     /* O and allocated */
                ) 
/* ------------------------------------------------------------------- */
{
    const char* func = "ReadPtsNeighs" ;
    int       weighted ;
    int       l ;
    int       ipt ;
    StatusET  err = STS_OK ;
    PtNeighsT *ptsneighsV ;
    int       nmax ;


    /* Indicator of weights presence */
    fscanf( Fnei, "%d", &weighted ) ;

    /* Allocate structure of all points' neighbours */
    if ( ( ptsneighsV = GenAlloc( NbPts, sizeof( PtNeighsT ), 
				  0, func, "ptsneighsV" ) )
           == NULL )
    {
        fprintf( stderr, "Cannot allocate list of neighbours of %d points\n" ,
                 NbPts ) ;
        err = STS_E_MEMORY ;
        return err ;
    }
    NeighDataP->PtsNeighsV = ptsneighsV ;

    /* For each point, default nb of neighbours is zero */
    for ( ipt = 0 ; ipt < NbPts ; ipt ++ )
    {
        ptsneighsV[ ipt ].NbNeigh = 0 ;
    }

    /* Read neighbours of each point */
    nmax = 0 ;
    for( l = 0 ; ( err == STS_OK ) && ( ! feof( Fnei ) ) ; l ++ )
    {
        int     ipt ;

        if ( fscanf( Fnei, "%d", &ipt ) == 1 )
        {
            int nbv ;   /* Declared nb of neighbours */
            int nv ;    /* Effective nb of neighbours */

            if ( fscanf( Fnei, "%d", &nbv ) == 1 )
            {
                NeighT   *neighsV ;
                int      iv ;

                if ( ( neighsV = GenAlloc( nbv, sizeof( NeighT ),
					   0, func, "neighsV" ) ) == NULL )
                {
                    fprintf( stderr, "Can't allocate %d neighb. for pt %d\n" ,
                             nbv, ipt ) ;
                    err = STS_E_MEMORY ;
                    return  err ;
                }

                ptsneighsV[ ipt - 1 ].NeighsV = neighsV ;

                for ( iv = 0, nv = 0 ; 
                      ( iv < nbv ) && ( ! feof( Fnei ) ) ; 
                      iv ++ )
                {
                    int     iptv ;

                    if ( fscanf( Fnei, "%d", &iptv ) == 1 )
                    {
                        if ( ( 1 <= iptv ) && ( iptv <= NbPts ) )
                        {
                            neighsV[ nv ].Index = iptv - 1 ;
                            nv ++ ;
                        }
                    }
                    else
                    {
                        fprintf( stderr, 
                                 "Error in neighb. file l.%d : neighbor %d\n",
                                 l, iv ) ;
                        err = STS_E_FILE ;
                    }
                } /* end  for ( iv ... ) */

                if ( weighted )
                {
                    for ( iv = 0, nv = 0 ; 
                         ( iv < nbv ) && ( ! feof( Fnei ) ) ; 
                         iv ++ )
                    {
                        float   weight ;

                        if ( fscanf( Fnei, "%g", &weight) == 1 )
                        {
                            if ( weight != 0.0 )
                            {
                                neighsV[ nv ].Weight = weight ;
                                nv ++ ;
                            }
                        }
                        else
                        {
                          fprintf( stderr, 
                                  "Error in neighb. file l.%d : weight %d\n",
                                  l, iv ) ;
                            err = STS_E_FILE ;
                        }
                    } /* end  for ( iv ... ) */
                } /* end if weighted */
                else
                {
                    for ( iv = 0 ; iv < nv ; iv ++ )
                        neighsV[ iv ].Weight = 1.0 ;
                }

                ptsneighsV[ ipt - 1 ].NbNeigh = nv ;

                if ( nv > nmax )     nmax = nv ;

            } /* end if ( fscanf( ... , &nbv ) == 1 ) */
            else 
                 {}
        } /* end if ( fscanf( ..., &ipt ) == 1 ) */
        else
            {}
    } /* end for( l = 0 ; ( ! feof( Fnei ) ) ; l ++ ) */

    *MaxNeiP = nmax ;  

    return err ;

}   /* end of ReadPtsNeighs() */


/* ------------------------------------------------------------------- */
static int      ReadImageNeigh
                (
                    FILE        *Fnei,          /* I/O */
                    NeighDataT  *NeighDataP     /* O and allocated */
                ) 
/* ------------------------------------------------------------------- */
{
    StatusET  err = STS_OK ;
    int       dlmin, dlmax, dcmin, dcmax ;
    int       dl, dc ;
    int       tnb ;     /* window size */
    int       nv ;      /* nb of non-zero weights */
    INeighT   *tneighV;  /* to be allocated */

    fscanf( Fnei, "%d %d %d %d", &dlmin, &dlmax, &dcmin, &dcmax ) ;
    tnb = ( dlmax - dlmin + 1 ) * ( dcmax - dcmin + 1 ) ;

    if ( ( tneighV = GenAlloc( tnb, sizeof( INeighT ),
			       0, "ReadImageNeigh", "tneighV" ) ) == NULL )
        return STS_E_MEMORY ;

    NeighDataP->Image.NeighsV = tneighV ;

    nv = 0 ;
    for ( dl = dlmin ; dl <= dlmax ; dl ++ )
    {
        for ( dc = dcmin ; dc <= dcmax ; dc ++ )
        {
            float   weight ;

            if ( fscanf( Fnei, "%g", &weight ) == 1 )
            {
                if ( weight != 0.0 )
                {
                    tneighV[ nv ].Dl = dl ;
                    tneighV[ nv ].Dc = dc ;
                    tneighV[ nv ].Weight = weight ;
                    nv ++ ;
                }
                else /* else this neighbour has no weight : skip it */
                { 
                }
            }
            else
            {
                fprintf( stderr, "Neighbors file error (dl = %d, dc = %d)\n",
                         dl, dc ) ;
                err = STS_E_FILE ;
            }
        } /* for dc */
    } /* for dl */

    NeighDataP->Image.NbNeigh = nv ;

    return err ;

}   /* end of ReadImageNeighs() */


/* ------------------------------------------------------------------- */
static int  SetImageNeigh
            ( 
                NeighET     NeighSpec,          /* I */
                char*       NeiDescS,           /* O [LEN_LINE+1] */
                NeighDataT* NeighDataP          /* O and allocated */
            ) 
/* ------------------------------------------------------------------- */
{
    INeighT   *neighV;  /* to be allocated */

    switch( NeighSpec )
    {
    case NEIGH_FOUR :
        if ( ( neighV = GenAlloc( 4, sizeof( INeighT ),
				  0, "SetImageNeigh", "neighV" ) ) == NULL )
        {
            fprintf( stderr, "Could not allocate %d image neighbours\n", 
                     4 ) ;
            return STS_E_MEMORY ;
        }
        NeighDataP->Image.NeighsV = neighV ;
        NeighDataP->Image.NbNeigh = 4 ;

        neighV[ 0 ].Dl = -1 ;
        neighV[ 0 ].Dc = 0 ;
        neighV[ 0 ].Weight = 1.0 ;

        neighV[ 1 ].Dl = 0 ;
        neighV[ 1 ].Dc = -1 ;
        neighV[ 1 ].Weight = 1.0 ;

        neighV[ 2 ].Dl = 0 ;
        neighV[ 2 ].Dc = +1 ;
        neighV[ 2 ].Weight = 1.0 ;

        neighV[ 3 ].Dl = +1 ;
        neighV[ 3 ].Dc = 0 ;
        neighV[ 3 ].Weight = 1.0 ;

        strncpy( NeiDescS , 
                 "  Default 1st-order neighbors (horizontal and vertical)\n",
                 LEN_LINE ) ;
        break ;

    default :
        fprintf( stderr, "Unknown neighborhood type (%d)\n", NeighSpec ) ;
        return STS_E_FUNCARG ;
    }

    return STS_OK ;
}   /* end of SetNeigh() */


/* ------------------------------------------------------------------- */
static int SaveResults
        ( const int          Npt,                   /* I */
          const int          Nd,                    /* I */
          const float*       ClassifM,              /* I */
          const SpatialT*    SpatialP,              /* I */
          const NemParaT*    NemParaP,              /* I */
          const StatModelT*  ModelP,                /* I */
	  const CriterT*     CriterP                /* I */ /*V1.03-d*/
        ) 
/* ------------------------------------------------------------------- */
{
    FILE*       fout ;
    FILE*       fmf ;
    char        mfname[ LEN_FILENAME ] ;
    int         ipt ;           /* point counter : 0..Npt */
    int         k ;             /* class counter : 0..nk-1 */
    int         d ;             /* dimension counter : 0..Nd-1 */
    int         nk = ModelP->Spec.K ;
    StatusET    err = STS_OK ;


    /* SaveResults algorithm :
       - save classification matrix, in hard or fuzzy partition format
       - save means of each class
    */

    if ( strcmp( NemParaP->OutBaseName, "-" ) != 0 )  /*V1.04-b*/
      {
	if ( ( fout = fopen( NemParaP->OutName, "w" ) ) == NULL )
	  {
	    fprintf( stderr, "Could not open file '%s' in write mode\n", 
		     NemParaP->OutName ) ;
	    return STS_E_FILEOUT ;
	  }
      }
    else
      fout = stdout ;

    if ( NemParaP->Format == FORMAT_HARD )
      {
        int*   kmaxesV;  /* classes having same maximum probabilites */

        if ( ( kmaxesV = GenAlloc( nk, sizeof( int ),
				   0, "SaveResults", "kmaxesV" ) ) == NULL ) 
            return STS_E_MEMORY ;

        /* For each record */
        for ( ipt = 0 ; ( ipt < Npt ) && ( err == STS_OK ) ; ipt ++ )
          {
            int   kmax ;

            /* Compute its MAP class */
            kmax = ComputeMAP( ClassifM, ipt, nk, NemParaP->TieRule, 
			       kmaxesV ) ;

            /* Save the record's class (add 1 to have k in 1..Nk) */
            if ( fprintf( fout, "%d ", kmax + 1 ) == EOF )
              {
                fprintf( stderr, "Cannot write Ci i = %d\n",
                            ipt + 1 ) ;
                err = STS_E_FILEOUT ;
              }

            if ( ( SpatialP->Type == TYPE_IMAGE ) &&
                 ( ((ipt + 1) % SpatialP->NeighData.Image.Nc) == 0 ) )
              {
                fprintf( fout, "\n" ) ;
              }
          } /* end for ( ipt ... ) */

	fprintf( fout, "\n" ) ; /*V1.07-b*/

        GenFree( kmaxesV ) ;

      } /* end if Format == HARD */
    else
      {
        for ( ipt = 0 ; ( ipt < Npt ) && ( err == STS_OK ) ; ipt ++ )
          {
            for ( k = 0 ; ( k < nk ) && ( err == STS_OK ) ; k ++ )
              {
                if ( fprintf( fout, " %5.3f ", 
                              ClassifM[ (ipt * nk) + k ] ) == EOF )
                  {
                    fprintf( stderr, "Cannot write Uik i = %d, k =%d\n",
                            ipt + 1, k + 1 ) ;
                    
                    err = STS_E_FILEOUT ;
                  }
              } /* end for ( k = 0 ... ) */
            fprintf( fout, "\n" ) ;
          } /* end for ( ipt ... ) */
      } /* end else Format != HARD */

    if ( fout != stdout )
      fclose( fout ) ;

    if ( strcmp( NemParaP->OutBaseName , "-" ) != 0 )  /*V1.04-b*/
      {
	strncpy( mfname, NemParaP->OutBaseName, LEN_FILENAME ) ;
	strncat( mfname, EXT_MFNAME, LEN_FILENAME ) ;
	if ( ( fmf = fopen( mfname, "w" ) ) == NULL )
	  {
	    fprintf( stderr, "Could not open file '%s' in write mode\n", 
		     mfname ) ;
	    return STS_E_FILEOUT ;
	  }
      }
    else
      fmf = stderr ;

    /*V1.03-d*/ /*V1.03-e*/ /*V1.05-j*/
    fprintf( fmf, 
	     "Criteria U=NEM, D=Hathaway, L=mixture, M=markov ps-like, Z=log pseudo-l, error\n\n" );
    fprintf( fmf, "  %g    %g    %g    %g   %g   %g\n\n", 
	     CriterP->U, CriterP->D, CriterP->L, CriterP->M, CriterP->Z,
	     CriterP->Errcur.Errorrate ) ; 

    fprintf( fmf, "Beta (%s)\n", BetaDesVC[ ModelP->Spec.BetaModel ] );
    fprintf( fmf, "  %6.4f\n", ModelP->Para.Beta ) ;

    /*V1.06-a*/
    switch( ModelP->Spec.ClassFamily )
      {
      case FAMILY_NORMAL:
	fprintf( fmf, "Mu (%d), Pk, and sigma (%d) of the %d classes\n\n", 
		 Nd, Nd, nk ) ;
	break;

      case FAMILY_LAPLACE:
	fprintf( fmf, "Mu (%d), Pk, and lambda (%d) of the %d classes\n\n", 
		 Nd, Nd, nk ) ;
	break;

      case FAMILY_BERNOULLI: /*V1.07-a*/
	fprintf( fmf, "Mu (%d), Pk, and disp (%d) of the %d classes\n\n", 
		 Nd, Nd, nk ) ;
	break;

      default:
	fprintf( fmf, "Mu (%d), Pk, and disp (%d) of the %d classes\n\n", 
		 Nd, Nd, nk ) ;
      }

    for ( k = 0 ; k < nk ; k ++ )
      {
        /* Save mean, proportion, volume, and normalized variance matrix */
        /* V1.03-b*/    /*V1.06-a*/
        for ( d = 0 ; d < Nd ; d ++ )
        {
            fprintf( fmf, " %10.3g ", ModelP->Para.Center_KD[ (k*Nd) + d ] ) ;
        }
        fprintf( fmf, "  %5.3g  ", ModelP->Para.Prop_K[ k ] ) ;
        for ( d = 0 ; d < Nd ; d ++ )
        {
	  switch( ModelP->Spec.ClassFamily )
	    {
	    case FAMILY_NORMAL:
	      fprintf( fmf, " %10g ", 
		       sqrt( ModelP->Para.Disp_KD[ (k*Nd) + d ] ) ) ;
	      break;

	    case FAMILY_LAPLACE:
	      fprintf( fmf, " %10g ", 
		       ModelP->Para.Disp_KD[ (k*Nd) + d ] ) ;
	      break;

	    case FAMILY_BERNOULLI: /*V1.07-a*/
	      fprintf( fmf, " %10g ", 
		       ModelP->Para.Disp_KD[ (k*Nd) + d ] ) ;
	      break;

	    default:
	      fprintf( fmf, " %10g ", 
		       ModelP->Para.Disp_KD[ (k*Nd) + d ] ) ;
	    }
        }
        fprintf( fmf, "\n" ) ;
      }

    if ( fmf != stderr )
      fclose( fmf ) ;

    return err ;

}  /* end of SaveResults() */

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


/* ------------------------------------------------------------------- */
static void FreeAllocatedData
(
 DataT*       DataP,       /* O and deallocated */
 SpatialT*    SpatialP,    /* O and deallocated */
 ModelParaT*  ModelParaP,  /* O and deallocated */  /*V1.06-a*/
 CriterT*     CriterP,     /* O and deallocated */  /*V1.06-d*/
 float*       ClassifM     /* O and deallocated */
)
/* ------------------------------------------------------------------- */
{
  int ipt ;  /* counter 0..Npt-1 to free each point neighbors */


  /* Free components of DataP */
  GenFree( DataP->PointsM ) ;     DataP->PointsM    = NULL ;
  GenFree( DataP->LabelV ) ;      DataP->LabelV     = NULL ;
  GenFree( DataP->SiteVisitV ) ;  DataP->SiteVisitV = NULL ;
  GenFree( DataP->SortPos_ND ) ;  DataP->SortPos_ND = NULL ;

  /* Free components of SpatialP */
  switch( SpatialP->Type )
    {
    case TYPE_SPATIAL: 
      /* deallocate each point's neighbors */
      for ( ipt = 0; ipt < DataP->NbPts ; ipt ++ )
	GenFree( SpatialP->NeighData.PtsNeighsV[ ipt ].NeighsV ) ;

      /* deallocate array of array of neighbors */
      GenFree( SpatialP->NeighData.PtsNeighsV ) ;

      break ;

    case TYPE_IMAGE:
      /* deallocate neighborhood window */
      GenFree( SpatialP->NeighData.Image.NeighsV ) ;
      break ;

    default: /* non spatial : no allocated spatial info */
      ;
    }

  /* Free components of ModelParaP */
  GenFree( ModelParaP->Center_KD ) ;
  GenFree( ModelParaP->Disp_KD ) ;
  GenFree( ModelParaP->Prop_K ) ;

  GenFree( ModelParaP->NbObs_K ) ;
  GenFree( ModelParaP->NbObs_KD ) ;
  GenFree( ModelParaP->Iner_KD ) ;


  /* Free components of CriterP */
  GenFree( CriterP->Errinfo.Refclas_N_Kr  ) ;
  GenFree( CriterP->Errinfo.Perm_Kmfac_Km ) ;
  GenFree( CriterP->Errcur.Agree_Km_Km    ) ;
  GenFree( CriterP->Errcur.Loclas_N_Kc    ) ;

  CriterP->Errinfo.Refclas_N_Kr  = NULL ;
  CriterP->Errinfo.Perm_Kmfac_Km = NULL ;
  CriterP->Errcur.Agree_Km_Km    = NULL ;
  CriterP->Errcur.Loclas_N_Kc    = NULL ;

  /* Free classification matrix */
  GenFree( ClassifM ) ;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
