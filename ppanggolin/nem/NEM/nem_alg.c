/*\

    NEM_ALG.C

    Programme NEM (Neighborhood EM) : routines de calcul algorithmique

    Van Mo DANG       Janvier 96


Vers-mod  Date         Who Description

1.03-a    07-JAN-1997  MD  Print criteria to screen, compute them
                           at each iteration only if logfile requested
1.03-b    07-JAN-1997  MD  Write to logfile only if not NULL (bug fix)
1.03-c    22-AUG-1997  MD  Removed ComputePartitionParallel() (outdated)
1.03-d    22-AUG-1997  MD  Add MakeParaFromLabeled() called by ClassifyByNem()
1.03-e    22-AUG-1997  MD  Fix error - Pk was not init in MakeRandomPara()
1.03-f    22-AUG-1997  MD  ClassifyByNem() now returns criteria
1.03-g    30-SEP-1997  MD  Criterion L renamed D, and add new criterion L
1.04-a    05-OCT-1997  MD  Heuristic estimation of beta
1.04-b    10-OCT-1997  MD  Convergence threshold as modifiable option
1.04-c    10-OCT-1997  MD  Sites visit order specified
1.04-d    13-OCT-1997  MD  Debug mode for beta heuristic
1.04-e    13-OCT-1997  MD  Compute classification error for 2 clusters
1.04-f    13-OCT-1997  MD  Choice of criterion to select best run
1.04-g    13-OCT-1997  MD  Save all criteria of best try in random init mode
1.04-h    01-DEC-1997  MD  Add ComputePartitionGEM (Monte-Carlo to get cih)
1.05-a    12-JAN-1998  MD  Add arg MissMode to ParaP_VkI() and fEstimNoise()
1.05-b    12-JAN-1998  MD  In MakeParaFromLabeled() Pk do not vary any more
1.05-c    12-JAN-1998  MD  In MakeParaFromLabeled() process missing data
1.05-d    12-JAN-1998  MD  In MakeRandomPara() process missing data
1.05-e    17-JAN-1998  MD  emptyk as arg to P*V*I and fEstimNoise
1.05-f    26-JAN-1998  MD  GenAlloc/GenFree instead of malloc/calloc/free
1.05-g    26-JAN-1998  MD  Changed Iterations appearance
1.05-h    27-JAN-1998  MD  Fix bug: clasiniM was freed *before* reading it
1.05-i    03-FEB-1998  MD  Initial vol = whole vol / K (was whole vol before)
1.05-j    03-FEB-1998  MD  Fix bug: L computed even for i/h with zero cih !
1.05-k    05-FEB-1997  MD  Add Markov fuzzy pseudo-like. (CritET and CriterT)
1.05-l    05-FEB-1997  MD  call RandomFloat() instead of drand48()
1.05-m    05-FEB-1997  MD  check that drawn center does not exist
1.06-a    03-JUN-1998  MD  Initial vol = whole vol / K (was whole vol before)
1.06-b    28-JUN-1998  MD  Reflect changes to nem_typ.h (NoiseModel)
1.06-c    28-JUN-1998  MD  fGetNeigh and fCompuDens fetched local to their use
1.06-d    29-JUN-1998  MD  MakeParaFromLabeled uses EstimPara to compute param.
1.06-e    30-JUN-1998  MD  Call MakeParaFromLabeled in INIT_SORT mode
1.06-f    30-JUN-1998  MD  Add ModelPreprocess 
1.06-g    01-JUL-1998  MD  Call MakeParaFromLabeled in INIT_FILE mode
1.06-h    17-JUL-1998  MD  Bug fix : CompSortValue now puts NaN at end
1.06-i    19-JUL-1998  MD  Change heuristic, save classif. just before L drops
1.06-j    03-AUG-1998  MD  Add INIT_MIXINI and INIT_MIXFIX
1.06-k    10-SEP-1998  MD  ComputePartitionNEM/GEM call ComputeLocalProba
1.06-l    10-SEP-1998  MD  Add SiteUpdate parallel or sequential
1.06-m    15-SEP-1998  MD  Add EstimBeta
1.06-n    21-SEP-1998  MD  Add TieRule arg to ComputeMAP
1.06-o    21-SEP-1998  MD  Change convergence test (test criterion also)
1.06-p    05-OCT-1998  MD  MakeParaFromLabeled draws random if empty class
                           and INIT_LABEL does not fail in that case
1.06-q    22-OCT-1998  MD  Bug : RandNemAlgo must save/restore best cen./disp.
1.06-r    03-NOV-1998  MD  Process case 0 < cumnum < EPSILON
1.06-t    01-DEV-1998  MD  pkfki and cinumv now double (was float)
1.08-a    20-JUI-2017  GG   Add param input by file rather than by arguments
\*/

#include "nem_typ.h"    /* DataT, ... */
#include "nem_nei.h"    /* FuncNeighNone, ... */
#include "nem_mod.h"    /* EstimPara, ... */
#include "nem_rnd.h"    /* RandomInteger */
#include "genmemo.h"    /* GenAlloc, GenFree */
#include <stdio.h>      /* printf, ... */
#include <string.h>     /* memcpy */
#include <math.h>       /* exp */
#include <time.h>       /* time */
#include <float.h>      /* FLT_MAX */

#define MINFLOAT FLT_MIN   /* alias value for compatibility issues */

#ifdef __TURBOC__
#include <alloc.h>      /* coreleft, ... */ 
#endif

#include "nem_alg.h"    /* prototypes */

#define NBCYCLES_STARTUP   1
#define NBCYCLES_RECORD    2


#define MAX_BETA           5    /* To check against divergence of beta */
#define MIN_BETA           -5   /* To check against divergence of beta */

#ifndef MAXFLOAT
 #ifdef FLT_MAX
   #define MAXFLOAT FLT_MAX
 #else
   #define MAXFLOAT 3.40282347e+38F
 #endif
#endif


typedef struct
{
    int     Index ;
    float   Value ;
} 
SortPtT ;

typedef struct
{
    int*        KmaxesV ;   /* values of ex-aequo MAP labels [nk] */
    double*     PkFkiM;     /* pk fk(xi) [npt,nk] */ /*V1.06-t*/
    float*      LogPkFkiM;  /* log( pk fk(xi) ) [npt,nk] */
    float*      CtmpM ;     /* last clas. for parallel update [npt,nk] */
    float*      ColdM ;     /* clas. of previous XEM iteration [npt,nk] */
    double*     CiNumV;     /* numerator of cik's for one i [nk] */ /*V1.06-t*/
    PtNeighsT   Neighs;     /* neigh. ind/wei, NeighsV [SpatialP->MaxNeighs] */
}
WorkingT ;


static float mknan()
{
   return atof( "nan" ) ;
}


/* ==================== LOCAL FUNCTION PROTOTYPING =================== */


/* Called by ClassifyByNem */

  static void ModelPreprocess
        ( 
	  const ModelSpecT*   SpecP,            /* I */ 
	  DataT*              DataP             /* I/O */
	) ; 

  static int ClassifyByNemOneBeta                              /*V1.04-a*/
        ( 
          const DataT         *DataP,           /* I */
          const NemParaT      *NemParaP,        /* I */
          const SpatialT      *SpatialP,        /* I */
          StatModelT          *StatModelP,      /* I/O */
          float               *ClassifM,        /* I/O */
	  CriterT             *CriterP          /* O */  /*V1.03-f*/
        ) ;

  static int ClassifyByNemHeuBeta                             /*V1.04-a*/
        ( 
          const DataT         *DataP,           /* I */
          NemParaT            *NemParaP,        /* I[M] */
          const SpatialT      *SpatialP,        /* I */
          StatModelT          *StatModelP,      /* I/O */
          float               *ClassifM,        /* I/O */
	  CriterT             *CriterP          /* O */  /*V1.03-f*/
        ) ;


/* Called by ModelPreprocess */

	static int CompSortValue( const void* elt1P, const void* elt2P ) ;



/* Called by ClassifyByNemHeuBeta */

  static int ClassifyByNemOneBeta                               /*V1.04-a*/
        ( 
          const DataT         *DataP,           /* I */
          const NemParaT      *NemParaP,        /* I */
          const SpatialT      *SpatialP,        /* I */
          StatModelT          *StatModelP,      /* I/O */
          float               *ClassifM,        /* I/O */
	  CriterT             *CriterP          /* O */  /*V1.03-f*/
        ) ;


/* Called by ClassifyByNemOneBeta */

    static void StartLogFile( const char* LogName, int Npt, FILE** FlogP ) ;

    static void InitPara
        (
         const DataT*       DataP,       /* I */
	 const SampleDesT*  DescP,       /* I */
	 const ModelSpecT*  SpecP,       /* I */
	 ModelParaT*        ParaP,       /* O */
         float*             C_NK         /* T */
        ) ;

    static int InitPartitionSort
        (
          const DataT*    DataP,          /* I */
          int             Nk,             /* I */
          int             SortedVar,      /* I : sorted variable : 0..dim-1 */
          float           *ClassifM       /* O (to allocate before call) */
        ) ;


    static int ComputePartitionFromPara
        ( 
	 const int          Needinit,    /* I : 1 if need init, 0 if normal */
	 const DataT*       DataP,       /* I */
	 const NemParaT*    NemParaP,    /* I */
	 const ModelSpecT*  SpecP,       /* I */
	 ModelParaT*        ParaP,       /* I (temporary modification) */
         const SpatialT*    SpatialP,    /* I */
         float*             C_NK,        /* O */
         CriterT*           CriterP,     /* O */
	 FILE*              Flog,        /* I/O */
	 WorkingT*          WorkingP     /* T */
	) ;


    static int NemAlgo
        (
         const DataT*       DataP,           /* I */
         const NemParaT*    NemParaP,        /* I */
         const SpatialT*    SpatialP,        /* I */
	 const ModelSpecT*  SpecP,     	     /* I */
	 const SampleDesT*  DescP,     	     /* I */
         FILE*              Flog,            /* I/O */
	 ModelParaT*        ParaP,     	     /* I/O */
         float*             CM,              /* I/O */
         WorkingT*          WorkP,           /* O */
         CriterT*           CriterP          /* O */
        ) ;

    static int MakeParaFromLabeled  /*V1.03-d*/
        ( 
            const DataT*        DataP,      /* I */
	    const float*        C_NK,       /* I */
	    const ModelSpecT*   SpecP,      /* I */
            const SampleDesT*   DescP,      /* I */
	    ModelParaT*         ParaP,      /* O */
	    int                 *misskP,    /* O */
	    int                 *missdP     /* O */
        ) ;

    static StatusET MakeRandomPara
        ( 
            const DataT*        DataP,      /* I */
	    const ModelSpecT*   SpecP,      /* I */
            const SampleDesT*   DescP,      /* I */
	    ModelParaT*         ParaP       /* O */
        ) ;

    static StatusET ComputePkFkiM
        (
            const DataT*        DataP,      /* I */
            const ModelSpecT*   SpecP,      /* I */
            const ModelParaT*   ParaP,      /* I */
            double*             PkFkiM,     /* O */
            float*              LogPkFkiM   /* O */
        ) ;

    static int ComputePartition
        (
          const ModelSpecT*   SpecP,      /* I */
          const ModelParaT*   ParaP,      /* I */
          const DataT*        DataP,      /* I */
          const SpatialT*     SpatialP,   /* I */
          const NemParaT*     NemParaP,   /* I */
          FILE*               Flog,       /* I/O */
          float*              CM,         /* I/O */
          WorkingT*           WorkP,      /* I:PkFkiM O:CtmpM,Neighs,CiNumV */
          CriterT*            CriterP     /* O */
        ) ;

    static void WriteLogClasses
        (
            FILE*               Flog,       /* I/O */
            int                 Nd,         /* I : dimension */
            const ModelSpecT*   SpecP,      /* I */
            const ModelParaT*   ParaP       /* I */
        ) ;


          /*V1.04-f*/
	  static float ChosenCrit( const CriterT* CriterP , CritET Which ) ; 


static StatusET RandNemAlgo
        (
         const DataT*       DataP,           /* I */
         const NemParaT*    NemParaP,        /* I */
         const SpatialT*    SpatialP,        /* I */
	 const ModelSpecT*  SpecP,     	     /* I */
	 const SampleDesT*  DescP,     	     /* I */
         FILE*              Flog,            /* I/O */
	 ModelParaT*        ParaP,     	     /* I/O */
         float*             CM,              /* I/O */
         WorkingT*          WorkP,           /* I/O */
         CriterT*           CriterP          /* O */
        ) ;


/* Called by InitPartitionSort */

    static int CompSortValue( const void* elt1P, const void* elt2P ) ;



/* Called by NemAlgo */

    static void WriteLogHeader
        (
            FILE*               Flog,       /* I/O */
            int                 NbEIters,   /* I */
            int                 Nd,         /* I : dimension */
            const ModelSpecT*   SpecP       /* I */
        ) ;

static int HasConverged       /* ret : TRUE if convergence test satisfied */
(
  const CvemET    CvTest,     /* I : which convergence test to use */
  const float     CvThres,    /* I : convergence threshold to use */
  const float*    ColdM,      /* I : previous classification matrix */
  const float*    CM,         /* I : current  classification matrix */
  const int       Npt,        /* I : number of objects */
  const int       Nk,         /* I : number of classes */
  const float     OldCrit,    /* I : previous criterion */ 
  const int       CritToDo,   /* I : 1 = criterion needs to be computed */
  const float     Beta,       /* I : value of beta to compute criterion */
  const SpatialT* SpatialP,   /* I : spatial info to compute criterion */
  const CritET    WhichCrit,  /* I : what criterion to use */
  CriterT*        CriterP,    /* I/O : current criterion */
  WorkingT*       WorkP       /* I: LogPkFkiM, T: Neighs */
) ;

#if 0
    static StatusET ComputePkFkiM
        (
            const DataT*        DataP,      /* I */
            const ModelSpecT*   SpecP,      /* I */
            const ModelParaT*   ParaP,      /* I */
            double*             PkFkiM,     /* O */
            float*              LogPkFkiM   /* O */
        ) ;

    static int ComputePartition
        (
          const ModelSpecT*   SpecP,      /* I */
          const ModelParaT*   ParaP,      /* I */
          const DataT*        DataP,      /* I */
          const SpatialT*     SpatialP,   /* I */
          const NemParaT*     NemParaP,   /* I */
          FILE*               Flog,       /* I/O */
          float*              CM,         /* I/O */
          WorkingT*           WorkP,      /* I:PkFkiM O:CtmpM,Neighs,CiNumV */
          CriterT*            CriterP     /* O */
        ) ;

    static void WriteLogClasses
        (
            FILE*               Flog,       /* I/O */
            int                 Nd,         /* I : dimension */
            const ModelSpecT*   SpecP,      /* I */
            const ModelParaT*   ParaP       /* I */
        ) ;

#endif


static int  /* 0 = OK,  1 = beta diverges to infinity */
EstimBeta
(
  const BetaET       BetaModel,  /* I : beta estimation mode */
  const BtaPsGradT*  BtaPsGradP, /* I : gradient ascent parameters */
  const SpatialT*    SpatialP,   /* I : neighborhood system */
  const float*       C_NK,       /* I : classification matrix (Npt, NbK)*/
  const int          Npt,        /* I : number of points */
  const int          Nbk,        /* I : number of classes */
  float*             BetaP,      /* I/O : previous and new beta */
  PtNeighsT*         NeighsP     /* W : to store a site's neighbors */
) ;


    static StatusET ComputeCrit   /*V1.03-a*/
     (
        int                 Npt,         /* I */
        int                 Nk,          /* I */
        float               Beta,        /* I : spatial coefficient */
        const float*        CM,          /* I : class. matrix [Npt*Nk] */
        const SpatialT*     SpatialP,    /* I : neighborhood system */
        WorkingT*           WorkP,       /* I :PkFkiM, O:Neighs */
        CriterT*            CriterP      /* O : computed criteria */ 
     ) ;

    static void CalcError                /*V1.04-e*/
      ( 
       const float*     CM,       /* I : found classification */
       const int        N,        /* I : number of points */
       const int        Harden,   /* I : 1 use hardened CM, 0 use CM itself */
       const ErrinfoT*  ErrinfoP, /* I : info to compute error */
       ErrcurT*         ErrcurP   /* O : computed error and agreement */
      ) ;


/* Called by ComputePartition */

  static int ComputePartitionGEM
      (
          const ModelSpecT*   SpecP,      /* I */
          const ModelParaT*   ParaP,      /* I */
          const DataT*        DataP,      /* I */
          const SpatialT*     SpatialP,   /* I */
          const NemParaT*     NemParaP,   /* I */
          FILE*               Flog,       /* I/O */
          float*              CM,         /* I/O */
          WorkingT*           WorkP,      /* I:PkFkiM O:CtmpM,Neighs,CiNumV */
          CriterT*            CriterP     /* O */
      ) ;

  static int ComputePartitionNEM
      (
          const ModelSpecT*   SpecP,      /* I */
          const ModelParaT*   ParaP,      /* I */
          const DataT*        DataP,      /* I */
          const SpatialT*     SpatialP,   /* I */
          const NemParaT*     NemParaP,   /* I */
          FILE*               Flog,       /* I/O */
          float*              CM,         /* I/O */
          WorkingT*           WorkP,      /* I:PkFkiM O:CtmpM,Neighs,CiNumV */
          CriterT*            CriterP     /* O */
      ) ;


/* Called by ComputePartitionNEM */

    static int ComputeLocalProba
    (
     const int            Ipt,        /* I: current site */
     const int            Nk,         /* I: nb of classes */
     const ModelParaT*    ParaP,      /* I: parameters (use Beta) */
     const NeighDataT*    NeighDataP, /* I: neighborhood system */
     GetNeighFT*          FGetNeigh,  /* I: function to fetch neighbors */
     const double*        PkfkiM,     /* I: pk fk(xi) */
     const float*         Cin_NK,     /* I: class of neighbors */
     float*               Cout_K,     /* O: class of site Ipt */
     PtNeighsT*           NeighsP,    /* W: to store fetched neighbors */
     double*              Cinum_K     /* W: to store cik's numerators */
    );

    static void WriteLogCrit
     (
        FILE*               Flog,        /* I/O */
        const int           Npt,         /* I */
        const int           Nk,          /* I */
        const float         Beta,        /* I : spatial coefficient */
        const float*        CM,          /* I : class. matrix [Npt*Nk] */
        const SpatialT*     SpatialP,    /* I : neighborhood system */
        WorkingT*           WorkP,       /* I :PkFkiM, O:Neighs */
        CriterT*            CriterP      /* O : computed criteria */ 
     ) ;


/* Called by ComputePartitionGEM */

    static int Multinomial(int km, const float *tk) ;

#if 0

    static int ComputeLocalProba
    (
     const int            Ipt,        /* I: current site */
     const int            Nk,         /* I: nb of classes */
     const ModelParaT*    ParaP,      /* I: parameters (use Beta) */
     const NeighDataT*    NeighDataP, /* I: neighborhood system */
     GetNeighFT*          FGetNeigh,  /* I: function to fetch neighbors */
     const double*        PkfkiM,     /* I: pk fk(xi) */
     const float*         Cin_NK,     /* I: class of neighbors */
     float*               Cout_K,     /* O: class of site Ipt */
     PtNeighsT*           NeighsP,    /* W: to store fetched neighbors */
     double*              Cinum_K     /* W: to store cik's numerators */
    );

    static void WriteLogCrit
     (
        FILE*               Flog,        /* I/O */
        const int           Npt,         /* I */
        const int           Nk,          /* I */
        const float         Beta,        /* I : spatial coefficient */
        const float*        CM,          /* I : class. matrix [Npt*Nk] */
        const SpatialT*     SpatialP,    /* I : neighborhood system */
        WorkingT*           WorkP,       /* I :PkFkiM, O:Neighs */
        CriterT*            CriterP      /* O : computed criteria */ 
     ) ;


#endif


/* Called by ComputeLocalProba */

    static float SumNeighsOfClass
     (
        int             K ,     /* I : class to search, 0..Nk-1 */
        int             Nbn ,   /* I : number of neighbours */
        int             Nk,     /* I : number of classes */
        const NeighT*   NeiV ,  /* I : neighbours, [Nbn] */
        const float*    CM      /* I : classification matrix [Npt*Nk] */
     ) ;


/* Called by WriteLogCrit */

    static StatusET ComputeCrit
     (
        int                 Npt,         /* I */
        int                 Nk,          /* I */
        float               Beta,        /* I : spatial coefficient */
        const float*        CM,          /* I : class. matrix [Npt*Nk] */
        const SpatialT*     SpatialP,    /* I : neighborhood system */
        WorkingT*           WorkP,       /* I :PkFkiM, O:Neighs */
        CriterT*            CriterP      /* O : computed criteria */ 
     ) ;


/* Called by ComputeCrit */

#if 0
    static float SumNeighsOfClass
     (
        int             K ,     /* I : class to search, 0..Nk-1 */
        int             Nbn ,   /* I : number of neighbours */
        int             Nk,     /* I : number of classes */
        const NeighT*   NeiV ,  /* I : neighbours, [Nbn] */
        const float*    CM      /* I : classification matrix [Npt*Nk] */
     ) ;

static void CalcError                /*V1.04-e*/
      ( 
       const float*     CM,       /* I : found classification */
       const int        N,        /* I : number of points */
       const int        Harden,   /* I : 1 use hardened CM, 0 use CM itself */
       const ErrinfoT*  ErrinfoP, /* I : info to compute error */
       ErrcurT*         ErrcurP   /* O : computed error and agreement */
      ) ;

#endif


/* Called by  */



/* ==================== GLOBAL FUNCTION DEFINITION =================== */

/* ------------------------------------------------------------------- */
int ClassifyByNem                                             /*V1.04-a*/
        ( 
          const NemParaT      *NemParaP,        /* I */
          const SpatialT      *SpatialP,        /* I */
          DataT               *DataP,           /* I/O */
          StatModelT          *StatModelP,      /* I/O */
          float               *ClassifM,        /* I/O */
	  CriterT             *CriterP          /* O */  /*V1.03-f*/
        ) 
/* ------------------------------------------------------------------- */
{
  StatusET   sts = STS_OK ;


  /*V1.06-f*/
  ModelPreprocess( & StatModelP->Spec, DataP ) ; 

  switch( StatModelP->Spec.BetaModel )
    {
    case BETA_FIX:
    case BETA_PSGRAD:
      sts = ClassifyByNemOneBeta( DataP, NemParaP, SpatialP, 
				    StatModelP, ClassifM, CriterP ) ;
      break;

    case BETA_HEUD:
    case BETA_HEUL:
      sts = ClassifyByNemHeuBeta( DataP, (NemParaT*) NemParaP, SpatialP, 
				  StatModelP, ClassifM, CriterP ) ;
      break ;

    default:
      fprintf( out_stderr , "Error : unkwnown beta estimation mode %d\n", 
	       StatModelP->Spec.BetaModel ) ;
      sts = STS_E_FUNCARG ;
    }

  return sts ;
}




/* ------------------------------------------------------------------- */
int ComputeMAP                     /* ret : MAP label in 0..Nk-1 */
        ( 
          const float* ClassifM,   /* I */
          int          Ipt,        /* I */
          int          Nk,         /* I */
	  TieET        TieRule,    /* I */
          int*         kmaxesV     /* I/O [Nk] */
        ) 
/* ------------------------------------------------------------------- */
{
    int   k ;
    int   nequal ;
    int   kmax = 0;
    float ukmax = ClassifM[ Ipt * Nk + 0 ] ;


    /* ... First pass : find out the first maximum uik */
    for ( k = 1 ; k < Nk ; k ++ ) {

        float uik = ClassifM[ (Ipt * Nk) + k ] ;        
        if ( uik > ukmax )
          {
            ukmax = uik ;
            kmax = k ;
          }
    } /* end for ( k = 1 ... ) */

    if ( TieRule == TIE_RANDOM ) {
      /* ... Second pass : check if there are other k / uik = ukmax */
      kmaxesV[ 0 ] = kmax ;
      for ( k = kmax + 1, nequal = 0 ; k < Nk ; k ++ ) {

        float uik = ClassifM[ (Ipt * Nk) + k ] ;        
        if ( uik == ukmax )
          {
            nequal ++ ;
            kmaxesV[ nequal ] = k ;
          }
      }

      /* ... If more than one k / uik = ukmax, choose randomly one of them */
      if ( nequal > 0 ) {
        /* Choose randomly one of the k */
        int  imax ;   /* index of chosen k : 0 .. nequal */
        
        imax = RandomInteger( 0, nequal ) ;
        kmax = kmaxesV[ imax ] ;
      }
      else {
        kmax = kmaxesV[ 0 ] ;
      }
    } /* else TieRule == TIE_FIRST : keep first chosen maximum */

    return kmax ;

}  /* end of ComputeMAP() */


/* ------------------------------------------------------------------- */
void LabelToClassVector
 ( 
  const int Nk,    /* I: number of classes */
  const int Label, /* I: class of interest 0..Nk-1 */
  float* Cout_K    /* O: classification vector [Nk] */
 )
/* ------------------------------------------------------------------- */
{
  int k;

  for ( k = 0 ; k < Nk ; k ++ ) 
    Cout_K[ k ] = 0.0 ;
  
  if ( ( 0 <= Label ) && ( Label < Nk ) )
    Cout_K[ Label ] = 1.0 ;
}



/* ==================== LOCAL FUNCTION DEFINITION =================== */


/* ------------------------------------------------------------------- */
static void ModelPreprocess
        ( 
	  const ModelSpecT*   SpecP,            /* I */ 
	  DataT*              DataP             /* I/O */
	)
/* ------------------------------------------------------------------- */
{
  const char* func = "ModelPreprocess" ;

  SortPtT*    tabptV ;
  int         npt = DataP->NbPts ;
  int         nd  = DataP->NbVars ;
  int         d ;
  int         ipt ;


  switch( SpecP->ClassFamily ) 
    {
    case FAMILY_LAPLACE:  	/* Compute sorted index */
    case FAMILY_BERNOULLI:  	/* Compute sorted index */
      {	
	/* Allocate globally the sorted index */
	DataP->SortPos_ND = GenAlloc( npt * nd, sizeof( int ), 
				      1, func, "SortPos_ND" ) ;

	/* Allocate the local sorting table to be sorted */
	tabptV = GenAlloc( npt, sizeof( SortPtT ), 1, func, "tabptV" ) ;

	/* For each variable */
	for( d = 0 ; d < nd ; d ++ )
	  {
	    /* Fill the local sorting table with data's selected variable */
	    for ( ipt = 0 ; ipt < npt ; ipt ++ )
	      {
		tabptV[ ipt ].Index = ipt ;
		tabptV[ ipt ].Value = DataP->PointsM[ ( ipt * nd ) + d ] ;
	      }

	    /* Sort table by ascending variable value */
	    qsort( (void *) tabptV, npt, sizeof( SortPtT ), CompSortValue ) ;

	    /* Save the sorted indices */
	    for ( ipt = 0 ; ipt < npt ; ipt ++ )
	      DataP->SortPos_ND[ ( ipt * nd ) + d ] = tabptV[ ipt ].Index ;
	  }
      
	/* Free the local sorting table */
	GenFree( tabptV ) ;
      }
      break ;

    default:
      DataP->SortPos_ND = NULL ;
    }
}



/* ------------------------------------------------------------------- */
static int ClassifyByNemHeuBeta                             /*V1.04-a*/
        ( 
          const DataT         *DataP,           /* I */
          NemParaT            *NemParaP,        /* I[M] */
          const SpatialT      *SpatialP,        /* I */
          StatModelT          *StatModelP,      /* I/O */
          float               *ClassifM,        /* I/O */
	  CriterT             *CriterP          /* O */  /*V1.03-f*/
        ) 
/* ------------------------------------------------------------------- */
{
  const char* func = "ClassifyByNemHeuBeta" ;

  int    npt = DataP->NbPts ;
  int    nk = StatModelP->Spec.K ;
  int    nbtamax ;  /* number of tested beta values */
  int    nbta ;     /* number of valid beta values */
  float* btaV ;     /* valid beta values [1..nbtamax] */
  float* criV ;     /* criteria value for each valid beta [1..nbtamax] */
  float* clasiniM ; /* initial classification matrix */
  float* clasbestM ;/* "best" obtained classification matrix */
  float  btatest ;  /* currently tested beta */
  float  Dmin ;     /* minimum criterion found so far */
  int    Dincreas ;  /* TRUE if increase of criterion already occurred */
  int    Ddrop ; /* TRUE if drop of criterion slope detected */
  float  prevSlope ;/* previously computed criterion slope */
  float  thisSlope ;/* currently computed criterion slope */
  float  DdropThres ;/* threshold of criterion slope change to detect */
  float  LlossThres ;
  float  Lmax ;
  int    Lfound ;
  int    stop ;
  float  btaEst ;   /* estimated beta */
  char   namecri[ LEN_FILENAME + 1 ] ; /* name of file to save criteria */
  FILE*  fcri ;     /* pointer to file to save criteria */


  /* Initialize to NaN to shut gcc warnings (and catch errors) */
  thisSlope = mknan() ;
  prevSlope = mknan() ;
  Lmax      = mknan() ;
  btaEst    = mknan() ;

  nbtamax = (int) ( NemParaP->BtaHeuMax / NemParaP->BtaHeuStep ) + 1 ;

  btaV = GenAlloc( nbtamax + 1 , sizeof( float ), 0, func, "btaV" ) ;
  criV = GenAlloc( nbtamax + 1 , sizeof( float ), 0, func, "criV"  ) ;
  clasiniM = GenAlloc( npt * nk, sizeof( float ), 0, func, "clasiniM" ) ;
  clasbestM = GenAlloc( npt * nk, sizeof( float ), 0, func, "clasbestM" ) ;

  if ( ( btaV == NULL ) || ( criV == NULL ) || ( clasiniM == NULL ) ) 
    return STS_E_MEMORY ;

  /* Save initial classification matrix */
  memcpy( clasiniM, ClassifM, 
	  npt * nk * sizeof( float ) ) ;

  /* For tests */  /*V1.04-d*/
  if ( NemParaP->Debug )
    {
      strncpy( namecri, NemParaP->OutBaseName, LEN_FILENAME ) ;
      strncat( namecri, ".bta", LEN_FILENAME ) ;
      fcri = fopen( namecri , "w" ) ;
    }
  else
    fcri = NULL ;


  
  /* While bta less than btamax and stopping condition not met */

  Dincreas = FALSE;  /* Initialize values for Hathaway heuristic */
  Dmin = 0.0; 
  Ddrop = FALSE ;
  DdropThres = - NemParaP->BtaHeuDDrop * DataP->NbPts ;

  Lfound = FALSE ;   /* Initialize values for likelihood heuristic */
  LlossThres= NemParaP->BtaHeuLLoss * DataP->NbPts;

  fprintf( out_stderr, "* * Starting heuristic * *\n" ) ; if(out_stderr) fflush( out_stderr );
  fprintf( out_stderr, "    Parameters : " ) ; if(out_stderr) fflush( out_stderr );

  if ( StatModelP->Spec.BetaModel == BETA_HEUD )
    fprintf( out_stderr, "Drop < %3.1f  or  (Dmax - D) < %3.1f (Dmax - Dmin)\n" ,
	     DdropThres, NemParaP->BtaHeuDLoss ) ; if(out_stderr) fflush( out_stderr );

  else
    fprintf( out_stderr, "L < Lmax - %3.1f \n" , LlossThres ) ; if(out_stderr) fflush( out_stderr );


  for ( btatest = 0.0, nbta = 0, stop = FALSE ;
	( btatest <= NemParaP->BtaHeuMax ) && ( ! stop ) ;
	btatest += NemParaP->BtaHeuStep )
    {
      StatusET sts ;

      fprintf( out_stderr, "\n * * Testing beta = %5.2f * * \n" , btatest ) ; if(out_stderr) fflush( out_stderr );


      StatModelP->Para.Beta = btatest ;

      /* Reload initial classification matrix */
      memcpy( ClassifM, clasiniM, npt * nk * sizeof( float ) ) ;

      sts = ClassifyByNemOneBeta( DataP, NemParaP, SpatialP, 
				    StatModelP, ClassifM, CriterP ) ;
      if ( sts == STS_OK )
	{
	  nbta ++ ;

	  btaV[ nbta ] = btatest ;

	  if ( StatModelP->Spec.BetaModel == BETA_HEUD )
	    {
	      criV[ nbta ] = CriterP->D ;

	      if ( CriterP->D < Dmin )
		Dmin = CriterP->D ;

    	      if ( nbta >= 2 )
    		{
    		  prevSlope = thisSlope ;
    		  thisSlope = ( criV[ nbta ] - criV[ nbta - 1 ] ) /
    		    ( btaV[ nbta ] - btaV[ nbta - 1 ] ) ;
    		  
    		  if ( thisSlope >= 0.5 * DataP->NbPts ) /* positive slope */
    		    Dincreas = TRUE ;
    		}
    
    	      if ( nbta >= 3 )
    		{
    		  fprintf( out_stderr, "    * Drop : %5.1f (threshold %5.1f) *\n", 
    			       ( thisSlope - prevSlope ) , DdropThres ) ; if(out_stderr) fflush( out_stderr );

    
    		  if ( (! Ddrop) && (! Dincreas) ) {

		    if ( (thisSlope - prevSlope) < DdropThres ) {

    		      Ddrop = TRUE ;
		      if ( ! NemParaP->Debug )  stop = TRUE ;
    		      fprintf( out_stderr, "    * ---> Drop detected *\n" ) ; if(out_stderr) fflush( out_stderr );

    		      /* btaEst = btaV[ nbta - 2 ] ; V1.06-i*/
		      btaEst = btaV[ nbta - 1 ] ; 

    		    }
		    else /* no drop of slope detected, save this class. */ {

		      memcpy( clasbestM, ClassifM, 
			      npt * nk * sizeof( float ) ) ;

		    }

		  } /* Else drop or increase already detected */

    		}
	      else /* 1st or 2nd tested beta -> save classification */
		memcpy( clasbestM, ClassifM, 
			npt * nk * sizeof( float ) ) ;
		
	    }
	  else  /* else Likelihood heuristic */
	    {
	      criV[ nbta ] = CriterP->L ;
	      
    	      if ( nbta < 2 )  /* 1st beta */ {

		  Lmax = CriterP->L ;
		  fprintf( out_stderr, " * *  L threshold  :  %3.1f  * * \n" ,
			   Lmax - LlossThres ) ; if(out_stderr) fflush( out_stderr );

		  memcpy( clasbestM, ClassifM, npt * nk * sizeof( float ) ) ;
		}
	      else /* >= 2nd beta */ 
		{
		  if ( CriterP->L > Lmax )
		    {
		      Lmax = CriterP->L ;
		      fprintf( out_stderr, " * *  L threshold  :  %3.1f  * * \n" ,
			       Lmax - LlossThres ) ; if(out_stderr) fflush( out_stderr );

		    }

    		  if (! Lfound) {

		    if ( CriterP->L < Lmax - LlossThres ) {

    		      Lfound = TRUE ;
		      if ( ! NemParaP->Debug )  stop = TRUE ;
    		      fprintf( out_stderr, "    * ---> L loss detected *\n" ) ; if(out_stderr) fflush( out_stderr );

    		      /* btaEst = btaV[ nbta - 2 ] ; V1.06-i*/
		      btaEst = btaV[ nbta - 1 ] ; 
    		    }
		    else /* no drop yet detected -> save classification */ {
		      memcpy( clasbestM, ClassifM, 
			      npt * nk * sizeof( float ) ) ;
		    }
		  }
		}
	    }
	  /* end-else Likelihood heuristic */

	  if ( fcri != NULL )
	    fprintf( fcri, "%5.2f  %10.1f  %6.3f  %d\n", 
		     btatest , 
		     StatModelP->Spec.BetaModel == BETA_HEUD ? 
		     CriterP->D : CriterP->L,
		     CriterP->Errcur.Errorrate,
		     StatModelP->Spec.BetaModel == BETA_HEUD ? 
		     Ddrop : Lfound );
	}
      /* end-if ClassifyByNemOneBeta OK */
    }
  /* end-While bta less than btamax and stopping condition not met */

  if ( ( StatModelP->Spec.BetaModel == BETA_HEUD ) && ( ! Ddrop ) ) {
    int ibta ;
    int found ;
    float  DThres ; /* threshold of criterion loss to detect */
    
    DThres = criV[ 1 ] - ( criV[ 1 ] - Dmin ) * NemParaP->BtaHeuDLoss ;

    for ( ibta = 1, found = FALSE ; 
	  ( ibta <= nbta ) && ( ! found ) ; ibta ++ )
      found = criV[ ibta ] <= DThres ;

    if ( found )
      btaEst = btaV[ ibta - 2 ] ;
    else {

      fprintf( out_stderr , "Warning : heuristic failed to detect beta\n" ) ; if(out_stderr) fflush( out_stderr );

      btaEst = 0.0 ;
    }

    /* Reload initial classification matrix */
    memcpy( ClassifM, clasiniM, npt * nk * sizeof( float ) ) ;
  }
  else /* likelihood heuristic, or D heuristic and drop detected */ {

    /* restore saved "best" classification, and use INIT_FILE */ //*V1.06-i*/
    memcpy( ClassifM, clasbestM, npt * nk * sizeof( float ) ) ;
    NemParaP->InitMode = INIT_FILE ;
  }

  if ( ( StatModelP->Spec.BetaModel == BETA_HEUL ) && ( ! Lfound ) )
    {
      /* btaEst = btaV[ nbta - 1 ] ; V1.06-i*/
      btaEst = btaV[ nbta ] ;
      fprintf( out_stderr, "    * L loss not detected *\n" ) ; if(out_stderr) fflush( out_stderr );

    }

  if ( fcri != NULL )
    fclose( fcri ) ;

  fprintf( out_stderr , "\n * * *  Estimated beta : %3.2f * * *\n" , btaEst ) ; if(out_stderr) fflush( out_stderr );

  if ( StatModelP->Spec.BetaModel == BETA_HEUD )
    fprintf( out_stderr , " * * *   Using %s * * *\n\n" , 
	     Ddrop ? "drop detection" : "loss thresholding" ) ; if(out_stderr) fflush( out_stderr );



  GenFree( btaV ) ;     /*V1.05-h*/
  GenFree( criV ) ;
  GenFree( clasiniM ) ;
  GenFree( clasbestM ) ;
  
  StatModelP->Para.Beta = btaEst ;
  return ClassifyByNemOneBeta( DataP, NemParaP, SpatialP, 
				 StatModelP, ClassifM, CriterP ) ;
}



/* ------------------------------------------------------------------- */
static int ClassifyByNemOneBeta                             /*V1.04-a*/
        ( 
          const DataT         *DataP,           /* I */
          const NemParaT      *NemParaP,        /* I */
          const SpatialT      *SpatialP,        /* I */
          StatModelT          *StatModelP,      /* I/O */
          float               *ClassifM,        /* I/O */
	  CriterT             *CriterP          /* O */  /*V1.03-f*/
        ) 
/* ------------------------------------------------------------------- */
{
    const char* func        = "ClassifyByNemOneBeta" ;
    WorkingT    working     = {NULL};   /* working variables to allocate */
    FILE*       flog        = NULL ;
    StatusET    err         = STS_OK ;
    int         npt         = DataP->NbPts;
    int         nk          = StatModelP->Spec.K ;
    int         k ;
    int         missk ;  /* 0..K-1 : class with missing data (-1 if none) */
    int         missd ;  /* 0..D-1 : variable with missing data (-1 if none) */


    /* Allocate all working variables */
    {
        working.KmaxesV         = GenAlloc( nk, sizeof( int ), 0, func, "KmaxesV" ) ;
        working.CtmpM           = GenAlloc( npt * nk, sizeof( float ), 0, func, "CtmpM" ) ;
        working.ColdM           = GenAlloc( npt * nk, sizeof( float ), 0, func, "ColdM" ) ;
        working.CiNumV          = GenAlloc( nk, sizeof( double ), 0, func, "CiNumV" ) ;
        working.PkFkiM          = GenAlloc( npt * nk, sizeof( double ), 0, func, "PkFkiM" ) ;
        working.LogPkFkiM       = GenAlloc( npt * nk, sizeof( float ), 0, func, "LogPkFkiM" ) ;
        working.Neighs.NbNeigh  = SpatialP->MaxNeighs ;
        working.Neighs.NeighsV  = GenAlloc( SpatialP->MaxNeighs, sizeof( NeighT ), 
					    0, func, "NeighsV" ) ;
    }

#ifdef __TURBOC__
    fprintf( out_stderr, "\n Remaining memory after allocating : %lu bytes\n", 
             (unsigned long) coreleft() ); if(out_stderr) fflush( out_stderr );

#endif

    if ( ( working.KmaxesV == NULL ) || 
	 ( working.CtmpM == NULL ) || 
	 ( working.CiNumV == NULL ) || 
         ( working.ColdM == NULL ) || 
	 ( working.PkFkiM == NULL ) ||
	 ( working.LogPkFkiM == NULL ) ||
         ( ( working.Neighs.NeighsV == NULL ) && 
           ( SpatialP->MaxNeighs > 0 ) ) )
    {
        fprintf( out_stderr, "Could not allocate NEM working variables\n" ); if(out_stderr) fflush( out_stderr );

        return STS_E_MEMORY ;
    }

    if ( NemParaP->DoLog )
    {
        StartLogFile( NemParaP->LogName, DataP->NbPts, &flog ) ;
    }
    else flog = NULL ;


    switch( NemParaP->InitMode )
    {
    case INIT_SORT:    /* initial partition to be computed by sorting */
        fprintf( out_stderr, 
                 "Computing initial partition (sort variable %d) ...\n",
                  NemParaP->SortedVar + 1 ) ; if(out_stderr) fflush( out_stderr );


        InitPara( DataP, & StatModelP->Desc, & StatModelP->Spec, 
		  & StatModelP->Para, 
		  working.CtmpM ) ;

        if ( ( err = InitPartitionSort( DataP, nk, 
                                        NemParaP->SortedVar, 
                                        ClassifM ) ) != STS_OK )
	  return err ;
	
	/*V1.06-e*/
	err = MakeParaFromLabeled( DataP, ClassifM, & StatModelP->Spec, 
				   & StatModelP->Desc,
				   & StatModelP->Para, & missk, & missd ) ; 
        if ( err != STS_OK )
	  return err ;

        if ( flog != NULL )
            fprintf( flog, "Initialization by sorting variable %d:\n",
                     NemParaP->SortedVar + 1 ) ; if(flog) fflush( flog );


        err = NemAlgo( DataP, NemParaP, SpatialP, & StatModelP->Spec, 
		       & StatModelP->Desc,
		       flog, & StatModelP->Para, ClassifM, 
		       & working, CriterP ) ;
        break ;


    case INIT_FILE:    /* initial classification is given by a file */
        if ( flog != NULL )
            fprintf( flog, "Initialization with specified partition :\n" ) ; if(flog) fflush( flog );


	/*V1.06-g*/
        InitPara( DataP, & StatModelP->Desc, & StatModelP->Spec, 
		  & StatModelP->Para, 
		  working.CtmpM ) ;

	err = MakeParaFromLabeled( DataP, ClassifM, & StatModelP->Spec, 
				   & StatModelP->Desc,
				   & StatModelP->Para, & missk, & missd ) ; 
        if ( err != STS_OK )
	  return err ;

        err = NemAlgo( DataP, NemParaP, SpatialP, & StatModelP->Spec, 
		       & StatModelP->Desc,
		       flog, & StatModelP->Para, ClassifM, 
		       & working, CriterP ) ;
        break ;


    case INIT_LABEL:
        fprintf( out_stderr, 
		 "Initializing centers from partially labeled sample\n"); if(out_stderr) fflush( out_stderr );


        if ( flog != NULL ) {
	  fprintf( flog, "Initialization with partially labeled sample :\n" ) ; if(flog) fflush( flog );

	  fprintf( flog, "%4d ", 0 ) ; if(flog) fflush( out_stderr );

	}

	/*V1.03-d*/
        InitPara( DataP, & StatModelP->Desc, & StatModelP->Spec, 
		  & StatModelP->Para, 
		  working.CtmpM ) ;

        /* Compute initial parameters using observations with known labels */
	err = MakeParaFromLabeled( DataP, ClassifM, & StatModelP->Spec, 
				   & StatModelP->Desc,
				   & StatModelP->Para, & missk, & missd ) ; 
        if ( ( err != STS_OK ) && ( err != STS_W_EMPTYCLASS ) )
	  return err ;

	/* Proportions are set equal */
	for ( k = 0 ; k < nk ; k ++ )
	  StatModelP->Para.Prop_K[ k ] = 1.0 / nk ; /*V1.05-b*/

	ComputePartitionFromPara( 1, DataP, NemParaP, & StatModelP->Spec, 
				  & StatModelP->Para, SpatialP, 
				  ClassifM, CriterP, flog, & working ) ;

	/* Run NEM with this initial partition */
        err = NemAlgo( DataP, NemParaP, SpatialP, & StatModelP->Spec, 
		       & StatModelP->Desc,
		       flog, & StatModelP->Para, ClassifM, 
		       & working, CriterP ) ;
        break ;


    case INIT_PARAM_FILE:
        fprintf( out_stderr, 
		 "Initializing parameters from given value\n");
        if ( flog != NULL ) {
	  fprintf( flog, "Initializing parameters from given value :\n" ) ; if(out_stderr) fflush( flog );

	  fprintf( flog, "%4d ", 0 ) ; if(flog) fflush( flog );

	}

	ComputePartitionFromPara( 1, DataP, NemParaP, & StatModelP->Spec, 
				  & StatModelP->Para, SpatialP, 
				  ClassifM, CriterP, flog, & working ) ;

	/* Run NEM with this initial partition */
        err = NemAlgo( DataP, NemParaP, SpatialP, & StatModelP->Spec, 
		       & StatModelP->Desc,
		       flog, & StatModelP->Para, ClassifM, 
		       & working, CriterP ) ;
	break ;

    default:   /* Else randomly generate initial classifications */
        err = RandNemAlgo ( DataP, NemParaP, SpatialP, & StatModelP->Spec, 
			    & StatModelP->Desc,
			    flog, & StatModelP->Para, ClassifM, 
			    & working, CriterP ) ;
	  
    } /* end switch( NemParaP->InitMode ) */

    /* Free allocated structures */
    GenFree( working.KmaxesV ) ;
    GenFree( working.CtmpM ) ;
    GenFree( working.CiNumV ) ;
    GenFree( working.PkFkiM ) ;
    GenFree( working.LogPkFkiM ) ;
    GenFree( working.Neighs.NeighsV ) ;
    GenFree( working.ColdM ) ;

    if ( flog != NULL ) fclose( flog ) ;

    return err ;

}  /* end of ClassifyByNem() */






/* ------------------------------------------------------------------- */
static void InitPara
        (
         const DataT*       DataP,       /* I */
	 const SampleDesT*  DescP,       /* I */
	 const ModelSpecT*  SpecP,       /* I */
	 ModelParaT*        ParaP,       /* O */
         float*             C_NK         /* T */
        ) 
/*\

    This function initializes to NaN all class parameters, computes
    the volume of the whole sample, and minimum and maximum within
    each variable.

\*/
/* ------------------------------------------------------------------- */
{
  int             npt         = DataP->NbPts;
  int             nd          = DataP->NbVars ;
  int             nk          = SpecP->K ;
  int             ipt ;
  int             d, k ;
  int             emptyk ;  /* index of empty class, 0 if none 
			       (not used here)*/


  /* Compute minimum and maximum within each variable */
  /*V1.05-d*/
  for ( d = 0 ; d < nd ; d ++ )
    {
      DescP->MiniSam_D[ d ] = MAXFLOAT ;
      DescP->MaxiSam_D[ d ] = -MAXFLOAT ;
      for ( ipt = 0 ; ipt < npt ; ipt ++ )
	{
	  float x = DataP->PointsM[ ipt * nd + d ] ;
	  
	  if ( ! isnan( x ) )
	    {
	      if ( x < DescP->MiniSam_D[ d ] )
		DescP->MiniSam_D[ d ] = x ;
	      
	      if ( x > DescP->MaxiSam_D[ d ] ) 
		DescP->MaxiSam_D[ d ] = x ;
	    }
	}
    }
  
  /* Compute dispersion of whole sample 
     by assigning all data to first class */
  for ( ipt = 0 ; ipt < npt ; ipt ++ )
    {
      C_NK[ ipt * nk + 0 ] = 1.0 ;
      for ( k = 1 ; k < nk ; k ++ )
	{
	  C_NK[ ipt * nk + k ] = 0.0 ;
	}
    }

  EstimPara( C_NK, DataP, nk, MISSING_IGNORE, SpecP,
	     &emptyk, ParaP ) ;
  
  for ( d = 0 ; d < nd ; d ++ )
    {
      DescP->DispSam_D[ d ] = 
	ParaP->Disp_KD[ 0 * nd + d ] ;
    }

  /* Set to nan all class parameters */
  for ( k = 0 ; k < nk ; k ++ )
    {
      ParaP->Prop_K[ k ]    = mknan() ;
      ParaP->NbObs_K [ k ]  = mknan() ;
      for ( d = 0 ; d < nd ; d ++ )
	{
	  ParaP->Center_KD[ k * nd + d ] = mknan() ;
	  ParaP->Disp_KD[ k * nd + d ]   = mknan() ;
	  ParaP->NbObs_KD[ k * nd + d ]  = mknan() ;
	}
    }

}



/* ------------------------------------------------------------------- */

static float ChosenCrit( const CriterT* CriterP , CritET Which )  /*V1.04-f*/
{
  switch( Which )
    {
    case CRIT_U: return CriterP->U ;
    case CRIT_M: return CriterP->M ;  /*V1.05-k*/
    case CRIT_D: return CriterP->D ;
    case CRIT_L: return CriterP->L ;
    default:     return CriterP->D ;
    }
}


/* ------------------------------------------------------------------- */
/*V1.03-d*/
static int MakeParaFromLabeled
        ( 
            const DataT*        DataP,      /* I */
	          const float*        C_NK,       /* I */
	          const ModelSpecT*   SpecP,      /* I */
            const SampleDesT*   DescP,      /* I */
	          ModelParaT*         ParaP,      /* O */
	          int                 *misskP,    /* O */
	          int                 *missdP     /* O */
        ) 
/*\

     This function initializes the parameters given a classification.
     If a class and variable has no observation, the mean is drawn randomly.
     If a class and variable has < 3 observation, 
     the dispersion is set to ( whole sample dispersion / nk ).

     In the latter case, the function returns ( k * D + d ), where k and
     d are the indices of the empty class / variable (0..K-1, 0..D-1).
     If there is no empty class / variable, -1 is returned.

\*/
/* ------------------------------------------------------------------- */
{
    int             nd          = DataP->NbVars ;
    int             nk          = SpecP->K ;
    int             sts         = STS_OK ;
    int             d, k ;
    int             emptyk ; /* index of empty class, 0 if none */

    /*V1.06-d*/

   *misskP = -1 ;
   *missdP = -1 ;

    /* Estimate parameters using known labels */
    sts = EstimPara( C_NK, DataP, nk, MISSING_IGNORE, SpecP,
		     &emptyk, ParaP ) ;

    /* If some class is empty, produce an error message */
    if ( sts != STS_OK )
      {
	if ( sts == STS_W_EMPTYCLASS )
	  fprintf( out_stderr, "Class %d has no labeled observation\n", emptyk ) ; if(out_stderr) fflush( out_stderr );

	return sts ;
      }

    /* Check mean and dispersion of each class and variable */
    for ( k = 0 ; k < nk ; k ++ )
      for ( d = 0 ; d < nd ; d ++ )
	{
	  /* If there is no observation, draw mean randomly in min-max range */
	  if ( ParaP->NbObs_KD[ k * nd + d ] < EPSILON )
	    {
	      fprintf( out_stderr , 
		       "Warning: no data in class k=%d, variable=%d\n" ,
		       k + 1 , d + 1 ); if(out_stderr) fflush( out_stderr );

	      *misskP = k ;
	      *missdP = d ;

	      ParaP->Center_KD[ k * nd + d ] = 
		RandomFloat( DescP->MiniSam_D[ d ], 
			     DescP->MaxiSam_D[ d ] ) ;
	    }
	  /* (Else, this mean could be computed from labeled data) */

	  /* If this class/variable had less than 3 observations, 
	     set its dispersion to ( whole sample dispersion / nk ) */
	  if ( ParaP->NbObs_KD[ k * nd + d ] < 3 )
	    {
	      ParaP->Disp_KD[ k * nd + d ] = DescP->DispSam_D[ d ] / nk ;
	    }
	}

    return sts ;

}  /* end of MakeParaFromLabeled() */


/* ------------------------------------------------------------------- */
static StatusET  MakeRandomPara
        ( 
            const DataT*        DataP,      /* I */
	    const ModelSpecT*   SpecP,      /* I */
            const SampleDesT*   DescP,      /* I */
	    ModelParaT*         ParaP       /* O */
        ) 
/* ------------------------------------------------------------------- */
{
    int             npt         = DataP->NbPts;
    int             nd          = DataP->NbVars ;
    int             nk          = SpecP->K ;
    int             k ;
    int             d ;


    /* For each class and variable */
    for ( k = 0 ; k < nk ; k ++ )
      for ( d = 0 ; d < nd ; d ++ )
	{
	  /* set its dispersion to ( whole sample dispersion / nk ) 
	   */
	  ParaP->Disp_KD[ k * nd + d ] = DescP->DispSam_D[ d ] / nk ;
	}

    /* Proportions are set equal */
    for ( k = 0 ; k < nk ; k ++ )
      ParaP->Prop_K[ k ] = 1.0 / nk ; /*V1.03-e*/


    /* For each class */ 
    for ( k = 0 ; k < nk ; k ++ )
    {
        int     ipt = 0 ;  /* = 0 to shut gcc warning */
	int     again ;
	int     ndraw ;

	/* Its mean is picked randomly from the data */

	/*V1.05-m*/
	/* Draw an object again while identical to a previously 
	   drawn center */
	for ( ndraw = 0, again = 1 ; again && ( ndraw < 100 ) ; ndraw ++ )
	  {
	    int h ;

	    ipt = RandomInteger( 0 , npt - 1 ) ;

	    /* Compare with a previous class center until one is identical */
	    for ( h = 0, again = 0 ; ( h < k ) && ( ! again ) ; h ++ )
	      {
		int different = 0 ; /* assume no difference until a different
				       coordinate is found */
		for ( d = 0 ; d < nd ; d ++ )
		  {
		    if ( ! isnan( DataP->PointsM[ ipt * nd + d ] )  )
		      {
			if ( ParaP->Center_KD[ h*nd + d ] != 
			     DataP->PointsM[ ipt*nd + d ] )
			  different = 1 ;
		      }
		    else
		      different = 1 ;  /* if this point has a missing
					  coordinate, assume it's different */
		  }

		if ( ! different )
		  again = 1 ;
	      }
	  }


	/* Get mean coordinates from drawn object non-NaN coordinates ;
	   or by uniform sampling if drawn object's coordinate is NaN */
        for ( d = 0 ; d < nd ; d ++ )
	  {
	    /* Use drawn point's dth coordinate if not NaN */
	    if( ! isnan( DataP->PointsM[ ipt * nd + d ] )  )  /*V1.05-d*/
	      {
		ParaP->Center_KD[ k*nd + d ] = DataP->PointsM[ ipt*nd + d ] ;
	      }
	    else  /* else draw randomly in min-max range of variable d */
	      {
		ParaP->Center_KD[ k*nd + d ] = 
		  RandomFloat( DescP->MiniSam_D[ d ], 
			       DescP->MaxiSam_D[ d ] ) ;
	      }
	  }
    }

    return STS_OK ;

}  /* end of MakeRandomPara() */



/* ------------------------------------------------------------------- */
static void StartLogFile( const char* LogName, int Npt, FILE** FlogP ) 
/* ------------------------------------------------------------------- */
{
  float mult; 

        if ( ( (*FlogP) = fopen( LogName, "w" ) ) == NULL )
        { 
            // FIX: setbuf requires a FILE*, not a FILE**
            setbuf(*FlogP, NULL);
            fprintf( out_stderr, "Could not open file '%s' in write mode\n", 
                     LogName ) ; if(out_stderr) fflush( out_stderr );

        }
        else
        {
            time_t timer = time( NULL ) ;

            fprintf( (*FlogP), "NEM log file  -  %s\n",
                     asctime( localtime( &timer ) ) ) ; if(*(FlogP)) fflush( *(FlogP) );


	    mult = exp( -((int) ( log(Npt/1000.)/log(10) )) * log( 10 ) ) ;
	    fprintf( (*FlogP), "  Criteria are multiplied by %f\n\n", mult ) ; if(*(FlogP)) fflush( *(FlogP) );

        }
}  /* end of StartLogFile() */


/* ------------------------------------------------------------------- */
static int  InitPartitionSort
        ( 
	 const DataT*    DataP,      /* I */
	 int             Nk,         /* I */
	 int             SortedVar,  /* I : sorted variable : 0..dim-1 */
	 float           *ClassifM   /* O (to allocate before call) */
	) 
/* ------------------------------------------------------------------- */
{
    int         ipt ;
    SortPtT*    tabptV ;
    int         npt = DataP->NbPts ;
    int         nd  = DataP->NbVars ;


    /* Allocate the local sorting table to be sorted */
    if ( ( tabptV = GenAlloc( npt, sizeof( SortPtT ),
			      0, "InitPartitionSort", "tabptV" ) ) == NULL )
      return STS_E_MEMORY ;


    /* Fill the local sorting table with data's selected variable */
    for ( ipt = 0 ; ipt < npt ; ipt ++ )
    {
        tabptV[ ipt ].Index = ipt ;
        tabptV[ ipt ].Value = DataP->PointsM[ ( ipt * nd ) + SortedVar ] ;
    }

    /* Sort table by ascending variable value */
    qsort( (void *) tabptV, npt, sizeof( SortPtT ), CompSortValue ) ;

    /* Use table to partition the data into Nk clusters of equal size */
    for ( ipt = 0 ; ipt < npt ; ipt ++ )
    {
        int index   = tabptV[ ipt ].Index ; /* true index of point */
        int ki      = ( ipt * Nk ) / npt ;  /* class to affect to point */
        int otherk ;                        /* other classes */

        for ( otherk = 0 ; otherk < Nk ; otherk ++ )
        {
            ClassifM[ (index * Nk) + otherk ] = 0.0 ;
        }

        ClassifM[ (index * Nk) + ki ] = 1.0 ;
    }

    /* Free the local sorting table */
    GenFree( tabptV ) ;

    return STS_OK ;

}   /* end of InitPartitionSort() */



static int CompSortValue( const void* elt1P, const void* elt2P )
{
    const SortPtT*  p1 = elt1P ;
    const SortPtT*  p2 = elt2P ;

    /*V1.06-h*/
    if ( isnan( p1->Value ) )           return 1 ;  /* 1 ">" 2 NaN at end */
    else if ( isnan( p2->Value ) )      return -1 ; /* 2 ">" 1 NaN at end */
    else if ( p1->Value < p2->Value )   return -1 ;
    else if ( p1->Value > p2->Value )   return 1 ;
    else                                return 0 ;
}




/* ------------------------------------------------------------------- */
static StatusET RandNemAlgo
        (
         const DataT*       DataP,           /* I */
         const NemParaT*    NemParaP,        /* I */
         const SpatialT*    SpatialP,        /* I */
	 const ModelSpecT*  SpecP,     	     /* I */
	 const SampleDesT*  DescP,     	     /* I */
         FILE*              Flog,            /* I/O */
	 ModelParaT*        ParaP,     	     /* I/O */
         float*             ClassifM,        /* I/O */
         WorkingT*          WorkP,           /* I/O */
         CriterT*           CriterP          /* O */
        ) 
/* ------------------------------------------------------------------- */
{
  const char* func      = "RandNemAlgo" ;
  StatusET  err         = STS_OK ;
  int       npt         = DataP->NbPts;
  int       nd          = DataP->NbVars;
  int       nk          = SpecP->K ;

  int       irandom ;
  int       nbsucc ;
  int       bestTry = -1 ;  /* stays -1 if no success 
			       (to shut gcc warning) */
  CriterT   bestCritS ;/*V1.04-g*/
  float*    bestClM ;  /* to allocate */
  float*    bestCen_KD ;  /* to allocate */ /*V1.06-q*/
  float*    bestDis_KD ;  /* to allocate */
  float     initbeta ; /* specified initial value of beta */
  float     bestbeta ; /* beta of best classification */


  InitPara( DataP, DescP, SpecP, 
	    ParaP, WorkP->CtmpM ) ;

  if ( ( bestClM = GenAlloc( npt * nk, sizeof(float), 
			     0, func, "bestClM" ) ) == NULL )
    return STS_E_MEMORY ;

  if ( ( bestCen_KD = GenAlloc( nk * nd, sizeof(float), 
			     0, func, "bestCen_KD" ) ) == NULL )
    return STS_E_MEMORY ;

  if ( ( bestDis_KD = GenAlloc( nk * nd, sizeof(float), 
			     0, func, "bestDis_KD" ) ) == NULL )
    return STS_E_MEMORY ;

  initbeta = ParaP->Beta ;

  /* For each random start */
  for ( irandom = 0, nbsucc = 0 ; 
      irandom < NemParaP->NbRandomInits ; 
      irandom ++ )
  {
    fprintf( out_stderr, 
             "\nRandom initial partition  %d\n", irandom + 1 ) ; if(out_stderr) fflush( out_stderr );


    if ( Flog != NULL )
    {
       fprintf( Flog, "\nRandom initialization %d :\n", irandom + 1 ) ; if(out_stderr) fflush( out_stderr );

       fprintf( Flog, "%4d ", 0 ) ; if(out_stderr) fflush( out_stderr );

    }

    /* Generate randomly initial parameters */
    err = MakeRandomPara( DataP , SpecP, DescP,
			  ParaP );
    if ( err != STS_OK )
	return err ;

    if (SpecP->BetaModel == BETA_PSGRAD ) {
      if ( NemParaP->BtaPsGrad.RandInit ) /* random initial beta */
	ParaP->Beta = RandomFloat( 0.0, 2.0 );
      else   /* Specified initial beta */
	ParaP->Beta = initbeta ;
    }

    /* Initialize classification with zeros to remove eventual NaN/Inf from
       previous runs */
    { 
      int ipt, k ;
      for ( ipt = 0 ; ipt < npt ; ipt ++ ) {
	for ( k = 1 ; k < nk ; k ++ ) {
	  ClassifM[ ipt * nk + k ] = 0.0 ;
	}
      }
    }

    ComputePartitionFromPara( 1, DataP, NemParaP, SpecP, 
			      ParaP, SpatialP, 
			      ClassifM, CriterP, Flog, WorkP ) ;

    /* Run NEM with this initial partition */
    err = NemAlgo( DataP, NemParaP, SpatialP, SpecP, DescP,
		   Flog, ParaP, ClassifM, 
		   WorkP, CriterP ) ;

    /* If NEM was successful */
    if ( err == STS_OK )
    {
        /* Count this success */
        nbsucc += 1 ;

        /* If it is the first success */
        if ( nbsucc == 1 )
        {
            /* Memorize the results as the currently best ones */
            memcpy( bestClM, ClassifM, npt * nk * sizeof(float) ) ;
	    bestCritS = *CriterP ;
            bestTry   = irandom ;
	    bestbeta  = ParaP->Beta ;
        }
        else    /* Else it is the 2nd or more success */     {
	  /* If the criterion is better than the current one,
	     Memorize the results */
	  if ( ChosenCrit( CriterP , NemParaP->Crit ) > 
	       ChosenCrit( & bestCritS , NemParaP->Crit ) )  {

	    bestCritS = *CriterP ;
	    bestTry   = irandom ;
	    bestbeta  = ParaP->Beta ;
	    memcpy( bestClM, ClassifM, npt * nk * sizeof(float) ) ;
	    memcpy( bestCen_KD, ParaP->Center_KD, nk * nd * sizeof(float) );
	    memcpy( bestDis_KD, ParaP->Disp_KD, nk * nd * sizeof(float) );
	  }
        }
    } /* end if NEM was successful */
  } /* end for each random start */

  /* If at least one success, return OK, 
     and restore best classification results */
  if ( nbsucc > 0 ) {
    int emptyk ; /* index of empty class, 0 if none */

    err = STS_OK ;
    memcpy( ClassifM, bestClM, npt * nk * sizeof(float) ) ;

    /* restore estimated centers/disp, 
       used if missing data and replace mode */
    memcpy( ParaP->Center_KD, bestCen_KD, nk * nd * sizeof(float) );
    memcpy( ParaP->Disp_KD,   bestDis_KD, nk * nd * sizeof(float) );
    EstimPara( ClassifM, DataP, nk, NemParaP->MissMode, SpecP,
	       &emptyk, ParaP ) ;
    ParaP->Beta = bestbeta ;

    *CriterP = bestCritS ;

    fprintf( out_stderr, 
             "Best start was %d (%s = %g)\n", bestTry+1, 
	     CritStrVC[ NemParaP->Crit ] , 
	     ChosenCrit( CriterP , NemParaP->Crit ) ) ; if(out_stderr) fflush( out_stderr );


    if ( CriterP->Errinfo.Kr != 0 )
      fprintf( out_stderr, "Error of best start = %5.1f\n", 
	       100.0 * CriterP->Errcur.Errorrate ) ; if(out_stderr) fflush( out_stderr );


    if ( Flog != NULL )
       fprintf( Flog, "Best start was %d (U = %g)\n", 
                bestTry+1, ChosenCrit( CriterP , NemParaP->Crit ) ) ; if(out_stderr) fflush( out_stderr );

  }
  /* otherwise, return last error status */

  GenFree( bestDis_KD ) ;
  GenFree( bestCen_KD ) ;
  GenFree( bestClM ) ;

  return err ;

} /* end RandNemAlgo() */


/* ------------------------------------------------------------------- */
static int NemAlgo
        (
         const DataT*       DataP,           /* I */
         const NemParaT*    NemParaP,        /* I */
         const SpatialT*    SpatialP,        /* I */
	 const ModelSpecT*  SpecP,     	     /* I */
	 const SampleDesT*  DescP,     	     /* I */
         FILE*              Flog,            /* I/O */
	 ModelParaT*        ParaP,     	     /* I/O */
         float*             CM,              /* I/O */
         WorkingT*          WorkP,           /* I/O */
         CriterT*           CriterP          /* O */
        ) 
/* ------------------------------------------------------------------- */
{
    StatusET        err         = STS_OK ;
    StatusET        sts         = STS_OK ;
    int             npt         = DataP->NbPts;
    int             nk          = SpecP->K ;

    int             ipt ;       /* 0..npt-1 : index of current point */
    int             ik ;        /* 0..nk-1  : index of current class */
    int             iter ;      /* 0..NbIters-1 : current iteration */
    int             converged ; /* TRUE if convergence reached */
    float           oldcrit ;   /* criterion at previous iteration */

    int             emptyk ;    /* 1..nk = index of empty class (0 = none) */


    /* Initialize old partition and criterion to extreme values
       (to make sure we skip the first convergence test)
    */
    for( ipt = 0 ; ipt < npt ; ipt ++ ) {
        for ( ik = 0 ; ik < nk ; ik ++ )
            WorkP->ColdM[ ( ipt * nk ) + ik ] = 0.0 ;
    }
    oldcrit = -MAXFLOAT ;

    WriteLogHeader( Flog, NemParaP->NbEIters, DataP->NbVars, SpecP ) ;
    fprintf( out_stderr, "  Iterations : " ) ; if(out_stderr) fflush( out_stderr );

    fprintf( out_stderr, "%4d ", 0 ) ; if(out_stderr) fflush( out_stderr );
              /*V1.05-g*/

    /* For each iteration of NEM (until : 
       convergence, or iteration count reached, or empty class) */
    for ( iter = 1, converged = FALSE ;         /*V1.05-g*/
          ( iter <= NemParaP->NbIters ) && 
          ( ! converged ) &&
          ( err == STS_OK ) ; 
          iter ++ )
    {
        fprintf( out_stderr, "\b\b\b\b\b" ) ; if(out_stderr) fflush( out_stderr );
       /*V1.05-g*/
        fprintf( out_stderr, "%4d ", iter ) ; if(out_stderr) fflush( out_stderr );
       /*V1.05-g*/

        if ( Flog != NULL ) fprintf( Flog, "%4d ", iter ) ; if(out_stderr) fflush( out_stderr );


        /* Save last computed partition as old partition */
        memcpy( WorkP->ColdM, CM, npt * nk * sizeof( float ) ) ;
	oldcrit = ChosenCrit( CriterP , NemParaP->Crit ) ;

        /* M-Step (Maximization) : compute current model parameters
           from old partition */
	if ( NemParaP->ParamFileMode != PARAM_FILE_FIX )
	  err = EstimPara( CM, DataP, nk,  NemParaP->MissMode, SpecP,
			   &emptyk, ParaP ) ;

	EstimBeta( SpecP->BetaModel, &NemParaP->BtaPsGrad, SpatialP, CM,
		   npt, nk,
		   & ParaP->Beta, & WorkP->Neighs ) ;

        /* If no empty class */
        if ( err == STS_OK )
        {
            /* E-Step (Expectation) : compute current partition from
               current model parameters and old partition */

	    sts = ComputePartitionFromPara( 0, DataP, NemParaP, SpecP, ParaP, 
				      SpatialP, 
				      CM, CriterP, Flog, WorkP ) ;

            /* Test convergence */
	    converged = HasConverged( NemParaP->CvTest, NemParaP->CvThres, 
				      WorkP->ColdM, CM, npt, nk, oldcrit, 
				      (Flog==NULL), ParaP->Beta, SpatialP,
				      NemParaP->Crit,
				      CriterP, WorkP ) ;  /*V1.06-o*/
        if ( sts == STS_E_INFINITE){
            err=sts;
            fprintf( out_stderr, "FATAL ERROR, NEM criteria reach infinite value\n"); if(out_stderr) fflush( out_stderr );

            if ( Flog != NULL )
               fprintf( Flog, "FATAL ERROR, NEM criteria reach infinite value\n"); if(Flog) fflush( Flog );

        }

        }
        else if ( err == STS_W_EMPTYCLASS ) /* Else -> class emptyk is empty */
        {
            fprintf( out_stderr, "Class %d empty at iteration %d\n", 
                     emptyk, iter ) ; if(out_stderr) fflush( out_stderr );

            if ( Flog != NULL )
               fprintf( Flog, " Class %d empty at iteration %d\n", 
                        emptyk, iter ) ; if(Flog) fflush( Flog );

        }
         
    } /* end for ( iter = 1, converged = FALSE ... ) */

    iter = iter - 1 ;  /*V1.05-g*/

    /* Compute and display value of criteria */  /*V1.03-a*/
    if ( iter == 0 )       /* model parameters estimation not yet done */
    {
        EstimPara( CM, DataP, nk, NemParaP->MissMode, SpecP, 
		   &emptyk, ParaP ) ;
        ComputePkFkiM( DataP, SpecP, ParaP, 
		       WorkP->PkFkiM, WorkP->LogPkFkiM ) ;
    }
    ComputeCrit( npt, nk, ParaP->Beta, CM, SpatialP, WorkP, CriterP ) ;
    CalcError( CM, npt, 1, & CriterP->Errinfo, & CriterP->Errcur );

    fprintf( out_stderr, "\n" ) ; if(out_stderr) fflush( out_stderr );
  /* -> to end iterations count line */
    fprintf( out_stderr, 
	     "  criterion NEM = %6.3f / Ps-Like = %6.3f / Lmix = %6.3f\n" ,
	     CriterP->U , CriterP->M , CriterP->L ) ; if(out_stderr) fflush( out_stderr );
    /*V1.03-g*/
    if ( CriterP->Errinfo.Kr != 0 ) 
      fprintf( out_stderr, "  error = %5.3f\n", CriterP->Errcur.Errorrate ) ; if(out_stderr) fflush( out_stderr );


    /* If convergence test was requested and no error occurred */
    if ( ( NemParaP->CvTest != CVTEST_NONE ) && ( err == STS_OK ) )
    {
        /* Display message about convergence */
        if ( converged )
        {
            fprintf( out_stderr, "  NEM converged after %d iterations\n", iter ) ; if(out_stderr) fflush( out_stderr );

        }
        else 
        {
            fprintf( out_stderr, "  NEM did not converge after %d iterations\n", 
                     iter ) ; if(out_stderr) fflush( out_stderr );

        }
    }

    return err ;

}  /* end of NemAlgo() */


/* ------------------------------------------------------------------- */
static void WriteLogHeader
        (
            FILE*               Flog,       /* I/O */
            int                 NbEIters,   /* I */
            int                 Nd,         /* I : dimension */
            const ModelSpecT*   SpecP       /* I */
        ) 
/* ------------------------------------------------------------------- */
{
    int     it, k, d ;

    if ( Flog == NULL )     return ;  /*V1.03-b*/

    fprintf( Flog, "%4s  %5s %5s %5s", "It", "UM", "PM", "Er" ) ; if(Flog) fflush( Flog );

    for ( it = 0 ; it < NbEIters ; it ++ )
    {
        fprintf( Flog, " %3s%-2d %3s%-2d %3s%-2d", 
		 "UE", it+1, "PE", it+1, "Er", it+1 ) ; if(Flog) fflush( Flog );

    }

    fprintf( Flog, " " ) ; if(Flog) fflush( Flog );


    fprintf( Flog, " %5s", "Beta" ) ; if(Flog) fflush( Flog );


    fprintf( Flog, " " ) ; if(Flog) fflush( Flog );


    for ( k = 0 ; k < SpecP->K ; k ++ )
    {
        fprintf( Flog, " %3s%02d", "P", k + 1 ) ; if(Flog) fflush( Flog );

    }

    fprintf( Flog, " " ) ; if(Flog) fflush( Flog );


    for ( k = 0 ; k < SpecP->K ; k ++ )
    {
        for ( d = 0 ; d < Nd ; d ++ )
        {
            fprintf( Flog, " %3s%02d_%1d", "M", k + 1, d + 1 ) ; if(Flog) fflush( Flog );

        }
    }

    fprintf( Flog, " " ) ; if(Flog) fflush( Flog );


    for ( k = 0 ; k < SpecP->K ; k ++ )
    {
        for ( d = 0 ; d < Nd ; d ++ )
        {
            fprintf( Flog, " %3s%02d_%1d", "D", k + 1, d + 1 ) ; if(Flog) fflush( Flog );

        }
    }

    fprintf( Flog, " " ) ; if(Flog) fflush( Flog );


    for ( k = 0 ; k < SpecP->K ; k ++ )
    {
        for ( d = 0 ; d < Nd ; d ++ )
        {
            fprintf( Flog, " %3s%02d_%1d", "n", k + 1, d + 1 ) ; if(Flog) fflush( Flog );

        }
    }

    fprintf( Flog, "\n" ) ; if(Flog) fflush( Flog );


}   /* end of WriteLogHeader() */



/* ------------------------------------------------------------------- */
static int ComputePartitionFromPara
        ( 
	 const int          Needinit,    /* I : 1 if need init, 0 if normal */
	 const DataT*       DataP,       /* I */
	 const NemParaT*    NemParaP,    /* I */
	 const ModelSpecT*  SpecP,       /* I */
	 ModelParaT*        ParaP,       /* I (temporary modification) */
         const SpatialT*    SpatialP,    /* I */
         float*             C_NK,        /* I/O */
         CriterT*           CriterP,     /* O */
	 FILE*              Flog,        /* I/O */
	 WorkingT*          WorkingP     /* T */
	)
/* ------------------------------------------------------------------- */
{
StatusET err = STS_OK ;

  /* Compute density of each point relatively to each class */
  ComputePkFkiM( DataP, SpecP, ParaP, 
                 WorkingP->PkFkiM, WorkingP->LogPkFkiM ) ;
  if ( Needinit ) {
    /* Compute initial partition by "blind" segmentation */
    float beta = ParaP->Beta ;
    ParaP->Beta = 0.0 ;
   ComputePartition( SpecP, ParaP, DataP, SpatialP, NemParaP, 
		      NULL, C_NK, WorkingP, CriterP ) ;
    ParaP->Beta = beta ;
  }

  /* Compute partition from parameters and previous partition */
 err = ComputePartition( SpecP, ParaP, DataP, SpatialP, NemParaP, 
		    Flog, C_NK, WorkingP, CriterP ) ;

  /* Write parameters to file */
  WriteLogClasses( Flog, DataP->NbVars, SpecP, ParaP ) ;

  if ( ( Needinit ) && ( Flog != NULL ) )
    fprintf( Flog, "\n" ) ; if(Flog) fflush( Flog );


return(err);

}   /* end of ComputePartitionFromPara() */


/* ------------------------------------------------------------------- */
static void WriteLogClasses
        (
            FILE*               Flog,       /* I/O */
            int                 Nd,         /* I : dimension */
            const ModelSpecT*   SpecP,      /* I */
            const ModelParaT*   ParaP       /* I */
        ) 
/* ------------------------------------------------------------------- */
{
    int     k, d ;

    if ( Flog == NULL )     return ;  /*V1.03-b*/

    fprintf( Flog, " " ) ; if(Flog) fflush( Flog );


    fprintf( Flog, " %5.3f", ParaP->Beta ) ; if(Flog) fflush( Flog );



    fprintf( Flog, " " ) ; if(Flog) fflush( Flog );


    for ( k = 0 ; k < SpecP->K ; k ++ )
    {
        fprintf( Flog, " %5.3f", ParaP->Prop_K[ k ] ) ; if(Flog) fflush( Flog );

    }

    fprintf( Flog, " " ) ; if(Flog) fflush( Flog );


    for ( k = 0 ; k < SpecP->K ; k ++ )
    {
        for ( d = 0 ; d < Nd ; d ++ )
        {
            fprintf( Flog, " %7.3f", ParaP->Center_KD[ (k*Nd) + d ] ) ; if(Flog) fflush( Flog );

        }
    }

    fprintf( Flog, " " ) ; if(Flog) fflush( Flog );


    for ( k = 0 ; k < SpecP->K ; k ++ )
    {
        for ( d = 0 ; d < Nd ; d ++ )
        {
            fprintf( Flog, " %7.3f", ParaP->Disp_KD[ (k*Nd) + d ] ) ; if(Flog) fflush( Flog );

        }
    }

    fprintf( Flog, " " ) ; if(Flog) fflush( Flog );


    for ( k = 0 ; k < SpecP->K ; k ++ )
    {
        for ( d = 0 ; d < Nd ; d ++ )
        {
            fprintf( Flog, " %7.1f", ParaP->NbObs_KD[ (k*Nd) + d ] ) ; if(Flog) fflush( Flog );

        }
    }

    fprintf( Flog, "\n" ) ; if(Flog) fflush( Flog );


}   /* end of WriteLogClasses() */


/* ------------------------------------------------------------------- */
static int HasConverged       /* ret : TRUE if convergence test satisfied */
(
  const CvemET    CvTest,     /* I : which convergence test to use */
  const float     CvThres,    /* I : convergence threshold to use */
  const float*    ColdM,      /* I : previous classification matrix */
  const float*    CM,         /* I : current  classification matrix */
  const int       Npt,        /* I : number of objects */
  const int       Nk,         /* I : number of classes */
  const float     OldCrit,    /* I : previous criterion */ 
  const int       CritToDo,   /* I : 1 = criterion needs to be computed */
  const float     Beta,       /* I : value of beta to compute criterion */
  const SpatialT* SpatialP,   /* I : spatial info to compute criterion */
  const CritET    WhichCrit,  /* I : what criterion to use */
  CriterT*        CriterP,    /* I/O : current criterion */
  WorkingT*       WorkP       /* I: LogPkFkiM, T: Neighs */
) 
/* ------------------------------------------------------------------- */
{
  switch( CvTest ) {
    case CVTEST_CLAS: {
      int ipt, ik ;
      float maxdif = 0.0 ;

      for ( ipt = 0 ; ipt < Npt ; ipt ++ )
	for ( ik = 0 ; ik < Nk ; ik ++ ) {
	  int i = (ipt * Nk) + ik ;
	  float dif = CM[ i ] - ColdM[ i ] ;
	  if ( dif < 0 ) dif = -dif ;
        
	  if ( dif > maxdif ) maxdif = dif ;
	}

      return ( maxdif < CvThres ) ;
    }

    case CVTEST_CRIT: {
      float curcrit ;
      float critdif ;

      if ( CritToDo )
	ComputeCrit( Npt, Nk, Beta, CM, SpatialP, WorkP, CriterP ) ;

      curcrit = ChosenCrit( CriterP , WhichCrit ) ;
      if ( curcrit != 0 ) 
	critdif = fabs( ( curcrit - OldCrit ) / curcrit ) ;
      else
	critdif = MAXFLOAT ;

      return ( critdif < CvThres ) ;
    }
    

    case CVTEST_NONE: return FALSE ;
    default:          return FALSE ;
  }

}   /* end of HasConverged() */






/* ------------------------------------------------------------------- */
static int  /* 0 = OK,  1 = beta diverges to infinity, 2 = call bug */
EstimBeta
(
  const BetaET       BetaModel,  /* I : beta estimation mode */
  const BtaPsGradT*  BtaPsGradP, /* I : gradient ascent parameters */
  const SpatialT*    SpatialP,   /* I : neighborhood system */
  const float*       C_NK,       /* I : classification matrix (Npt, NbK)*/
  const int          Npt,        /* I : number of points */
  const int          Nbk,        /* I : number of classes */
  float*             BetaP,      /* I/O : previous and new beta */
  PtNeighsT*         NeighsP     /* W : to store a site's neighbors */
)
/* ------------------------------------------------------------------- */
{
  int    iter ;           /* current iteration : 0..NbIter-1 */
  int    conv ;           /* TRUE if convergence reached */
  int    ipt ;            /* current site  : 0..Npt-1 */
  int    k ;              /* current class : 0..Nbk-1 */
  float  con_ik ;         /* sum_{j\neighbor{i}} v_{ij} c_{jk} */
  float  exp_ik ;         /* exp{ beta * con_ik } */
  float  sumk_exp ;       /* sum_k exp_ik */
  float  sumk_con_exp ;   /* sum_k con_ik exp_ik */
  float  sumk_con2_exp ;  /* sum_k (con_ik)^2 exp_ik */
  float  sumk_cik_con ;   /* sum_k c_ik con_ik */
  float  crit ;           /* log pseudo-likelihood */
  float  grad ;           /* first deriv of log pseudo-likelihood */
  float  dsec ;           /* - second deriv of log pseudo-likelihood */

  GetNeighFT* fGetNeigh ; /* function fetching neighbors of a site */


  if ( ( BetaModel != BETA_PSGRAD ) || ( SpatialP->Type == TYPE_NONSPATIAL ) )
    return 0 ;

  /* Get function specific to chosen spatial model */
  if ( GetSpatialFunc( SpatialP->Type, 
			       & fGetNeigh ) != STS_OK )
    return 2 ;

  /* For each iteration of gradient ascent */
  for ( iter = 0,
	  conv = FALSE ;
	( iter < BtaPsGradP->NbIter ) &&
	  ( ! conv ) ;
	iter ++, conv = (fabs( grad ) < ( BtaPsGradP->ConvThres * Npt )) ) {

    /* Compute crit, grad and dsec by summation over all sites */
    for ( ipt = 0, crit = 0.0, grad = 0.0, dsec = 0.0 ;
	  ipt < Npt ; ipt ++ ) {

      /* Get site's neighbours */
      int nbn = fGetNeigh( ipt, &SpatialP->NeighData, NeighsP ) ;

      for ( k = 0,
	      sumk_exp      = 0.0,
	      sumk_con_exp  = 0.0,
	      sumk_con2_exp = 0.0,
	      sumk_cik_con  = 0.0 ;
	    k < Nbk ; k ++ )  {

	con_ik = SumNeighsOfClass( k, nbn, Nbk, NeighsP->NeighsV, C_NK );
	exp_ik = exp( (*BetaP) * con_ik ) ;

	sumk_exp        += exp_ik ;
	sumk_con_exp    += con_ik * exp_ik ;
	sumk_con2_exp   += con_ik * con_ik * exp_ik ;
	sumk_cik_con    += C_NK[ ( ipt * Nbk ) + k ] * con_ik ;
      }

      crit += (*BetaP) * sumk_cik_con - log( sumk_exp ) ;
      grad += sumk_cik_con - sumk_con_exp / sumk_exp ;

      /* might get NaN because of high sumk_exp*/
      dsec += (float) ( (double)  
	( sumk_con2_exp * sumk_exp - sumk_con_exp * sumk_con_exp )
	/ (double) ( sumk_exp * sumk_exp ) ) ; 
    }

    /* Update beta */
    if ( BtaPsGradP->Step <= 0.0 ) {
      dsec = dsec * 4 ;  /* reduce the step to avoid overstepping */
      if ( dsec < ( Npt / 10 ) ) 
	dsec = Npt / 10 ;

      (*BetaP) += grad / dsec ;
    }
    else {
      (*BetaP) += grad * ( BtaPsGradP->Step / Npt ) ;
    }

  } /* end  For each iteration of gradient ascent */

  if ( (*BetaP) > MAX_BETA ) {
    (*BetaP) = MAX_BETA ;
    return 1 ;
  }

  if ( (*BetaP) < MIN_BETA ) {
    (*BetaP) = MIN_BETA ;
    return 1 ;
  }

  if ( isnan(*BetaP) ) {
    (*BetaP) = 0 ;
    return 1 ;
  }


  return 0 ;
  /*???*/
}   /* end of EstimBeta() */


/* ------------------------------------------------------------------- */
static StatusET ComputePkFkiM
        (
            const DataT*        DataP,      /* I */
            const ModelSpecT*   SpecP,      /* I */
            const ModelParaT*   ParaP,      /* I */
            double*              PkFkiM,     /* O */
            float*              LogPkFkiM   /* O */
        ) 
/* ------------------------------------------------------------------- */
{
    StatusET        err = STS_OK ;
    int             i ;
    int             k ;
    int             npt   = DataP->NbPts;
    int             nk    = SpecP->K ;
    CompuDensFT*    fCompuDens ;      /*V1.06-c*/

    /* Get functions specific to chosen noise model */
    if ( ( err = GetDensityFunc( SpecP ,
                                 & fCompuDens ) ) != STS_OK )     /*V1.06-c*/
    {
        return err ;
    }


    /* Compute densities of each point relative to each class */
    for ( k = 0 ; k < nk ; k ++ )
    {
        double pk    = ParaP->Prop_K[ k ] ;
        float logpk ;

	if ( pk > EPSILON ) 
	  logpk = log( pk ) ;
	else
	  {
	    logpk = atof( "-Inf" ) ;
	    err = STS_W_EMPTYCLASS ;
	  }

        for ( i = 0 ; i < npt ; i ++ )
        {
            double fki ;
	    float  logfki ;

            fCompuDens( DataP->NbVars, k, ParaP, 
                        & (DataP->PointsM[ i * DataP->NbVars ]) ,
                        & fki , & logfki );

            PkFkiM   [ ( i * nk ) + k ] = pk * fki ;
            LogPkFkiM[ ( i * nk ) + k ] = logpk + logfki ;
        }
    }

    return err ;

}   /* end of ComputePkFkiM() */



/* ------------------------------------------------------------------- */
static int ComputePartition
      (
          const ModelSpecT*   SpecP,      /* I */
          const ModelParaT*   ParaP,      /* I */
          const DataT*        DataP,      /* I */
          const SpatialT*     SpatialP,   /* I */
          const NemParaT*     NemParaP,   /* I */
          FILE*               Flog,       /* I/O */
          float*              CM,         /* I/O */
          WorkingT*           WorkP,      /* I:PkFkiM O:CtmpM,Neighs,CiNumV */
          CriterT*            CriterP     /* O */
      ) 
/* ------------------------------------------------------------------- */
{
  StatusET sts ;

  if ( NemParaP->Algo == ALGO_GEM ) 
    {
      sts = ComputePartitionGEM( SpecP, ParaP, DataP, SpatialP, NemParaP, 
				 Flog, CM, WorkP, CriterP ) ;
    }
  else
    {
      sts = ComputePartitionNEM( SpecP, ParaP, DataP, SpatialP, NemParaP, 
				 Flog, CM, WorkP, CriterP ) ;
    }
  return sts;
}   /* end of ComputePartition() */




/* ------------------------------------------------------------------- */
static int ComputePartitionNEM
      (
          const ModelSpecT*   SpecP,      /* I */
          const ModelParaT*   ParaP,      /* I */
          const DataT*        DataP,      /* I */
          const SpatialT*     SpatialP,   /* I */
          const NemParaT*     NemParaP,   /* I */
          FILE*               Flog,       /* I/O */
          float*              CM,         /* I/O */
          WorkingT*           WorkP,      /* I:PkFkiM O:CtmpM,Neighs,CiNumV */
          CriterT*            CriterP     /* O */
      ) 
/* ------------------------------------------------------------------- */
{
    StatusET        err ;
    int             ivis ;  /*V1.04-c*/
    int             itere ;
    const double*   pkfkiM = WorkP->PkFkiM ;
    float*          ctmpM = WorkP->CtmpM;
    int             npt   = DataP->NbPts;
    int             nk    = SpecP->K ;
    GetNeighFT*     fGetNeigh ;     /*V1.06-c*/


    /* Get function specific to chosen spatial model */
    if ( ( err = GetSpatialFunc( SpatialP->Type, 
                                 & fGetNeigh ) ) != STS_OK )     /*V1.06-c*/
      return err ;

    /* Compute criterion after M-Step */ /*V1.03-a*/
    WriteLogCrit( Flog, npt, nk, ParaP->Beta, CM, SpatialP, WorkP, CriterP ) ;

    /* Iteratively compute new partition */
    for ( itere = 0 ; itere < NemParaP->NbEIters ; itere ++ ) {

        /* Save last computed value of partition in temporary buffer */
        memcpy( ctmpM, CM, npt * nk * sizeof( float ) ) ;

        /* For each point xi */
        for ( ivis = 0 ; ivis < npt ; ivis ++ ) {

	    int ipt = DataP->SiteVisitV[ ivis ] ;  /*V1.04-c*/

            /* If not in init_label mode or label is not known */
	    if ( ( NemParaP->InitMode   != INIT_LABEL ) || 
		 ( DataP->LabelV[ ipt ] == 0 ) ) {

		ComputeLocalProba( ipt, nk, ParaP, &(SpatialP->NeighData), 
				   fGetNeigh, pkfkiM, 
				   (NemParaP->SiteUpdate == UPDATE_SEQ)? 
				   CM : ctmpM, 
				   &CM[ipt*nk], 
				   &WorkP->Neighs, WorkP->CiNumV );

		/* Eventual C-step */
		if ( NemParaP->Algo == ALGO_NCEM ) {
		  int kmap   = ComputeMAP( CM , ipt , nk , NemParaP->TieRule, 
					   WorkP->KmaxesV ) ;

		  LabelToClassVector( nk, kmap, &CM[ ipt * nk ] ) ;
		}
	    }
            /* Else (label known) -> do not recompute its classification */

        } /* end for ( ivis = 0 ; ivis < npt ; ivis ++ ) */

	/* Compute criterion after this internal iteration of E-Step */
	WriteLogCrit( Flog, npt, nk, ParaP->Beta, CM, SpatialP, 
		      WorkP, CriterP ) ;

    } /* end for ( itere = 0 ; itere < NemParaP->NbEIters ; itere ++ ) */

if( CriterP->U == -INFINITY ){
    return(STS_E_INFINITE);
}else{
    return STS_OK ;
}
}   /* end of ComputePartitionNEM() */




/* ------------------------------------------------------------------- */
/* Compute fuzzy (or hard) partition with fixed parameters, using
   Monte-Carlo simulations */ /*1.04-h*/

static int ComputePartitionGEM /* +++ */
      (
          const ModelSpecT*   SpecP,      /* I */
          const ModelParaT*   ParaP,      /* I */
          const DataT*        DataP,      /* I */
          const SpatialT*     SpatialP,   /* I */
          const NemParaT*     NemParaP,   /* I */
          FILE*               Flog,       /* I/O */
          float*              CM,         /* I/O */
          WorkingT*           WorkP,      /* I:PkFkiM O:CtmpM,Neighs,CiNumV */
          CriterT*            CriterP     /* O */
      ) 
/* ------------------------------------------------------------------- */
{
    StatusET        err ;
    const double*   pkfkiM = WorkP->PkFkiM ;
    float*          ctmpM = WorkP->CtmpM;
    int             npt   = DataP->NbPts;
    int             nk    = SpecP->K ;
    GetNeighFT*     fGetNeigh ;      /*V1.06-c*/

    float*          z_nk ;     /* currently simulated partition */
    int*            occur_nk ; /* occurrence count for class h at site i */
    int             icycle ;   /* current Monte-Carlo cycle */
    int             kdraw ;    /* drawn class : 0..nk-1 */
    int             ivis ;


    if ( ( err = GetSpatialFunc( SpatialP->Type, 
                                 & fGetNeigh ) ) != STS_OK )     /*V1.06-c*/
        return err ;

    /* Compute criterion after M-Step */ /*V1.03-a*/
    WriteLogCrit( Flog, npt, nk, ParaP->Beta, CM, SpatialP, WorkP, CriterP ) ;

    /* Allocate all memory */
    if ( ( z_nk = GenAlloc( npt * nk, sizeof( float ), 
		     0, "ComputePartitionGEM", "z_nk" ) ) == NULL )
	return STS_E_MEMORY ;

    if ( ( occur_nk = GenAlloc( npt * nk, sizeof( int ),
			 0, "ComputePartitionGEM", "occur_nk") ) == NULL )
	return STS_E_MEMORY ;

    /* Initialize values */
    {
      int ipt, k ;

      /* Initialize hard classification from fuzzy classification, by MPM */
      for ( ipt = 0 ; ipt < npt ; ipt ++ )
	{
	  int kmap   = ComputeMAP( CM , ipt , nk , NemParaP->TieRule, 
				   WorkP->KmaxesV ) ;

	  LabelToClassVector( nk, kmap, &z_nk[ ipt * nk ] ) ;
	}

      /* Initialize to zero all occurrence counts */
      for ( ipt = 0 ; ipt < npt ; ipt ++ )
	for ( k = 0 ; k < nk ; k ++ )
	  occur_nk[ ( ipt * nk ) + k ] = 0 ;
    }


    /* Iteratively run Monte-Carlo simulation with startup and
       recording durations */
    for ( icycle = 0 ; 
	  icycle < NemParaP->NbEIters * (NBCYCLES_STARTUP + NBCYCLES_RECORD) ;
	  icycle ++ ) {

        /* Save last computed value of partition in temporary buffer */
        memcpy( ctmpM, z_nk, npt * nk * sizeof( float ) ) ;

	/* Run through each site */
        for ( ivis = 0 ; ivis < npt ; ivis ++ ) {
	    int ipt = DataP->SiteVisitV[ ivis ] ;  /*V1.04-c*/

            /* If not in init_label mode or label is not known */
	    if ( ( NemParaP->InitMode   != INIT_LABEL ) || 
		 ( DataP->LabelV[ ipt ] == 0 ) ) {

		ComputeLocalProba( ipt, nk, ParaP, &(SpatialP->NeighData), 
				   fGetNeigh, pkfkiM, 
				   (NemParaP->SiteUpdate == UPDATE_SEQ)? 
				   z_nk : ctmpM, 
				   &CM[ipt*nk], 
				   &WorkP->Neighs, WorkP->CiNumV );

		/* Update zi by multinomial draw from ti's k proportions */
		kdraw = Multinomial( nk , & CM[ ipt * nk ] ) - 1 ;
                    
		LabelToClassVector( nk, kdraw, &z_nk[ ipt * nk ] ) ;

		/* Update occurrence counter if in record phase */
		if ( icycle >= NemParaP->NbEIters * NBCYCLES_STARTUP )
		  occur_nk[ ( ipt * nk ) + kdraw ] += 1 ;
	    }
            /* Else (label known) -> do not recompute its classification */

        } /* end for ( ivis = 0 ; ivis < npt ; ivis ++ ) */

      } /* end for ( icycle = 0 ; ... ) */


    /* Compute frequency of classes at each site */
    if ( NemParaP->NbEIters * NBCYCLES_RECORD > 0 ) {
      int ipt, k ;
      int nbrecord = NemParaP->NbEIters * NBCYCLES_RECORD ;

      for ( ipt = 0 ; ipt < npt ; ipt ++ )
	{
	  for ( k = 0 ; k < nk ; k ++ )
	    CM[ ( ipt * nk ) + k ] = 
	      ( (float) occur_nk[ ( ipt * nk ) + k ] ) / nbrecord ;
	}
    }

    /* Compute criterion after this internal iteration of E-Step */
    WriteLogCrit( Flog, npt, nk, ParaP->Beta, CM, SpatialP, WorkP, CriterP ) ;

    /* Free all memory */
    GenFree( z_nk ) ;
    GenFree( occur_nk ) ;


    return STS_OK ;

}   /* end of ComputePartitionGEM() */



/*-------------------------------------------------------*/
static int ComputeLocalProba
 (
  const int            Ipt,        /* I: current site */
  const int            Nk,         /* I: nb of classes */
  const ModelParaT*    ParaP,      /* I: parameters (use Beta) */
  const NeighDataT*    NeighDataP, /* I: neighborhood system */
  GetNeighFT*          FGetNeigh,  /* I: function to fetch neighbors */
  const double*        PkfkiM,     /* I: pk fk(xi) */
  const float*         Cin_NK,     /* I: class of neighbors */
  float*               Cout_K,     /* O: class of site Ipt */
  PtNeighsT*           NeighsP,    /* W: to store fetched neighbors */
  double*              Cinum_K     /* W: to store cik's numerators */
 )
{
  static int first=TRUE; /* FALSE once a zero density point has occurred */

  int     k ;      /* current class : 0..Nk-1 */
  int     nbn ;    /* nb of neighbours of ipt */
  double  cumnum ; /* for point i, sum over k of pkfki*context */

  /* Get point's neighbours */
  nbn = FGetNeigh( Ipt, NeighDataP, NeighsP ) ;

  /* For each class k : 
     - compute the contribution of xi's neighbors
     (j), based on their last cjk
     - based on the neighbors contribution and on the density for 
     class k at xi, compute cik's numerator
     - increment the cumulated sum of ci's K numerators
  */
  for ( k = 0, cumnum = 0.0 ; k < Nk ; k ++ ) {

    /* i's class k context contrib. */
    float   context = SumNeighsOfClass( k, nbn, Nk, NeighsP->NeighsV, Cin_NK );

    Cinum_K[ k ] = PkfkiM[ ( Ipt * Nk ) + k ] * 
      exp( (double) ParaP->Beta * context ) ;

    cumnum = cumnum + Cinum_K[ k ] ;

  } /* end for ( k = 0, cumnum = 0.0 ; k < nk ; k ++ ) */

  /* cik = numerator of cik / sum of ci's K numerators */
  if ( cumnum > 0 ) {
    if ( cumnum > EPSILON ) {
      double   invZ = ( 1 / cumnum ) ;

      for ( k = 0 ; k < Nk ; k ++ )
	Cout_K[ k ] = invZ * Cinum_K[ k ] ;
    }
    else {    /*V1.06-r*/
      double   invZ = ( 1 / ( cumnum / EPSILON ) ) ;

      for ( k = 0 ; k < Nk ; k ++ )
	Cout_K[ k ] = invZ * ( Cinum_K[ k ] / EPSILON ) ;
    }
  }
  else {
    double   invZ = ( 1.0 / Nk ) ;

    for ( k = 0 ; k < Nk ; k ++ )
      Cout_K[ k ] = invZ ;

    if ( first ) {
      first = FALSE ;
      fprintf( out_stderr, "Warning : pt %d density = 0\n", Ipt ) ; if(out_stderr) fflush( out_stderr );

    }
  }

 return ( cumnum > 0.0 ) ;
}

 
/*-------------------------------------------------------*/
static void WriteLogCrit
     (
        FILE*               Flog,        /* I/O */
        const int           Npt,         /* I */
        const int           Nk,          /* I */
        const float         Beta,        /* I : spatial coefficient */
        const float*        CM,          /* I : class. matrix [Npt*Nk] */
        const SpatialT*     SpatialP,    /* I : neighborhood system */
        WorkingT*           WorkP,       /* I :PkFkiM, O:Neighs */
        CriterT*            CriterP      /* O : computed criteria */ 
     )
/*-------------------------------------------------------*/
{
  float mult;

  if ( Flog != NULL ) {

    ComputeCrit( Npt, Nk, Beta, CM, SpatialP, 
		 WorkP, CriterP ) ;

    mult = exp( - ((int) ( log(Npt/1000.)/log(10) )) * log( 10 ) ) ;

    fprintf( Flog, " %5.0f %5.0f %5.3f",
	     CriterP->U * mult, CriterP->M * mult,
	     CriterP->Errcur.Errorrate ) ; if(Flog) fflush( Flog );

  }
}



/*-------------------------------------------------------*/
static int Multinomial(int km, const float *tk) 
/*
 Entrees
     km           nombre de modalites
     tk           probabilites de chaque modalites

 Retour
     valeur retenue
*/

{
    int        k ;
    float      pk,x;

    x = RandomFloat( 0.0, 1.0 ) ;  /*V1.05-l*/
    pk = 0.0;
    for ( k = 0 ; k < km ; k ++ )
    {
	pk = pk + tk[ k ] ;
	if ( x <= pk ) return ( k + 1 ) ;
    }
    return km;
}



/* ------------------------------------------------------------------- */
static StatusET ComputeCrit
     (
        int                 Npt,         /* I */
        int                 Nk,          /* I */
        float               Beta,        /* I : spatial coefficient */
        const float*        CM,          /* I : class. matrix [Npt*Nk] */
        const SpatialT*     SpatialP,    /* I : neighborhood system */
        WorkingT*           WorkP,       /* I :LogPkFkiM, T:Neighs */
        CriterT*            CriterP      /* O : computed criteria */ 
     )
/* ------------------------------------------------------------------- */
{
    StatusET    err ;
    int         i ;
    int         k ;
    NeighT*     neiV = WorkP->Neighs.NeighsV ;
    GetNeighFT* fGetNeigh ;     /*V1.06-c*/

    if ( ( err = GetSpatialFunc( SpatialP->Type, 
                                 & fGetNeigh ) ) != STS_OK )    /*V1.06-c*/
    {
        return err ;
    }

    CriterP->D = 0.0 ; /* D = sum[i]sum[k] cik (log pkfki-log cik) */
    CriterP->G = 0.0 ; /* G = sum[i]sum[j]sum[k] cik cjk wij */
    CriterP->U = 0.0 ; /* U = D + 0.5 * beta * G */
    CriterP->M = 0.0 ; /* M = D + beta * G + Z */
    CriterP->L = 0.0 ; /* L = sum[i] log( sum[k] pkfki ) = sum[i] log fi */
    CriterP->Z = 0.0 ; /* Z =-sum[i] log(sum[k] exp(bta * sum[j~i] wij cjk)) */

    /* For each point i */
    for ( i = 0 ; i < Npt ; i ++ )
    {
        double fi ;  /* fi = sum[k] pkfki */
	float zi ;  /* zi = sum[k] exp(bta * sum[j~i] wij cjk) */

        /* Get point i's neighbors */
        int nbn = fGetNeigh( i, & SpatialP->NeighData, &(WorkP->Neighs) ) ;

        /* For each class k */
        for ( k = 0 , fi = 0.0, zi = 0.0 ; 
	      k < Nk ; 
	      k ++ )
        {
            float   cik      = CM[ i * Nk + k ] ;
	    float   pik      = SumNeighsOfClass( k, nbn, Nk, neiV, CM ) ;

            /* If point i has non-null membership in class k */
            if ( cik > MINFLOAT )
            {
                /* Add to criteria the contribution of i's membership in k */
                float logpkfki = WorkP->LogPkFkiM[ i * Nk + k ] ;
                float dik = cik * ( logpkfki - log( cik ) ) ;
                float gik = cik * pik ;

                CriterP->D = CriterP->D + dik ;
                CriterP->G = CriterP->G + gik ;
            }
            /* Else point i for class k has no contribution to criteria D/G*/

	    fi = fi + WorkP->PkFkiM[ i * Nk + k ] ;/*V1.05-j*/
	    zi = zi + exp( Beta * pik ) ;

        } /* end   For each class k */

	CriterP->L = CriterP->L + log( fi ) ;
	CriterP->Z = CriterP->Z - log( zi ) ;

    } /* end   For each point i */

    /* Criterion U is a combination of criteria D and G */
    CriterP->U = CriterP->D + 0.5 * Beta * CriterP->G ; /*V1.05-k*/
    CriterP->M = CriterP->D + Beta * CriterP->G + CriterP->Z ; /*V1.05-k*/

    CalcError( CM, Npt, 1, & CriterP->Errinfo, & CriterP->Errcur );

    return STS_OK ;

}   /* end of ComputeCrit() */



/* ------------------------------------------------------------------- */
static void CalcError                /*V1.04-e*/
      ( 
       const float*     Cla_N_Kc, /* I : found classification */
       const int        N,        /* I : number of points */
       const int        Harden,   /* I : 1 use hardened CM, 0 use CM itself */
       const ErrinfoT*  ErrinfoP, /* I : info to compute error */
       ErrcurT*         ErrcurP   /* O : computed error and agreement */
      )
/* ------------------------------------------------------------------- */
{
  int    Kc = ErrinfoP->Kc ;
  int    Kr = ErrinfoP->Kr ;
  int    Km = ErrinfoP->Km ;
  float* Loclas_N_Kc = ErrcurP->Loclas_N_Kc ;

  int*   kmaxes_Kc ;  /* used by ComputeMAP */
  int    kmap ;       /* MAP class for current object : 0..Kc-1 */
  int    in ;         /* cur object : 0..N-1 */
  int    ikc ;        /* cur class of current classification : 0..Km-1 */
  int    ikr ;        /* cur class of ref. classification : 0..Km-1 */

  float  bestagree ;  /* cur max agreement over checked permutations : 0..N */
  int    ikmfact ;    /* cur permutation : 0..Km!-1 */
  float  sumagree ;   /* summed diagonal of agreement matrix : 0..N */
  int    ikm ;        /* cur diagonal element of agreement matrix : 0..Km-1 */
  int    kperm ;      /* permuted ikm using permutation #ikmfact  : 0..Km-1 */


  if ( Kr == 0 ) {
    ErrcurP->Errorrate = mknan() ;
    return ;
  }

  kmaxes_Kc = GenAlloc( Kc, sizeof( int ), 0, "CalcError", "kmaxes_Kc" ) ;
  if ( kmaxes_Kc == NULL ) 
    return ;

  /* Copy classification locally, eventually harden */
  memcpy( Loclas_N_Kc, Cla_N_Kc, N * Kc * sizeof(float) ) ;
  if ( Harden ) {
    for ( in = 0 ; in < N ; in ++ ) {
      kmap = ComputeMAP( Loclas_N_Kc, in, Kc, ErrinfoP->TieRule, kmaxes_Kc ) ;
      LabelToClassVector( Kc, kmap, & Loclas_N_Kc[ in * Kc ] ) ;
    }
  }

  /* 
   *  Compute agreement matrix : for each class ikc of current partition
   *  and each class ikr of reference partition, their agreement is the
   *  scalar product of their [N] classification vectors.
   *  For added classes (Kc to Km-1 and Kr to Km-1) : agreement = 0.
   */
  for ( ikc = 0 ; ikc < Km ; ikc ++ ) {
    for ( ikr = 0 ; ikr < Km ; ikr ++ ) {
      if ( ( ikc < Kc ) && ( ikr < Kr ) ) {
	ErrcurP->Agree_Km_Km[ ikc * Km + ikr ] = 0.0 ;
	for ( in = 0 ; in < N ; in ++ )
	  ErrcurP->Agree_Km_Km[ ikc * Km + ikr ] += 
	    Loclas_N_Kc[ in * Kc + ikc ] *
	    ErrinfoP->Refclas_N_Kr[ in * Kr + ikr ] ;
      }
      else
	ErrcurP->Agree_Km_Km[ ikc * Km + ikr ] = 0.0 ;
    }
  }

  /* Find permutation of lines with maximum sum of diagonal */
  bestagree            = 0.0 ; 
  ErrcurP->Ibestpermut = 0 ;
  for ( ikmfact = 0 ; ikmfact < ErrinfoP->Kmfac ; ikmfact ++ ) {
    /* sum diagonal of permuted matrix and compare to best agreement */
    for ( ikm = 0, sumagree = 0.0 ; ikm < Km ; ikm ++ ) {
      kperm     = ErrinfoP->Perm_Kmfac_Km[ ( ikmfact * Km ) + ikm ] ;
      sumagree += ErrcurP->Agree_Km_Km[ kperm * Km + ikm ] ;
    }
    if ( sumagree > bestagree ) {
      bestagree            = sumagree ;
      ErrcurP->Ibestpermut = ikmfact ;
    }
  }

  ErrcurP->Errorrate = ( N - bestagree ) / N ;

  GenFree( kmaxes_Kc ) ;
}


/* ------------------------------------------------------------------- */
static float SumNeighsOfClass
     (
        int             K ,     /* I : class to search, 0..Nk-1 */
        int             Nbn ,   /* I : number of neighbours */
        int             Nk,     /* I : number of classes */
        const NeighT*   NeiV ,  /* I : neighbours, [Nbn] */
        const float*    CM      /* I : classification matrix [Npt*Nk] */
     ) 
/* ------------------------------------------------------------------- */
{
    float   Sumj_wij_cjk ; /* sum of neighbours'membership in class K */

    int     n ;            /* neighbour counter : 0..nbn-1 */
    float   sumj_wij ;     /* divide sumj_wij_cjk by this normalizing factor */

    for ( n = 0, Sumj_wij_cjk = 0.0, sumj_wij = 0.0 ; 
          n < Nbn ; 
          n ++ )
    {
        int     j   = NeiV[ n ].Index ;
        float   wij = NeiV[ n ].Weight ;
        float   cjk = CM[ (j * Nk) + K ] ;

        Sumj_wij_cjk  = Sumj_wij_cjk + ( wij * cjk ) ;
        sumj_wij      = sumj_wij + wij ;
    }

    if ( sumj_wij != 0.0 )
    {
        /* sumj_wij_cjk  = sumj_wij_cjk / sumj_wij ; +++*/
    }

    return Sumj_wij_cjk ;

}   /* end of SumNeighsOfClass() */


/* ~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
