#ifndef NEM_TYP_H
#define NEM_TYP_H

/*\
    NEM_TYP.H

    Programme NEM (Neighborhood EM) : tous les types echanges entre modules

    Van Mo DANG       Janvier 96


    Vers-mod  Date         Description

    V1.03-a   01-NOV-1996  Add log-likelihood criterion in CriterT
    V1.04-a   04-OCT-1997  Add BetaET for beta estimation mode
    V1.04-b   04-OCT-1997  Add BtaStep and BtaMode in NemParaT
    V1.04-c   10-OCT-1997  Add CvThres in NemParaT
    V1.04-d   10-OCT-1997  Add Seed in NemParaT
    V1.04-e   10-OCT-1997  Type OrderET and add VisitOrder in NemParaT
    V1.04-f   13-OCT-1997  Debug and RefName in NemParaT
    V1.04-g   13-OCT-1997  Error in CriterT
    V1.04-h   13-OCT-1997  CritET and Crit in NemParaT
    V1.05-a   12-JAN-1998  NbMissing in DataT, MisModeET, MissMode in NemParaT
    V1.05-b   12-JAN-1998  MissMode in parameters of EstimNoiseFT
    V1.05-c   16-JAN-1998  EstimNoiseFT : EmptyK_P as parameter, return status
    V1.05-d   05-FEB-1998  Add Markov fuzzy pseudo-like. (CritET and CriterT)
    V1.06-a   17-JUN-1998  NoiseModel -> StatModel, structure changed
    V1.06-b   23-JUN-1998  LEN_LINE <- nem_exe.c
    V1.06-c   30-JUN-1998  Add define EPSILON
    V1.06-d   03-AUG-1998  Add INIT_MIXFIX and INIT_MIXINI
    V1.06-e   10-SEP-1998  Add SiteUpdate in NemParaT and UpdET
    V1.06-f   10-SEP-1998  Add TieRule in NemParaT and TieET
    V1.06-g   15-SEP-1998  Add BtaPsGrad in NemParaT and BtaPsGradT struct
    V1.06-h   20-SEP-1998  Add ErrInfo in CriterT and ErrInfoT struct
    V1.06-i   21-SEP-1998  Change CvTest in CriterT and add CvemET
    V1.06-j   30-NOV-1998  Add EPSILON_INV
    V1.06-k   01-DEC-1998  FkP double* instead of *float in compudensft
    V1.07-a   26-FEB-1999  FAMILY_BERNOULLI added
    1.08-a    20-JUI-2017  GG   Add param input by file rather than by arguments 
\*/

/*
 *  Constant definitions
 */

#define    TRUE             1
#define    FALSE            0

#define    LEN_NOISEMODEL   20
#define    LEN_FILENAME     200
#define    LEN_LINE         30000  /*V1.06-b*/


#define    MAX_PTS          1000000L
#define    MAX_VARS         100

#define    CONV_THRES       0.001

#define    EXT_OUTNAMEHARD  ".cf"
#define    EXT_OUTNAMEFUZZY ".uf"
#define    EXT_MFNAME       ".mf"
#define    EXT_LOGNAME      ".log"
#define    EXT_INITPARAM    ".m"

#define    EPSILON          1e-20  /* to check for FP zero or equality */
#define    EPSILON_INV      1e20   /* multiply by this for small floats */

#define DEFAULT_ALGO         ALGO_NEM
#define DEFAULT_BETA         1.0
#define DEFAULT_BTAMODE      BETA_FIX           /*V1.04-b*/
#define DEFAULT_BTAHEUSTEP   0.01                /*V1.04-b*/
#define DEFAULT_BTAHEUMAX    0.5                /*V1.04-b*/
#define DEFAULT_BTAHEUDDROP  0.8                /*V1.04-b*/
#define DEFAULT_BTAHEUDLOSS  0.5                /*V1.04-b*/
#define DEFAULT_BTAHEULLOSS  0.02               /*V1.04-b*/
#define DEFAULT_BTAGRADNIT   1                  /*V1.06-f*/
#define DEFAULT_BTAGRADCVTH  0.001              /*V1.06-f*/
#define DEFAULT_BTAGRADSTEP  0.01                /*V1.06-f*/
#define DEFAULT_BTAGRADRAND  0              /*V1.06-f*/
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
#define DEFAULT_NBRANDINITS  50                 /*V1.06-h*/
#define DEFAULT_ORDER        ORDER_DIRECT       /*V1.04-f*/
#define DEFAULT_UPDATE       UPDATE_SEQ         /*V1.06-d*/
#define DEFAULT_TIE          TIE_RANDOM         /*V1.06-e*/

#include <stdio.h>   /* FILE */
FILE* out_stderr;

/*
 *  Enumerated types
 */

typedef enum { 
        STS_OK = 0,
        STS_I_DONE,
        STS_W_EMPTYCLASS,
        STS_E_ARG,
        STS_E_MEMORY, 
        STS_E_FILEIN,
        STS_E_FILEOUT,
        STS_E_FILE,
        STS_E_FUNCARG,
        STS_E_INFINITE
        } 
        StatusET ;


typedef enum { 
        ALGO_NEM ,    
        ALGO_NCEM ,
        ALGO_GEM ,
        ALGO_NB
        } 
        AlgoET;


typedef enum { 
        BETA_FIX ,        /* Use given fixed beta */
        BETA_PSGRAD ,     /* Estimate beta using pseudo-likelihood gradient */
        BETA_HEUD,        /* Estimate beta using Hathaway heuristic */  
        BETA_HEUL,        /* Estimate beta using likelihood heuristic */  
        BETA_NB
        } 
        BetaET;           /* Beta estimation mode */ /*V1.04-a*/


typedef enum { 
        CRIT_U ,    /* Use NEM criterion */
        CRIT_M ,    /* Use markovian fuzzy log pseudo-likelihood */ /*V1.05-d*/
        CRIT_D ,    /* Use Hathaway criterion */
        CRIT_L ,    /* Use Mixture log-likelihood */
        CRIT_NB  
        }
        CritET;     /* Which criterion to choose local max */ /*V1.04-h*/



typedef enum { 
        TYPE_SPATIAL, 
        TYPE_IMAGE, 
        TYPE_NONSPATIAL,
        TYPE_NB
        } 
        TypeET;


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

/*V1.06-a*/
#if 0
typedef enum {
	   MODEL_P_VkI,
	   MODEL_PkVkI,
	   MODEL_P_V_I,
	   MODEL_PkV_I,
	   MODEL_NB
	   }
	   ModelET ;
#endif

typedef enum { 
  FAMILY_NORMAL,
  FAMILY_LAPLACE,
  FAMILY_BERNOULLI,
  FAMILY_NB
} 
FamilyET;

typedef enum { 
  DISPER___,  /* 1 dispersion for all classes and variables (V.I) */
  DISPER_K_,  /* 1 dispersion for each class, same in all variables (Vk.I) */
  DISPER__D,  /* 1 dispersion for each variable, same in all classes (V.B) */
  DISPER_KD,  /* 1 dispersion for each class and variable (VkBk) */
  DISPER_NB
} 
DisperET;

typedef enum { 
  PROPOR__,  /* equal proportion for all classes = 1/K (P) */
  PROPOR_K,  /* 1 proportion for each class (Pk) */
  PROPOR_NB
} 
ProporET;


typedef enum {
        FORMAT_HARD,
        FORMAT_FUZZY,
        FORMAT_NB
        }
        FormET ;

typedef enum {
        INIT_SORT,
        INIT_RANDOM,
        INIT_PARAM_FILE,      /* EM mixture estimate at start */
        INIT_FILE,
        INIT_LABEL,
        INIT_NB
}
        InitET ;

typedef enum
{
  MISSING_REPLACE,  /* Replace missing statistics with expectation as in EM */
  MISSING_IGNORE,   /* Ignore missing statistics as in CEM */
  MISSING_NB
}
MissET ;   /*V1.05-a*/

typedef enum
{
  PARAM_FILE_FIX,  /* EM mixture parameters fixed at start an during all the clustering process */
  PARAM_FILE_INIT,   /* EM mixture parameters fixed at start */
  NO_PARAM_FILE,  /* EM mixture parameters not fixed at start */
}
ParamFileET ;   /*V1.08-a*/

typedef enum {
        NEIGH_FOUR,
        NEIGH_FILE,
        NEIGH_NB
        }
        NeighET ;

typedef enum {
        ORDER_DIRECT,
        ORDER_RANDOM,
        ORDER_NB
        }
        OrderET ;       /*V1.04-e*/


typedef enum {
        UPDATE_SEQ,
        UPDATE_PARA,
        UPDATE_NB
        }
        UpdET ;         /*V1.06-e*/


typedef enum {
        TIE_RANDOM,
        TIE_FIRST,
        TIE_NB
        }
        TieET ;         /*V1.06-f*/


typedef enum {
        CVTEST_NONE,
        CVTEST_CLAS,
	CVTEST_CRIT,
        CVTEST_NB
        }
        CvemET ;        /*V1.06-g*/



/*
 *  Structured types
 */

typedef struct
{
    int         NbPts ;       /* number of observation vectors */
    int         NbVars ;      /* number of variables */
    int         NbMiss ;      /* number of missing data 0..Npts*NbVars */
                              /*V1.05-a*/
    float       *PointsM ;    /* observations (NbPts,NbVars) to allocate */
    int         *LabelV ;     /* fixed labels (NbPts) to allocate: 0..k */
    int         *SiteVisitV ; /* site to visit (NbPts) to allocate: 0..Npts-1*/
    int         *SortPos_ND ; /* PointsM[ SortPos_ND[i*D+d]*D+d ] : +++*/
}
DataT ;         /* Matrix of observed data (each line = 1 point) */


/*V1.06-g*/
typedef struct
{
  int     NbIter ;     /* Max number of iterations of gradient ascent */
  float   ConvThres ;  /* Convergence if gradient <  threshold * N */
  float   Step ;       /* >0 : bta += grad*(step/N), 0: bta += grad/dsec */
  int     RandInit ;   /* 1 = random initial beta, 0 = specified by -b */
}
BtaPsGradT ;   /* parameters of beta pseudo-likelihood gradient estimation */


typedef struct
{
    AlgoET  Algo ;      /* Type of algorithm */
  /*V1.06-a*/
#if 0
     float   Beta ;      /* context weight : >= 0 */
     BetaET  BtaMode ;   /* type of beta estimation */  /*V1.04-b*/
#endif
    float   BtaHeuStep ;/* step of beta for heuristic estimation */ /*V1.04-b*/
    float   BtaHeuMax ; /* maximum beta for heuristic */
    float   BtaHeuDDrop ;/* drop of Hathaway slope for beta heuristic */
    float   BtaHeuDLoss ;/* proportion of Hathaway loss for beta heuristic */
    float   BtaHeuLLoss ;/* proportion of likelihood loss for beta heuristic */
    BtaPsGradT BtaPsGrad ; /* parameters of beta gradient estimation */
    CritET  Crit ;      /* criterion to choose local max */ /*V1.04-h*/
    float   CvThres ;   /* convergence threshold */    /*V1.04-c*/
    CvemET  CvTest ;    /* which convergence test to use */
    int     DoLog ;     /* TRUE if log file requested */
    int     NbIters ;   /* nb of iterations for NEM */
    int     NbEIters ;  /* nb of iterations for E-step */
    int     NbRandomInits ;  /* nb of random initializations */
    long    Seed ;      /* random generator seed */   /*V1.04-d*/
    FormET  Format ;    /* output file format (hard or fuzzy) */
    InitET  InitMode ;  /* initialization mode (histogram, random, file) */
    MissET  MissMode ;  /* how to process missing statistics */ /*V1.05-a*/
    ParamFileET ParamFileMode ; /* parameters used at begin or throughout the clustering process*/ /*V1.08-a*/
    int     SortedVar ; /* variable to be sorted : 0..NbVars */
    NeighET NeighSpec ; /* neighborhood specification */
    OrderET VisitOrder ;/* order of visit at E-step */ /*V1.04-e*/
    UpdET   SiteUpdate ;/* site update scheme at E-step */
    TieET   TieRule ;   /* rule for equal probabilities when computing MAP */
    int     Debug ;     /* TRUE if in debug mode */    /*V1.04-f*/
    char    OutBaseName[ LEN_FILENAME + 1 ] ; /* base name of output file */
    char    OutName[ LEN_FILENAME + 1 ] ;   /* name of output file */
    char    LogName[ LEN_FILENAME + 1 ] ;   /* name of log file ("" = no log) */
    char    StartName[ LEN_FILENAME + 1 ] ; /* name of initial partition file */
    char    NeighName[ LEN_FILENAME + 1 ] ; /* name of neighborhood file */
    char    LabelName[ LEN_FILENAME + 1 ] ; /* name of fixed labels file */
    char    RefName[ LEN_FILENAME + 1 ] ;   /* name of reference labels file *//*V1.04-f*/
    char    ParamName[ LEN_FILENAME + 1 ] ; /* name of initialization param file *//*V1.08-a*/
}
NemParaT ;      /* NEM running parameters */


typedef struct
{
    int     Dl ;    /* neighbour shift in line :     -2 .. 2 */
    int     Dc ;    /* neighbour shift in column :   -2 .. 2 */
    float   Weight ; /* neighbour weight >= 0 */
}
INeighT ;       /* one neighbour (in image configurations) */

typedef struct
{
    int     Nl ;        /* Image number of lines (height) : > 0 */
    int     Nc ;        /* Image number of columns (width) : > 0 */
    int     NbNeigh ;   /* nb of allocated neighbours : >= 0 */
    INeighT *NeighsV ;  /* to be allocated : pixel's neighbours */
}
ImageNeighT ;   /* neighbourhood system (in image configurations) */

typedef struct
{
    int         Index ; /* index of neighbour : 0 .. Npt-1 */
    float       Weight ;/* weight of neighbour (default : 1.0) */
}
NeighT ;        /* one neighbour (in non-image configurations) */

typedef struct
{
    int     NbNeigh ;   /* nb of allocated neighbours */
    NeighT  *NeighsV ;  /* to be allocated : point's neighbours */
}
PtNeighsT;      /* one point's neighbours (in non-image configuration) */

typedef union
{
    ImageNeighT Image ;
    PtNeighsT   *PtsNeighsV ;   /* to be allocated : all points' neighbours */
}
NeighDataT ;    /* generic neighbourhood system */

typedef struct
{
    NeighDataT  NeighData ; /* generic neighbourhood system */
    int         MaxNeighs ; /* maximum number of neighbours (>= 0) */
    TypeET      Type ;      /* type of spatial configuration */
}
SpatialT ;


/*V1.06-a*/

#if 0
typedef struct
{
    float   *Pk ;   /* proportions (K) */   /* to be allocated */
    float   *Vk ;   /* volumes (K) */       /* to be allocated */
    float   *Ck ;   /* shapes (d*d*K) */    /* to be allocated */
    float   *Mk ;   /* means (d*K) */       /* to be allocated */
}
NoiseParaT ;


typedef struct
{
    ModelET     ModelNum ;
    int         Nk ;
    NoiseParaT  NoisePara ;
}
NoiseModelT ;
#endif

typedef struct
{
  int       K ;                /* number of classes */
  FamilyET  ClassFamily ;
  DisperET  ClassDisper ;
  ProporET  ClassPropor ;
  BetaET    BetaModel ;
}
ModelSpecT ;  /* Model specification */

typedef struct
{
  float     Beta ;
  float*    Center_KD ;  /* Center in each class and variable (K*D) */
  float*    Disp_KD ;    /* Dispersion in each class and variable (K*D) */
  float*    Prop_K ;     /* Proportion of each class (in ]0,1[) */

  float*    NbObs_K ;    /* Nb of observations in each class (K) */
  float*    NbObs_KD ;   /* Nb of observations in each class/variable (K*D) */
  float*    Iner_KD ;    /* Inertia = sum_i cik * Dist(xid, mkd) (K*D) */
}
ModelParaT ;  /* Model parameters */

typedef struct
{
  float*     DispSam_D ; /* Dispersion of whole sample in each variable (D) */
  float*     MiniSam_D ; /* Minimum of whole sample in each variable (D) */
  float*     MaxiSam_D ; /* Maximum of whole sample in each variable (D) */
}
SampleDesT ;  /* Sample description */

typedef struct
{
  ModelSpecT   Spec ;
  ModelParaT   Para ;
  SampleDesT   Desc ;
}
StatModelT ;  /* Model description */


typedef struct
{
  int      Kc ;            /* # of classes in computed classification */
  int      Kr ;            /* # of classes in reference classification */
  int      Km ;            /* max( Kr, Kc ) */
  int      Kmfac ;         /* Km! */
  TieET    TieRule ;       /* same value as in NemPara */
  float*   Refclas_N_Kr ;  /* reference classification : 0/1 */
  int*     Perm_Kmfac_Km ; /* all permutation of Km classes */
}
ErrinfoT ;    /* Classification error information */


typedef struct
{
  float*   Agree_Km_Km ;   /* #common elements between found and ref classes */
  float*   Loclas_N_Kc ;   /* local copy of computed classification */
  int      Ibestpermut ;   /* index of best agreement permutation 0..Km!-1 */
  float    Errorrate ;     /* # misclassified objects of best permut / N */
}
ErrcurT ;    /* Currently computed classification error */


typedef struct
{
  float    D ; /* hathaway crit. D = sum[i]sum[k] cik (log pkfki-log cik) */
  float    G ; /* geog. cohesion G = sum[i]sum[j]sum[k] cik cjk wij */
  float    U ; /* NEM maximized criterion U = D + 0.5 * beta * G */
  float    M ; /* markovian fuzzy class. like. M = D + beta * G - Z */
  float    L ; /* mixture likelihood crit. L = sum[i] log sum[k] pkfki ) */
  float    Z ; /* log pseudo-l. Z =-sum[i]log(sum[k]e(bta*sum[j~i]wij cjk)) */
  ErrinfoT Errinfo ; /* information to compute error */   /*V1.06-h*/
  ErrcurT  Errcur ;  /* current error rate */   /*V1.06-h*/
} /*V1.05-d*/
CriterT ;       /*V1.03-a*/




typedef StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
EstimNoiseFT          
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  StatModelT    *StatModelP /* O : estimated parameters */
 ) ;

typedef int CompuDensFT         /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,      /* I : point dimension */
            int                Ik,      /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP,   /* I : noise parameters *//*V1.06-a*/
            const float*       XV,      /* I : point (dim d) */
            double*            FkP,     /* O : density for class Ik */
            float*             LogFkP   /* O : log of density */
        ) ;


typedef int GetNeighFT         /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) ;





static const char *TypeDesC[ TYPE_NB ] = { "Spatial", "Image", "NoSpatial" } ;
static const char *AlgoDesC[ ALGO_NB ] = { "NEM" , "NCEM (C-step)" , 
             "GEM (Monte-Carlo at E-step)" } ;
static const char *BetaDesVC[ BETA_NB ] = { "fixed", 
              "pseudo-likelihood gradient",
              "heuristic Hathaway crit",
              "heuristic mixture likelihood"} ;
static const char *FamilyDesVC[ FAMILY_NB ] = { "Normal", "Laplace", 
            "Bernoulli" } ;
static const char *DisperDesVC[ DISPER_NB ] = { "S__", "SK_", "S_D", "S_KD" } ;
static const char *ProporDesVC[ PROPOR_NB ] = { "P_", "Pk" } ;
static const char *InitDesC[ INIT_NB ] = { "Sort_variables" ,
             "Random_centers" , 
             "Param_file" ,
             "Full_labels" , 
             "Partial_labels" } ;    /*V1.03-f*/

static const char *TypeStrC[ TYPE_NB ] = { "S", "I", "N" } ; /*V1.06-a*/

static const char *AlgoStrVC[ ALGO_NB ] = { "nem" , "ncem" , "gem" } ;
static const char *BtaStrVC[ BETA_NB ] = { "fix" , "psgrad", 
                       "heu_d" , "heu_l" } ;
static const char *CritStrVC[ CRIT_NB ] = { "U" , "M", "D" , "L" } ;
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


static const char *OrderDesC[ ORDER_NB ] = { "Direct_order" , "Random_order" };



#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
