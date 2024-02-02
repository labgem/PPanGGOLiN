/*\

    NEM_MOD.C

    Programme NEM (Neighborhood EM) : routines calcul parametres classes

    Van Mo DANG       Janvier 96


    Vers-mod  Date         Description

    1.03-a    01-NOV-1996  Add function SetIdMatrix + calls by ParaPk.. 
                           to output variance matrix in file don.mf
    1.04-a    08-JAN-1998  Fixed error in DensPkVkI log(vk) needs * Nd
    1.05-a    12-JAN-1998  MissMode in ParaP*V*I
    1.05-b    12-JAN-1998  Process missing data in DensPkVkI
    1.05-c    12-JAN-1998  Add ParaPkV*I_Missing() (empty)
    1.05-d    15-JAN-1998  Add CommonGaussDiagMissing() <= ParaPkV*I_Missing()
    1.05-e    26-JAN-1998  Call to GenAlloc
    1.06-a    29-JUN-1998  Reflect change to nem_mod.h
    1.06-b    29-JUN-1998  Reflect change to nem_typ.h
    1.06-c    29-JUN-1998  Compute Para.NbObs_KD
    1.06-d    02-JUL-1998  Bug: CommonGaus/replace, sum(x-m)^2 neq sum(x^2)-m^2
    1.06-e    02-JUL-1998  Bug: in VkI_Missing/replace, div inertia by nk*D
    1.06-f    02-JUL-1998  Add EstimParaLaplace and called subfunctions
    1.06-g    01-DEC-1998  FkP double* instead of *float in dens...
    1.07-a    26-FEB-1999  Add FAMILY_BERNOULLI in GetDensityFunc and EstimPara
    1.07-b    26-FEB-1999  Add DensBernoulli
    1.07-c    03-MAR-1999  Fix bug DensBernoulli: disp==0 may give nonzero dens
\*/

#include "genmemo.h"    /* GenAlloc */
#include "nem_typ.h"    /* DataT, ... */
#include "nem_mod.h"    /* ParaP_V_I, ... */
#include <stdio.h>      /* printf, ... */
#include <stdlib.h>     /* malloc, ... */
#include <string.h>     /* memcpy, ... */
#include <math.h>       /* exp, ... */
#include <float.h>      /* FLT_MAX */

#define TWO_PI    2 * 3.14159 

/*V1.05-d*/
#define _IJ       ( ( i * D ) + j )    /* access a (.,D) matrix by (i,j) */
#define _HJ       ( ( h * D ) + j )    /* access a (.,D) matrix by (h,j) */
#define _IH       ( ( i * K ) + h )    /* access a (.,K) matrix by (i,h) */
#define sqr(x)    ((x)*(x))            /* macro for x^2 */

#ifndef MAXFLOAT
 #ifdef FLT_MAX
   #define MAXFLOAT FLT_MAX
 #else
   #define MAXFLOAT 3.40282347e+38F
 #endif
#endif



/* ==================== LOCAL FUNCTION PROTOTYPING =================== */


/* Indirect call by GetDensityFunc() */

static int DensNormalDiag      /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,         /* I : point dimension */
            int                Ik,         /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP ,     /* I : model parameters */
            const float*       XV,         /* I : point (dim d) */
            double*            FkP,        /* O : density for class Ik */
            float*             LogFkP      /* O : log of density */
        ) ;


static int DensLaplaceDiag      /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,         /* I : point dimension */
            int                Ik,         /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP ,     /* I : model parameters */
            const float*       XV,         /* I : point (dim d) */
            double*            FkP,        /* O : density for class Ik */
            float*             LogFkP      /* O : log of density */
        ) ;


static int DensBernoulli        /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,         /* I : point dimension */
            int                Ik,         /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP ,     /* I : model parameters */
            const float*       XV,         /* I : point (dim d) */
            double*            FkP,        /* O : density for class Ik */
            float*             LogFkP      /* O : log of density */
        ) ;


/* Called by EstimPara() */

static StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY */ /*V1.06-b*/
EstimParaNormal 
(
  const float   *C_NK,      /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           Nk,         /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/
  const ModelSpecT* SpecP,  /* I : model specification */ /*V1.06-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 */ /*V1.05-b*/
  ModelParaT    *ParaP      /* O : estimated parameters */
) ;




static StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY */ /*V1.06-b*/
EstimParaLaplace
(
  const float   *C_NK,      /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           Nk,         /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/
  const ModelSpecT* SpecP,  /* I : model specification */ /*V1.06-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 */ /*V1.05-b*/
  ModelParaT    *ParaP      /* O : estimated parameters */
) ;




/* Called by EstimParaNormal() */

/*V1.05-d*/
static StatusET          /* Return status : OK, EMPTY or MEMORY */
CommonGaussDiag
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  int           N,             /* I : nb of observation vectors */
  int           D,             /* I : nb of variables */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  int           K,             /* I : number of clusters */
  MissET        Miss,          /* I : how to treat missing data */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */

  float*        OldNewCen_KD,  /* I/O : old then updated means (K,D) */

  int*          EmptyK_P,      /* O : index of empty class (0 or 1 ..K) */
  float*        N_K,           /* O : size of a class (K) */
  float*        N_KD,          /* O : size of a class and variable (K,D) */
  float*        Iner_KD        /* O : inertia of a class/variable (K,D) */
 ) ; 


static void InerToDisp
(
 DisperET      DispType, /* I : dispersion model */
 int           N,        /* I : nb of observation vectors */
 int           Nk,       /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD,  /* O : dispersion in each class/variable */
 StatusET*     StsP      /* [O] : STS_E_FUNCARG if unknown DispType */
) ;



/* Called by InerToDisp() */

static void InerToDisp__
(
 int           N,        /* I : nb of observation vectors */
 int           Nk,       /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
) ;

static void InerToDispK_
(
 int           N,        /* I : nb of observation vectors */
 int           K,        /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
) ;


static void InerToDisp_D
(
 int           N,        /* I : nb of observation vectors */
 int           K,        /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
) ;


static void InerToDispKD
(
 int           N,        /* I : nb of observation vectors */
 int           K,        /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
) ;



/* Called by EstimParaLaplace() */


/*V1.06-f*/
static StatusET          /* Return status : OK, EMPTY or MEMORY */
CommonLaplaceDiag
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  int           N,             /* I : nb of observation vectors */
  int           D,             /* I : nb of variables */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  int           K,             /* I : number of clusters */
  MissET        Miss,          /* I : how to treat missing data */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */

  float*        OldNewCen_KD,  /* I/O : old then updated means (K,D) */

  int*          EmptyK_P,      /* O : index of empty class (0 or 1 ..K) */
  float*        N_K,           /* O : size of a class (K) */
  float*        N_KD,          /* O : size of a class and variable (K,D) */
  float*        Iner_KD        /* O : inertia of a class/variable (K,D) */
 ) ; 


#if 0 /* already declared above */

static void InerToDisp
(
 DisperET      DispType, /* I : dispersion model */
 int           N,        /* I : nb of observation vectors */
 int           Nk,       /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
) ;

#endif


/* Called by CommonLaplaceDiag() */

static void EstimSizes
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  float*        N_K,           /* O : size of a class (K) */
  float*        N_KD           /* O : size of a class and variable (K,D) */  
 ) ;


static StatusET          /* Return status : OK, EMPTY */
EstimLaplaceCenters
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const MissET  Miss,          /* I : how to treat missing data */
  const float*  N_K,           /* I : size of a class (K) */
  const float*  N_KD,          /* I : size of a class and variable (K,D) */  
  const float*  OldCen_KD,     /* I : old centers */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */

  float*        OldNewCen_KD,  /* I/O : eventually updated centers (K,D) */
  int*          EmptyK_P       /* O : index of empty class (0 or 1 ..K) */
 ) ;


static void EstimLaplaceIner
 ( 
  const float*  X_ND,          /* I : data matrix (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const MissET  Miss,          /* I : how to treat missing data */
  const float*  N_K,           /* I : size of a class (K) */
  const float*  N_KD,          /* I : size of a class and variable (K,D) */  
  const float*  OldCen_KD,     /* I : old centers */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */
  const float*  NewCen_KD,     /* I : new centers */
  float*        Iner_KD        /* O : inertia of a class/variable (K,D) */  
 ) ;



/* Called by EstimLaplaceCenters() */

static void ComputeMedian
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const int     H,             /* I : current cluster */
  const int     J,             /* I : current variable */
  const float   totwei,        /* I : total weight */
  int*          ImedP,         /* O : index of median of observed values */
  float*        CumweiP,       /* O : cumulated weights until median obs. */
  float*        MedvalP        /* O : median value (eventually midway) */
 ) ;



static float FindMinInerLaplaceEM  /* ret: minimizer of expected inertium */
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const int     H,             /* I : current cluster */
  const int     J,             /* I : current variable */
  const int     Imed,          /* I : median position of observed data */
  const float   Cumwei,        /* I : cumulated weight until Xmed */
  const float   Cen0,          /* I : previous center */
  const float   Disp0,         /* I : previous dispersion */
  const float   Nhj,           /* I : number of observations */
  const float   Nmis           /* I : number of missing data */
 ) ;



/* Called by FindMinInerLaplaceEM() */

static float DerivInerDir      /* ret: +/- R'(Y), R expected inertia */ 
 (
  const float   Y,             /* I : point at which to compute R' */
  const float   Weidif,        /* I : "after+med-bef" or "bef+med-after" */
  const float   Intwei,        /* I : cumulated intermediate weight */
  const float   Nmis,          /* I : number of missing data */
  const float   Cen0,          /* I : previous center */
  const float   Disp0          /* I : previous dispersion */
 ) ;




/* ==================== GLOBAL FUNCTION DEFINITION =================== */


/* ------------------------------------------------------------------- */
int GetDensityFunc  /* STS_OK or STS_E_FUNCARG */
        (
            const ModelSpecT  *SpecP,           /* I */
            CompuDensFT**     CompuDensFP       /* O */
        )
/* ------------------------------------------------------------------- */
{

    switch( SpecP->ClassFamily )
    {
        case FAMILY_NORMAL:
	  *CompuDensFP = DensNormalDiag ;
	  return STS_OK ;

        case FAMILY_LAPLACE:
	  *CompuDensFP = DensLaplaceDiag ;
	  return STS_OK ;

        case FAMILY_BERNOULLI:  /*V1.07-a*/
	  *CompuDensFP = DensBernoulli ;
	  return STS_OK ;

        default :
	  *CompuDensFP = NULL ;
	  fprintf( stderr, "GetDensityFunc bad arg : family = %d\n",
                   SpecP->ClassFamily ) ;
	  return STS_E_FUNCARG ;
    }

  /*???*/
}   /* end of GetDensityFunc() */




/* ------------------------------------------------------------------- */
StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY */ /*V1.06-b*/
EstimPara 
(
  const float   *C_NK,      /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           Nk,         /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/
  const ModelSpecT* SpecP,  /* I : model specification */ /*V1.06-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 */ /*V1.05-b*/
  ModelParaT    *ParaP      /* I/O : previous and new estimated parameters */
)
/* ------------------------------------------------------------------- */
{
  StatusET    sts ;  /* return status */

  int         k ;

  /* Family dependent estimation method 
   */
  switch( SpecP->ClassFamily ) {
  case FAMILY_NORMAL:
    sts = EstimParaNormal( C_NK, DataP, Nk, MissMode, SpecP, 
			   EmptyK_P, ParaP ) ;
    break ;

  case FAMILY_LAPLACE:
    sts = EstimParaLaplace( C_NK, DataP, Nk, MissMode, SpecP, 
			    EmptyK_P, ParaP ) ;      
    break;

  case FAMILY_BERNOULLI: /*+++DANGER: MissMode forced to MISSING_IGNORE+++*/
    sts = EstimParaLaplace( C_NK, DataP, Nk, MISSING_IGNORE, SpecP, 
			    EmptyK_P, ParaP ) ;  /*V1.07-a*/
    break;

  default:
    sts = STS_E_FUNCARG ;
  }
  
  /* Estimate proportions */
  if ( SpecP->ClassPropor == PROPOR_K ) {

    for ( k = 0; k < Nk ; k ++ )
      ParaP->Prop_K[ k ] = ParaP->NbObs_K[ k ] / DataP->NbPts ;
  }
  else {
    
    for ( k = 0; k < Nk ; k ++ )
      ParaP->Prop_K[ k ] = 1.0 / Nk ;
  }

  return sts ;
  /*???*/
}   /* end of EstimPara() */




/* ==================== LOCAL FUNCTION DEFINITION =================== */





/* ------------------------------------------------------------------- */
static int DensNormalDiag      /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,         /* I : point dimension */
            int                Ik,         /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP ,     /* I : noise parameters */
            const float*       XV,         /* I : point (dim d) */
            double*            FkP,        /* O : density for class Ik */
            float*             LogFkP      /* O : log of density */
        )
/* ------------------------------------------------------------------- */
{
  int     d ;      /* 0..Nd-1 : current variable */
  float   dk ;     /* sum_d [ log(2pi skd2) + (xd-mkd)^2/skd2 ] */
  int     nbobs ;  /* 0..Nd-1 : nb of observed variables */ /*V1.05-b*/
  int     nuldisp; /* TRUE if one of the dispersions is zero */

  /* Gaussian mixture density, diagonal models :

     Complete data :
     fk(x) = [ prod_d (2*pi * vkd)^(-0.5) ] * exp( -0.5 * dist_x_mk )
     where  dist_x_mk = sum[d=1,Nd]( (x(d) - mk(d))^2 ) / vkd

     Incomplete data :
     fk(x) = [ prod_{d in oi} (2*pi * vkd)^(-0.5) ] * exp( -0.5 * dist_x_mk )
     where   dist_x_mk = sum[d in oi]  (x(d) - mk(d))^2 / vk
     =>
     log fk(x) = - 0.5 sum_{d in oi}[ log(2*pi * vkd) + (xd-mkd)^2/vkd ]

  */

  for ( d = 0, 
          dk = 0.0, nbobs = 0, nuldisp = 0 ; 
        d < Nd ; 
        d ++ )
    {
      if ( ! isnan( XV[ d ] ) )   /*V1.05-b*/
        {
          float dif = XV[ d ] - ParaP->Center_KD[ (Ik * Nd) + d ] ;

          if ( ParaP->Disp_KD[ (Ik * Nd) + d ] > EPSILON )
            dk = dk 
	      + log( TWO_PI * ParaP->Disp_KD[ (Ik * Nd) + d ] ) 
	      + (dif * dif) / ParaP->Disp_KD[ (Ik * Nd) + d ] ;
          else
            nuldisp = 1 ;

          nbobs ++ ;
        }
  }

  if ( ! nuldisp )
    {
      /*V1.04-a*//*V1.05-b*/
      *LogFkP = -0.5 * dk ;
      *FkP = exp( *LogFkP ) ;
      return 0 ;
    }
  else
    {
      *LogFkP = - MAXFLOAT ;
      *FkP = 0.0 ;
      return -1 ;
    }


}   /* end of DensNormalDiag() */



/* ------------------------------------------------------------------- */
static int DensLaplaceDiag      /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,         /* I : point dimension */
            int                Ik,         /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP ,     /* I : noise parameters */
            const float*       XV,         /* I : point (dim d) */
            double*            FkP,        /* O : density for class Ik */
            float*             LogFkP      /* O : log of density */
        )
/* ------------------------------------------------------------------- */
{
  int     d ;      /* 0..Nd-1 : current variable */
  float   dk ;     /* sum_d [ log(2 lkd) + |xd-mkd| / lkd ] */
  int     nbobs ;  /* 0..Nd-1 : nb of observed variables */ /*V1.05-b*/
  int     nuldisp; /* TRUE if one of the dispersions is zero */

  /* Laplace mixture density, diagonal models :

     Complete data :
     fk(x) = [ prod_d (2*vkd)^(-1) ] * exp( - dist_x_mk )
     where  dist_x_mk = sum[d=1,Nd] | x(d) - mk(d) | / vkd

     Incomplete data :
     fk(x) = [ prod_{d in oi} (2*vkd)^(-1) ] * exp( - dist_x_mk )
     where   dist_x_mk = sum[d in oi]  | x(d) - mk(d) | / vkd
     =>
     log fk(x) = - sum_{d in oi}[ log(2*vkd) + |xd-mkd| / vkd ]

  */

  for ( d = 0, 
          dk = 0.0, nbobs = 0, nuldisp = 0 ; 
        d < Nd ; 
        d ++ )
    {
      if ( ! isnan( XV[ d ] ) )   /*V1.05-b*/
        {
          float dif = XV[ d ] - ParaP->Center_KD[ (Ik * Nd) + d ] ;

          if ( ParaP->Disp_KD[ (Ik * Nd) + d ] > EPSILON )
            dk = dk 
      	+ log( 2 * ParaP->Disp_KD[ (Ik * Nd) + d ] ) 
      	+ fabs( dif ) / ParaP->Disp_KD[ (Ik * Nd) + d ] ;
          else
            nuldisp = 1 ;

          nbobs ++ ;
        }
  }

  if ( ! nuldisp )
    {
      /*V1.04-a*//*V1.05-b*/
      *LogFkP = - dk ;
      *FkP = exp( *LogFkP ) ;
      return 0 ;
    }
  else
    {
      *LogFkP = - MAXFLOAT ;
      *FkP = 0.0 ;
      return -1 ;
    }

}   /* end of DensLaplaceDiag() */


/* ------------------------------------------------------------------- */
static int DensBernoulli      /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,         /* I : point dimension */
            int                Ik,         /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP ,     /* I : noise parameters */
            const float*       XV,         /* I : point (dim d) */
            double*            FkP,        /* O : density for class Ik */
            float*             LogFkP      /* O : log of density */
        )
/* ------------------------------------------------------------------- */
{
  int     d ;      /* 0..Nd-1 : current variable */
  float   dk ;     /* sum_d [ -log(1-vkd) + log{(1-vkd)/vkd} |xd-mkd| ] */
  int     nbobs ;  /* 0..Nd-1 : nb of observed variables */ /*V1.05-b*/
  int     nuldens; /* TRUE if the probability is zero */

  /* Bernoulli density, diagonal models :

     Complete data :
     fk(x) = prod_d v_kd^|x_d - m_kd| (1-v_kd)^( 1 - |x_d - m_kd| )

     Incomplete data :
     fk(x) = prod_{d in oi} v_kd^|x_d - m_kd| (1-v_kd)^( 1 - |x_d - m_kd| )
     =>
     log fk(x) = - sum_{d in oi} [ -log(1-vkd) + |xd-mkd| * log((1-vkd)/vkd) ]

     if there is a d such that vkd = 0 and |xd-mkd| != 0, prob=0.

  */

  for ( d = 0, 
          dk = 0.0, nbobs = 0, nuldens = 0 ; 
        d < Nd ; 
        d ++ ) {

      if ( ! isnan( XV[ d ] ) )   /*V1.05-b*/
        {
	  float disp = ParaP->Disp_KD[ (Ik * Nd) + d ] ;
          int   absdif = 
	    abs( (int) ( XV[ d ] - ParaP->Center_KD[ (Ik * Nd) + d ] ) ) ;

          if ( disp > EPSILON )
            dk = dk + absdif * log( ( 1 - disp ) / disp ) - log( 1 - disp ) ;
          else  /* null dispersion */ {

	    if ( absdif != 0 ) /* prob(xid != center) = 0 */ {
	      nuldens = 1 ;
	    }
	    else /* prob(xid = center) = 1 => no change to dk */ {
	      dk = dk + 0.0 - 0.0 ;
	    }
	  }

          nbobs ++ ;
        }
  }

  if ( ! nuldens )
    {
      /*V1.04-a*//*V1.05-b*/
      *LogFkP = - dk ;
      *FkP = exp( *LogFkP ) ;
      return 0 ;
    }
  else
    {
      *LogFkP = - MAXFLOAT ;
      *FkP = 0.0 ;
      return -1 ;
    }

}   /* end of DensBernoulli() */


/* ------------------------------------------------------------------- */
static StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY */ /*V1.06-b*/
EstimParaNormal 
(
  const float   *C_NK,      /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           Nk,         /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/
  const ModelSpecT* SpecP,  /* I : model specification */ /*V1.06-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 */ /*V1.05-b*/
  ModelParaT    *ParaP      /* O : estimated parameters */
)
/* ------------------------------------------------------------------- */
{
  int      N = DataP->NbPts ;
  int      D = DataP->NbVars ;

  StatusET sts ;


  /* Common computation for any diagonal gaussian model with missing data */
  sts = CommonGaussDiag
    ( DataP->PointsM, N, D, C_NK, Nk, MissMode, ParaP->Disp_KD, 
      ParaP->Center_KD, 
      EmptyK_P, ParaP->NbObs_K, ParaP->NbObs_KD, ParaP->Iner_KD ) ;

  /* Dispersion is computed depending on dispersion model */
  InerToDisp( SpecP->ClassDisper, N, Nk, D, ParaP->NbObs_K, ParaP->NbObs_KD, 
	      ParaP->Iner_KD, MissMode, 
	      ParaP->Disp_KD, & sts ) ;

  return sts ;

}   /* end of EstimParaNormal() */



/* ------------------------------------------------------------------- */
static StatusET          /* Return status : OK, EMPTY or MEMORY */
CommonGaussDiag
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  int           N,             /* I : nb of observation vectors */
  int           D,             /* I : nb of variables */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  int           K,             /* I : number of clusters */
  MissET        Miss,          /* I : how to treat missing data */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */

  float*        OldNewCen_KD,  /* I/O : old then updated means (K,D) */

  int*          EmptyK_P,      /* O : index of empty class (0 or 1 ..K) */
  float*        N_K,           /* O : size of a class (K) */
  float*        N_KD,          /* O : size of a class and variable (K,D) */
  float*        Iner_KD        /* O : inertia of a class/variable (K,D) */
 ) 
/*\

   This function does the common computation that is needed at the
   M-step, in all the diagonal Gaussian models with missing data.

   The routine takes as input a data matrix Xij_nd (possibly
   containing NaN values), a classification matrix Cih_nk, and the
   kind of missing data processing Miss (REPLACE or IGNORE) ; and also
   the old volumes OldVh_k, covariance matrices OldChjj_kdd, and means
   OldNewMhj_kd of the classes, to be used in the REPLACE mode.

   It produces as output the reestimated mean vectors in OldNewMhj_kd;
   reestimated proportions in NewP_k; the last encountered empty class
   in EmptyK_P (or 0 if no empty class); the {sum[i] cih} in Nh_k, the
   {sum[i] cih rij} in Nhj_kd (number of present observations in class h
   and variable j); and the "missing data dispersion" of each
   class/variable in Disphj_kd (i.e. the missing data equivalent of
   {sum[i] Cih (xij - mhj)^2}).

\*/
/* ------------------------------------------------------------------- */
{
  const char* func = "CommonGaussDiag" ;

  float*   sumdata_KD ;    /* sum[i] cih rij xij   (K,D), local alloc */
  float*   sumsquare_KD ;  /* sum[i] cih rij xij^2 (K,D), local alloc */
  float*   oldmean_KD ;    /* saved previous means (K,D), local alloc */

  int      h ;             /* current class    0..K-1 */
  int      j ;             /* current variable 0..D-1 */
  int      i ;             /* current object   0..N-1 */
  StatusET sts ;           /* return status */

  float    obsinerhj ;


  sts = STS_OK ;

  /* Allocate local structures */
  sumdata_KD   = GenAlloc( K * D, sizeof(float), 1, func, "sumdata" ) ;
  sumsquare_KD = GenAlloc( K * D, sizeof(float), 1, func, "sumsquare" );
  oldmean_KD   = GenAlloc( K * D, sizeof(float), 1, func, "oldmean" ) ;

  /* Save old mean value */
  memcpy( oldmean_KD , OldNewCen_KD , K * D * sizeof( float ) ) ;

  /* Set no empty class by default */
  (*EmptyK_P) = 0 ;

  /* For each class h */
  for ( h = 0 ; h < K ; h ++ ) {

      /* For each variable j */
      for ( j = 0 ; j < D ; j ++ ) {

	  /* Initialize the sum[i] quantities to 0 */
	  N_K[ h ]             = 0.0 ;
	  N_KD[ _HJ ]          = 0.0 ;
	  sumdata_KD[ _HJ ]    = 0.0 ;
	  sumsquare_KD[ _HJ ]  = 0.0 ;

	  /* For each object i */
	  for ( i = 0 ; i < N ; i ++ ) {

	      float cih = C_NK[ _IH ] ;
	      float xij = X_ND[ _IJ ] ;

	      /* Increment class size */
	      N_K[ h ] += cih ;

	      /* If this object's j variable is observed */
	      if ( ! isnan( xij ) )
		{
		  /* Increment all sum[i] of observed quantities */
		  N_KD[ _HJ ]          += cih ;
		  sumdata_KD[ _HJ ]    += cih * xij ;
		  sumsquare_KD[ _HJ ]  += cih * xij * xij ;
		}
	      /* Else xij is missing, don't change sum[i] of observed data */
	  }
	  /* End  For each object i (all sums[i] are now done) */

	  /* If this class size > 0 */
	  if ( N_K[ h ] > 0 ) {
	      /* 
	       * Compute means Mhj and dispersions Dhj from computed sums,
	       * depending on REPLACE or IGNORE missing mode 
	       */
	      if ( Miss == MISSING_REPLACE ) {
		  /* 
		   * new_mhj = 
		   *  ( sum[i] cih rij xij + ( nh - nhj ) old_mhj ) / nh 
		   *
		   * disp_hj = 
		   *  "observed dispersion"_hj +
		   *  ( nh - nhj ) [ (new_mhj - old_mhj)^2 + old_vh old_chjj ]
		   */

		  OldNewCen_KD[ _HJ ] = 
		    ( sumdata_KD[ _HJ ] + 
		      ( N_K[ h ] - N_KD[ _HJ ] ) * oldmean_KD[ _HJ ] )
		    / N_K[ h ] ;

		  /* V1.06-d
		   * "observed dispersion"_hj = 
		   *      sum[i] cih rij (xij - new_mhj)^2 =
		   *      sum[i] cih rij xij^2 - 
		   *       new_mhj * ( 2 * sum[i] cih rij xij - nhj new_mhj )
		   */
		  obsinerhj = sumsquare_KD[ _HJ ] - 
		    OldNewCen_KD[ _HJ ] * 
		    (2*sumdata_KD[ _HJ ] - N_KD[ _HJ ] * OldNewCen_KD[ _HJ ]) ;

		  Iner_KD[ _HJ ] = 
		    obsinerhj +
		    ( N_K[ h ] - N_KD[ _HJ ] ) *
		    ( sqr( OldNewCen_KD[ _HJ ] - oldmean_KD[ _HJ ] ) +
		      OldDisp_KD[ _HJ ] ) ;
	      }
	      else  { /* assume Miss == MISSING_IGNORE */

		  /* 
		   * new_mhj = 
		   *  . if  nhj <> 0,  ( sum[i] cih rij xij ) / nhj  
		   *  . if  nhj == 0,  old_mhj
		   *
		   * disp_hj = 
		   *  "observed dispersion"_hj 
		   */

		  if ( N_KD[ _HJ ] > 0 )
		    OldNewCen_KD[ _HJ ] = sumdata_KD[ _HJ ] / N_KD[ _HJ ];
		  else
		    OldNewCen_KD[ _HJ ] = oldmean_KD[ _HJ ] ;

		  /*
		   * "observed dispersion"_hj = 
		   *      sum[i] cih rij (xij - new_mhj)^2 =
		   *      sum[i] cih rij xij^2 - nhj new_mhj^2
		   */
		  obsinerhj = sumsquare_KD[ _HJ ] - 
		    N_KD[ _HJ ] * sqr( OldNewCen_KD[ _HJ ] ) ;

		  Iner_KD[ _HJ ] = obsinerhj ;
	      }

	  }
	  else /* then this class size == 0 */ {
	      /* signal empty class and which one */
	      sts = STS_W_EMPTYCLASS ;
	      (*EmptyK_P) = h + 1 ;
	  }

      }
      /* end  For each variable j */

  }
  /* end  For each class h */


  /* Free local structures */
  GenFree( oldmean_KD ) ;
  GenFree( sumsquare_KD ) ;
  GenFree( sumdata_KD  ) ;

  return sts ;

}   /* end of CommonGaussDiag() */



/* ------------------------------------------------------------------- */
static void InerToDisp
(
 DisperET      DispType, /* I : dispersion model */
 int           N,        /* I : nb of observation vectors */
 int           Nk,       /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD,  /* O : dispersion in each class/variable */
 StatusET*     StsP      /* [O] : STS_E_FUNCARG if unknown DispType */
)
/* ------------------------------------------------------------------- */
{
  switch( DispType ) {
  case DISPER___:  /* Common dispersion for all classes and variables */
    InerToDisp__( N, Nk, D, NbObs_K, NbObs_KD, Iner_KD, MissMode, 
		  Disp_KD ) ;
    break ;

  case DISPER_K_:  /* In each class, common dispersion for all variables */
    InerToDispK_( N, Nk, D, NbObs_K, NbObs_KD, Iner_KD, MissMode, 
		  Disp_KD ) ;
    break ;

  case DISPER__D:  /* In each class, common dispersion for all variables */
    InerToDisp_D( N, Nk, D, NbObs_K, NbObs_KD, Iner_KD, MissMode, 
		  Disp_KD ) ;
    break ;

  case DISPER_KD:  /* In each class, common dispersion for all variables */
    InerToDispKD( N, Nk, D, NbObs_K, NbObs_KD, Iner_KD, MissMode, 
		  Disp_KD ) ;
    break ;

  default:
    (*StsP) = STS_E_FUNCARG ;
  }
}   /* end of InerToDisp() */


/* ------------------------------------------------------------------- */
static void InerToDisp__
(
 int           N,        /* I : nb of observation vectors */
 int           Nk,       /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
)
/* ------------------------------------------------------------------- */
{
  int      k ;         /* current class    0..K-1 */
  int      d ;         /* current variable 0..D-1 */

  float    vol ;       /* common volume */
  float    nobs ;      /* number of effective observations : ]0,N*D] */

  /* Compute common volume : 
     sum dispersions over the classes and variables 
     then divide by all/observed data size */
  
  vol = 0.0 ;
  nobs = 0.0 ;

  for ( k = 0; k < Nk ; k ++ )
    {
      /* If not empty class */
      if ( NbObs_K[ k ] > 0 )
	{
	  for ( d = 0 ; d < D ; d ++ )
	    {
	      vol  += Iner_KD[ k * D + d ] ;
	      nobs += NbObs_KD[ k * D + d ] ; 
	    }
	}
    }

  if ( MissMode == MISSING_REPLACE )
    vol /= ( N * D ) ;
  else
    vol /= nobs ;  /* nobs must be = N * D - DataP->NbMiss */

  /* For each class k and variable d */
  for ( k = 0 ; k < Nk ; k ++ )
    for ( d = 0 ; d < D ; d ++ )
    {
      /* Set its volume to the common volume */
      Disp_KD[ k * D + d ] = vol ;
    }

}   /* end of InerToDisp__() */


/* ------------------------------------------------------------------- */
static void InerToDispK_
(
 int           N,        /* I : nb of observation vectors */
 int           K,        /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
)
/* ------------------------------------------------------------------- */
{
  int      k ;         /* current class    0..K-1 */
  int      d ;         /* current variable 0..D-1 */

  float    sumd_nkd ;    /* size of obs class */
  float    sumd_inerkd ; /* obs inertia in class k */
  float    dispk ;       /* dispersion in class k */


  /* For each class k */
  for ( k = 0; k < K ; k ++ ) {

      /* If not empty class */
      if ( NbObs_K[ k ] > 0 ) {

	  /* Reestimate its volume: sum inertia over the variables 
	     then divide by all/observed data size in this class */

	  sumd_nkd    = 0.0 ;
	  sumd_inerkd = 0.0 ;  /* use first component as sum first */

	  for ( d = 0 ; d < D ; d ++ )
	    {
	      sumd_nkd     += NbObs_KD[ k * D + d ] ;
	      sumd_inerkd  += Iner_KD[ k * D + d ] ;
	    }

	  if ( MissMode == MISSING_REPLACE )
	    dispk = sumd_inerkd /= ( D * NbObs_K[ k ] ) ;  /*V1.06-e*/
	  else
	    dispk = sumd_inerkd / sumd_nkd ;

	  /* Copy same dispersion into all variables of this class */
	  for ( d = 0 ; d < D ; d ++ )
	    {
	      Disp_KD[ k * D + d ] = dispk ;
	    }
      }
      /* (Else class k is empty : keep old dispersion */ 

    }

  /*??*/

}   /* end of InerToDispK_() */

/* ------------------------------------------------------------------- */
static void InerToDisp_D
(
 int           N,        /* I : nb of observation vectors */
 int           K,        /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
)
/* ------------------------------------------------------------------- */
{
  int      k ;         /* current class    0..K-1 */
  int      d ;         /* current variable 0..D-1 */

  float    sumk_nkd ;    /* obs size of variable d */
  float    sumk_inerkd ; /* obs inertia in variable d */
  float    dispd ;       /* dispersion of variable d */

  /* For each variable d */
  for ( d = 0; d < D ; d ++ ) {

    /* Reestimate its volume: sum inertia over the classes 
       then divide by all/observed data size in this variable */

    sumk_nkd    = 0.0 ;
    sumk_inerkd = 0.0 ;

    for ( k = 0 ; k < K ; k ++ ) {
      sumk_nkd    += NbObs_KD[ k * D + d ] ;
      sumk_inerkd += Iner_KD[ k * D + d ] ;
    }

    if ( MissMode == MISSING_REPLACE )
      dispd = sumk_inerkd / N ;
    else
      dispd = sumk_inerkd / sumk_nkd ;

    /* Copy same dispersion of this variable into all classes */
    for ( k = 0 ; k < K ; k ++ ) {
      Disp_KD[ k * D + d ] = dispd ;
    }
  }

  /*???*/

}   /* end of InerToDisp_D() */



/* ------------------------------------------------------------------- */
static void InerToDispKD
(
 int           N,        /* I : nb of observation vectors */
 int           K,        /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
)
/* ------------------------------------------------------------------- */
{
  int      k ;         /* current class    0..K-1 */
  int      d ;         /* current variable 0..D-1 */

  /* For each class k */
  for ( k = 0 ; k < K ; k ++ ) {

    /* For each variable d */
    for ( d = 0; d < D ; d ++ ) {

      /* Reestimate its volume: divide inertia by all/observed data 
	 size in this variable/ class if not empty */

      if ( MissMode == MISSING_REPLACE ) {
	if ( NbObs_K[ k ] > EPSILON )
	  Disp_KD[ k * D + d ] = Iner_KD[ k * D + d ] / NbObs_K[ k ] ;
	/* else empty class : do not update dispersion */
      }
      else /* IGNORE */ {
	if ( NbObs_KD[ k * D + d ] > EPSILON )
	  Disp_KD[ k * D + d ] = Iner_KD[ k * D + d ] / NbObs_KD[ k * D + d ] ;
      }
    }
  }

  /*???*/

}   /* end of InerToDispKD() */




/* ------------------------------------------------------------------- */
static StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY */ /*V1.06-b*/
EstimParaLaplace
(
  const float   *C_NK,      /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           Nk,         /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/
  const ModelSpecT* SpecP,  /* I : model specification */ /*V1.06-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 */ /*V1.05-b*/
  ModelParaT    *ParaP      /* O : estimated parameters */
)
/* ------------------------------------------------------------------- */
{
  int      N = DataP->NbPts ;
  int      D = DataP->NbVars ;

  StatusET sts ;

  /* Common computation for any diagonal gaussian model with missing data */
  sts = CommonLaplaceDiag
    ( DataP->PointsM, DataP->SortPos_ND, N, D, C_NK, Nk, MissMode, 
      ParaP->Disp_KD, 
      ParaP->Center_KD, 
      EmptyK_P, ParaP->NbObs_K, ParaP->NbObs_KD, ParaP->Iner_KD ) ;

  /* Dispersion is computed depending on dispersion model */
  InerToDisp( SpecP->ClassDisper, N, Nk, D, ParaP->NbObs_K, ParaP->NbObs_KD, 
	      ParaP->Iner_KD, MissMode, 
	      ParaP->Disp_KD, & sts ) ;

  return sts ;

  /*???*/
}   /* end of EstimParaLaplace() */


/*V1.06-f*/
/* ------------------------------------------------------------------- */
static StatusET          /* Return status : OK, EMPTY or MEMORY */
CommonLaplaceDiag
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  int           N,             /* I : nb of observation vectors */
  int           D,             /* I : nb of variables */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  int           K,             /* I : number of clusters */
  MissET        Miss,          /* I : how to treat missing data */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */

  float*        OldNewCen_KD,  /* I/O : old then updated means (K,D) */

  int*          EmptyK_P,      /* O : index of empty class (0 or 1 ..K) */
  float*        N_K,           /* O : size of a class (K) */
  float*        N_KD,          /* O : size of a class and variable (K,D) */
  float*        Iner_KD        /* O : inertia of a class/variable (K,D) */
 )
/* ------------------------------------------------------------------- */
{
  const char* func = "CommonLaplaceDiag" ;

  float*   oldcent_KD ;    /* saved previous centers (K,D), local alloc */

  StatusET sts ;           /* return status */


  /* Allocate local structures */
  oldcent_KD   = GenAlloc( K * D, sizeof(float), 1, func, "oldcent_KD" ) ;

  /* Save old mean value */
  memcpy( oldcent_KD , OldNewCen_KD , K * D * sizeof( float ) ) ;

  EstimSizes( X_ND, C_NK, N, D, K, 
	      N_K, N_KD ) ;

  sts = EstimLaplaceCenters( X_ND, Sort_ND, C_NK, N, D, K, Miss, 
			     N_K, N_KD, oldcent_KD, OldDisp_KD, 
			     OldNewCen_KD, EmptyK_P ) ;

  EstimLaplaceIner( X_ND, C_NK, N, D, K, Miss, 
		    N_K, N_KD, oldcent_KD, OldDisp_KD, OldNewCen_KD, 
		    Iner_KD ) ;

  /* Free local structures */
  GenFree( oldcent_KD ) ;

  return sts ;

  /*???*/
}   /* end of CommonLaplaceDiag() */


/* ------------------------------------------------------------------- */
static void EstimSizes
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  float*        N_K,           /* O : size of a class (K) */
  float*        N_KD           /* O : size of a class and variable (K,D) */  
 )
/* ------------------------------------------------------------------- */
{
  int      h ;             /* current class    0..K-1 */
  int      j ;             /* current variable 0..D-1 */
  int      i ;             /* current object   0..N-1 */

  
  /* For each class h and variable j */
  for ( h = 0 ; h < K ; h ++ ) {

      for ( j = 0 ; j < D ; j ++ ) {

	  /* Initialize the sum[i] quantities to 0 */
	  N_K[ h ]             = 0.0 ;
	  N_KD[ _HJ ]          = 0.0 ;

	  /* For each object i, increment class size and
	     eventually observed data size if xij not nan */
	  for ( i = 0 ; i < N ; i ++ ) {

	      float cih = C_NK[ _IH ] ;
	      float xij = X_ND[ _IJ ] ;

	      N_K[ h ] += cih ;

	      if ( ! isnan( xij ) ) {
		  N_KD[ _HJ ]          += cih ;
	      }
	  } /* End  For each object i (all sums[i] are now done) */
      }
  }  /* End   for each class h and variable j */

}   /* end of EstimSizes() */



/* ------------------------------------------------------------------- */
static StatusET          /* Return status : OK, EMPTY */
EstimLaplaceCenters
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const MissET  Miss,          /* I : how to treat missing data */
  const float*  N_K,           /* I : size of a class (K) */
  const float*  N_KD,          /* I : size of a class and variable (K,D) */  
  const float*  OldCen_KD,     /* I : old centers */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */

  float*        NewCen_KD,     /* O : updated centers (K,D) */
  int*          EmptyK_P       /* O : index of empty class (0 or 1 ..K) */
 )
/* ------------------------------------------------------------------- */
{
  StatusET sts ;           /* return status */

  int      h ;             /* current class    0..K-1 */
  int      j ;             /* current variable 0..D-1 */

  int      imed ;          /* index of median of observed values */
  float    medval ;        /* median value (eventually midway) */
  float    cumwei ;        /* cumulated weights until median observation */


  /* Set no empty class by default */
  (*EmptyK_P) = 0 ;
  sts         = STS_OK ;


  /* For each class h and variable j */
  for ( h = 0 ; h < K ; h ++ ) {

    for ( j = 0 ; j < D ; j ++ ) {

      /* If this class size > 0 */
      if ( N_K[ h ] > EPSILON ) {
	
	/* First, find median position of observed values 
	 */
	ComputeMedian( X_ND, Sort_ND, C_NK, N, D, K, h, j, N_KD[ _HJ ],
		       & imed, & cumwei, & medval ) ;

	/* 
	 * Compute center from observed median, depending on presence
	 * of missing data, REPLACE or IGNORE missing mode 
	 */

	if ( N_KD[ h * D + j ] == N_K[ h ] ) /* no missing data */ {
	  NewCen_KD[ h * D + j ] = medval ;
	}
	else if ( N_KD[ h * D + j ] >= EPSILON ) /* some missing data */ {
	   if ( Miss == MISSING_REPLACE ) {

	     NewCen_KD[ h * D + j ] = 
	       FindMinInerLaplaceEM( X_ND, Sort_ND, C_NK, N, D, K, h, j, 
				     imed, cumwei, OldCen_KD[ h * D + j ], 
				     OldDisp_KD[ h * D + j ],
				     N_KD[ h * D + j ], 
				     N_K[ h ] - N_KD[ h * D + j ] ) ;

	     /* 
		NewCen_KD[ h * D + j ] =  
		( N_KD[ h * D + j ] * medval + 
		( N_K[ h ] - N_KD[ h * D + j ] ) * OldCen_KD[ h * D + j ] )
		/ N_K[ h ] ;
		*/
	   }
	   else  /* assume Miss == MISSING_IGNORE */ {
	     NewCen_KD[ h * D + j ] = medval ;
	   }
	}
	else /* N_KD[ h * D + j ] < EPSILON : all data missing */ {
	  NewCen_KD[ h * D + j ] = OldCen_KD[ h * D + j ] ;
	}

      }
      else /* then this class size == 0 */ {
	NewCen_KD[ h * D + j ] = OldCen_KD[ h * D + j ] ;
	/* signal empty class and which one */
	sts = STS_W_EMPTYCLASS ;
	(*EmptyK_P) = h + 1 ;
      }

    }
  }
  /* end  For each class h and variable j */

  return sts ;
  /*???*/
}   /* end of EstimLaplaceCenters() */



/* ------------------------------------------------------------------- */
static void ComputeMedian
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const int     H,             /* I : current cluster */
  const int     J,             /* I : current variable */
  const float   totwei,        /* I : total weight */
  int*          ImedP,         /* O : index of median of observed values */
  float*        CumweiP,       /* O : cumulated weights until median obs. */
  float*        MedvalP        /* O : median value (eventually midway) */
 )
/* ------------------------------------------------------------------- */
{
  float    halfwei = totwei / 2 ;
  int      i ;             /* current index in ascending order  0..N-1 */
  int      pos ;           /* actual position of current object 0..N-1 */
  int      posmed ;        /* actual position of median 0..N-1 */
  int      posnext ;       /* actual position of next to median 0..N-1 */

  
  /* Scan the data in ascending order (skip NaN values), and find position 
     where the cumulated weights >= half total weight */
  for ( i = 0, 
	  (*CumweiP) = 0.0 ;
	( i < N ) &&
	  (*CumweiP) < halfwei ;
	i ++ ) {

    pos = Sort_ND[ i * D + J ] ;
    if ( !isnan( X_ND[ pos * D + J ] ) ) 
      (*CumweiP) += C_NK[ pos * K + H ] ;

  }  /* At this point, (*CumweiP) >= halfwei */

  (*ImedP) = i - 1 ;
  posmed = Sort_ND[ (*ImedP) * D + J ] ;

  /* Median value depends if = or > */
  if ( (*CumweiP) > halfwei + EPSILON ) /* > : median value is here*/ {

    (*MedvalP) = X_ND[ posmed * D + J ] ;
  }
  else  /* = : median value is midway to next non nan observation */ {

    for ( i = (*ImedP) + 1 ;
	  isnan( X_ND[ Sort_ND[ i * D + J ] * D + J ] ) ||
	    ( C_NK[ Sort_ND[ i * D + J ] * K + H ] < EPSILON ) ;
	  i ++ ) {
    }
    posnext = Sort_ND[ i * D + J ] ;
    (*MedvalP) = 0.5 * ( X_ND[ posmed * D + J ] + X_ND[ posnext * D + J ] ) ;
  }

}   /* end of ComputeMedian() */




/* ------------------------------------------------------------------- */
static float FindMinInerLaplaceEM  /* ret: minimizer of expected inertium */
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const int     H,             /* I : current cluster */
  const int     J,             /* I : current variable */
  const int     Imed,          /* I : median position of observed data */
  const float   Cumwei,        /* I : cumulated weight until Xmed */
  const float   Cen0,          /* I : previous center */
  const float   Disp0,         /* I : previous dispersion */
  const float   Nhj,           /* I : number of observations */
  const float   Nmis           /* I : number of missing data */
 )
/* ------------------------------------------------------------------- */
{
  float   xmed ;       /* observation at median position */
  float   weidif ;     /* (before + med) - after or (med + after) - before */
  float   weibefore ;  /* cumulated weight before median position */
  int     direction ;  /* 1 if xmed < cen0, -1 otherwise */

  float   a ;       /* lower bound of current interval */
  float   b ;       /* upper bound of current interval */
  float   intwei ;  /* cumulated intermediate weight */
  float   d ;       /* current derivative R'(.) */
  float   res ;     /* value to return (minimizes expected inertia) */
  int     i ;       /* current index in ascending order  Ifirst..Ilast */

  int     sts ;     /* 0 : R' < 0 until now
		       1 : R'(a) >= 0
		       2 : R'(a) < 0 and R'(b) >= 0 */

  xmed = X_ND[ Sort_ND[ Imed * D + J ] * D + J ] ;

  if ( xmed < Cen0 - EPSILON ) {
    direction =  1 ;
    weidif    =  Cumwei - ( Nhj - Cumwei ) ;
  }
  else if ( xmed >  Cen0 + EPSILON ) {
    direction =  -1 ;
    weibefore =  Cumwei - C_NK[ Sort_ND[ Imed * K + H ] * K + H ] ;
    weidif    =  (Nhj - weibefore) - weibefore ;
  }
  else   /* xmed == Cen0, no search necessary */ {
    return Cen0 ;
  }


  /* Start at median position */
  a      = xmed ;
  intwei = 0.0 ;
  sts    = 0 ; 

  if ( (d = DerivInerDir( a, weidif, intwei, Nmis, Cen0, Disp0 )) >= 0 ) 
    /* R'(a) >= 0,  minimum at a */ {
    sts = 1 ;  
  }
  else  /* R'(a) < 0, look at b */ {

    for ( i = Imed + direction ;

	  ( sts == 0 ) && 
	    ( (X_ND[ Sort_ND[ i * D + J ] * D + J ] - Cen0) * direction < 0 ) 
	    && ( i >= 0 ) && ( i < N ) ;      /* "safety belt" */

	  i = i + direction ) 
      /* For each intermediate point */ {

      if ( ( ! isnan( X_ND[ Sort_ND[ i * D + J ] * D + J ] ) )
	   && ( C_NK[ Sort_ND[ i * D + J ] * K + H ] > EPSILON ) ) {

	b = X_ND[ Sort_ND[ i * D + J ] * D + J ] ;
	if ( (d = DerivInerDir( b, weidif, intwei, Nmis, Cen0, Disp0 )) >= 0 ) 
	  /* R'(b) >= 0, and R'(a) < 0 */ {

	  sts = 2 ;
	}
	else /* R'(b) < 0 */ {

	  /* Start checking new interval */
	  a = b ;
	  intwei += 2 * C_NK[ Sort_ND[ i * D + J ] * K + H ] ;
	  if ( (d = DerivInerDir( a, weidif, intwei, Nmis, Cen0, Disp0 )) >= 0)
	    /* R'(a) >= 0,  minimum at a */ {

	    sts = 1 ;
	  }
	}
      }
    }    /* end for each intermediate point */
  }

  if ( sts == 0 ) 
    /* R'(b) < 0, R'(a) < 0 => check for b = Cen0 */ {
    b = Cen0 ;
    if ( (d = DerivInerDir( b, weidif, intwei, Nmis, Cen0, Disp0 )) >= 0 ) 
      /* R'(b) >= 0, and R'(a) < 0 */ {

      sts = 2 ;
    }
  }

  switch( sts ) {
  case 0:
    /* R'(a) < 0, R'(b) < 0 => Cen0 is the minimizer */
    res = Cen0 ;
    break ;

  case 1:
    /* R'(a) >= 0 => a is the minimizer */
    res = a ;
    break ;

  default:
    /* R'(a) < 0, R'(b) >= 0 => 
       y in [a,b] minimizing R is solution of :
       | y - Cen0 | = - Disp0 log( 1 - ( weidif - intwei ) / Nmis )
       */
    res = Cen0 + direction * Disp0 * log( 1 - ( weidif + intwei ) / Nmis ) ;
  }

  return res ;
  /*???*/
}   /* end of FindMinInerLaplaceEM() */


/* ------------------------------------------------------------------- */
static float DerivInerDir      /* ret: +/- R'(Y), R expected inertia */ 
 (
  const float   Y,             /* I : point at which to compute R' */
  const float   Weidif,        /* I : "after+med-bef" or "bef+med-after" */
  const float   Intwei,        /* I : cumulated intermediate weight */
  const float   Nmis,          /* I : number of missing data */
  const float   Cen0,          /* I : previous center */
  const float   Disp0          /* I : previous dispersion */
 )
/*\
     This function returns 
       R'(Y)   if Xmed < Cen0
     - R'(Y)   if Xmed > Cen0
     Testing the positiveness of this function thus tests R'(Y) > 0
     in the first case and R'(Y) < 0 in the second case.
\*/
/* ------------------------------------------------------------------- */
{

  float d = Weidif + Intwei 
    - Nmis * ( 1 - exp( - fabs( Y - Cen0 ) / Disp0 ) ) ;

  return d ;

  /*???*/
}   /* end of DerivIner() */




/* ------------------------------------------------------------------- */
static void EstimLaplaceIner
 ( 
  const float*  X_ND,          /* I : data matrix (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const MissET  Miss,          /* I : how to treat missing data */
  const float*  N_K,           /* I : size of a class (K) */
  const float*  N_KD,          /* I : size of a class and variable (K,D) */  
  const float*  OldCen_KD,     /* I : old centers */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */
  const float*  NewCen_KD,     /* I : new centers */
  float*        Iner_KD        /* O : inertia of a class/variable (K,D) */  
 )
/* ------------------------------------------------------------------- */
{
  int      h ;             /* current class    0..K-1 */
  int      j ;             /* current variable 0..D-1 */
  int      i ;             /* current object   0..N-1 */

  
  /* For each class h and variable j */
  for ( h = 0 ; h < K ; h ++ ) {

      for ( j = 0 ; j < D ; j ++ ) {

	  /* Initialize the inertium to 0 */
	  Iner_KD[ _HJ ] = 0.0 ;

	  /* For each non nan object i, increment class inertium */
	  for ( i = 0 ; i < N ; i ++ ) {

	      float cih = C_NK[ _IH ] ;
	      float xij = X_ND[ _IJ ] ;

	      if ( ! isnan( xij ) ) {
		Iner_KD[ _HJ ] += cih * fabs( xij - NewCen_KD[ _HJ ] ) ;
	      }

	  } /* End  For each non-nan object i (all sums[i] are now done) */

	  
	  if ( ( N_KD[ _HJ ] < N_K[ h ] ) &&  ( Miss == MISSING_REPLACE ) )
	       /* If some missing data and replace mode */ {

	    Iner_KD[ _HJ ] += 
	      ( N_K[ h ] - N_KD[ _HJ ] ) *
	      ( fabs( OldCen_KD[ _HJ ] - NewCen_KD[ _HJ ] ) +
		OldDisp_KD[ _HJ ] * 
		exp( - fabs( OldCen_KD[ _HJ ] - NewCen_KD[ _HJ ] ) / 
		     OldDisp_KD[ _HJ ] ) ) ; 
	  }

      }
  }  /* End   for each class h and variable j */
  
  /*???*/
}   /* end of EstimLaplaceIner() */




/*+++*/

#if 0


/* ------------------------------------------------------------------- */
StatusET                    /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaP_V_I
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{
  int      k ;
  StatusET sts ;


  sts = ParaPkV_I( Cih_nk, DataP, K, MissMode, EmptyK_P, NoiseParaP ) ;

  for ( k = 0 ; k < K ; k ++ )
    {
        NoiseParaP->Pk[ k ] = 1.0 / K ;
    }
  return sts ;
}

/* ------------------------------------------------------------------- */
StatusET ParaPkV_I          /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{

    if ( DataP->NbMiss == 0 )
      return ParaPkV_I_Complete( Cih_nk, DataP, K, MissMode, 
				 EmptyK_P, NoiseParaP ) ;
    else
      return ParaPkV_I_Missing( Cih_nk, DataP, K, MissMode, 
				EmptyK_P, NoiseParaP ) ;

}   /* end of ParaPkV_I() */


/* ------------------------------------------------------------------- */
static StatusET             /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaPkV_I_Complete
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{
    int         h ;         /* class counter : 0 .. K - 1 */
    float       sumh_trTh ; /* sum over h of class h's dispersion matrix's trace */
    float       V ;         /* common volume = sum[h] trTh / (N * Nd) */
    int         nd  = DataP->NbVars ;   /* Nd : nb of variables */
    int         npt = DataP->NbPts ;    /* N : nb of points */

    /* For each class */
    for ( h = 0 , (*EmptyK_P) = 0, sumh_trTh = 0.0 ; 
          h < K ; 
          h ++ )
    {
        int     ipt ;
        int     d ;
        float   sumi_uih ;
        float   sumi_uih_yi2 ;

        /* 'fuzzy' cardinal : sizeh = sum[i]( uih ) */
        /* proportion :       ph    = 1 / K */
        /* mean :             mh(d) = sum[i]( uih * yi(d) ) / sizeh */
        /* trace of its dispersion matrix : 
           trTh = sum[i]( uih * sum[d] yi(d)™ ) - sizeh * sum[d] mh(d)™ */
        /* volume :           Vh = V = sum[h] trTh / (N * Nd) */

        /* Initialize means to zero */
        for ( d = 0 ; d < nd ; d ++ )
        {
            NoiseParaP->Mk[ ( h * nd ) + d ] = 0 ;
        }

        /* Do cumulated sums over all points */
        for ( ipt = 0, sumi_uih = 0.0, sumi_uih_yi2 = 0.0 ; 
              ipt < npt ; 
              ipt ++ )
        {
            float   sumd_uih_yid2 ;
            float   uih                             = Cih_nk[ ( ipt * K ) + h ] ;

            for ( d = 0, sumd_uih_yid2 = 0.0 ; 
                  d < nd ; 
                  d ++ )
            {
                float   yid   = DataP->PointsM[ (ipt * nd) + d ] ;

                NoiseParaP->Mk[ ( h * nd ) + d ]    += uih * yid ;
                sumd_uih_yid2                       += uih * yid * yid ;
            }

            sumi_uih        += uih ;
            sumi_uih_yi2    += sumd_uih_yid2 ;
        }

        /* Normalize the sums by the cardinal of the class */
        if ( sumi_uih != 0.0 )
        {
            float   sumd_mhd2 ;
            float   trTh ;
            float   invsizeh = 1 / sumi_uih ;

            for ( d = 0, sumd_mhd2 = 0.0 ; 
                  d < nd ; 
                  d ++ )
            {
                float   mhd         = NoiseParaP->Mk[ (h * nd) + d ] ;

                mhd                 *= invsizeh ;
                sumd_mhd2           += mhd * mhd ;

                NoiseParaP->Mk[ (h * nd) + d ]  = mhd ;
            }

            trTh                    = sumi_uih_yi2 - sumi_uih * sumd_mhd2 ;
            sumh_trTh               += trTh ;
        }
        else
        {
            (*EmptyK_P) = h + 1 ;
        }

        NoiseParaP->Pk[ h ] = sumi_uih / npt ;

        /* Set normalized variance matrix to identity */
        SetIdMatrix( nd , & NoiseParaP->Ck[ h*nd*nd ] ) ;  /*V1.03-a*/

    } /* end for ( h = 0 , (*EmptyK_P) = 0 ; h < K ; h ++ ) */

    V = sumh_trTh / ( npt * nd ) ;

    for ( h = 0 ; h < K ; h ++ )
    {
        NoiseParaP->Vk[ h ] = V ;
    }

    if ( (*EmptyK_P) == 0 )
      return STS_OK ;
    else
      return STS_W_EMPTYCLASS ;

}   /* end of ParaPkV_I_Complete() */



/* ------------------------------------------------------------------- */
static StatusET             /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaPkV_I_Missing
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{
  const char* func = "ParaPkV_I_Missing" ;

  int      N = DataP->NbPts ;
  int      D = DataP->NbVars ;

  int      h ;         /* current class    0..K-1 */
  int      j ;         /* current variable 0..D-1 */

  float*   nh_k ;      /* size / class (K) local alloc (unused) */
  float*   nhj_kd ;    /* size / obs class & var (K,D) local alloc (unused)*/
  float*   disphj_kd ; /* dispersion / class & var (K,D) local alloc */
  float    vol ;       /* common volume */

  StatusET sts ;


  /* Allocate local data */
  nh_k      = GenAlloc( K,     sizeof( float ), 1, func, "nh_k" ) ;
  nhj_kd    = GenAlloc( K * D, sizeof( float ), 1, func, "nhj_kd" ) ;
  disphj_kd = GenAlloc( K * D, sizeof( float ), 1, func, "disphj_kd" ) ;

  /* Common computation for any diagonal gaussian model with missing data */
  sts = CommonGaussDiagMissing
    ( DataP->PointsM, DataP->NbPts, DataP->NbVars, Cih_nk, K, MissMode,
      NoiseParaP->Vk, NoiseParaP->Ck, 
      NoiseParaP->Mk, NoiseParaP->Pk, EmptyK_P, nh_k, nhj_kd, disphj_kd ) ;

  /* Compute common volume : 
     sum dispersions over the classes and variables 
     then divide by all/observed data size */
  
  vol = 0.0 ;

  for ( h = 0; h < K ; h ++ )
    {
      /* If not empty class */
      if ( nh_k[ h ] > 0 )
	{
	  for ( j = 0 ; j < D ; j ++ )
	    {
	      vol += disphj_kd[ _HJ ] ;
	    }
	}
    }

  if ( MissMode == MISSING_REPLACE )
    vol /= ( N * D ) ;
  else
    vol /= ( N * D - DataP->NbMiss ) ;

  /* For each class h */
  for ( h = 0; h < K ; h ++ )
    {
      /* Set its volume to the common volume */
      NoiseParaP->Vk[ h ] = vol ;

      /* Set its normalized covariance to identity */
      SetIdMatrix( D, & NoiseParaP->Ck[ h * D * D ] ) ;
    }
      
  /* Free local data */
  GenFree( disphj_kd ) ;
  GenFree( nhj_kd ) ;
  GenFree( nh_k ) ;

  return sts ;
}   /* end of ParaPkV_I_Missing() */



/* ------------------------------------------------------------------- */
StatusET             /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaP_VkI
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{
  StatusET sts ;
  int k ;

  sts = ParaPkVkI( Cih_nk, DataP, K, MissMode, EmptyK_P, NoiseParaP ) ;

  for ( k = 0 ; k < K ; k ++ )
    {
        NoiseParaP->Pk[ k ] = 1.0 / K ;
    }
  return sts ;
}


/* ------------------------------------------------------------------- */
StatusET             /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaPkVkI
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{

    if ( DataP->NbMiss == 0 )
      return ParaPkVkI_Complete( Cih_nk, DataP, K, MissMode, 
				 EmptyK_P, NoiseParaP ) ;
    else
      return ParaPkVkI_Missing( Cih_nk, DataP, K, MissMode, 
				EmptyK_P, NoiseParaP ) ;

}   /* end of ParaPkV_I() */



/* ------------------------------------------------------------------- */
static StatusET             /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaPkVkI_Complete
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{
    int         k ;

    /* For each class */
    for ( k = 0 , (*EmptyK_P) = 0 ; k < K ; k ++ )
    {
        int     ipt ;                  
        int     id ;
        int     nd  = DataP->NbVars ;
        int     npt = DataP->NbPts ;    /* = n nb points */
        float   sizek = 0.0 ;            /* = nk fuzzy cardinal [0,n] */
        float   trWk = 0.0 ;

        /* Compute its 'fuzzy' cardinal : nk = sum[i=1,n]( uik ) */
        /* Compute its proportion : pk = nk / n */
        /* Compute its mean (mk) : mk = sum[i=1,n]( uik * Y(i,:) ) / nk */
        /* Compute its intra-class inertia matrix : 
           Wk = sum[i=1,n]( uik * ( Y(i,:) - mk ) * ( Y(i,:) - mk )' ) / nk */
        /* Compute its volume : Vk = tr(Wk) / d */

        /* Initialize means to zero */
        for ( id = 0 ; id < nd ; id ++ )
        {
            NoiseParaP->Mk[ (k * nd) + id ] = 0 ;
        }

        /* Do cumulated sums over all points */
        for ( ipt = 0 ; ipt < npt ; ipt ++ )
        {
            float   uik = Cih_nk[ (ipt * K) + k ] ;

            sizek += uik ;
            for ( id = 0 ; id < nd ; id ++ )
            {
                float   yid = DataP->PointsM[ (ipt * nd) + id ] ;
                float   mkd = NoiseParaP->Mk[ (k * nd) + id ] ;

                mkd         = mkd + uik * yid ;
                trWk        = trWk + uik * yid * yid ;

                NoiseParaP->Mk[ (k * nd) + id ] = mkd ;
            }
        }

        /* Normalize the sums by the cardinal of the class */
        if ( sizek != 0.0 )
        {
            trWk /= sizek ;

            for ( id = 0 ; id < nd ; id ++ )
            {
                float   mkd = NoiseParaP->Mk[ (k * nd) + id ] ;

                mkd         = mkd / sizek ;
                trWk        = trWk - ( mkd * mkd ) ;

                NoiseParaP->Mk[ (k * nd) + id ] = mkd ;
            }
            NoiseParaP->Vk[ k ] = trWk / nd ;
        }
        else
        {
            (*EmptyK_P) = k + 1 ;
        }

        NoiseParaP->Pk[ k ] = sizek / npt ;

        /* Set normalized variance matrix to identity */
        SetIdMatrix( nd , & NoiseParaP->Ck[ k*nd*nd ] ) ;  /*V1.03-a*/

    } /* end for ( k = 0 , (*EmptyK_P) = 0 ; k < K ; k ++ ) */

    if ( (*EmptyK_P) == 0 )
      return STS_OK ;
    else
      return STS_W_EMPTYCLASS ;

}   /* end of ParaPkVkI_Complete() */


/* ------------------------------------------------------------------- */
static StatusET             /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaPkVkI_Missing
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{
  const char* func = "ParaPkVkI_Missing" ;

  int      N = DataP->NbPts ;
  int      D = DataP->NbVars ;

  int      h ;         /* current class    0..K-1 */
  int      j ;         /* current variable 0..D-1 */

  float*   nh_k ;      /* size / class (K) local alloc */
  float*   nhj_kd ;    /* size / obs class and variable (K,D) local alloc */
  float*   disphj_kd ; /* dispersion / class and variable (K,D) local alloc */

  float    sumj_nhj ;  /* size of obs class */

  StatusET sts ;


  /* Allocate local data */
  nh_k      = GenAlloc( K,     sizeof( float ), 1, func, "nh_k" ) ;
  nhj_kd    = GenAlloc( K * D, sizeof( float ), 1, func, "nhj_kd" ) ;
  disphj_kd = GenAlloc( K * D, sizeof( float ), 1, func, "disphj_kd" ) ;

  /* Common computation for any diagonal gaussian model with missing data */
  sts = CommonGaussDiagMissing
    ( DataP->PointsM, DataP->NbPts, DataP->NbVars, Cih_nk, K, MissMode,
      NoiseParaP->Vk, NoiseParaP->Ck, 
      NoiseParaP->Mk, NoiseParaP->Pk, EmptyK_P, nh_k, nhj_kd, disphj_kd ) ;

  /* For each class h */
  for ( h = 0; h < K ; h ++ )
    {
      /* If not empty class */
      if ( nh_k[ h ] > 0 )
	{
	  /* Reestimate its volume: sum dispersions over the variables 
	     then divide by all/observed data size */
	  sumj_nhj = 0 ;
	  NoiseParaP->Vk[ h ] = 0.0 ;  /* use as sum first */

	  for ( j = 0 ; j < D ; j ++ )
	    {
	      sumj_nhj            += nhj_kd[ _HJ ] ;
	      NoiseParaP->Vk[ h ] += disphj_kd[ _HJ ] ;
	    }

	  if ( MissMode == MISSING_REPLACE )
	    NoiseParaP->Vk[ h ] /= nh_k[ h ] ;
	  else
	    NoiseParaP->Vk[ h ] /= sumj_nhj ;
	}
      else /* class h is empty */
	{
	  NoiseParaP->Vk[ h ] = 0 ;
	}

      /* Set its normalized covariance to identity */
      SetIdMatrix( D, & NoiseParaP->Ck[ h * D * D ] ) ;
    }

  /* Free local data */
  GenFree( disphj_kd ) ;
  GenFree( nhj_kd ) ;
  GenFree( nh_k ) ;

  return sts ;
}



/*V1.05-d*/
/* ------------------------------------------------------------------- */
static StatusET                /* Returns : OK, W_EMPTYCLASS or E_MEMORY */
CommonGaussDiagMissing
 (
  const float*  Xij_nd,        /* I : data matrix (N,D) */
  int           N,             /* I : nb of observation vectors */
  int           D,             /* I : nb of variables */
  const float*  Cih_nk,        /* I : classification matrix (N,K) */
  int           K,             /* I : number of clusters */
  MissET        Miss,          /* I : how to treat missing data */
  const float*  OldVh_k,       /* I : old volumes (K) */
  const float*  OldChjj_kdd,   /* I : old normalized covariances (K,D,D) */

  float*        OldNewMhj_kd,  /* I/O : old then updated means (K,D) */
  float*        NewPh_k,       /* O : updated proportions (K) */
  int*          EmptyK_P,      /* O : index of empty class (0 or 1 ..K) */
  float*        Nh_k,          /* O : size of a class (K) */
  float*        Nhj_kd,        /* O : size of a class and variable (K,D) */
  float*        Disphj_kd      /* O : dispersion of a class/variable (K,D) */
 )
/*\

   This function does the common computation that is needed at the
   M-step, in all the diagonal Gaussian models with missing data.

   The routine takes as input a data matrix Xij_nd (possibly
   containing NaN values), a classification matrix Cih_nk, and the
   kind of missing data processing Miss (REPLACE or IGNORE) ; and also
   the old volumes OldVh_k, covariance matrices OldChjj_kdd, and means
   OldNewMhj_kd of the classes, to be used in the REPLACE mode.

   It produces as output the reestimated mean vectors in OldNewMhj_kd;
   reestimated proportions in NewP_k; the last encountered empty class
   in EmptyK_P (or 0 if no empty class); the {sum[i] cih} in Nh_k, the
   {sum[i] cih rij} in Nhj_kd (number of present observations in class h
   and variable j); and the "missing data dispersion" of each
   class/variable in Disphj_kd (i.e. the missing data equivalent of
   {sum[i] Cih (xij - mhj)^2}).

\*/
/* ------------------------------------------------------------------- */
{
  const char* func = "CommonGaussDiagMissing" ;

  float*   sumdata_hj_kd ;    /* sum[i] cih rij xij   (K,D), local alloc */
  float*   sumsquare_hj_kd ;  /* sum[i] cih rij xij^2 (K,D), local alloc */
  float*   oldmean_hj_kd ;    /* saved previous means (K,D), local alloc */
  int      h ;                /* current class    0..K-1 */
  int      j ;                /* current variable 0..D-1 */
  int      i ;                /* current object   0..N-1 */
  StatusET sts ;              /* return status */


  sts = STS_OK ;

  /* Allocate local structures */
  sumdata_hj_kd   = GenAlloc( K * D, sizeof(float), 1, func, "sumdata" ) ;
  sumsquare_hj_kd = GenAlloc( K * D, sizeof(float), 1, func, "sumsquare" );
  oldmean_hj_kd   = GenAlloc( K * D, sizeof(float), 1, func, "oldmean" ) ;

  /* Save old mean value */
  memcpy( oldmean_hj_kd , OldNewMhj_kd , K * D * sizeof( float ) ) ;

  /* Set no empty class by default */
  (*EmptyK_P) = 0 ;

  /* For each class h */
  for ( h = 0 ; h < K ; h ++ )
    {
      /* For each variable j */
      for ( j = 0 ; j < D ; j ++ )
	{
	  /* Initialize the sum[i] quantities to 0 */
	  Nh_k[ h ]               = 0.0 ;
	  Nhj_kd[ _HJ ]           = 0.0 ;
	  sumdata_hj_kd[ _HJ ]    = 0.0 ;
	  sumsquare_hj_kd[ _HJ ]  = 0.0 ;
	  OldNewMhj_kd[ _HJ ]     = 0.0 ; /* use Mhj as sum first */

	  /* For each object i */
	  for ( i = 0 ; i < N ; i ++ )
	    {
	      float cih = Cih_nk[ _IH ] ;
	      float xij = Xij_nd[ _IJ ] ;

	      /* Increment class size */
	      Nh_k[ h ] += cih ;

	      /* If this object's j variable is observed */
	      if ( ! isnan( xij ) )
		{
		  /* Increment all sum[i] of observed quantities */
		  Nhj_kd[ _HJ ]           += cih ;
		  sumdata_hj_kd[ _HJ ]    += cih * xij ;
		  sumsquare_hj_kd[ _HJ ]  += cih * xij * xij ;
		}
	      /* Else xij is missing, don't change sum[i] of observed data */
	    }
	  /* End  For each object i (all sums[i] are now done) */

	  /* If this class size > 0 */
	  if ( Nh_k[ h ] > 0 )
	    {
	      /* 
	       * Compute means Mhj and dispersions Dhj from computed sums,
	       * depending on REPLACE or IGNORE missing mode 
	       */
	      if ( Miss == MISSING_REPLACE )
		{
		  /* 
		   * new_mhj = 
		   *  ( sum[i] cih rij xij + ( nh - nhj ) old_mhj ) / nh 
		   *
		   * disp_hj = 
		   *  "observed dispersion"_hj +
		   *  ( nh - nhj ) [ (new_mhj - old_mhj)^2 + old_vh old_chjj ]
		   */
		  float obsdisphj ;

		  OldNewMhj_kd[ _HJ ] = 
		    ( sumdata_hj_kd[ _HJ ] + 
		      ( Nh_k[ h ] - Nhj_kd[ _HJ ] ) * oldmean_hj_kd[ _HJ ] )
		    / Nh_k[ h ] ;

		  /*
		   * "observed dispersion"_hj = 
		   *      sum[i] cih rij xij^2 - nhj new_mhj^2 
		   */
		  obsdisphj = sumsquare_hj_kd[ _HJ ] - 
		    Nhj_kd[ _HJ ] * sqr( OldNewMhj_kd[ _HJ ] ) ;

		  Disphj_kd[ _HJ ] = 
		    obsdisphj +
		    ( Nh_k[ h ] - Nhj_kd[ _HJ ] ) *
		    ( sqr( OldNewMhj_kd[ _HJ ] - oldmean_hj_kd[ _HJ ] ) +
		      OldVh_k[ h ] * OldChjj_kdd[ h * D * D + j * D + j ] ) ;
		}
	      else  /* assume Miss == MISSING_IGNORE */
		{
		  float obsdisphj ;

		  /* 
		   * new_mhj = 
		   *  . if  nhj <> 0,  ( sum[i] cih rij xij ) / nhj  
		   *  . if  nhj == 0,  old_mhj
		   *
		   * disp_hj = 
		   *  "observed dispersion"_hj 
		   */

		  if ( Nhj_kd[ _HJ ] > 0 )
		    OldNewMhj_kd[ _HJ ] = sumdata_hj_kd[ _HJ ] / Nhj_kd[ _HJ ];
		  else
		    OldNewMhj_kd[ _HJ ] = oldmean_hj_kd[ _HJ ] ;

		  /*
		   * "observed dispersion"_hj = 
		   *      sum[i] cih rij xij^2 - nhj new_mhj^2 
		   */
		  obsdisphj = sumsquare_hj_kd[ _HJ ] - 
		    Nhj_kd[ _HJ ] * sqr( OldNewMhj_kd[ _HJ ] ) ;

		  Disphj_kd[ _HJ ] = obsdisphj ;
		}

	    } 
	  else /* then this class size == 0 */
	    {
	      /* signal empty class and which one */
	      sts = STS_W_EMPTYCLASS ;
	      (*EmptyK_P) = h + 1 ;
	    }
	}
      /* end  For each variable j */

      /* Update proportions as  ph = nh / n */
      NewPh_k[ h ] = Nh_k[ h ] / N ;
    }
  /* end  For each class h */


  /* Free local structures */
  GenFree( oldmean_hj_kd ) ;
  GenFree( sumsquare_hj_kd ) ;
  GenFree( sumdata_hj_kd  ) ;

  return sts ;

}  /*  end of CommonGaussDiagMissing() */





/* ------------------------------------------------------------------- */
void SetIdMatrix ( int Nd , float* M )
/* Set matrix M to identity */ /*V1.03-a*/
/* ------------------------------------------------------------------- */
{
  int l , c ;  /* line and column counters : 0..Nd-1 */

  for ( l = 0 ; l < Nd ; l ++ )
    {
      for ( c = 0 ; c < Nd ; c ++ )
        {
          if ( c == l )
            M[ (l * Nd) + c ] = 1.0 ;
          else
            M[ (l * Nd) + c ] = 0.0 ;
        }
    }
}


#endif

/* ------------------------------------------------------------------- */


/* ~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
