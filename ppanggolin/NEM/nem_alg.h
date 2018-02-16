#include "nem_typ.h"    /* DataT, ... */

/*
  V1.06-a   17-JUN-1998  NoiseModel -> StatModel
  V1.06-b   30-JUN-1998  ClassifyByNem const NemParaP, no const DataP (sort)
  V1.06-c   20-SEP-1998  Add LabelToClassVector
  V1.06-d   21-SEP-1998  Add TieRule arg to ComputeMAP
*/

    int ClassifyByNem
        ( 
          const NemParaT      *NemParaP,        /* I */
          const SpatialT      *SpatialP,        /* I */
          DataT               *DataP,           /* I/O */
          StatModelT          *StatModelP,      /* I/O */
          float               *ClassifM,        /* I/O */
	  CriterT             *CriterP          /* O */  /*V1.03-f*/
        ) ;

    int ComputeMAP
        ( 
          const float* ClassifM,   /* I */
          int          Ipt,        /* I */
          int          Nk,         /* I */
	  TieET        TieRule,    /* I */ /*V1.06-d*/
          int*         kmaxesV     /* T [Nk] */
        ) ;


void LabelToClassVector
 ( 
  const int Nk,    /* I: number of classes */
  const int Label, /* I: class of interest 0..Nk-1 */
  float* Cout_K    /* O: classification vector [Nk] */
 ) ;
