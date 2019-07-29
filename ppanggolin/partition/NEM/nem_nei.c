/*\
    nem_nei.c

    Fonctions de voisinage.
\*/
#include    "nem_nei.h"     /* prototypes */
#include    <stdio.h>       /* fprintf */
#include    <string.h>      /* memcpy */


static int      GetNeighNone                 /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) ;


static int      GetNeighImage                /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) ;


static int      GetNeighIrreg                /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) ;



/* ------------------------------------------------------------------- */
int GetSpatialFunc  /* STS_OK or STS_E_FUNCARG */
    (
     TypeET          SpatialType,    /* I */
     GetNeighFT**    GetNeighFP      /* O */
    ) 
/* ------------------------------------------------------------------- */
{
    switch( SpatialType )
    {
        case TYPE_NONSPATIAL :
             *GetNeighFP = GetNeighNone ;
        return STS_OK ;

        case TYPE_SPATIAL :
             *GetNeighFP = GetNeighIrreg ;
        return STS_OK ;

        case TYPE_IMAGE :
             *GetNeighFP = GetNeighImage ;
        return STS_OK ;

        default :
             *GetNeighFP = NULL ;
                fprintf( stderr, "GetSpatialFuncs bad arg : Type = %d\n",
                         SpatialType ) ;
        return STS_E_FUNCARG ;
    }
}   /* end of GetSpatialFuncs() */




/* ------------------------------------------------------------------- */
int         GetNeighNone                 /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) 
/* ------------------------------------------------------------------- */
{
    return 0 ;
}


/* ------------------------------------------------------------------- */
int         GetNeighImage                /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) 
/* ------------------------------------------------------------------- */
{
    int     in ;
    int     rnn ;       /* real number of neighbours */
    int     nbn     = NeighDataP->Image.NbNeigh ;
    INeighT *neiV   = NeighDataP->Image.NeighsV ;
    int     nl      = NeighDataP->Image.Nl ;
    int     nc      = NeighDataP->Image.Nc ;
    int     line    = Ipt / nc ;
    int     col     = Ipt % nc ;

    if ( nbn > PtNeighsP->NbNeigh )
    {
        nbn = PtNeighsP->NbNeigh ;
    }

    for ( in = 0, rnn = 0 ; in < nbn ; in ++ )
    {
        int l = line + neiV[ in ].Dl ;
        int c = col + neiV[ in ].Dc ;

        if ( ( 0 <= l ) && ( l < nl ) && ( 0 <= c ) && ( c < nc ) )
        {
            PtNeighsP->NeighsV[ rnn ].Index  = l * nc + c ;
            PtNeighsP->NeighsV[ rnn ].Weight = neiV[ in ].Weight ;
            rnn ++ ;
        }
    }

    return rnn ;
}

/* ------------------------------------------------------------------- */
int         GetNeighIrreg                /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) 
/* ------------------------------------------------------------------- */
{
    PtNeighsT*  ptnP = &(NeighDataP->PtsNeighsV[ Ipt ]) ;
    int         nbn = ptnP->NbNeigh ;

    if ( nbn > PtNeighsP->NbNeigh )
    {
        nbn = PtNeighsP->NbNeigh ;
    }

    memcpy( PtNeighsP->NeighsV, ptnP->NeighsV, nbn * sizeof( NeighT ) );

    return nbn ;        
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
