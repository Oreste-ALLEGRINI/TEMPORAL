#ifndef DAMEVENT_H
#define DAMEVENT_H

#include <stdint.h>
#include "quadratorcts.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _TuileEvt{
    int32_t tuile;      // Numero de tuile
    int32_t htimestamp; // Timestamps (numero de frame valide)
    uint32_t ltimestamp; // Timestamps (frame<<24 | timestamp)
    int32_t temper;     // Temperature
    int32_t photons;    // Nombre de photon
    int32_t nbDie;      // Nombre de die allume
    int32_t energie;   // Energie en eV
    int32_t tmax;       // Timestamp relatif dernier die allume
    int32_t valide;     // Valide pour le calcul de vektor

    int32_t mX;         // Valeur entiere de X * LRS_MULT_VEKTOR
    int32_t mY;         // Valeur entiere de Y * LRS_MULT_VEKTOR
    int32_t mZ;         // Valeur entiere de Z * LRS_MULT_VEKTOR
    int32_t mT;         // !!!
}TuileEvt;

typedef struct _TuileImg{
    int16_t g_pix[DIM_PIX_X][DIM_PIX_Y]; // Grille des photons par pixels
    int16_t g_die[DIM_DIE_X][DIM_DIE_Y]; // Grille des timestamp par die
}TuileImg;

typedef struct _TuileDam{
    int numero;
    TuileEvt tuile;
    TuileImg image;
}TuileDam;

typedef struct _PaireDam{
    int numero;
    TuileEvt tuiles[DIM_TUILE];
    TuileImg images[DIM_TUILE];
}PaireDam;

void videQuadratorDam(TuileDam *tuileDam);

#ifdef __cplusplus
}
#endif



#endif
