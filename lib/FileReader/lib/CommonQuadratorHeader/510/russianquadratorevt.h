#ifndef RUSSIANQUADRATOREVT_H
#define RUSSIANQUADRATOREVT_H

#include "quadratorcts.h"
#include "quadratorevt.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct _TuileEvt_r{
  int32_t tuile;      // Numero de tuile
  int32_t photons;    // Nombre de photon
  int32_t nbDie;      // Nombre de die allume
  int32_t energie;   // Energie en eV

  int32_t temper;

  int32_t mX;         // Valeur entiere de X * LRS_MULT_VEKTOR
  int32_t mY;         // Valeur entiere de Y * LRS_MULT_VEKTOR
  int32_t mZ;         // Valeur entiere de Z * LRS_MULT_VEKTOR
}TuileEvt_r;


typedef struct _PaireDam_r{
  int numero;
  TuileEvt_r tuiles[DIM_TUILE];
}PaireDam_r;

//Transforme un évènement PaireDam en PaireDam_r, qui pourra etre enregistré.
inline void pDam2pDam_l(PaireDam *pDam, PaireDam_r *pDam_r){
  //TODO
  pDam_r->numero=pDam->numero;
  for(int i=0;i<DIM_TUILE;i++){

    pDam_r->tuiles[i].tuile=pDam->tuiles[i].tuile;
    pDam_r->tuiles[i].photons=pDam->tuiles[i].photons;
    pDam_r->tuiles[i].nbDie=pDam->tuiles[i].nbDie;
    pDam_r->tuiles[i].energie=pDam->tuiles[i].energie;
    pDam_r->tuiles[i].temper=pDam->tuiles[i].temper;

    pDam_r->tuiles[i].mX=pDam->tuiles[i].mX;
    pDam_r->tuiles[i].mY=pDam->tuiles[i].mY;
    pDam_r->tuiles[i].mZ=pDam->tuiles[i].mZ;
  }
}

//Transforme un évènement PaireDam en PaireDam_r, qui pourra etre enregistré.
inline void pDam_l2pDam(PaireDam_r *pDam_r, PaireDam *pDam){
  //TODO
  pDam->numero=pDam_r->numero;
  memset(pDam, 0, sizeof(PaireDam));
  for(int i=0;i<DIM_TUILE;i++){

    pDam->tuiles[i].tuile=pDam_r->tuiles[i].tuile;
    pDam->tuiles[i].photons=pDam_r->tuiles[i].photons;
    pDam->tuiles[i].nbDie=pDam_r->tuiles[i].nbDie;
    pDam->tuiles[i].energie=pDam_r->tuiles[i].energie;
    pDam->tuiles[i].temper=pDam_r->tuiles[i].temper;

    pDam->tuiles[i].mX=pDam_r->tuiles[i].mX;
    pDam->tuiles[i].mY=pDam_r->tuiles[i].mY;
    pDam->tuiles[i].mZ=pDam_r->tuiles[i].mZ;
  }
}

#ifdef __cplusplus
}
#endif


#endif // RUSSIANQUADRATOREVT_H

