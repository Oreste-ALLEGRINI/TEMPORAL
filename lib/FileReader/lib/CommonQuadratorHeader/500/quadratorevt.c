/*
 * Gestion des evenements a la sauce damavan
 */

#include "quadratorevt.h"

#include <string.h>
#include <math.h>

void videQuadratorDam(TuileDam *tuileDam)
{   tuileDam->numero=0;
    memset((void *)&tuileDam->tuile,0,sizeof(tuileDam->tuile));
    memset((void *)&(tuileDam->image.g_pix[0][0]),0,sizeof(tuileDam->image.g_pix));
    memset((void *)&(tuileDam->image.g_die[0][0]),0xFF,sizeof(tuileDam->image.g_die));
    tuileDam->tuile.energie=0;
    tuileDam->tuile.htimestamp = 0;
    tuileDam->tuile.ltimestamp = 0;
    tuileDam->tuile.mT = 0;
    tuileDam->tuile.mX = 0;
    tuileDam->tuile.mY = 0;
    tuileDam->tuile.mZ = 0;
    tuileDam->tuile.nbDie = 0;
    tuileDam->tuile.photons=0;
    tuileDam->tuile.temper = 0;
    tuileDam->tuile.tmax = 0;
    tuileDam->tuile.tuile = 0;
    tuileDam->tuile.valide = 0;
}
