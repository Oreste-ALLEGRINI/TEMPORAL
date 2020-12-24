#ifndef KRISTAL_H
#define KRISTAL_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _kristal{
    int32_t version;
    int32_t matiere;   // Matiere  MAT_LYSO ou MAT_CEBR
    int32_t indice;    // Indice du cristal multiplie par LRS_MULT_INDICE pour entier
    int32_t energie;   // Energie de mesure du photopic (511000 eV) en eV
    int32_t photopik;  // Position du photo pic 511 keV en photons
    int32_t epaisseur; // Epaisseur du cristal en mm
}Kristal;

#ifdef __cplusplus
}
#endif

#endif // KRISTAL_H
