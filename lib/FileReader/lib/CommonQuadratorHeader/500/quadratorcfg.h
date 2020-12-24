#ifndef KONFIGURATIONPET
#define KONFIGURATIONPET

#include <stdint.h>

#include "quadratorcts.h"
#include "kristal.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _Ponderetemps{
    int32_t version;
    int32_t ptemps[16];
}Ponderetemps;

typedef struct _Configset{
    int32_t p_tri; // Trigger mode
    int32_t p_x7f; // 0x7f
    int32_t p_val; // validation
    int32_t p_int; // integration
    int32_t n_logic; //neighborLogic (1 activée, 0 sinon)
}Configset;

typedef struct _Kompress{
    int32_t masque; // masque =| MSQ_PTEMPS 
    int32_t e_min;  // energie min en photon
    int32_t e_max;  //         max 
    int32_t d_min;  // nombre de die
}Kompress;

typedef struct _VektorKonf{
    int32_t dt;             // = 30;
    Kristal kristal;
}VektorKonf;

typedef struct _Ctuile{
    int32_t p_til;
    Configset c_set;
    Kompress c_kom;
    Ponderetemps ptemps;
    VektorKonf c_vek;
}Ctuile;


// Taille de .kqua  552
typedef struct _QuadatorKonf{
    int32_t version;
    int32_t p_dur;
    int32_t p_coi;
    int32_t p_smi; // Somme energie Minimun
    int32_t p_sma; // Somme energie Maximun
    int32_t p_mod; // PET 2Tuiles Compton
    int32_t n_tuile; // Nombre de tuile
    Ctuile ctuile[DIM_TUILE];
} QuadatorKonf;

#ifdef __cplusplus
}
#endif

#endif // KONFIGURATIONPET

