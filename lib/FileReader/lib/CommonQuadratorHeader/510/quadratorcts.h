#ifndef CONSTANTESPET_H
#define CONSTANTESPET_H

#ifndef ENERGY0
#define ENERGY0 511 // E0 = 511 KeV
#endif

#ifndef DIM_DIE_X
#define DIM_DIE_X 4
#endif

#ifndef DIM_DIE_Y
#define DIM_DIE_Y 4
#endif

#ifndef DIM_PIX_X
#define DIM_PIX_X (2*DIM_DIE_X)
#endif

#ifndef DIM_PIX_Y
#define DIM_PIX_Y (2*DIM_DIE_Y)
#endif

#ifndef DIM_TUILE
#define DIM_TUILE 4
#endif

#ifndef DIM_TUILE_PET
#define DIM_TUILE_PET 2
#endif

/**** ****/
#ifndef DIM_Z
#define DIM_Z 41
#endif

#ifndef DIM_XBAR
#define DIM_XBAR 15
#endif

#ifndef DIM_YBAR
#define DIM_YBAR 15
#endif

#ifndef DIM_SATVAL
#define DIM_SATVAL 1059
#endif

#define LRS_DESAT 2880.0
#define LRS_ERR_TEMPER -100000

#define LRS_VERSION 510

// Mode 1 tuile
#define LRS_MODE_1TUILE 0x00000001
#define LRS_MODE_SPECTRO 0x00000002

// Mode 2 tuiles
#define LRS_MODE_2TUILE 0x00000010
#define LRS_MODE_1XPET  0x00000020
#define LRS_MODE_1XCOMP 0x00000040
#define LRS_MODE_1XSAND 0x00000080

// Mode 3 tuiles
#define LRS_MODE_3TUILE 0x00000100
#define LRS_MODE_SYNCHRO 0x00000200

// Mode 4 tuiles
#define LRS_MODE_4TUILE 0x00001000
#define LRS_MODE_2XPET  0x00002000
#define LRS_MODE_2XCOMP 0x00004000
#define LRS_MODE_2XSAND 0x00008000

#define MSQ_PTEMPS  0x00000001
#define MSQ_PPHOTON 0x00000002
#define MSQ_SATUR   0x00000004
#define MSQ_CENTRE  0x00000008
#define MSQ_VOISIN  0x00000010
#define MSQ_VEKTOR  0x00000020
#define MSQ_PILOTE  0x00000040
#define MSQ_ACTIVE  0x00000100
#define MSQ_PATOR   0x00001000
#define MSQ_ALL_IMG 0x00010000

#define LRS_TIMESTAMP ((5.0e-9/256.0)*1.0e12)

#define LRS_BORNE_CRT_INF 459900
#define LRS_BORNE_CRT_SUP 562100

#define MAT_LYSO 0x0001
#define MAT_CEBR 0x0002

#define LRS_MULT_INDICE 1000.0
#define LRS_MULT_VEKTOR 1000.0


#define LRS_HISTO_FREQ 100
// 5÷(327,68×10⁻⁶)
#define LRS_CINQS_NFRAME 15250

#endif // CONSTANTESPET_H

