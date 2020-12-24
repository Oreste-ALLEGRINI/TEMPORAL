#ifndef TCPPROTOCOL_H
#define TCPPROTOCOL_H

#define TIMERINTERVAL 60000
#define TIMERTO 300000

#define LRS_STATUS_ALUM     0x00000001 // Allumage
#define LRS_STATUS_CONF     0x00000002 // Configuration
#define LRS_STATUS_POWU     0x00000003 // Power Up
#define LRS_STATUS_VBIA     0x00000004 // Vbias scan
#define LRS_STATUS_CSET     0x00000005 // Config set
#define LRS_STATUS_CAPA     0x00000006 // Capture activate
#define LRS_STATUS_CAPT     0x00000007 // Capture --sec
#define LRS_STATUS_CAPD     0x00000008 // Capture deactivate
#define LRS_STATUS_POWD     0x00000009 // Power down
#define LRS_STATUS_PHOT     0x00000010 // Savoir si une photo est possible ou non (si LRS_STATUS_ERR_MASQ = erreur lors de la prise de photo si LRS_STATUS_INT_MASQ = il n'y a pas d'appareil photos)
#define LRS_STATUS_VERS     0x00000011 // Disconect
#define LRS_STATUS_CZT      0x00000012 // Disconect
#define LRS_STATUS_TELE     0x00000013 // Disconect
#define LRS_STATUS_COPB     0x00000014

#define LRS_STATUS_ERR_MASQ 0x80000000 //Erreur logicielle
#define LRS_STATUS_INT_MASQ 0x40000000 //Action interdite (capture sans avoir fait de configset par exemple (surtout du debug))
#define LRS_STATUS_PENDING_MASQ 0x20000000 //Action interdite (capture sans avoir fait de configset par exemple (surtout du debug))

#define LRS_ETAT_ALUM     0x00000001 // Allumage
#define LRS_ETAT_POWU     0x00000002 // Power Up
#define LRS_ETAT_CONF     0x00000004 // Configuration
#define LRS_ETAT_VBIA     0x00000008 // Vbias scan
#define LRS_ETAT_CSET     0x00000010 // Config set
#define LRS_ETAT_CAPA     0x00000020 // Capture activate
#define LRS_ETAT_CAPT     0x00000040 // Capture --sec

#define LRS_TEMPO_TEMPE 30000

const quint32 magik = 0x34363631; // nombre magique pour les petites trames

const quint32 SALU=0x53414c55; // Trame Salutation
const quint32 ULAS=0x554c4153;
const quint32 salu=0x73616c75;
const quint32 ulas=0x756c6173;

const quint32 CONF=0x434f4e46; // commande Configuration
const quint32 FNOC=0x464e4f43;
const quint32 conf=0x636f6e66;
const quint32 fnoc=0x666e6f63;

const quint32 EVEN=0x4556454e; // Trame evenement
const quint32 NEVE=0x4e455645;
const quint32 even=0x6576656e;
const quint32 neve=0x6e657665;

const quint32 IMAG=0x47414d49; // Trame image
const quint32 GAMI=0x494d4147;
const quint32 imag=0x67616d69;
const quint32 gami=0x696d6167;

const quint32 QUIT=0x51554954; // Trame Quitter
const quint32 TIUQ=0x54495551;
const quint32 quit=0x71756974;
const quint32 tiuq=0x74697571;

const quint32 INFO=0x4f464e49; // Trame information
const quint32 ONFI=0x494e464f;
const quint32 info=0x6f666e69;
const quint32 ofni=0x696e666f;

const quint32 STAT=0x54415453; // Trame status
const quint32 TATS=0x53544154;
const quint32 xstat=0x74617473;
const quint32 tats=0x73746174;

const quint32 ALLU=0x554c4c41; // commande Allumage
const quint32 ULLA=0x414c4c55;
const quint32 allu=0x756c6c61;
const quint32 ulla=0x616c6c75;

const quint32 CAPT=0x54504143; // commande Acquistion
const quint32 TPAC=0x43415054;
const quint32 capt=0x74706163;
const quint32 tpac=0x63617074;

const quint32 FCOM=0x4d4f4346; // Trame fin compile
const quint32 MOCF=0x46434f4d;
const quint32 fcom=0x6d6f6366;
const quint32 mocf=0x66636f6d;

const quint32 ETEI=0x49455445; // commande Eteindre
const quint32 IETE=0x45544549;
const quint32 etei=0x69657465;
const quint32 iete=0x65746569;

const quint32 CHRO =0x4f524843; // Trame chrono
const quint32 ORHC =0x4348524f;
const quint32 chro =0x6f726863;
const quint32 ohrc =0x6372686f;

const quint32 TEMP=0x504d4554; // Trame temperature
const quint32 PMET=0x54454d50;
const quint32 temp=0x706d6574;
const quint32 pmet=0x74656d70;

const quint32 DETE=0x45544544; // commande Demarre temperature
const quint32 ETED=0x44455445;
const quint32 dete=0x65746564;
const quint32 eted=0x64657465;

const quint32 ARTE=0x45545241; // commande Arrete temperature
const quint32 ETRA=0x41525445;
const quint32 arte=0x65747261;
const quint32 etra=0x61727465;


const quint32 ERCO=0x4552434f; // commande Erreur Client connecté
const quint32 OCRE=0x4f435245;
const quint32 erco=0x6572636f;
const quint32 ocre=0x6f637265;

const quint32 phot=0x746f6870; // Prendre une photo
const quint32 tohp=0x70686f74;
const quint32 PHOT=0x544f4850;
const quint32 TOHP=0x50484f54;

const quint32 PING=0x50494e47; //Ping du serveur
const quint32 GNIP=0x47e44950;
const quint32 ping=0x70696e67; //Ping du Clt
const quint32 gnip=0x676e6970;

const quint32 PONG=0x504f4e47; //Reponse du serveur
const quint32 GNOP=0x474e4f50;
const quint32 pong=0x706f6e67; //Reponse du Clt
const quint32 gnop=0x676e6f70;

const quint32 SPEC=0x53504543;// Envoi du spectre CZT
const quint32 CEPS=0x43455053;
const quint32 spec=0x73706563;// Réponse du client (pas d'implémentation)
const quint32 ceps=0x63657073;

const quint32 TELE=0x54454c45;// Envoi du spectre CZT
const quint32 ELET=0x454c4554;
const quint32 tele=0x74656c65;// Réponse du client (pas d'implémentation)
const quint32 elet=0x656c6574;

#endif /* end of include guard TCPPROTOCOL_H */

