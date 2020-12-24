#ifndef DAMEVENT_H
#define DAMEVENT_H

#include <stdint.h>
#include "quadratorcts.h"
#include <string.h>
#include <QDataStream>
#include <math.h>
#include <stdlib.h>
#include "json.hpp"
#ifndef NO_YAML
#include <yaml-cpp/yaml.h>
#endif
#include <QDebug>


using json=nlohmann::json;


typedef struct _TuileEvt_l{
  int32_t tuile       = 0;// Numero de tuile
  int32_t photons     = 0;// Nombre de photon
  uint32_t htimestamp = 0;// Timestamps (numero de frame valide)
  uint32_t ltimestamp = 0;// Timestamps (frame<<24 | timestamp)
  int32_t nbDie       = 0;// Nombre de die allume
  int32_t energie     = 0;// Energie en eV

  int32_t temper      = 0;//temperature

  int32_t mX          = 0;// Valeur entiere de X * LRS_MULT_VEKTOR
  int32_t mY          = 0;// Valeur entiere de Y * LRS_MULT_VEKTOR
  int32_t mZ          = 0;// Valeur entiere de Z * LRS_MULT_VEKTOR
}TuileEvt_l;


typedef struct _PaireDam_l{
  int        numero = 0;
  TuileEvt_l tuiles[DIM_TUILE];
}PaireDam_l;

typedef struct _TuileDam_l{
  int             numero = 0;
  TuileEvt_l      tuile;
}TuileDam_l;

//Transforme un évènement PaireDam en PaireDam_r, qui pourra etre enregistré.

typedef struct _TuileEvt{
  int32_t  tuile      = 0;// Numero de tuile
  uint32_t htimestamp = 0;// Timestamps (numero de frame valide)
  uint32_t ltimestamp = 0;// Timestamps (frame<<24 | timestamp)
  int32_t  temper     = 0;// Temperature
  int32_t  photons    = 0;// Nombre de photon
  int32_t  nbDie      = 0;// Nombre de die allume
  int32_t  energie    = 0;// Energie en eV
  int32_t  tmax       = 0;// Timestamp relatif dernier die allume
  int32_t  valide     = 0;// Valide pour le calcul de vektor

  int32_t  mX         = 0;// Valeur entiere de X * LRS_MULT_VEKTOR
  int32_t  mY         = 0;// Valeur entiere de Y * LRS_MULT_VEKTOR
  int32_t  mZ         = 0;// Valeur entiere de Z * LRS_MULT_VEKTOR
  int32_t  mT         = 0;// !!!
  _TuileEvt():tuile(0),htimestamp(0),ltimestamp(0),temper(0),photons(0),nbDie(0),energie(0),tmax(0),valide(0),mX(0),mY(0),mZ(0),mT(0) {}
}TuileEvt;

typedef struct _TuileImg{
  int16_t g_pix[DIM_PIX_X][DIM_PIX_Y] = {}; // Grille des photons par pixels
  int16_t g_die[DIM_DIE_X][DIM_DIE_Y] = {}; // Grille des timestamp par die
}TuileImg;

typedef struct _TuileDam{
  int      numero = 0;
  TuileEvt tuile;
  TuileImg image;
}TuileDam;


typedef struct _PaireDam{
  int      numero = 0;
  TuileEvt tuiles[DIM_TUILE];
  TuileImg images[DIM_TUILE];
}PaireDam;


inline void videQuadratorDam(TuileDam *tuileDam)
{
  tuileDam->numero=0;
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

template<class PaireEvent = PaireDam>
inline json paireDam2json(PaireEvent evt){
  json j;
  j["eventNumber"] = evt.numero;
  j["tiles"][0]["tileNumber"] = evt.tuiles[0].tuile;
  j["tiles"][0]["temperature"] = evt.tuiles[0].temper;
  j["tiles"][0]["htimestamp"] = evt.tuiles[0].htimestamp;
  j["tiles"][0]["ltimestamp"] = evt.tuiles[0].ltimestamp;
  j["tiles"][0]["X"] = evt.tuiles[0].mX;
  j["tiles"][0]["Y"] = evt.tuiles[0].mY;
  j["tiles"][0]["Z"] = evt.tuiles[0].mZ;
  j["tiles"][0]["photons"] = evt.tuiles[0].photons;
  j["tiles"][0]["energie"] = evt.tuiles[0].energie;
  j["tiles"][0]["nbDie"] = evt.tuiles[0].nbDie;
  j["tiles"][1]["tileNumber"] = evt.tuiles[1].tuile;
  j["tiles"][1]["temperature"] = evt.tuiles[1].temper;
  j["tiles"][1]["htimestamp"] = evt.tuiles[1].htimestamp;
  j["tiles"][1]["ltimestamp"] = evt.tuiles[1].ltimestamp;
  j["tiles"][1]["X"] = evt.tuiles[1].mX;
  j["tiles"][1]["Y"] = evt.tuiles[1].mY;
  j["tiles"][1]["Z"] = evt.tuiles[1].mZ;
  j["tiles"][1]["photons"] = evt.tuiles[1].photons;
  j["tiles"][1]["energie"] = evt.tuiles[1].energie;
  j["tiles"][1]["nbDie"] = evt.tuiles[1].nbDie;
  return j;
}


template<class PaireEvent = PaireDam>
inline PaireEvent json2PaireDam(json evt){
  PaireEvent dam;
  try{
    evt.at("eventNumber").get_to(dam.numero);
    evt.at("tiles").at(0).at("tileNumber").get_to(dam.tuiles[0].tuile);
    evt.at("tiles").at(0).at("temperature").get_to(dam.tuiles[0].temper);
    evt.at("tiles").at(0).at("htimestamp").get_to(dam.tuiles[0].htimestamp);
    evt.at("tiles").at(0).at("ltimestamp").get_to(dam.tuiles[0].ltimestamp);
    evt.at("tiles").at(0).at("X").get_to(dam.tuiles[0].mX);
    evt.at("tiles").at(0).at("Y").get_to(dam.tuiles[0].mY);
    evt.at("tiles").at(0).at("Z").get_to(dam.tuiles[0].mZ);
    evt.at("tiles").at(0).at("photons").get_to(dam.tuiles[0].photons);
    evt.at("tiles").at(0).at("energie").get_to(dam.tuiles[0].energie);
    evt.at("tiles").at(0).at("nbDie").get_to(dam.tuiles[0].nbDie);
    evt.at("tiles").at(1).at("tileNumber").get_to(dam.tuiles[1].tuile);
    evt.at("tiles").at(1).at("temperature").get_to(dam.tuiles[1].temper);
    evt.at("tiles").at(1).at("htimestamp").get_to(dam.tuiles[1].htimestamp);
    evt.at("tiles").at(1).at("ltimestamp").get_to(dam.tuiles[1].ltimestamp);
    evt.at("tiles").at(1).at("X").get_to(dam.tuiles[1].mX);
    evt.at("tiles").at(1).at("Y").get_to(dam.tuiles[1].mY);
    evt.at("tiles").at(1).at("Z").get_to(dam.tuiles[1].mZ);
    evt.at("tiles").at(1).at("photons").get_to(dam.tuiles[1].photons);
    evt.at("tiles").at(1).at("energie").get_to(dam.tuiles[1].energie);
    evt.at("tiles").at(1).at("nbDie").get_to(dam.tuiles[1].nbDie);
  }catch(nlohmann::detail::out_of_range e){
    dam.numero = -1;//Mise en erreur du numéro de séquence, car nous n'avons pas de retour booléen sur celle ci
    return dam;
  }
  return dam;

}

inline int tuileDam2PaireDam(TuileDam *tuileDam, PaireDam *paireDam, int numtuile){
  paireDam->numero = tuileDam->numero;
  for(int i = 0; i< 4; i++){
    if(paireDam->tuiles[i].tuile == numtuile){
      paireDam->tuiles[i] = tuileDam->tuile;
      paireDam->images[i] = tuileDam->image;
      return 1;
    }
  }
  return 0;
}

inline int tuileDam_l2PaireDam_l(TuileDam_l *tuileDam, PaireDam_l *paireDam, int numtuile){
  paireDam->numero = tuileDam->numero;
  for(int i = 0; i< 4; i++){
    if(paireDam->tuiles[i].tuile == numtuile){
      paireDam->tuiles[i] = tuileDam->tuile;
      return 1;
    }
  }
  return 0;
}

inline int paireDam2TuileDam(PaireDam *paireDam, TuileDam *toReturn, int numtuile){
  toReturn->numero = paireDam->numero;
  for(int i = 0; i < 4; i++){
    if(paireDam->tuiles[i].tuile == numtuile){
      toReturn->tuile = paireDam->tuiles[i];
      toReturn->image = paireDam->images[i];
      return 1;
    }
  }

  return 0;
}
inline int paireDam_l2TuileDam_l(PaireDam_l *paireDam, TuileDam_l *toReturn, int numtuile){
  toReturn->numero = paireDam->numero;
  for(int i = 0; i < 4; i++){
    if(paireDam->tuiles[i].tuile == numtuile){
      toReturn->tuile = paireDam->tuiles[i];
      return 1;
    }
  }
  return 0;
}

inline void tDam2tDam_l(TuileDam *tDam, TuileDam_l *tDam_l){
  //TODO
  tDam_l->numero=tDam->numero;
  tDam_l->tuile.tuile=tDam->tuile.tuile;
  tDam_l->tuile.photons=tDam->tuile.photons;
  tDam_l->tuile.nbDie=tDam->tuile.nbDie;
  tDam_l->tuile.energie=tDam->tuile.energie;
  tDam_l->tuile.temper=tDam->tuile.temper;
  tDam_l->tuile.htimestamp=tDam->tuile.htimestamp;
  tDam_l->tuile.ltimestamp=tDam->tuile.ltimestamp;

  tDam_l->tuile.mX=tDam->tuile.mX;
  tDam_l->tuile.mY=tDam->tuile.mY;
  tDam_l->tuile.mZ=tDam->tuile.mZ;
}

//Transforme un évènement PaireDam en PaireDam_r, qui pourra etre enregistré.
inline void tDam_l2tDam(TuileDam_l *tDam_l, TuileDam *tDam){
  //TODO
  tDam->numero=tDam_l->numero;
  tDam->tuile.tuile=tDam_l->tuile.tuile;
  tDam->tuile.photons=tDam_l->tuile.photons;
  tDam->tuile.nbDie=tDam_l->tuile.nbDie;
  tDam->tuile.energie=tDam_l->tuile.energie;
  tDam->tuile.temper=tDam_l->tuile.temper;
  tDam->tuile.htimestamp=tDam_l->tuile.htimestamp;
  tDam->tuile.ltimestamp=tDam_l->tuile.ltimestamp;

  tDam->tuile.mY=tDam_l->tuile.mY;
  tDam->tuile.mZ=tDam_l->tuile.mZ;
  tDam->tuile.mX=tDam_l->tuile.mX;
}


//Transforme un évènement PaireDam en PaireDam_r, qui pourra etre enregistré.
inline void pDam2pDam_l(PaireDam *pDam, PaireDam_l *pDam_l){
  //TODO
  pDam_l->numero=pDam->numero;
  for(int i=0;i<DIM_TUILE;i++){
    pDam_l->tuiles[i].tuile=pDam->tuiles[i].tuile;
    pDam_l->tuiles[i].photons=pDam->tuiles[i].photons;
    pDam_l->tuiles[i].nbDie=pDam->tuiles[i].nbDie;
    pDam_l->tuiles[i].energie=pDam->tuiles[i].energie;
    pDam_l->tuiles[i].temper=pDam->tuiles[i].temper;
    pDam_l->tuiles[i].htimestamp=pDam->tuiles[i].htimestamp;
    pDam_l->tuiles[i].ltimestamp=pDam->tuiles[i].ltimestamp;

    pDam_l->tuiles[i].mX=pDam->tuiles[i].mX;
    pDam_l->tuiles[i].mY=pDam->tuiles[i].mY;
    pDam_l->tuiles[i].mZ=pDam->tuiles[i].mZ;
  }
}

//Transforme un évènement PaireDam en PaireDam_r, qui pourra etre enregistré.
inline void pDam_l2pDam(PaireDam_l *pDam_l, PaireDam *pDam){
  //TODO
  memset(pDam, 0, sizeof(PaireDam));
  pDam->numero=pDam_l->numero;
  for(int i=0;i<DIM_TUILE;i++){
    pDam->tuiles[i].tuile=pDam_l->tuiles[i].tuile;
    pDam->tuiles[i].photons=pDam_l->tuiles[i].photons;
    pDam->tuiles[i].nbDie=pDam_l->tuiles[i].nbDie;
    pDam->tuiles[i].energie=pDam_l->tuiles[i].energie;
    pDam->tuiles[i].temper=pDam_l->tuiles[i].temper;
    pDam->tuiles[i].htimestamp=pDam_l->tuiles[i].htimestamp;
    pDam->tuiles[i].ltimestamp=pDam_l->tuiles[i].ltimestamp;

    pDam->tuiles[i].mX=pDam_l->tuiles[i].mX;
    pDam->tuiles[i].mY=pDam_l->tuiles[i].mY;
    pDam->tuiles[i].mZ=pDam_l->tuiles[i].mZ;
  }
}

inline QDataStream &operator<<(QDataStream& out, const PaireDam_l& dam){
    out.writeBytes(reinterpret_cast<const char*>(&dam), sizeof (PaireDam_l));
    return out;
}

inline QDataStream &operator>>(QDataStream& in, PaireDam_l& dam){
    QByteArray param;
    in >> param;
    std::memcpy(&dam, param.data(), sizeof (PaireDam_l));
    return in;
}


inline QDataStream &operator<<(QDataStream& out, const TuileDam_l& dam){
    out.writeBytes(reinterpret_cast<const char*>(&dam), sizeof (TuileDam_l));
    return out;
}

inline QDataStream &operator>>(QDataStream& in, TuileDam_l& dam){
    QByteArray param;
    in >> param;
    std::memcpy(&dam, param.data(), sizeof (TuileDam_l));
    return in;
}

inline QDataStream &operator<<(QDataStream& out, const TuileDam& dam){
    out.writeBytes(reinterpret_cast<const char*>(&dam), sizeof (TuileDam));
    return out;
}

inline QDataStream &operator>>(QDataStream& in, TuileDam& dam){
    QByteArray param;
    in >> param;
    std::memcpy(&dam, param.data(), sizeof (TuileDam));
    return in;
}

inline QDataStream &operator<<(QDataStream& out, const PaireDam& dam){
    out.writeBytes(reinterpret_cast<const char*>(&dam), sizeof (PaireDam));
    return out;
}

inline QDataStream &operator>>(QDataStream& in, PaireDam& dam){
    QByteArray param;
    in >> param;
    std::memcpy(&dam, param.data(), sizeof (PaireDam));
    return in;
}

// Numero de tuile :  xx
// Nombre de photon : xx
// timestamp = ....
// nbDie
// energie ??,
// temper
//     i1        i2        i3         i4
//   +=========+=========+=========+=========+
// x0|-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|
//   |----ts---|----x----|----x----|----x----| j1
// x1|-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|
//   +=========+=========+=========+=========+
// . |-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|
// . |----x----|----x----|----x----|----x----| j2
// . |-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|
// . +=========+=========+=========+=========+
//   |-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|
//   |----x----|----x----|----x----|----x----| j3
//   |-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|
//   +=========+=========+=========+=========+
//   |-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|
//   |----x----|----x----|----x----|----x----| j4
// x7|-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|-XXX-XXX-|
//   +=========+=========+=========+=========+
//      y0  y1   ...                ...  y7
// XXX nb de photon
inline std::string prettyPrint(TuileDam &even)
{
  std::string result;
  int width   = 9;
  char corner = '+';
  char hline  = '=';
  char vline  = '|';
  char sep    = '.';
  result = "\n/**********************\\";
  result += "\n*Numérop d'événement:" + std::to_string(even.numero);
  result += "\n*Numéro de tuile:" + std::to_string(even.tuile.tuile);
  long long int ts = 0;
  ts = even.tuile.ltimestamp << 24;
  ts |= even.tuile.htimestamp;
  result+= "\n*TimeStamp: " + std::to_string(ts);
  result+= "\n*Temperature: " + std::to_string(even.tuile.temper);
  result+= "\n*Photons: " + std::to_string(even.tuile.photons);
  result+= "\n*nbDie: " + std::to_string(even.tuile.nbDie);
  result+= "\n*Energie: " + std::to_string(even.tuile.energie);
  result+= "\n*tMax: " + std::to_string(even.tuile.tmax);
  result+= "\n";

  // ligne horyzontale +===+==+
  auto hLineSep = [&]() {
    for(int i=0; i<4;i++){
      result.push_back(corner);
      result.append(width,hline);
    }
    result.push_back(corner);
    result.push_back('\n');
  };
  // ligne horizontale de pixel
  auto hLinePixel = [&](int k){
    for(int i=0;i<4;i++){
      result.push_back(vline);
      std::string tmp  = std::to_string(even.image.g_pix[2*i][k]);
      std::string tmp2 = std::to_string(even.image.g_pix[2*i+1][k]);
      int tmpSize  = tmp.size();
      int tmp2Size = tmp2.size();
      result.append(((width)/2-tmpSize),sep);
      result.append(tmp);
      result.push_back(sep);
      result.append((width)/2-tmp2Size,sep);
      result.append(tmp2);
    }
    result.push_back(vline);
    result.push_back('\n');
  };

  // ligne horyzontale timestamp
  auto hLineDie = [&](int k){
    for(int i=0;i<4;i++){
      std::string tmp = std::to_string(even.image.g_die[i][k]);
      int tmpSize = tmp.size();
      result.push_back(vline);
      result.append((width-tmpSize)/2 + (tmpSize+1)%2,sep);
      result.append(tmp);
      result.append((width-tmpSize)/2,sep);
    }
    result.push_back(vline);
    result.push_back('\n');
  };

  //parcours des dies
  for(int j=0;j<4;j++){
    hLineSep();
    hLinePixel(2*j);
    hLineDie(j);
    hLinePixel(2*j+1);
  }
  hLineSep();
  return result;
}


inline QDebug operator<< (QDebug dbg, TuileDam &dam){
  QDebugStateSaver stateSaver(dbg);
  dbg << prettyPrint(dam).c_str();
  return dbg;
}

inline std::ostream & operator<< (std::ostream &out, TuileDam &dam){
  out << prettyPrint(dam);
  return out;
}

inline std::string prettyPrint(PaireDam &even)
{
  std::string result;
  int width   = 9;
  char corner = '+';
  char hline  = '=';
  char vline  = '|';
  char sep    = '.';
  result = "\n/**********************\\";
  result += "\n*Numérop d'événement:" + std::to_string(even.numero);
  for(int i =0; i< 2; i++){
    result += "\n*Numéro de tuiles[i]:" + std::to_string(even.tuiles[i].tuile);
    long long int ts = 0;
    ts = even.tuiles[i].ltimestamp << 24;
    ts |= even.tuiles[i].htimestamp;
    result+= "\n*TimeStamp: " + std::to_string(ts);
    result+= "\n*Temperature: " + std::to_string(even.tuiles[i].temper);
    result+= "\n*Photons: " + std::to_string(even.tuiles[i].photons);
    result+= "\n*nbDie: " + std::to_string(even.tuiles[i].nbDie);
    result+= "\n*Energie: " + std::to_string(even.tuiles[i].energie);
    result+= "\n*tMax: " + std::to_string(even.tuiles[i].tmax);
    result+= "\n*mX: " + std::to_string(even.tuiles[i].mX);
    result+= "\n*mY: " + std::to_string(even.tuiles[i].mY);
    result+= "\n*mZ: " + std::to_string(even.tuiles[i].mZ);
    result+= "\n";

    // ligne horyzontale +===+==+
    auto hLineSep = [&]() {
      for(int j=0; j<4;j++){
        result.push_back(corner);
        result.append(width,hline);
      }
      result.push_back(corner);
      result.push_back('\n');
    };
    // ligne horizontale de pixel
    auto hLinePixel = [&](int k){
      for(int j=0;j<4;j++){
        result.push_back(vline);
        std::string tmp  = std::to_string(even.images[i].g_pix[2*j][k]);
        std::string tmp2 = std::to_string(even.images[i].g_pix[2*j+1][k]);
        int tmpSize  = tmp.size();
        int tmp2Size = tmp2.size();
        result.append(((width)/2-tmpSize),sep);
        result.append(tmp);
        result.push_back(sep);
        result.append((width)/2-tmp2Size,sep);
        result.append(tmp2);
      }
      result.push_back(vline);
      result.push_back('\n');
    };

    // ligne horyzontale timestamp
    auto hLineDie = [&](int k){
      for(int j=j;j<4;i++){
        std::string tmp = std::to_string(even.images[i].g_die[j][k]);
        int tmpSize = tmp.size();
        result.push_back(vline);
        result.append((width-tmpSize)/2 + (tmpSize+1)%2,sep);
        result.append(tmp);
        result.append((width-tmpSize)/2,sep);
      }
      result.push_back(vline);
      result.push_back('\n');
    };

    //parcours des dies
    for(int j=0;j<4;j++){
      hLineSep();
      hLinePixel(2*j);
      hLineDie(j);
      hLinePixel(2*j+1);
    }
    hLineSep();
  }
  return result;
}


inline QDebug operator<< (QDebug dbg, PaireDam &dam){
  QDebugStateSaver stateSaver(dbg);
  dbg << prettyPrint(dam).c_str();
  return dbg;
}

inline std::ostream & operator<< (std::ostream &out, PaireDam &dam){
  out << prettyPrint(dam);
  return out;
}

#ifndef NO_YAML
namespace YAML {
template<>
  struct convert<TuileEvt> {
  static Node encode(const TuileEvt& te) {
      Node node;
      node["tuile"]	 = te.tuile;
      node["htimestamp"] = te.htimestamp;
      node["ltimestamp"] = te.ltimestamp;
      node["temper"]	 = te.temper;
      node["photons"]	 = te.photons;
      node["nbDie"]	 = te.nbDie;
      node["energie"]	 = te.energie;
      node["tmax"]	 = te.tmax;
      node["valide"]	 = te.valide;
      node["mX"]	 = te.mX;
      node["mY"]	 = te.mY;
      node["mZ"]	 = te.mZ;
      node["mT"]	 = te.mT;
      return node;
    }
  static bool decode(const Node& node, TuileEvt& te) {
    if(!node.IsMap() ) {
      return false;
    }
    te.tuile	  = node["tuile"].as<int32_t>();
    te.htimestamp = node["htimestamp"].as<uint32_t>();
    te.ltimestamp = node["ltimestamp"].as<uint32_t>();
    te.temper	  = node["temper"].as<int32_t>();
    te.photons	  = node["photons"].as<int32_t>();
    te.nbDie	  = node["nbDie"].as<int32_t>();
    te.energie	  = node["energie"].as<int32_t>();
    te.tmax	  = node["tmax"].as<int32_t>();
    te.valide	  = node["valide"].as<int32_t>();
    te.mX	  = node["mX"].as<int32_t>();
    te.mY	  = node["mY"].as<int32_t>();
    te.mZ	  = node["mZ"].as<int32_t>();
    te.mT	  = node["mT"].as<int32_t>();
    return true;
  }
 };

 template<>
  struct convert<TuileImg> {
  static Node encode(const TuileImg& ti) {
      Node node;
      std::array<int16_t,DIM_PIX_X> tmp1;
      for (int i = 0; i< DIM_PIX_Y;i++){
	memcpy(tmp1.data(),ti.g_pix[i], DIM_PIX_X*sizeof(int16_t));
	node.push_back(tmp1);
      }
      std::array<int16_t,DIM_DIE_X> tmp2;
      for (int i = 0; i< DIM_DIE_Y;i++){
	memcpy(tmp2.data(),ti.g_pix[i], DIM_DIE_X*sizeof(int16_t));
	node.push_back(tmp2);
      }

      return node;
    }
  static bool decode(const Node& node, TuileImg& ti) {
    if(!node.IsMap()) {
      return false;
    }
    for (int i = 0; i< DIM_PIX_Y;i++)
      memcpy(ti.g_pix[i],node["g_pix"][i].as<std::array<int16_t,DIM_PIX_X>>().data(), DIM_PIX_X*sizeof(int16_t));

    for (int i = 0; i< DIM_DIE_Y;i++)
      memcpy(ti.g_die[i],node["g_die"][i].as<std::array<int16_t,DIM_DIE_X>>().data(), DIM_DIE_X*sizeof(int16_t));

    return true;
  }
 };

  template<>
    struct convert<TuileDam> {

    static Node encode(const TuileDam& pd) {
      Node node;
      node.push_back(pd.numero);
      node.push_back(pd.tuile);
      node.push_back(pd.image);
      return node;
    }

    static bool decode(const Node& node, TuileDam& pd) {
      if(!node.IsMap() ) {
	return false;
      }

      pd.numero = node["numero"].as<int>();
      pd.tuile	= node["TuileEvt"].as<TuileEvt>();
      pd.image	= node["TuileImg"].as<TuileImg>();
      return true;
    }
  };

  // Emitter
  inline Emitter& operator << (YAML::Emitter& out, const TuileImg& ti) {
    out << YAML::BeginMap;
    out << YAML::Key << "g_pix" << YAML::Value;

    out << YAML::BeginSeq;
    for(int i = 0; i< DIM_PIX_X; i++){
      out << YAML::Flow;
      out << YAML::BeginSeq;
      for(int j = 0;j< DIM_PIX_Y;j++)
	out << ti.g_pix[i][j] ;
      out << YAML::EndSeq;
    }
    out << YAML::EndSeq;
    out << YAML::Key << "g_die" << YAML::Value;
    out << YAML::BeginSeq;
    for(int i = 0; i< DIM_DIE_X; i++){
      out << YAML::Flow;
      out << YAML::BeginSeq;
      for(int j = 0;j< DIM_DIE_Y;j++)
	out << ti.g_die[i][j] ;
      out << YAML::EndSeq;
    }
    out << YAML::EndSeq;
    out << YAML::EndMap;
    return out;
  }

  inline Emitter& operator << (YAML::Emitter& out, const TuileEvt& te) {
    out  << YAML::Node(te) ;
    return out;
  }

  inline Emitter& operator << (YAML::Emitter& out, const TuileDam& td) {
    out   << YAML::BeginMap << YAML::Key << "numero" << YAML::Value << td.numero << YAML::Key << "TuileEvt" << YAML::Value << td.tuile << YAML::Key << "TuileImg" << YAML::Value << td.image << YAML::EndMap;
    return out;
  }

}


#endif
#endif
