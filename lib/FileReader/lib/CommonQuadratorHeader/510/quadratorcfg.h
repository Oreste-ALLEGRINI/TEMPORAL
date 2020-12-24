#ifndef KONFIGURATIONPET
#define KONFIGURATIONPET

#include <cstring>
#include <QDataStream>
#include "quadratorcts.h"
#include "kristal.h"
#include <yaml-cpp/yaml.h>

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
  int32_t p_coi;
  int32_t p_smi; // Somme energie Minimun
  int32_t p_sma; // Somme energie Maximun
  int32_t p_mod; // PET 2Tuiles Compton
  int32_t n_tuile; // Nombre de tuile
  Ctuile ctuile[DIM_TUILE];
} QuadatorKonf;

typedef struct _QuadatorKonf QuadratorKonf;


inline QDataStream &operator<<(QDataStream& out, const QuadatorKonf& konf){
    out.writeBytes(reinterpret_cast<const char*>(&konf), sizeof (QuadatorKonf));
    return out;
}

inline QDataStream &operator>>(QDataStream& in, QuadatorKonf& konf){
    QByteArray param;
    in >> param;
    memcpy(&konf, param.data(), sizeof (QuadatorKonf));
    return in;
}


namespace YAML {
  template<>
    /// Pondere temps => yaml
    struct convert<Ponderetemps> {
    static Node encode(const Ponderetemps& pt) {
      Node node;
      node["version"] = pt.version;
      std::array<int32_t,16> tmp;
      memcpy(tmp.data(),pt.ptemps,16*sizeof(int32_t)),
	node["ptemps"]  = tmp;
      return node;
    }
    static bool decode(const Node& node, Ponderetemps& pt) {
      if(!node.IsMap() ) {
	return false;
      }
      pt.version = node["version"].as<int32_t>();
      memcpy(pt.ptemps,
	     node["ptemps"].as<std::array<int32_t,16>>().data(),
	     16*sizeof(int32_t));
      return true;
    }
  };

  /// Config Set => yaml
  template<>
    struct convert<Configset> {
    static Node encode(const Configset& cs) {
      Node node;
      node["p_tri"]   = cs.p_tri;
      node["p_x7f"]   = cs.p_x7f;
      node["p_val"]   = cs.p_val;
      node["p_int"]   = cs.p_int;
      node["n_logic"] = cs.n_logic;
      return node;
    }
    static bool decode(const Node& node, Configset& cs) {
      if(!node.IsMap() )
	return false;
      cs.p_tri	= node["p_tri"].as<int32_t>();
      cs.p_x7f	= node["p_x7f"].as<int32_t>();
      cs.p_val	= node["p_val"].as<int32_t>();
      cs.p_int	= node["p_int"].as<int32_t>();
      cs.n_logic = node["n_logic"].as<int32_t>();
      return true;
    }
  };

  /// Kompress => yaml
  template<>
    struct convert<Kompress> {
    static Node encode(const Kompress& kmp) {
      Node node;
      node["masque"] = kmp.masque;
      node["e_min"]  = kmp.e_min;
      node["e_max"]  = kmp.e_max;
      node["d_min"]  = kmp.d_min;
      return node;
    }
    static bool decode(const Node& node, Kompress& kmp) {
      if(!node.IsMap() )
	return false;
      kmp.masque = node["masque"].as<int32_t>();
      kmp.e_min	= node["e_min"].as<int32_t>();
      kmp.e_max	= node["e_max"].as<int32_t>();
      kmp.d_min	= node["d_min"].as<int32_t>();
      return true;
    }
  };

  /// VektorKonf => to yaml
  template<>
    struct convert<VektorKonf> {
    static Node encode(const VektorKonf& vk) {
      Node node;
      node["dt"]      = vk.dt;
      node["kristal"] = vk.kristal;
      return node;
    }
    static bool decode(const Node& node, VektorKonf& vk) {
      if(!node.IsMap() )
	return false;
      vk.dt	 = node["dt"].as<int32_t>();
      vk.kristal = node["kristal"].as<Kristal>();
      return true;
    }
  };

  /// Ctuile => to yaml
  template<>
    struct convert<Ctuile> {
    static Node encode(const Ctuile& ct) {
      Node node;
      node["p_til"]  = ct.p_til ;
      node["c_set"]  = ct.c_set ;
      node["c_kom"]  = ct.c_kom ;
      node["ptemps"] = ct.ptemps ;
      node["c_vek"]  = ct.c_vek ;
      return node;
    }
    static bool decode(const Node& node, Ctuile& ct) {
      if(!node.IsMap() )
	return false;
      ct.p_til  = node["p_til"].as<int32_t>();
      ct.c_set  = node["c_set"].as<Configset>();
      ct.c_kom  = node["c_kom"].as<Kompress>();
      ct.ptemps = node["ptemps"].as<Ponderetemps>();
      ct.c_vek  = node["c_vek"].as<VektorKonf>();
      return true;
    }
  };

  /// QuadatorKonf => to yaml
  template<>
    struct convert<QuadratorKonf> {
    static Node encode(const QuadratorKonf& qk) {
      Node node;
      node["version"] = qk.version;
      node["p_coi"]   = qk.p_coi;
      node["p_smi"]   = qk.p_smi;
      node["p_sma"]   = qk.p_sma;
      node["p_mod"]   = qk.p_mod;
      node["n_tuile"] = qk.n_tuile;
      std::array<Ctuile,DIM_TUILE> tmp;
      memcpy(tmp.data(),qk.ctuile,DIM_TUILE*sizeof(Ctuile));
      node["ctuile"]  = tmp;
      return node;
    }
    static bool decode(const Node& node, QuadratorKonf& qk) {
      if(!node.IsMap() )
	return false;
      qk.version = node["version"].as<int32_t>();
      qk.p_coi	 = node["p_coi"].as<int32_t>();
      qk.p_smi	 = node["p_smi"].as<int32_t>();
      qk.p_sma	 = node["p_sma"].as<int32_t>();
      qk.p_mod	 = node["p_mod"].as<int32_t>();
      qk.n_tuile = node["n_tuile"].as<int32_t>();
      memcpy(qk.ctuile,
	     node["ctuile"].as<std::array<Ctuile,DIM_TUILE>>().data(),
	     sizeof(Ctuile)*DIM_TUILE);
      return true;
    }
  };
}

#endif // KONFIGURATIONPET
