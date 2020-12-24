#ifndef KRISTAL_H
#define KRISTAL_H
#include <yaml-cpp/yaml.h>


typedef struct _kristal{
  int32_t version;
  int32_t matiere;   // Matiere  MAT_LYSO ou MAT_CEBR
  int32_t indice;    // Indice du cristal multiplie par LRS_MULT_INDICE pour entier
  int32_t energie;   // Energie de mesure du photopic (511000 eV) en eV
  int32_t photopik;  // Position du photo pic 511 keV en photons
  int32_t epaisseur; // Epaisseur du cristal en mm
}Kristal;


namespace YAML {
  template<>
    struct convert<Kristal> {
    static Node encode(const Kristal& krt) {
      Node node;
      node["version"]   = krt.version ;
      node["matiere"]   = krt.matiere ;
      node["indice"]    = krt.indice ;
      node["energie"]   = krt.energie ;
      node["photopik"]  = krt.photopik ;
      node["epaisseur"] = krt.epaisseur ;
      return node;
    }
    static bool decode(const Node& node, Kristal& krt) {
      if(!node.IsMap() )
	return false;
      krt.version = node["version"].as<int32_t>();
      krt.matiere = node["matiere"].as<int32_t>();
      krt.indice = node["indice"].as<int32_t>();
      krt.energie = node["energie"].as<int32_t>();
      krt.photopik = node["photopik"].as<int32_t>();
      krt.epaisseur = node["epaisseur"].as<int32_t>();
      return true;
    }
  };
}



#endif // KRISTAL_H
