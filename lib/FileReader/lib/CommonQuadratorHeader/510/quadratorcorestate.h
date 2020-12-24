#ifndef QUADRATORCORESTATE_H
#define QUADRATORCORESTATE_H

namespace QuadratorCoreState{
  /// \brief Enumération représentant les différents état de la caméra. Ces états peuvent être combinés pour certain (ex : QuadratorCore::CONNECTED et QuadratorCore::POWER_UP.
  enum CoreState{
    NOT_CONNECTED      = 0x00000001,///< Le core n'est pas connecté à une caméra. C'est l'état de départ.
    ON_CONNECTION      = 0x00000002,///< Le core est en train de se connecter à une caméra.
    CONNECTED          = 0x00000004,///< Le core est connecté. Ce flag est laissé à 1 tant que le core est connecté
    POWERED_UP         = 0x00000008,///< Le TEK est bien alimenté. Ce flag est persistant.
    ON_CONFIGURATION   = 0x00000010,///< Une configuration est en train d'être envoyée à la caméra
    CONFIGURED         = 0x00000020,///< La caméra est configurée, Ce flag est persistant.
    ON_VBIAS           = 0x00000040,///< Le VBIAS est en cours
    VBIAS_DONE         = 0x00000080,///< Le VBIAS est effectué. Ce flag est persistant.
    ON_ACQUISITION     = 0x00000100,///< La caméra est en acquisition.
    ON_DESACTIVATION   = 0x00000200,///< La caméra est en train d'être désactivée.
    READY_TO_SHUTDOWN  = 0x00000400 ///< La caméra est prête à être éteinte. Peut être passé dans un status plutôt (à voir).
  };
}
#endif /* end of include guard QUADRATORCORESTATE_H */

