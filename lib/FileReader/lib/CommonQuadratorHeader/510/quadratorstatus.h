#ifndef QUADRATORSTATUS_H
#define QUADRATORSTATUS_H

#include <QString>
#include <QTcpSocket>
#include <qabstractsocket.h>
#include <QDataStream>
#include <QCoreApplication>
#include <cstring>

static const QString s_infoMessage[16] = {
  QCoreApplication::translate("QuadratorStatus", "Connecté avec succès à %1"),
  QCoreApplication::translate("QuadratorStatus", "Mise en Alimentation de la caméra Ok"),
  QCoreApplication::translate("QuadratorStatus", "VBias effectué avec succès"),
  QCoreApplication::translate("QuadratorStatus", "Les paramètres de configuration ont été envoyés avec succès"),
  QCoreApplication::translate("QuadratorStatus", "La configuration à été effectuée avec succès"),
  QCoreApplication::translate("QuadratorStatus", "Pré acquisition Ok, début de l'acquisition"),
  QCoreApplication::translate("QuadratorStatus", "Acquisition effectués avec succès en %1,%2 secondes"),
  QCoreApplication::translate("QuadratorStatus", "Photo Prise avec succès"),
  QCoreApplication::translate("QuadratorStatus", "Désactivation effectuée avec succès"),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", "Client déconnecté"),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", "Distance calculée par le télémètre : ")
};
static const QString s_pendingMessage[11] = {
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", "D'autre actions sont en cours sur le serveur, désactivation en attente..."),
  QCoreApplication::translate("QuadratorStatus", ""),
  QCoreApplication::translate("QuadratorStatus", "")
};
static const QString s_errorMessage[11] = {
  QCoreApplication::translate("QuadratorStatus", "La connection n'a pas pu aboutir :\n"),
  QCoreApplication::translate("QuadratorStatus", "Erreur lors de l'alimentation de la caméra."),
  QCoreApplication::translate("QuadratorStatus", "Erreur lors du VBias."),
  QCoreApplication::translate("QuadratorStatus", "Les paramètres de configuration ne sont pas correct pour cette caméra, ou la configuration à été corrompue."),
  QCoreApplication::translate("QuadratorStatus", "Erreur lors de la configuration de la caméra."),
  QCoreApplication::translate("QuadratorStatus", "Erreur lors de la préparation de l'acquisition."),
  QCoreApplication::translate("QuadratorStatus", "Erreur pendant l'acquisition."),
  QCoreApplication::translate("QuadratorStatus", "La photo n'a pas pu être prise."),
  QCoreApplication::translate("QuadratorStatus", "Erreur lors de la désactivation de la caméra."),
  QCoreApplication::translate("QuadratorStatus", "Erreur de connection :\n"),
  QCoreApplication::translate("QuadratorStatus", "")
};
static const QString s_forbiddenMessage[11] = {
  QCoreApplication::translate("QuadratorStatus", "Cette action n'est pas autorisée pour le moment."),
  QCoreApplication::translate("QuadratorStatus", "Cette action n'est pas autorisée pour le moment."),
  QCoreApplication::translate("QuadratorStatus", "Cette action n'est pas autorisée pour le moment."),
  QCoreApplication::translate("QuadratorStatus", "Cette action n'est pas autorisée pour le moment."),
  QCoreApplication::translate("QuadratorStatus", "Cette action n'est pas autorisée pour le moment."),
  QCoreApplication::translate("QuadratorStatus", "Cette action n'est pas autorisée pour le moment."),
  QCoreApplication::translate("QuadratorStatus", "Cette action n'est pas autorisée pour le moment."),
  QCoreApplication::translate("QuadratorStatus", "Cette action n'est pas autorisée pour le moment."),
  QCoreApplication::translate("QuadratorStatus", "Cette action n'est pas autorisée pour le moment."),
  QCoreApplication::translate("QuadratorStatus", "Cette action n'est pas autorisée pour le moment."),
  QCoreApplication::translate("QuadratorStatus", "Cette action n'est pas autorisée pour le moment.")
};

static const QString s_connectionError[3] = {
  QCoreApplication::translate("QuadratorStatus", "Version du client trop ancienne."),
  QCoreApplication::translate("QuadratorStatus", "Version du serveur trop ancienne."),
  QCoreApplication::translate("QuadratorStatus", "Un client est déjà connecté.")
};

static const QString s_vektorError[2] = {
  QCoreApplication::translate("QuadratorStatus", "Vektor : Les bases de données n'ont pas été trouvées."),
  QCoreApplication::translate("QuadratorStatus", "Vektor : La version de Vektor n'est pas compatible avec les structures de données utilisée.")
};

class QuadratorStatus
{
public:

  /// \brief Enumération listant les différent type de status.
  enum StatusType{
    INFO,///< Le status est une information. Les informations sont envoyée à la fin de l'action concernée. Ces status ont surtout pour but d'afficher des informations à un utilisateur.
    PENDING,///< Le status informe qu'un action est en attente car d'autres actions sont en cours. Cette fonctionnalité est seulement implémentée pour la désactivation (QuadratorStatus::DESACTIVATION), les autres action emettent le stype suivant.
    FORBBIDEN,///< Permet d'informer que l'action demandée n'est pas possible actuellement, car d'autres action sont en cours.
    AUTHORIZED,///< Status informant la possibilité d'effectuer l'action concernée. Cette fonctionnalité est implémentée seulement pour la prise de photo (QuadratorStatus::PHOTO).
    ERROR///< Status informant qu'une erreur a été rencontrée, pendant l'action correspondante.
  };

  /// \brief Enumération listant les différent Status pouvant être envoyé.
  enum StatusName{
    /// \brief Information sur la connection.
    ///
    /// Plusieurs cas sont possible, selon le type :
    /// - QuadratorStatus::INFO : La connection s'est bien déroulée. On peut récupérer via quadratorStatus::connectionAddress l'adresse de connection
    /// - QuadratorStatus::ERROR : La connection n'est pas possible avec ce serveur. Le type d'erreur peut être récupéré via quadratorStatus::connectionError. Cette méthode sort un QuadratorStatus::ConnectionErrorType.
    CONNECTION = 0,

    /// \brief Information concernant le power UP.
    ///
    /// Plusieurs cas sont possible selon le type :
    /// - QuadratorStatus::INFO : Le power up s'est bien passé et est terminé.
    /// - QuadratorStatus::ERROR : Le power up n'a pas pu aboutir.
    POWER_UP = 1,

    /// \brief Information concernant le VBIAS
    ///
    /// Plusieurs cas sont possible selon le type :
    /// - QuadratorStatus::INFO : Le VBIAS s'est correctement déroulé
    /// - QuadratorStatus::ERROR : Le VBIAS n'a pas pu aboutire à cause d'une erreur.
    /// - QuadratorStatus::FORBIDDEN : Une action est déjà en cours, le VBIAS est donc annulé.
    VBIAS = 2,

    /// \brief Information concernant les paramètres de configuration
    ///
    /// Plusieurs cas sont possible selon le type :
    /// - QuadratorStatus::INFO : Les paramêtres de configuration sont bon.
    /// - QuadratorStatus::ERROR : les paramêtres de configuration ne sont pas bon, ou ont été corrompu, ou la structure n'est pas de la bonne version.
    CONFIGURATION = 3,

    /// \brief Information concernant la configuration
    ///
    /// Plusieurs cas sont possible selon le type :
    /// - QuadratorStatus::INFO : La configuration s'est correctement déroulée
    /// - QuadratorStatus::ERROR : La configuration n'a pas pu aboutire à cause d'une erreur.
    /// - QuadratorStatus::FORBIDDEN : Une action est déjà en cours, la configuration est donc annulée.
    CONFIG_SET = 4,

    /// \brief Information concernant la préparation de l'acquisition.
    ///
    /// Pour l'instant cette fonctionnalité n'est pas implémentée, mais permettra de connaitre le début réel de l'acquisition.
    PRE_ACQUISITION = 5,

    /// \brief Information concernant l'acquisition
    ///
    /// Plusieurs cas sont possible selon le type :
    /// - QuadratorStatus::INFO : L'acquisition est terminée sans erreur. On peut retrouver le temps total d'acquisition grâce à la méthode QuadratorStatus::time.
    /// - QuadratorStatus:: ERREUR : L'acquisition n'a pas pu aboutir. Il est possible tout de même de recevoir des événement avant une erreur. Dans ce cas la il n'est pas possible de récupérer le temps total d'acquisition.
    /// - QuadratorStatus::FORBIDDEN : Une action est déjà en cours ou le vbias et la configuration n'ont pas été faites, l'acquisition est donc annulée.
    ACQUISITION = 6,

    /// \brief Information concernant la prise de photo
    ///
    /// Plusieurs cas sont possible selon le type :
    /// - QuadratorStatus::AUTHORIZED : Cela indique que la prise de photos est possible avec cette caméra. Ce status est envoyé juste après la connection, donc s'il nest pas reçus, il faut partir du principe que la photo est impossible.
    /// - QuadratorStatus::ERROR : Une erreur est survenu lors de la prise de photo.
    PHOTO = 7,

    /// \brief Information concernant la desactivation
    ///
    /// Plusieurs cas sont possible selon le type :
    /// - QuadratorStatus::INFO : La désactivation s'est erminée sans erreur.
    /// - QuadratorStatus::ERROR : Une erreur est survenu lors de la désactivation.
    /// - QuadratorStatus::PENDING : La désactivation est en attente de la fin d'une action en cours.
    DESACTIVATION = 8,

    /// \brief Information concernant une erreur du socket TCP
    ///
    /// Seul le type QuadratorStatus::ERROR est associé à cette information. Cela signifie qu'il y a eu une erreur de connection. Il est possible de récupérer l'erreur via la méthode QuadratorStatus::connectionError.
    SOCKET_ERROR = 9,

    /// \brief Information concernant la déconnection de la caméra
    ///
    /// Seul le type QuadratorStatus::IFNO est envoyé avec ce status. Il permet d'informer l'utilisateur de la déconnection de la caméra.
    DISCONNECTION = 10,

    /// \brief Information concernant Vektor
    ///
    /// Seule le type QuadratorStatus::ERROR est implémenté avec ce status. Il permet d'e connaître d'informer l'interface d'une erreur concernant la vectorisation. On peut connaitre le type d'erreur via la méthode QuadratorStatus::vektorError. Cette méthode retourne un QuadratorStatus::VektorErrorType.
    VEKTOR = 11,

    /// \brief Permet de créer une information customisée selon le type que l'on veut.
    ///
    /// Il faut renseigner l'attribut message afin de pouvoir afficher une information avec QuadratorStatus::what comme les autre status. Ce status est utile pour des ui plus simple.
    CUSTOM = 12,

    UNKNOWN = 13,

    /// \brief Permet de créer une information sur le spectre CZT
    ///
    /// Seul le type QuadratorStatus::AUTHORIZED est implémenté. Permet de savoir que le spectre CZT est présent.
    CZT = 14,

    /// \brief Permet de créer un statut pour le télémètre
    ///
    /// Les type QuadratorStatus::AUTHORIZED et QuadratorStatus::INFO sont implémenté. Le status INFO permet de récupérer la distance de l'objet situé en face de l'acquisition.
    TELEMETER = 15
  };

  /// \brief Enumération listant les différent type d'erreur de connection (côté protocole serveur/client)
  enum ConnectionErrorType{
    VERSION_ERROR_CLIENT_TOO_LOW,///< La version du client est trop faible et donc icompatible avec la version du serveur de la caméra.
    VERSION_ERROR_SERVER_TOO_LOW,///< La version du serveur de la caméra est trop faible et incompatible avec le client.
    USER_ALREADY_CONNECTED,///< Un utilisateur est déjà connecté à la caméra.
    NOT_CONNECTED ///< L'utilisateur n'est pas connecté, impossible de faire l'action demandée
  };

  /// \brief Enumération listant les erreur possible avec Vektoe
  enum VektorErrorType{
    NO_FILE,///< Les bases de données ne sont pas trouvé à leur emplacement par default.
    WRONG_HEADER_VERSION///< La version de Vektor ne fonctionne pas avec les mêmes structure de données.
  };

  /// \brief Constructeur de la classe.
  ///
  /// \param t_name Le nom du status. Voir QuadratorStatus::StatusName.
  /// \param t_type Le type de status. Voir QuadratorStatus::StatusType.
  QuadratorStatus(QuadratorStatus::StatusName t_name = UNKNOWN, QuadratorStatus::StatusType t_type = QuadratorStatus::INFO): m_type(t_type), m_name(t_name)
  {
      //Création du message d'erreur en fonction d'un tableau de message (cas général)
  }

  /// \brief Méthode retournant une phrase d'explication du status. utile pour les interface utilisateur manuelle.
  QString what(){
    if(m_name == CUSTOM){
      return message;
    }else if(m_type == INFO){
      if(m_name == CONNECTION){
        return s_infoMessage[m_name].arg(m_connectionAddress);
      }else if(m_name == ACQUISITION){
        QString info = s_infoMessage[m_name];
        return info.arg(m_time[0]).arg(m_time[1]);
      }else if(m_name == TELEMETER){
        QString info = s_infoMessage[m_name];
        return info.arg(m_telemeterDistance);
      }
      else{
        return s_infoMessage[m_name];
      }
    }else if(m_type == PENDING){
      return s_pendingMessage[m_name];
    }else if(m_type == ERROR){
      if(m_name == SOCKET_ERROR){
        QString error = s_errorMessage[m_name];
        return error.append(message);
      }else if(m_name == CONNECTION){
        QString error = s_errorMessage[m_name];
        return error.append(s_connectionError[m_connectionError]);
      }else if(m_name == VEKTOR){
        return s_vektorError[m_vektorError];
      }else{
        return s_errorMessage[m_name];
      }
    }else{
      return s_forbiddenMessage[m_name];
    }
  }

  /// \brief un message utillisé pour un status custom ou pour garder en mémoire la phrase d'erreur de socket.
  QString message; //TODO: passer en private avec la liste des messages dans un tableau dans le C

  /// \brief Getter de l'attribut m_disconnectionError, qui contient l'erreur de Socket.
  QAbstractSocket::SocketError disconnectionError() const{ return m_disconnectionError; }
  /// \brief Setter de l'attribut m_disconnectionError, qui contient l'erreur de Socket.
  void setDisconnectionError(const QAbstractSocket::SocketError err) { m_disconnectionError = err; }

  /// \brief Getter de l'attribut m_connectionError, qui contient l'erreur de connection.
  ConnectionErrorType connectionError() const{ return m_connectionError; }
  /// \brief Setter de l'attribut m_connectionError, qui contient l'erreur de connection.
  void setConnectionError(const ConnectionErrorType &connectionError){ m_connectionError = connectionError; }

  int telemeterDistance() const{ return m_telemeterDistance; }

  void setTelemeterDistance(const int &telemeterDistance){ m_telemeterDistance = telemeterDistance; }

  /// \brief Getter de l'attribut m_connectionAddress, qui contient l'adresse de la caméra connectée.
  QString connectionAddress() const { return m_connectionAddress; }
  /// \brief Setter de l'attribut m_connectionAddress, qui contient l'adresse de la caméra connectée.
  void setConnectionAddress(const QString &connectionAddress){ m_connectionAddress = connectionAddress; }

  /// \brief Getter de l'attribut m_time, qui contient le temps total d'acquisition. C'est un tableau à deux dimension contenant le nombre de seconde dans le premier élément et le nombre de milliseconde dans le second.
  const int *time() const{ return m_time; }
  /// \brief Setter de l'attribut m_time, qui contient le temps total d'acquisition. C'est un tableau à deux dimension contenant le nombre de seconde dans le premier élément et le nombre de milliseconde dans le second.
  void setTime(const int sec, const int usec){ m_time[0] = sec; m_time[1] = usec; }

  /// \brief Getter de l'attribut m_name, qui contient le nom du status. Cet attribut est initialisé dans le constructeur.
  StatusName name() const{ return m_name; }
  /// \brief Setter de l'attribut m_name, qui contient le nom du status. Cet attribut est initialisé dans le constructeur.
  void setName(const StatusName &name){ m_name = name; }

  /// \brief Getter de l'attribut m_type, qui contient le type du status. Cet attribut est initialisé dans le constructeur.
  StatusType type() const{ return m_type; }
  /// \brief Setter de l'attribut m_type, qui contient le type du status. Cet attribut est initialisé dans le constructeur.
  void setType(const StatusType &type){ m_type = type; }


  /// \brief Getter de l'attribut m_vektorError, qui contient le l'erreur de Vektor.
  VektorErrorType vektorError() const{ return m_vektorError; }
  /// \brief Setter de l'attribut m_vektorError, qui contient le l'erreur de Vektor.
  void setVektorError(const VektorErrorType &vektorError){ m_vektorError = vektorError; }

private:
  /// \brief contient le nom du status. Cet attribut est initialisé dans le constructeur.
  StatusType m_type;
  /// \brief contient le type du status. Cet attribut est initialisé dans le constructeur.
  StatusName m_name;

  union{
    QAbstractSocket::SocketError m_disconnectionError;///< Situationnel, sert en cas d'erreur de SOCKET
    ConnectionErrorType m_connectionError;///< Situationnel, seulement lors d'une erreur de connection type erreur de version ou client déjà connecté.
    VektorErrorType m_vektorError;///< Situationnel, seulement lors d'une erreur de Vektor.
    int m_time[2];///< Situationnel, permet de connaitre le temps total d'acquisition.
    int m_telemeterDistance;
  };
  QString m_connectionAddress;///< Situationnel, permet de connaitre l'adresse de connection.

};

inline QDataStream &operator<<(QDataStream &out, const QuadratorStatus &myObj)
{
  out << myObj.name();
  out << myObj.type();
  out << myObj.message;
  const int* time = myObj.time();
  out << time[0] << time[1];
  out << myObj.connectionAddress();
  return out;
}

inline QDataStream &operator>>(QDataStream &in, QuadratorStatus &myObj)
{
  int param;
  in >> param;
  myObj.setName(static_cast<QuadratorStatus::StatusName>(param));
  in >> param;
  myObj.setType(static_cast<QuadratorStatus::StatusType>(param));
  in >> myObj.message;
  int time[2];
  in >> time[0] >> time[1];
  myObj.setTime(time[0], time[1]);
  QString connectionAddress;
  in >> connectionAddress;
  myObj.setConnectionAddress(connectionAddress);
  return in;
}

#endif // QUADRATORSTATUS_H
