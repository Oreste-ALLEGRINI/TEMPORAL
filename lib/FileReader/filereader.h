#ifndef FILEREADER_H
#define FILEREADER_H
/// \file filereader.h
/// \brief fichier de la classe FileReader, qui contient toute les primitive de lecture des fichiers, ainsi que les fonction de manipulation des fichiers.
///
/// Afin d'utiliser cette librairie, seul ce header est à inclure. Cette libriaire a absolument besoin de la librairie common-quadrator-header.

#include "filereader_global.h"
#include <QString>
#include <quadratorevt.h>

/// \brief FileReader est la classe de manipulation des fichiers type damavan. 
///
/// Une fois la classe instanciée, elle peut ouvrir deux fichier maximum, un en écriture et un en lecture. Les fichiers en question peuvent être de 5 types différent :
/// - *.lqua -> Correspond au événements brut sans coïncidences temporelles. Ce type de fichier est l'écriture binaire d'une structure TuileDam. Il contient donc les informations relatives à l'événement ainsi que les cartes de photons et de temps.
/// - *.cqua -> Correspond au événements appairé en coïncidence. Ce type de fichier est l'écriture binaire d'une structure PaireDam. Il contient donc les informations relatives aux événements en coïncidence ainsi que les cartes de photons et de temps.
/// - *.lql -> Correspond au événements brut sans coïncidences temporelles, allégé. Ce type de fichier est l'écriture binaire d'une structure TuileDam_l. Il contient donc les informations relatives à l'événement sans les cartes de photons et de temps.
/// - *.cql -> Correspond au événements appairé en coïncidence, allégé. Ce type de fichier est l'écriture binaire d'une structure PaireDam_L. Il contient donc les informations relatives aux événements en coïncidence sans les cartes de photons et de temps.
/// - *.json -> Correspond, comme les fichiers *cql, au événements appairé en coïncidence, allégé. Ce type de fichier est l'écriture json d'une structure PaireDam_L. Il contient donc les informations relatives aux événements en coïncidence sans les cartes de photons et de temps.
///
/// l'instance de cette classe peut ouvrir deux fichiers maximums, un en écriture, et un en lecture. L'ouverture des fichiers peut être faite en mode automatique, c'est à dire que le type de fichiers est automatiquement trouvé grâce au suffixe de celui-ci.
/// Selon le type de fichier, certaine lecture et ecriture des structure damavan peut être interdite, comme par exemple la lecture d'une structure PaireDam dans un fichier *.cql, en effet le *cql ne contient pas les information necessaires (pas de carte de photons/temps).
///
/// Les méthodes de cette clase retourne des FileReaderResult, qui permettent de comprndre le type d'erreur rencontré lors de l'ouverture/écriture/lecture d'un fichier.
///
/// Voici un exemple d'utilisation de la classe :
///
/// \code{.cpp}
/// FileReader freader;
/// freader.openReadFile("monfichier", FileReader::CQUA);//Ouverture d'un fichier cqua en lecture
/// freader.openWriteFile("monfichier2", FileReader::LQUA);//Ouverture d'un fichier lqua en écriture
/// PaireDam dam;
/// freader.getNextEvent(&dam); //Récupération d'un PaireDam
/// freader.writeNextEvent(&dam); //Ecriture du contenu du PaireDam (écrit entre 1 à 4 TuileDam)
/// freader.closeFiles();
/// \endcode
class FileReader
{

public:
    /// \brief Enumération présentant les différents type de fichiers que l'on peut avoir, ainsi que deux type alternatif utilisé pour trouver le type d'un fichier de façon automatique.
    FILEREADERSHARED_EXPORT enum TypeFichier{
        CQUA,///< Type de fichier d'un *.cqua
        CQL,///< Type de fichier d'un *.cql
        LQUA,///< Type de fichier d'un *.lqua
        LQL,///< Type de fichier d'un *.lql
        JSON,///< Type de fichier d'un *.json
        CSV,///< Type de fichier d'un *.csv
        AUTO,///< Permet le recherche automatique du type de fichier à l'ouverture de celui ci
        FROMREADFILE///< Permet d'ouvrir un fichier en écriture avec le même type que celui en lecture
    };

    /// \brief Enumération représentant tout les retours des méthode de la classe. mis à part OK et END_OF_FILE, tous sont des codes d'erreurs.
    FILEREADERSHARED_EXPORT enum FileReaderResult{
        OK =                         1,///< la méthode à fonctionné sans erreur
        END_OF_FILE =                0,///< La fin du fichier en lecture est atteinte
        OPENNING_ERROR =            -1,///< le fichier en lecture ou en écriture n'a pas pu être ouvert
        FILE_ALREADY_OPENNED =      -2,///< Un fichier en lecture ou ecriture est déjà ouvert. Il faut donc utiliser closeFiles() afin de fermer les fichiers précédemment ouvert
        READ_ERROR =                -3,///< Une erreur est survenue lors de la lecture d'un fichier
        WRITE_ERROR =               -4,///< Une erreur est survenue lors de l'écriture d'un fichier
        CONVERSION_NOT_ALLOWED =    -5,///< La méthode ne permet pas de lire ou d'ecrire un type de structure avec le type de fichier ouvert. Pour passer outre, il faut songer a manuellement faire la conversion vers un type de fichier autorisé
        AUTO_DETECTION_FAILED =     -6,///< L'auto détection du type de fichier n'a pas pu fonctionner (erreur ou absence de suffixe dans le nom du fichier
        FILE_CORRUPTED =            -7,///< Le fichier est corrompu
        READ_FILE_NOT_OPENNED =     -8,///< Le fichier en lecture n'est pas ouvert. Ce code d'erreur apparait lors d'une tentative d'ouverture de fichier avec le flag FromReadFile quand le fichier en lecture n'est pas ouvert
        FILE_NOT_OPENNED =          -9,///< Le ficheir de lecture ou d'écriture n'est pas ouvert
        INVALID_FILE_TYPE =         -10///< Le type de fichier est invalide.
    };

    /// \brief constructeur de la classe
    FILEREADERSHARED_EXPORT FileReader();
    /// \brief Destructeur de la classe
    FILEREADERSHARED_EXPORT ~FileReader();
    /// \brief Méthode permettant l'ouverture d'un fichier en lecture
    ///
    /// @param entrant le nom du fichier à ouvrir
    /// @param type le type de fichier à ouvrir. Peut être soit un type de fichier, soit FileReader::AUTO. par default le flag est FileReader:AUTO.
    ///
    /// @return
    /// - FileReader::OK si le fichier est ouvert
    /// - FileReader::OPENNING_ERROR si le fichier n'existe pas, ou que le fichier ne peut pas être ouvert
    /// - FileReader::INVALID_FILE_TYPE si le flag FROMREADFILE est utilisé
    /// - FileReader::FILE_ALREADY_OPENNED si un fichier en lecture est déjà ouvert
    /// - FileReader::AUTO_DETECTION_FAILED si le flag FileReader::AUTO etait renseigné et que le suffixe ne correspond pas a un type de fichier connu
    /// - FileReader::FILE_CORRUPTED si le fichier n'est pas de la bonne taille par rapport à un fichier du type trouvé/renseigné. Cette protection ne fonctionne pas pour les fichier json.
    FILEREADERSHARED_EXPORT FileReaderResult openReadFile(QString entrant, TypeFichier type = FileReader::AUTO);

    /// \brief Méthode permettant l'ouverture d'un fichier en écriture
    ///
    /// @param sortant le nom du fichier à ouvrir
    /// @param type le type de fichier à ouvrir. Peut être soit un type de fichier, soit FileReader::AUTO, soit FileReader::FROMREADFILE. par default le flag est FileReader:AUTO.
    ///
    /// @return
    /// - FileReader::OK si le fichier est ouvert
    /// - FileReader::OPENNING_ERROR si le fichier n'existe pas, ou que le fichier ne peut pas être ouvert
    /// - FileReader::FILE_ALREADY_OPENNED si un fichier en écriture est déjà ouvert
    /// - FileReader::READ_FILE_NOT_OPENNED si le flag FROMREADFILE est utilisé et que le fichier en lecture n'a pas été préalablement ouvert
    /// - FileReader::AUTO_DETECTION_FAILED si le flag FileReader::AUTO etait renseigné et que le suffixe ne correspond pas a un type de fichier connu
    FILEREADERSHARED_EXPORT FileReaderResult openWriteFile(QString sortant, TypeFichier type = FileReader::AUTO);

    /// \brief Méthode permettant la récupération d'un PaireDam
    ///
    /// @param event le pointeur sur l'événement à lire.
    ///
    /// @return
    /// - FileReader::OK si un événement à été correctement lu est placé dans event
    /// - FileReader::FILE_NOT_OPENNED si le fichier de lecture n'a pas été ouvert
    /// - FileReader::CONVERSION_NOT_ALLOWED si le type de fichier est un *.cql/*.lql/*.json
    /// - FileReader::READ_ERROR si'il y a eu un problème de lecture dans le fichier
    /// - FileReader::END_OF_FILE lorsque la fin du fichier à été atteinte. Notez qu'event n'est pas correct dans ces cas là.
    /// - FileReader::FILE_CORRUPTED si le numéro de séquence n'est pas correct.
    FILEREADERSHARED_EXPORT FileReaderResult getNextEvent(PaireDam *event, bool forceConversion = false);

    /// \brief Méthode permettant la récupération d'un PaireDam_l
    ///
    /// @param event le pointeur sur l'événement à lire.
    ///
    /// @return
    /// - FileReader::OK si un événement à été correctement lu est placé dans event
    /// - FileReader::FILE_NOT_OPENNED si le fichier de lecture n'a pas été ouvert
    /// - FileReader::READ_ERROR si'il y a eu un problème de lecture dans le fichier
    /// - FileReader::END_OF_FILE lorsque la fin du fichier à été atteinte. Notez quevent n'est pas correct dans ces cas là.
    /// - FileReader::FILE_CORRUPTED si le numéro de séquence n'est pas correct.
    FILEREADERSHARED_EXPORT FileReaderResult getNextEvent(PaireDam_l *event);

    /// \brief Méthode permettant la récupération d'un TuileDam
    ///
    /// @param event le pointeur sur l'événement à lire.
    ///
    /// @return
    /// - FileReader::OK si un événement à été correctement lu est placé dans event
    /// - FileReader::FILE_NOT_OPENNED si le fichier de lecture n'a pas été ouvert
    /// - FileReader::CONVERSION_NOT_ALLOWED si le type de fichier est un *.cql/*.lql/*.json
    /// - FileReader::READ_ERROR si'il y a eu un problème de lecture dans le fichier
    /// - FileReader::END_OF_FILE lorsque la fin du fichier à été atteinte. Notez qu'event n'est pas correct dans ces cas là.
    /// - FileReader::FILE_CORRUPTED si le numéro de séquence n'est pas correct.
    FILEREADERSHARED_EXPORT FileReaderResult getNextEvent(TuileDam *event);

    /// \brief Méthode permettant la récupération d'un TuileDam_l
    ///
    /// @param event le pointeur sur l'événement à lire.
    ///
    /// @return
    /// - FileReader::OK si un événement à été correctement lu est placé dans event
    /// - FileReader::FILE_NOT_OPENNED si le fichier de lecture n'a pas été ouvert
    /// - FileReader::READ_ERROR si'il y a eu un problème de lecture dans le fichier
    /// - FileReader::END_OF_FILE lorsque la fin du fichier à été atteinte. Notez qu'event n'est pas correct dans ces cas là.
    /// - FileReader::FILE_CORRUPTED si le numéro de séquence n'est pas correct.
    FILEREADERSHARED_EXPORT FileReaderResult getNextEvent(TuileDam_l *event);

    /// \brief Méthode permettant l'écriture d'un PaireDam
    ///
    /// @param event le pointeur sur l'événement à écrire.
    ///
    /// @return
    /// - FileReader::OK si un événement à été correctement écrit dans le fichier de sortie
    /// - FileReader::FILE_NOT_OPENNED si le fichier de écriture n'a pas été ouvert
    /// - FileReader::WRITE_ERROR si'il y a eu un problème d'écriture dans le fichier
    FILEREADERSHARED_EXPORT FileReaderResult writeNextEvent(PaireDam *toWrite);

    /// \brief Méthode permettant l'écriture d'un PaireDam_l
    ///
    /// @param event le pointeur sur l'événement à écrire.
    ///
    /// @return
    /// - FileReader::OK si un événement à été correctement écrit dans le fichier de sortie
    /// - FileReader::FILE_NOT_OPENNED si le fichier de écriture n'a pas été ouvert
    /// - FileReader::CONVERSION_NOT_ALLOWED si le type de fichier est un *.cqua/*.lqua
    /// - FileReader::WRITE_ERROR si'il y a eu un problème d'écriture dans le fichier
    FILEREADERSHARED_EXPORT FileReaderResult writeNextEvent(PaireDam_l *toWrite);

    /// \brief Méthode permettant l'écriture d'un TuileDam
    ///
    /// @param event le pointeur sur l'événement à écrire.
    ///
    /// @return
    /// - FileReader::OK si un événement à été correctement écrit dans le fichier de sortie
    /// - FileReader::FILE_NOT_OPENNED si le fichier de écriture n'a pas été ouvert
    /// - FileReader::CONVERSION_NOT_ALLOWED si le type de fichier est un *.cqua/*.json
    /// - FileReader::WRITE_ERROR si'il y a eu un problème d'écriture dans le fichier
    FILEREADERSHARED_EXPORT FileReaderResult writeNextEvent(TuileDam *toWrite);

    /// \brief Méthode permettant l'écriture d'un TuileDam_l
    ///
    /// @param event le pointeur sur l'événement à écrire.
    ///
    /// @return
    /// - FileReader::OK si un événement à été correctement écrit dans le fichier de sortie
    /// - FileReader::FILE_NOT_OPENNED si le fichier de écriture n'a pas été ouvert
    /// - FileReader::CONVERSION_NOT_ALLOWED si le type de fichier est un *.cqua/*.cql/*.lqua/*.json
    /// - FileReader::WRITE_ERROR si'il y a eu un problème d'écriture dans le fichier
    FILEREADERSHARED_EXPORT FileReaderResult writeNextEvent(TuileDam_l *toWrite);

    ///\brief Méthode permettant de connaitre l'état d'ouverture du fichier de lecture
    ///
    /// @return true si le fichier est ouvert, false sinon
    FILEREADERSHARED_EXPORT bool isReadFileOpen();

    ///\brief Méthode permettant de connaitre l'état d'ouverture du fichier de lecture
    ///
    /// @return true si le fichier est ouvert, false sinon
    FILEREADERSHARED_EXPORT bool isWriteFileOpen();

    ///\brief Méthode permettant la fermeture des fichier en lecture et écriture
    FILEREADERSHARED_EXPORT void closeFiles();

    /// \brief retourne de type du fichier en lecture
    FILEREADERSHARED_EXPORT TypeFichier getType() const;



private:
    /// \internal
    /// \brief pointeur sur le fichier entrant
    FILE *entrant = NULL;

    /// \internal
    /// \brief pointeur sur le fichier sortant
    FILE *sortant = NULL;

    ///\internal
    /// \brief type du fichier entrant
    TypeFichier typeEntrant;

    /// \internal
    /// \brief type du fichier sortant
    TypeFichier typeSortant;

    ///\internal
    /// \brief permet de tester le numéro de séquence des événement du fichier entrant. Si celui-ci n'est pas séquentiel, on jette une erreur, car cel veut dire que le fichier est corrompu
    int numseqEntrant =0;

    ///\internal
    /// \brief peremt d'écrire le fichier avec un numéro de séquence. Permet d'e s'assurer la lecture du fichier 'être sur d'avoir un numéro sequentiel.
    int numseqSortant = 0;

    /// \internal
    /// \brief PaireDam intermédiaire, permet de pouvoir garder en mémoire un PaireDam lorsque l'on lit un fichier cqua/cql/json et qu'on veut une sortie en TuileDam ou TuileDam_l
    PaireDam intermediaire;

    /// \internal
    /// \brief indice de tuile du PaireDam intermédiaire. Permet, avec une boucle for, de tester les quatre numéro de tuile afin de pouvoir donner à l'utilisateur tout les TuileDam d'un PaireDam
    int indiceIntermediaire = 0;

    /// \internal
    /// \brief Méthode Permettant la lecture d'un fichier json
    ///
    /// @param event un événement à lire qui peut être du type PaireDam ou PaireDam_l
    ///
    /// @return
    /// - FileReader::OK si le json à bien été lu
    /// - FileReader::FILE_CORRUPTED Si le fichier n'est pas en json, s'il n'a pas tout les champs demandé ou si le numéro de séquence n'est pas correct.
    /// - FileReader::READ_ERROR Si le fichier n'a pas u être lu correctement.
    /// - FileReader::END_OF_FILE Si la fin du fichier à été atteint.
    template<class ToRead>
    FileReaderResult writeNextEventToJSON(ToRead *event);

    /// \internal
    /// \brief Méthode Permettant l'écriture d'un fichier json
    ///
    /// @param event un événement à écrire qui peut être du type PaireDam ou PaireDam_l
    ///
    /// @return
    /// - FileReader::OK si le fichier à bien été lu
    /// - FileReader::WRITE_ERROR Si le fichier n'a pas pu être lu correctement.
    template<class PaireEvent>
    FileReaderResult getNextEventFromJSON(PaireEvent *event);

    /// \internal
    /// \brief Méthode Permettant l'écriture d'un fichier csv
    ///
    /// @param event un événement à écrire qui peut être du type PaireDam ou PaireDam_l
    ///
    /// @return
    /// - FileReader::OK si le fichier à bien été lu
    /// - FileReader::WRITE_ERROR Si le fichier n'a pas pu être lu correctement.
    template<class PaireEvent>
    FileReaderResult writeNextEventToCsv(PaireEvent *event);

    /// \internal
    /// \brief Méthode Permettant lla lecture d'un fichier csv
    ///
    /// @param event un événement à lire qui peut être du type PaireDam ou PaireDam_l
    ///
    /// @return
    /// - FileReader::OK si le csv à bien été lu
    /// - FileReader::FILE_CORRUPTED Si le fichier n'est pas en csv, s'il n'a pas tout les champs demandé ou si le numéro de séquence n'est pas correct.
    /// - FileReader::READ_ERROR Si le fichier n'a pas u être lu correctement.
    /// - FileReader::END_OF_FILE Si la fin du fichier à été atteint.
    template<class PaireEvent>
    FileReaderResult getNextEventFromCsv(PaireEvent *event);

    /// \internal
    /// \brief Méthode Permettant la lecture d'un fichier binaire
    ///
    /// @param event un événement à lire qui peut être du type PaireDam, PaireDam_l, TuileDam ou TuileDam_l
    ///
    /// @return
    /// - FileReader::OK si le fichier à bien été lu
    /// - FileReader::FILE_CORRUPTED Si le numéro de séquence n'est pas correct.
    /// - FileReader::READ_ERROR Si le fichier n'a pas pu être lu correctement.
    /// - FileReader::END_OF_FILE Si la fin du fichier à été atteint.
    template<class PaireEvent>
    FileReaderResult getNextEventFromFile(PaireEvent *event);

    /// \internal
    /// \brief Méthode Permettant la l'ecriture d'un fichier binaire
    ///
    /// @param event un événement à écrire qui peut être du type PaireDam, PaireDam_l, TuileDam ou TuileDam_l
    ///
    /// @return
    /// - FileReader::OK si le fichier à bien été lu
    /// - FileReader::WRITE_ERROR Si le fichier n'a pas pu être lu correctement.
    template<class PaireEvent>
    FileReaderResult writeNextEventToFile(PaireEvent *event);
};


#endif // FILEREADER_H
