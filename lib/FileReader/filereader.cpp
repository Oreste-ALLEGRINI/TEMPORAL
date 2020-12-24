#include "filereader.h"
#include "quadratorevt.h"
#include <QFileInfo>
FileReader::FileReader()
{

}

FileReader::~FileReader()
{
    closeFiles();
}

FileReader::FileReaderResult FileReader::openReadFile(QString entrant, FileReader::TypeFichier type)
{
    //Trouver le type de fichier
    if(this->entrant == NULL){
        typeEntrant = type;
        switch(typeEntrant){
        case FROMREADFILE:
            return INVALID_FILE_TYPE;
        case AUTO:
            if( entrant.endsWith(".cqua"))
                typeEntrant = CQUA;
            else if( entrant.endsWith(".lqua"))
                typeEntrant = LQUA;
            else if( entrant.endsWith(".json"))
                typeEntrant = JSON;
            else if( entrant.endsWith( ".lql"))
                typeEntrant = LQL;
            else if( entrant.endsWith(".cql"))
                typeEntrant = CQL;
            else if( entrant.endsWith(".csv"))
                typeEntrant = CSV;
            else return FileReader::AUTO_DETECTION_FAILED;
            break;
        default:
            break;
        }
        this->entrant = fopen(entrant.toStdString().data(), "rb");
        if(this->entrant == NULL){
            return OPENNING_ERROR;
        }
        if(typeEntrant != JSON && typeEntrant != CSV){
            fseek(this->entrant, 0L, SEEK_END);
            int sz = ftell(this->entrant);
            rewind(this->entrant);
            switch(typeEntrant){
            case CQUA:
                if(sz ==0 || sz%sizeof(PaireDam) != 0) return FILE_CORRUPTED;
                break;
            case LQUA:
                if(sz ==0 || sz%sizeof(TuileDam) != 0) return FILE_CORRUPTED;
                break;
            case CQL:
                if(sz ==0 || sz%sizeof(PaireDam_l) != 0) return FILE_CORRUPTED;
                break;
            case LQL:
                if(sz ==0 || sz%sizeof(TuileDam_l) != 0) return FILE_CORRUPTED;
                break;
            }
        }
    }else return FileReader::FILE_ALREADY_OPENNED;

    return OK;
}


FileReader::FileReaderResult FileReader::openWriteFile(QString sortant, TypeFichier type)
{
    //Trouver le type de fichier
    if(this->sortant == NULL){
        typeSortant = type;
        if(type == FROMREADFILE){
            if(this->entrant != NULL){
                typeSortant = typeEntrant;
            }else return READ_FILE_NOT_OPENNED; //PAS DE BREAK POUR POUVOIR METTRE LE BON SUFFIXE
        }


        QString suffix;
        switch(typeSortant){
        case AUTO:
            if( sortant.endsWith(".cqua"))
                typeSortant = CQUA;
            else if( sortant.endsWith(".lqua"))
                typeSortant = LQUA;
            else if( sortant.endsWith(".json"))
                typeSortant = JSON;
            else if( sortant.endsWith( ".lql"))
                typeSortant = LQL;
            else if( sortant.endsWith(".cql"))
                typeSortant = CQL;
            else if( sortant.endsWith(".csv"))
                typeSortant = CSV;
            else return FileReader::AUTO_DETECTION_FAILED;
            break;
        case CQUA:
            suffix = ".cqua";
            break;
        case LQUA:
            suffix = ".lqua";
            break;
        case CQL:
            suffix = ".cql";
            break;
        case CSV:
            suffix = ".csv";
            break;
        case LQL:
            suffix = ".lql";
            break;
        case JSON:
            suffix = ".json";
            break;
        default:
            break;
        }
        if(!suffix.isEmpty()){
            QFileInfo sortantfile(sortant);
            sortant = sortantfile.absolutePath() + "/" + sortantfile.baseName() + suffix;
        }


        this->sortant = fopen(sortant.toStdString().data(), "wb");
        if(this->sortant == NULL){
            return OPENNING_ERROR;
        }
    }else return FileReader::FILE_ALREADY_OPENNED;

    return OK;

}


FileReader::FileReaderResult FileReader::getNextEvent(PaireDam* event, bool forceConversion){
    if(entrant == NULL) return FILE_NOT_OPENNED;

    if(typeEntrant == CQUA) return getNextEventFromFile<PaireDam>(event);

    else if(typeEntrant == LQUA){
        TuileDam dam;
        FileReaderResult res = getNextEventFromFile<TuileDam>(&dam);
        if(res == OK){
            event->tuiles[0].tuile = dam.tuile.tuile;
            tuileDam2PaireDam(&dam, event, dam.tuile.tuile);
            return res;
        }
        return res;
    }

    else if(forceConversion && typeEntrant == CQL){
      PaireDam_l dam;
      FileReaderResult res = getNextEventFromFile<PaireDam_l>(&dam);
      if(res == OK){
        pDam_l2pDam(&dam, event);
      }
      return res;
    }

    else if(forceConversion && typeEntrant == JSON){
      return getNextEventFromJSON(event);
    }

    else if(forceConversion && typeEntrant == LQL){
      TuileDam_l dam;
      FileReaderResult res = getNextEventFromFile<TuileDam_l>(&dam);
      event->tuiles[0].tuile = dam.tuile.tuile;
      TuileDam damdam;
      tDam_l2tDam(&dam, &damdam);
      tuileDam2PaireDam(&damdam, event, dam.tuile.tuile);
      return res;
    }

    else if(forceConversion && typeEntrant == CSV){
      return getNextEventFromCsv(event);
    }
    return CONVERSION_NOT_ALLOWED;
}



FileReader::FileReaderResult FileReader::getNextEvent(PaireDam_l* event){
    if(entrant == NULL) return FILE_NOT_OPENNED;

    if(typeEntrant == CQUA){
        PaireDam dam;
        FileReaderResult res = getNextEventFromFile<PaireDam>(&dam);
        if(res == OK){
            pDam2pDam_l(&dam, event);
            return OK;
        }else return res;
    }
    else if(typeEntrant == CQL) return getNextEventFromFile<PaireDam_l>(event);

    else if(typeEntrant == LQUA){
        TuileDam dam;
        FileReaderResult res = getNextEventFromFile<TuileDam>(&dam);
        if(res == OK){
            PaireDam paire;
            paire.tuiles[0].tuile = dam.tuile.tuile;
            tuileDam2PaireDam(&dam, &paire, dam.tuile.tuile);
            pDam2pDam_l(&paire, event);
            return OK;
        }else return res;
    }

    else if(typeEntrant == LQL){
        TuileDam_l dam;
        FileReaderResult res = getNextEventFromFile<TuileDam_l>(&dam);
        if(res == OK){
            event->tuiles[0].tuile = dam.tuile.tuile;
            tuileDam_l2PaireDam_l(&dam, event, dam.tuile.tuile);
            return OK;
        }else return res;
    }

    else if(typeEntrant == CSV){
      return getNextEventFromCsv(event);
    }

    else if(typeEntrant == JSON)return getNextEventFromJSON(event);
    return CONVERSION_NOT_ALLOWED;
}



FileReader::FileReaderResult FileReader::getNextEvent(TuileDam *event){
    if(entrant == NULL) return FILE_NOT_OPENNED;

    if(typeEntrant == CQUA){
        for(int i = 0; i < 2; i ++){
            for(; indiceIntermediaire > 0 && indiceIntermediaire <= 4; indiceIntermediaire++){
                if(paireDam2TuileDam(&intermediaire, event, indiceIntermediaire)){
                    indiceIntermediaire ++;
                    return OK;
                }
            }
            FileReaderResult res = getNextEventFromFile<PaireDam>(&intermediaire);
            if(res != OK ) return res;
            indiceIntermediaire = 1;
        }return FILE_CORRUPTED;
    }

    else if(typeEntrant == LQUA)return getNextEventFromFile<TuileDam>(event);
    return CONVERSION_NOT_ALLOWED;
}


FileReader::FileReaderResult FileReader::getNextEvent(TuileDam_l *event){
    if(entrant == NULL) return FILE_NOT_OPENNED;

    if(typeEntrant == CQUA){
        for(int i = 0; i < 2; i ++){
            for(; indiceIntermediaire > 0 && indiceIntermediaire <= 4; indiceIntermediaire++){
                TuileDam dam;
                if(paireDam2TuileDam(&intermediaire, &dam, indiceIntermediaire)){
                    tDam2tDam_l(&dam, event);
                    indiceIntermediaire ++;
                    return OK;
                }
            }
            FileReaderResult res = getNextEventFromFile<PaireDam>(&intermediaire);
            if(res != OK) return res;
            indiceIntermediaire = 1;
        }return FILE_CORRUPTED;
    }

    else if(typeEntrant == CQL){
        for(int i = 0; i < 2; i ++){
            for(; indiceIntermediaire > 0 && indiceIntermediaire <= 4; indiceIntermediaire++){
                TuileDam dam;
                if(paireDam2TuileDam(&intermediaire, &dam, indiceIntermediaire)){
                    tDam2tDam_l(&dam, event);
                    indiceIntermediaire ++;
                    return OK;
                }
            }
            PaireDam_l dam;
            FileReaderResult res = getNextEventFromFile<PaireDam_l>(&dam);
            if(res != OK) return res;
            pDam_l2pDam(&dam, &intermediaire);
            indiceIntermediaire = 1;
        }return FILE_CORRUPTED;
    }

    else if(typeEntrant == LQUA){
        TuileDam dam;
        FileReaderResult res = getNextEventFromFile<TuileDam>(&dam);
        if(res == OK){
            tDam2tDam_l(&dam, event);
            return OK;
        }else return res;
    }

    else if(typeEntrant == LQL)return getNextEventFromFile<TuileDam_l>(event);

    else if(typeEntrant == JSON){

        for(int i = 0; i < 2; i ++){
            for(; indiceIntermediaire > 0 && indiceIntermediaire <= 4; indiceIntermediaire++){
                TuileDam dam;
                if(paireDam2TuileDam(&intermediaire, &dam, indiceIntermediaire)){
                    tDam2tDam_l(&dam, event);
                    indiceIntermediaire ++;
                    return OK;
                }
            }
            FileReaderResult res = getNextEventFromJSON<PaireDam>(&intermediaire);
            if(res != OK) return res;
            indiceIntermediaire = 1;
        }return FILE_CORRUPTED;
    }

    else if(typeEntrant == CSV){

        for(int i = 0; i < 2; i ++){
            for(; indiceIntermediaire > 0 && indiceIntermediaire <= 4; indiceIntermediaire++){
                TuileDam dam;
                if(paireDam2TuileDam(&intermediaire, &dam, indiceIntermediaire)){
                    tDam2tDam_l(&dam, event);
                    indiceIntermediaire ++;
                    return OK;
                }
            }
            FileReaderResult res = getNextEventFromCsv<PaireDam>(&intermediaire);
            if(res != OK) return res;
            indiceIntermediaire = 1;
        }return FILE_CORRUPTED;
    }
    return CONVERSION_NOT_ALLOWED;
}

FileReader::FileReaderResult FileReader::writeNextEvent(PaireDam *toWrite){
    if(sortant==NULL) return FILE_NOT_OPENNED;

    if(typeSortant == CQUA) return writeNextEventToFile<PaireDam>(toWrite);

    else if(typeSortant == CQL){
        PaireDam_l dam;
        pDam2pDam_l(toWrite, &dam);
        return writeNextEventToFile<PaireDam_l>(&dam);
    }

    else if(typeSortant == LQUA){
        FileReaderResult flag;
        for(int i =0; i < 4; i++){
            TuileDam dam;
            if(paireDam2TuileDam(toWrite, &dam, i+1)){
                flag = writeNextEventToFile<TuileDam>(&dam);
                if(flag != OK) return flag;
            }
        }
        return flag;
    }

    else if(typeSortant == LQL){
        FileReaderResult flag;
        for(int i =0; i < 4; i++){
            TuileDam dam;
            if(paireDam2TuileDam(toWrite, &dam, i+1)){
                TuileDam_l dam_l;
                tDam2tDam_l(&dam, &dam_l);
                flag = writeNextEventToFile<TuileDam_l>(&dam_l);
                if(flag != OK) return flag;
            }
        }
        return flag;
    }

    else if(typeSortant == JSON){
        return writeNextEventToJSON<PaireDam>(toWrite);
    }

    else if(typeSortant == CSV){
        return writeNextEventToCsv<PaireDam>(toWrite);
    }
    return CONVERSION_NOT_ALLOWED;
}


FileReader::FileReaderResult FileReader::writeNextEvent(PaireDam_l *toWrite){
    if(sortant==NULL) return FILE_NOT_OPENNED;

    if(typeSortant == CQL) return writeNextEventToFile<PaireDam_l>(toWrite);

    else if(typeSortant == JSON)
        return writeNextEventToJSON<PaireDam_l>(toWrite);

    else if(typeSortant == CSV)
        return writeNextEventToCsv<PaireDam_l>(toWrite);

    else if(typeSortant == LQL){
        FileReaderResult flag;
        for(int i =0; i < 4; i++){
            TuileDam_l dam_l;
            if(paireDam_l2TuileDam_l(toWrite, &dam_l, i+1)){
                flag = writeNextEventToFile<TuileDam_l>(&dam_l);
                if(flag != OK) return flag;
            }
        }
        return flag;
    }
    return CONVERSION_NOT_ALLOWED;
}


FileReader::FileReaderResult FileReader::writeNextEvent(TuileDam *toWrite){
    if(sortant==NULL) return FILE_NOT_OPENNED;

    if(typeSortant == LQUA) return writeNextEventToFile<TuileDam>(toWrite);
    else if(typeSortant == LQL){
        TuileDam_l dam_l;
        tDam2tDam_l(toWrite, &dam_l);
        return writeNextEventToFile<TuileDam_l>(&dam_l);
    }
    return CONVERSION_NOT_ALLOWED;
}

bool FileReader::isReadFileOpen()
{
    return entrant != NULL;
}

bool FileReader::isWriteFileOpen()
{
    return sortant != NULL;
}



FileReader::FileReaderResult FileReader::writeNextEvent(TuileDam_l *toWrite){
    if(sortant==NULL) return FILE_NOT_OPENNED;

    if(typeSortant == LQL)return writeNextEventToFile<TuileDam_l>(toWrite);
    return CONVERSION_NOT_ALLOWED;
}

void FileReader::closeFiles()
{
    if(entrant != NULL){
        ::fclose(entrant);
        this->entrant = NULL;
        numseqEntrant = 0;
        indiceIntermediaire = 0;
    }
    if (sortant != NULL){
        ::fclose(sortant);
        this->sortant = NULL;
        numseqSortant = 0;
    }
}

FileReader::TypeFichier FileReader::getType() const
{
    return typeEntrant;
}


template<class PaireEvent>
FileReader::FileReaderResult FileReader::writeNextEventToJSON(PaireEvent *event)
{
    event->numero = numseqSortant++;
    if(fprintf(sortant, "%s\n", paireDam2json(*event).dump().c_str()) > 0)
        return OK;
    return WRITE_ERROR;
}


template<class PaireEvent>
FileReader::FileReaderResult FileReader::getNextEventFromJSON(PaireEvent *event)
{
    char string[1000] = {0};
    if(::fgets(string, 1000, entrant) != NULL){
        try{
            *event = json2PaireDam<PaireEvent>(json(json::parse(string)));
        }catch(nlohmann::detail::parse_error e){
            return FILE_CORRUPTED;
        }

        if(numseqEntrant++ == event->numero){
            return OK;
        }else return FILE_CORRUPTED;
    }
    else if(feof(entrant)){
        return END_OF_FILE;
    }
    return READ_ERROR;
}

template<class PaireEvent>
FileReader::FileReaderResult FileReader::getNextEventFromCsv(PaireEvent *event)
{
  //On jette le header
  if(numseqEntrant == 0){
    char string[1000] = {0};
    ::fgets(string, 1000, entrant);
  }
  for(int i = 0 ; i < 4; i++){
    char string[1000] = {0};
    if(::fgets(string, 1000, entrant) != NULL){
      QString line(string);
      QStringList list = line.split(";");
      if(list.size() == 11){
        event->numero = list.front().toInt();
        list.pop_front();
        if(numseqEntrant != event->numero){
          return FILE_CORRUPTED;
        }
        event->tuiles[i].tuile = list.front().toInt();
        list.pop_front();
        event->tuiles[i].ltimestamp = list.front().toInt();
        list.pop_front();
        event->tuiles[i].htimestamp = list.front().toInt();
        list.pop_front();
        event->tuiles[i].temper = list.front().toInt();
        list.pop_front();
        event->tuiles[i].photons = list.front().toInt();
        list.pop_front();
        event->tuiles[i].nbDie = list.front().toInt();
        list.pop_front();
        event->tuiles[i].energie = list.front().toInt();
        list.pop_front();
        event->tuiles[i].mX = list.front().toInt();
        list.pop_front();
        event->tuiles[i].mY = list.front().toInt();
        list.pop_front();
        event->tuiles[i].mZ = list.front().toInt();
        list.pop_front();
      }else return FILE_CORRUPTED;
    }
    else if(feof(entrant)){
      return END_OF_FILE;
    }else return READ_ERROR;
  }
  numseqEntrant++;
  return OK;
}

template<class PaireEvent>
FileReader::FileReaderResult FileReader::writeNextEventToCsv(PaireEvent *event)
{
  if(numseqSortant == 0){
    fprintf(sortant, "Event Number;Tile Number;TimeStamp;Frame Number;Temperature;Photons;Die Number;Energy;mX;mY;mZ\n");
  }
  event->numero = numseqSortant++;
  
  std::string result;
  for(int i =0; i< 4; i++){
    result += std::to_string(event->numero);
    result += ";" + std::to_string(event->tuiles[i].tuile);
    result+= ";" + std::to_string(event->tuiles[i].ltimestamp);
    result+= ";" + std::to_string(event->tuiles[i].htimestamp);
    result+= ";" + std::to_string(event->tuiles[i].temper);
    result+= ";" + std::to_string(event->tuiles[i].photons);
    result+= ";" + std::to_string(event->tuiles[i].nbDie);
    result+= ";" + std::to_string(event->tuiles[i].energie);
    result+= ";" + std::to_string(event->tuiles[i].mX);
    result+= ";" + std::to_string(event->tuiles[i].mY);
    result+= ";" + std::to_string(event->tuiles[i].mZ);
    result+= "\n";
  }
  if(fprintf(sortant, "%s", result.c_str()) > 0)
    return OK;
  return WRITE_ERROR;
}

template<class PaireEvent>
/** @author: Oreste ALLEGRINI @date: 4/11/2020 @brief: supression de la verification de l'incrementation @note: Le code natif reste en commentaire, la nouvelle fonction est active.*/
/** FileReader::FileReaderResult FileReader::getNextEventFromFile(PaireEvent *event)
{

    if(::fread((void *)event,sizeof(PaireEvent),1,entrant) == 1){
        if(numseqEntrant++ == event->numero){
            return OK;
        }else return FILE_CORRUPTED;

    }else if(feof(entrant)){
        return END_OF_FILE;
    }

    return READ_ERROR;
}*/
FileReader::FileReaderResult FileReader::getNextEventFromFile(PaireEvent *event)
{

    if(::fread((void *)event,sizeof(PaireEvent),1,entrant) == 1){
        return OK;
    }
    else if(feof(entrant)){
        return END_OF_FILE;
    }
}
/** Fin de modification */

template<class PaireEvent>
FileReader::FileReaderResult FileReader::writeNextEventToFile(PaireEvent *event)
{
    event->numero = numseqSortant++;
    if(::fwrite((void*)event, 1,sizeof(PaireEvent), sortant)!= 1)
        return OK;
    return WRITE_ERROR;
}
