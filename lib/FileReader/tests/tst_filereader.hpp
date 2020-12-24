#ifndef TST_FILEREADER
#define TST_FILEREADER

#include <QString>
#include <QtTest>
#include "filereader.h"
#include "quadratorevt.h"



inline char *toString(const FileReader::FileReaderResult &p)
{
    return QTest::toString(QString::number(p));
}

inline bool operator==(const PaireDam& lhs, const PaireDam& rhs)
{
    bool flag = true;
    flag = lhs.numero == rhs.numero;
    for(int i = 0; i <4; i++){
        flag &= lhs.tuiles[i].energie    == rhs.tuiles[i].energie   ;
        flag &= lhs.tuiles[i].htimestamp == rhs.tuiles[i].htimestamp;
        flag &= lhs.tuiles[i].ltimestamp == rhs.tuiles[i].ltimestamp;
        flag &= lhs.tuiles[i].mT         == rhs.tuiles[i].mT        ;
        flag &= lhs.tuiles[i].mX         == rhs.tuiles[i].mX        ;
        flag &= lhs.tuiles[i].mY         == rhs.tuiles[i].mY        ;
        flag &= lhs.tuiles[i].mZ         == rhs.tuiles[i].mZ        ;
        flag &= lhs.tuiles[i].nbDie      == rhs.tuiles[i].nbDie     ;
        flag &= lhs.tuiles[i].photons    == rhs.tuiles[i].photons   ;
        flag &= lhs.tuiles[i].temper     == rhs.tuiles[i].temper    ;
        flag &= lhs.tuiles[i].tmax       == rhs.tuiles[i].tmax      ;
        flag &= lhs.tuiles[i].tuile      == rhs.tuiles[i].tuile     ;
        flag &= lhs.tuiles[i].valide     == rhs.tuiles[i].valide    ;
        for(int l =0; l < 8; l++){
            for(int m = 0; m < 8; m++){
                flag &= lhs.images[i].g_pix[l][m]       == rhs.images[i].g_pix[l][m];
                flag &= lhs.images[i].g_die[l/2][m/2]   == rhs.images[i].g_die[l/2][m/2];
            }
        }
    }
    return flag;
}

inline bool operator==(const TuileDam& lhs, const TuileDam& rhs)
{
    bool flag = true;
    flag = lhs.numero == rhs.numero;
    flag &= lhs.tuile.energie    == rhs.tuile.energie   ;
    flag &= lhs.tuile.htimestamp == rhs.tuile.htimestamp;
    flag &= lhs.tuile.ltimestamp == rhs.tuile.ltimestamp;
    flag &= lhs.tuile.mT         == rhs.tuile.mT        ;
    flag &= lhs.tuile.mX         == rhs.tuile.mX        ;
    flag &= lhs.tuile.mY         == rhs.tuile.mY        ;
    flag &= lhs.tuile.mZ         == rhs.tuile.mZ        ;
    flag &= lhs.tuile.nbDie      == rhs.tuile.nbDie     ;
    flag &= lhs.tuile.photons    == rhs.tuile.photons   ;
    flag &= lhs.tuile.temper     == rhs.tuile.temper    ;
    flag &= lhs.tuile.tmax       == rhs.tuile.tmax      ;
    flag &= lhs.tuile.tuile      == rhs.tuile.tuile     ;
    flag &= lhs.tuile.valide     == rhs.tuile.valide    ;
    for(int l =0; l < 8; l++){
        for(int m = 0; m < 8; m++){
            flag &= lhs.image.g_pix[l][m]       == rhs.image.g_pix[l][m];
            flag &= lhs.image.g_die[l/2][m/2]   == rhs.image.g_die[l/2][m/2];
        }
    }
    return flag;
}

inline bool operator==(const PaireDam_l& lhs, const PaireDam_l& rhs)
{
    bool flag = true;
    flag = lhs.numero == rhs.numero;
    for(int i = 0; i <4; i++){
        flag &= lhs.tuiles[i].energie    == rhs.tuiles[i].energie   ;
        flag &= lhs.tuiles[i].htimestamp == rhs.tuiles[i].htimestamp;
        flag &= lhs.tuiles[i].ltimestamp == rhs.tuiles[i].ltimestamp;
        flag &= lhs.tuiles[i].mX         == rhs.tuiles[i].mX        ;
        flag &= lhs.tuiles[i].mY         == rhs.tuiles[i].mY        ;
        flag &= lhs.tuiles[i].mZ         == rhs.tuiles[i].mZ        ;
        flag &= lhs.tuiles[i].nbDie      == rhs.tuiles[i].nbDie     ;
        flag &= lhs.tuiles[i].photons    == rhs.tuiles[i].photons   ;
        flag &= lhs.tuiles[i].temper     == rhs.tuiles[i].temper    ;
        flag &= lhs.tuiles[i].tuile      == rhs.tuiles[i].tuile     ;
    }
    return flag;
}

inline bool operator==(const TuileDam_l& lhs, const TuileDam_l& rhs)
{
    bool flag = true;
    flag = lhs.numero == rhs.numero;
    flag &= lhs.tuile.energie    == rhs.tuile.energie   ;
    flag &= lhs.tuile.htimestamp == rhs.tuile.htimestamp;
    flag &= lhs.tuile.ltimestamp == rhs.tuile.ltimestamp;
    flag &= lhs.tuile.mX         == rhs.tuile.mX        ;
    flag &= lhs.tuile.mY         == rhs.tuile.mY        ;
    flag &= lhs.tuile.mZ         == rhs.tuile.mZ        ;
    flag &= lhs.tuile.nbDie      == rhs.tuile.nbDie     ;
    flag &= lhs.tuile.photons    == rhs.tuile.photons   ;
    flag &= lhs.tuile.temper     == rhs.tuile.temper    ;
    flag &= lhs.tuile.tuile      == rhs.tuile.tuile     ;
    return flag;
}

class FileReaderTest : public QObject
{
    Q_OBJECT

public:
    FileReaderTest(){}
    FileReader freader;
    PaireDam testPaireDam[10];
    PaireDam_l testPaireDam_l[10];
    TuileDam testTuileDam[20];
    TuileDam_l testTuileDam_l[20];

private Q_SLOTS:
    void initTestCase()
    {
        //Initialisation des tableau d'événement test
        for(int i =0; i < 10; i++){
            testPaireDam[i].numero = i;
            testPaireDam_l[i].numero = i;
            for(int j = 0; j < 2; j ++){
                testTuileDam[i*2+j].numero = i*2 + j;
                testTuileDam_l[i*2+j].numero = i*2 + j;

                testPaireDam[i].tuiles[j].tuile = i%3 + j%2 + 1;
                testPaireDam[i].tuiles[j].energie = i*1000 + j*100;
                testPaireDam[i].tuiles[j].htimestamp = i * 1500 + j* 150;
                testPaireDam[i].tuiles[j].ltimestamp = i*2000 + j* 200;
                testPaireDam[i].tuiles[j].mT = i*2500 + j* 250;
                testPaireDam[i].tuiles[j].mX = i*500 + j* 50;
                testPaireDam[i].tuiles[j].mY = i*250 + j* 25;
                testPaireDam[i].tuiles[j].mZ = i*100 + j* 10;
                testPaireDam[i].tuiles[j].nbDie = i*4%16 + 1;
                testPaireDam[i].tuiles[j].photons = i*150 + j*15;
                testPaireDam[i].tuiles[j].temper = i*10 + j*5;
                testPaireDam[i].tuiles[j].tmax = i*5 + j;
                testPaireDam[i].tuiles[j].valide = 1;

                testPaireDam_l[i].tuiles[j].tuile =  i%3 + j%2 + 1;
                testPaireDam_l[i].tuiles[j].energie = i*1000 + j*100;
                testPaireDam_l[i].tuiles[j].htimestamp = i * 1500 + j* 150;
                testPaireDam_l[i].tuiles[j].ltimestamp = i*2000 + j* 200;
                testPaireDam_l[i].tuiles[j].mX = i*500 + j* 50;
                testPaireDam_l[i].tuiles[j].mY = i*250 + j* 25;
                testPaireDam_l[i].tuiles[j].mZ = i*100 + j* 10;
                testPaireDam_l[i].tuiles[j].nbDie = i*4%16 + 1;
                testPaireDam_l[i].tuiles[j].photons = i*150 + j*15;
                testPaireDam_l[i].tuiles[j].temper = i*10 + j*5;

                testTuileDam[i*2 + j].tuile.tuile = i%3 + j%2 + 1;
                testTuileDam[i*2 + j].tuile.energie = i*1000 + j*100;
                testTuileDam[i*2 + j].tuile.htimestamp = i * 1500 + j* 150;
                testTuileDam[i*2 + j].tuile.ltimestamp = i*2000 + j* 200;
                testTuileDam[i*2 + j].tuile.mT = i*2500 + j* 250;
                testTuileDam[i*2 + j].tuile.mX = i*500 + j* 50;
                testTuileDam[i*2 + j].tuile.mY = i*250 + j* 25;
                testTuileDam[i*2 + j].tuile.mZ = i*100 + j* 10;
                testTuileDam[i*2 + j].tuile.nbDie = i*4%16 + 1;
                testTuileDam[i*2 + j].tuile.photons = i*150 + j*15;
                testTuileDam[i*2 + j].tuile.temper = i*10 + j*5;
                testTuileDam[i*2 + j].tuile.tmax = i*5 + j;
                testTuileDam[i*2 + j].tuile.valide = 1;

                testTuileDam_l[i*2 + j].tuile.tuile =  i%3 + j%2 + 1;
                testTuileDam_l[i*2 + j].tuile.energie = i*1000 + j*100;
                testTuileDam_l[i*2 + j].tuile.htimestamp = i * 1500 + j* 150;
                testTuileDam_l[i*2 + j].tuile.ltimestamp = i*2000 + j* 200;
                testTuileDam_l[i*2 + j].tuile.mX = i*500 + j* 50;
                testTuileDam_l[i*2 + j].tuile.mY = i*250 + j* 25;
                testTuileDam_l[i*2 + j].tuile.mZ = i*100 + j* 10;
                testTuileDam_l[i*2 + j].tuile.nbDie = i*4%16 + 1;
                testTuileDam_l[i*2 + j].tuile.photons = i*150 + j*15;
                testTuileDam_l[i*2 + j].tuile.temper = i*10 + j*5;

                for(int k = 0; k < 8; k ++){
                    for(int l = 0; l < 8; l++){
                        testPaireDam[i].images[j].g_pix[k][l] = i*10 + 2*j + l + k;
                        testTuileDam[i*2+ j].image.g_pix[k][l] = i*10 + 2*j + l + k;
                        if(k%2 == 0 && l%2 == 0){
                            testPaireDam[i].images[j].g_die[k/2][l/2] = i*4 + j + l + k;
                            testTuileDam[i*2+j].image.g_die[k/2][l/2] = i*4 + j + l + k;
                        }
                    }
                }
            }
        }

    }

    void cleanupTestCase()
    {
        //Suppression des fichier générés.
    }

    //Comme les tests d'ouverture sont simple, on les fait tous dans un même test
    void openReadSaneFilesAuto()
    {
        //Ouverture d’un fichier *.cqua sain en autodetection
        QCOMPARE(freader.openReadFile(QString(QString(SOURCES_DIR)) + "/testFiles/test.cqua"), FileReader::OK);
        freader.closeFiles();

        //Ouverture d’un fichier *.lqua sain en autodetection
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua"), FileReader::OK);
        freader.closeFiles();

        //Ouverture d’un fichier *.json sain en autodetection
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.json"), FileReader::OK);
        freader.closeFiles();

        //Ouverture d’un fichier *.lql sain en autodetection
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lql"), FileReader::OK);
        freader.closeFiles();

        //Ouverture d’un fichier *.csv sain en autodetection
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.csv"), FileReader::OK);
        freader.closeFiles();

        //Ouverture d’un fichier *.cql sain en autodetection
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cql"), FileReader::OK);
        freader.closeFiles();
    }

    //Comme les tests d'ouverture sont simple, on les fait tous dans un même test
    void openReadSaneFilesManuel()
    {
        //Ouverture d’un fichier *.cqua sain en autodetection
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cqua", FileReader::CQUA), FileReader::OK);
        freader.closeFiles();

        //Ouverture d’un fichier *.lqua sain en autodetection
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua", FileReader::LQUA), FileReader::OK);
        freader.closeFiles();

        //Ouverture d’un fichier *.json sain en autodetection
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.json", FileReader::JSON), FileReader::OK);
        freader.closeFiles();

        //Ouverture d’un fichier *.lql sain en autodetection
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lql", FileReader::LQL), FileReader::OK);
        freader.closeFiles();

        //Ouverture d’un fichier *.cql sain en autodetection
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cql", FileReader::CQL), FileReader::OK);
        freader.closeFiles();

        //Ouverture d’un fichier *.csv sain en autodetection
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.csv", FileReader::CSV), FileReader::OK);
        freader.closeFiles();
    }

    void openReadFilesFromReadFile()
    {
        //Ouverture d’un fichier *.cqua sain en FROME_READ_FILE
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cqua", FileReader::FROMREADFILE), FileReader::INVALID_FILE_TYPE);
        freader.closeFiles();

        //Ouverture d’un fichier *.lqua sain en FROME_READ_FILE
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua", FileReader::FROMREADFILE), FileReader::INVALID_FILE_TYPE);
        freader.closeFiles();

        //Ouverture d’un fichier *.json sain en FROME_READ_FILE
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.json", FileReader::FROMREADFILE), FileReader::INVALID_FILE_TYPE);
        freader.closeFiles();

        //Ouverture d’un fichier *.lql sain en FROME_READ_FILE
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lql", FileReader::FROMREADFILE), FileReader::INVALID_FILE_TYPE);
        freader.closeFiles();

        //Ouverture d’un fichier *.cql sain en FROME_READ_FILE
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cql", FileReader::FROMREADFILE), FileReader::INVALID_FILE_TYPE);
        freader.closeFiles();

        //Ouverture d’un fichier *.csv sain en FROME_READ_FILE
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.csv", FileReader::FROMREADFILE), FileReader::INVALID_FILE_TYPE);
        freader.closeFiles();
    }

    void openReadFilesDontExist()
    {
        //Ouverture d’un fichier *.cqua inexistant
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/inexistant.cqua"), FileReader::OPENNING_ERROR);

        //Ouverture d’un fichier *.lqua inexistant
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/inexistant.lqua"), FileReader::OPENNING_ERROR);

        //Ouverture d’un fichier *.json inexistant
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/inexistant.json"), FileReader::OPENNING_ERROR);

        //Ouverture d’un fichier *.lql inexistant
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/inexistant.lql"), FileReader::OPENNING_ERROR);

        //Ouverture d’un fichier *.cql inexistant
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/inexistant.cql"), FileReader::OPENNING_ERROR);

        //Ouverture d’un fichier *.csv inexistant
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/inexistant.csv"), FileReader::OPENNING_ERROR);
    }

    void openReadSizeCorruptedFile()
    {
        //Ouverture d’un fichier *.cqua avec une mauvaise taille -> utilisation d'un fichier différent : *.lqua
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua", FileReader::CQUA), FileReader::FILE_CORRUPTED);
        freader.closeFiles();

        //Ouverture d’un fichier *.lqua avec une mauvaise taille -> utilisation d'un fichier différent : *.cqua
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cqua", FileReader::LQUA), FileReader::FILE_CORRUPTED);
        freader.closeFiles();

        //Ouverture d’un fichier *.json avec une mauvaise taille -> utilisation d'un fichier différent : *.lqua
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua", FileReader::JSON), FileReader::OK); // En effet on ne peut pas connaitre la taille d'un fichier json
        freader.closeFiles();

        //Ouverture d’un fichier *.csv avec une mauvaise taille -> utilisation d'un fichier différent : *.lqua
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua", FileReader::JSON), FileReader::OK); // En effet on ne peut pas connaitre la taille d'un fichier json
        freader.closeFiles();

        //Ouverture d’un fichier *.lql avec une mauvaise taille -> utilisation d'un fichier différent : *.lqua
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua", FileReader::LQL), FileReader::FILE_CORRUPTED);
        freader.closeFiles();

        //Ouverture d’un fichier *.cql avec une mauvaise taille -> utilisation d'un fichier différent : *.lqua
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua", FileReader::CQL), FileReader::FILE_CORRUPTED);
        freader.closeFiles();
    }

    void openReadFileAlreadyOpenned()
    {
        //Ouverture d’un fichier alors qu'un fichier est déjà ouvert
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cqua");
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cqua"), FileReader::FILE_ALREADY_OPENNED);
        freader.closeFiles();

    }

    void openReadFileWithoutSuffixAuto()
    {
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/sans_suffixe"), FileReader::AUTO_DETECTION_FAILED);
    }

    void openWriteSaneFilesAuto()
    {
        //Ouverture d’un fichier *.cqua sain en autodetection
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua"), FileReader::OK);
        freader.closeFiles();
        QFile file(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua");
        file.remove();

        //Ouverture d’un fichier *.lqua sain en autodetection
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lqua"), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lqua");
        file.remove();

        //Ouverture d’un fichier *.json sain en autodetection
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.json"), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.json");
        file.remove();

        //Ouverture d’un fichier *.lql sain en autodetection
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql"), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql");
        file.remove();

        //Ouverture d’un fichier *.cql sain en autodetection
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cql"), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cql");
        file.remove();

        //Ouverture d’un fichier *.csv sain en autodetection
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.csv"), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.csv");
        file.remove();
    }

    void openWriteSaneFilesManuel()
    {
        //Ouverture d’un fichier *.cqua sain en autodetection
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::CQUA), FileReader::OK);
        freader.closeFiles();
        QFile file(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua");
        QCOMPARE(file.exists(), true);
        file.remove();

        //Ouverture d’un fichier *.lqua sain en autodetection
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::LQUA), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lqua");
        QCOMPARE(file.exists(), true);
        file.remove();

        //Ouverture d’un fichier *.json sain en autodetection
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::JSON), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.json");
        QCOMPARE(file.exists(), true);
        file.remove();

        //Ouverture d’un fichier *.lql sain en autodetection
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::LQL), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql");
        QCOMPARE(file.exists(), true);
        file.remove();

        //Ouverture d’un fichier *.cql sain en autodetection
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::CQL), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cql");
        QCOMPARE(file.exists(), true);
        file.remove();

        //Ouverture d’un fichier *.csv sain en autodetection
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::CSV), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.csv");
        QCOMPARE(file.exists(), true);
        file.remove();

    }

    void openWriteFilesFromReadFile()
    {
        //Ouverture d’un fichier *.cqua sain en autodetection
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cqua");
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::FROMREADFILE), FileReader::OK);
        freader.closeFiles();
        QFile file(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua");
        QCOMPARE(file.exists(), true);
        file.remove();

        //Ouverture d’un fichier *.lqua sain en autodetection
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua");
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::FROMREADFILE), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lqua");
        QCOMPARE(file.exists(), true);
        file.remove();

        //Ouverture d’un fichier *.json sain en autodetection
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.json");
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::FROMREADFILE), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.json");
        QCOMPARE(file.exists(), true);
        file.remove();

        //Ouverture d’un fichier *.lql sain en autodetection
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lql");
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::FROMREADFILE), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql");
        QCOMPARE(file.exists(), true);
        file.remove();

        //Ouverture d’un fichier *.cql sain en autodetection
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cql");
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::FROMREADFILE), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cql");
        QCOMPARE(file.exists(), true);
        file.remove();

        //Ouverture d’un fichier *.csv sain en autodetection
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.csv");
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::FROMREADFILE), FileReader::OK);
        freader.closeFiles();
        file.setFileName(QString(SOURCES_DIR) + "/testFiles/test_ecriture.csv");
        QCOMPARE(file.exists(), true);
        file.remove();
    }

    void openWriteFilesDontExist()
    {
        //Ouverture d’un fichier *.cqua inexistant
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/inexistant/test_ecriture.cqua"), FileReader::OPENNING_ERROR);

        //Ouverture d’un fichier *.lqua inexistant
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/inexistant/test_ecriture.lqua"), FileReader::OPENNING_ERROR);

        //Ouverture d’un fichier *.json inexistant
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/inexistant/test_ecriture.json"), FileReader::OPENNING_ERROR);

        //Ouverture d’un fichier *.lql inexistant
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/inexistant/test_ecriture.lql"), FileReader::OPENNING_ERROR);

        //Ouverture d’un fichier *.cql inexistant
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/inexistant/test_ecriture.cql"), FileReader::OPENNING_ERROR);

        //Ouverture d’un fichier *.csv inexistant
        QCOMPARE(freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/inexistant/test_ecriture.csv"), FileReader::OPENNING_ERROR);
    }

    void openWriteFileFromReadFileNotOpenned()
    {
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture", FileReader::FROMREADFILE), FileReader::READ_FILE_NOT_OPENNED);
        freader.closeFiles();
    }

    void openWriteFileAlreadyOpenned()
    {
        //Ouverture d’un fichier alors qu'un fichier est déjà ouvert
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua");
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua"), FileReader::FILE_ALREADY_OPENNED);
        QFile file(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua");
        file.remove();
        freader.closeFiles();

    }

    void openWriteFileWithoutSuffixAuto()
    {
        QCOMPARE(freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/sans_suffixe"), FileReader::AUTO_DETECTION_FAILED);
    }

    void readCquaFiles()
    {
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cqua");
        for(int i = 0; i < 10; i++){
            PaireDam dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            QCOMPARE(dam, testPaireDam[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam()), FileReader::END_OF_FILE);
        freader.closeFiles();

        //On force la lecture
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cqua");
        for(int i = 0; i < 10; i++){
            PaireDam dam;
            QCOMPARE(freader.getNextEvent(&dam, true), FileReader::OK);
            QCOMPARE(dam, testPaireDam[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cqua");
        for(int i = 0; i < 10; i++){
            PaireDam_l dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            QCOMPARE(dam, testPaireDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cqua");
        for(int i = 0; i < 20; i++){
            TuileDam dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            //ATTENTION LE PASSAGE DE PAIRE DAM VERS TUILE DAM NE MODIFIE PAS LE NUMERO D'EVENEMENT, ON DOIT DONC LE METTRE MANUELLEMENT AFIN DE GARANTIR LA R2USSITE DU TEST
            dam.numero = i;
            QCOMPARE(dam, testTuileDam[i]);
        }
        QCOMPARE(freader.getNextEvent(new TuileDam()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cqua");
        for(int i = 0; i < 20; i++){
            TuileDam_l dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            //ATTENTION LE PASSAGE DE PAIRE DAM VERS TUILE DAM NE MODIFIE PAS LE NUMERO D'EVENEMENT, ON DOIT DONC LE METTRE MANUELLEMENT AFIN DE GARANTIR LA R2USSITE DU TEST
            dam.numero = i;
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new TuileDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();
    }

    void readCqlFiles()
    {
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cql");
        QCOMPARE(freader.getNextEvent(new PaireDam()), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();

        //On force la lecture
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cql");
        for(int i = 0; i < 10; i++){
            PaireDam dam;
            QCOMPARE(freader.getNextEvent(&dam, true), FileReader::OK);
            PaireDam_l test;
            pDam2pDam_l(&dam, &test);
            QCOMPARE(test, testPaireDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cql");
        for(int i = 0; i < 10; i++){
            PaireDam_l dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            QCOMPARE(dam, testPaireDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cql");
        QCOMPARE(freader.getNextEvent(new TuileDam()), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.cql");
        for(int i = 0; i < 20; i++){
            TuileDam_l dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            //ATTENTION LE PASSAGE DE PAIRE DAM VERS TUILE DAM NE MODIFIE PAS LE NUMERO D'EVENEMENT, ON DOIT DONC LE METTRE MANUELLEMENT AFIN DE GARANTIR LA R2USSITE DU TEST
            dam.numero = i;
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new TuileDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();
    }

    void readLquaFiles()
    {
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua");
        for(int i = 0; i < 20; i++){
            PaireDam dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            TuileDam tuile;
            paireDam2TuileDam(&dam, &tuile, dam.tuiles[0].tuile);
            QCOMPARE(tuile, testTuileDam[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam()), FileReader::END_OF_FILE);
        freader.closeFiles();

        //On force la lecture
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua");
        for(int i = 0; i < 20; i++){
            PaireDam dam;
            QCOMPARE(freader.getNextEvent(&dam, true), FileReader::OK);
            TuileDam test;
            paireDam2TuileDam(&dam, &test, dam.tuiles[0].tuile);
            QCOMPARE(test, testTuileDam[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua");
        for(int i = 0; i < 20; i++){
            PaireDam_l dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            TuileDam_l tuile;
            paireDam_l2TuileDam_l(&dam, &tuile, dam.tuiles[0].tuile);
            QCOMPARE(tuile, testTuileDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua");
        for(int i = 0; i < 20; i++){
            TuileDam dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            //ATTENTION LE PASSAGE DE PAIRE DAM VERS TUILE DAM NE MODIFIE PAS LE NUMERO D'EVENEMENT, ON DOIT DONC LE METTRE MANUELLEMENT AFIN DE GARANTIR LA R2USSITE DU TEST
            dam.numero = i;
            QCOMPARE(dam, testTuileDam[i]);
        }
        QCOMPARE(freader.getNextEvent(new TuileDam()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lqua");
        for(int i = 0; i < 20; i++){
            TuileDam_l dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            //ATTENTION LE PASSAGE DE PAIRE DAM VERS TUILE DAM NE MODIFIE PAS LE NUMERO D'EVENEMENT, ON DOIT DONC LE METTRE MANUELLEMENT AFIN DE GARANTIR LA R2USSITE DU TEST
            dam.numero = i;
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new TuileDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();
    }

    void readLqlFiles()
    {

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lql");
        QCOMPARE(freader.getNextEvent(new PaireDam()), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();

        //On force la lecture
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lql");
        for(int i = 0; i < 20; i++){
            PaireDam dam;
            QCOMPARE(freader.getNextEvent(&dam, true), FileReader::OK);
            PaireDam_l dam_l;
            pDam2pDam_l(&dam, &dam_l);
            TuileDam_l test;
            paireDam_l2TuileDam_l(&dam_l, &test, dam.tuiles[0].tuile);
            QCOMPARE(test, testTuileDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lql");
        for(int i = 0; i < 20; i++){
            PaireDam_l dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            TuileDam_l tuile;
            paireDam_l2TuileDam_l(&dam, &tuile, dam.tuiles[0].tuile);
            QCOMPARE(tuile, testTuileDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lql");
        QCOMPARE(freader.getNextEvent(new TuileDam()), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.lql");
        for(int i = 0; i < 20; i++){
            TuileDam_l dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            //ATTENTION LE PASSAGE DE PAIRE DAM VERS TUILE DAM NE MODIFIE PAS LE NUMERO D'EVENEMENT, ON DOIT DONC LE METTRE MANUELLEMENT AFIN DE GARANTIR LA R2USSITE DU TEST
            dam.numero = i;
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new TuileDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();
    }

    void readJsonFiles()
    {
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.json");
        QCOMPARE(freader.getNextEvent(new PaireDam()), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();

        //On force la lecture
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.json");
        for(int i = 0; i < 10; i++){
            PaireDam dam;
            QCOMPARE(freader.getNextEvent(&dam, true), FileReader::OK);
            PaireDam_l test;
            pDam2pDam_l(&dam, &test);
            QCOMPARE(test, testPaireDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.json");
        for(int i = 0; i < 10; i++){
            PaireDam_l dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            QCOMPARE(dam, testPaireDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.json");
        QCOMPARE(freader.getNextEvent(new TuileDam()), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.json");
        for(int i = 0; i < 20; i++){
            TuileDam_l dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            //ATTENTION LE PASSAGE DE PAIRE DAM VERS TUILE DAM NE MODIFIE PAS LE NUMERO D'EVENEMENT, ON DOIT DONC LE METTRE MANUELLEMENT AFIN DE GARANTIR LA R2USSITE DU TEST
            dam.numero = i;
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new TuileDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();
    }

    void readCsvFiles()
    {
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.csv");
        QCOMPARE(freader.getNextEvent(new PaireDam()), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();

        //On force la lecture
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.csv");
        for(int i = 0; i < 10; i++){
            PaireDam dam;
            QCOMPARE(freader.getNextEvent(&dam, true), FileReader::OK);
            PaireDam_l test;
            pDam2pDam_l(&dam, &test);
            QCOMPARE(test, testPaireDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.csv");
        for(int i = 0; i < 10; i++){
            PaireDam_l dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            QCOMPARE(dam, testPaireDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new PaireDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.csv");
        QCOMPARE(freader.getNextEvent(new TuileDam()), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test.csv");
        for(int i = 0; i < 20; i++){
            TuileDam_l dam;
            QCOMPARE(freader.getNextEvent(&dam), FileReader::OK);
            //ATTENTION LE PASSAGE DE PAIRE DAM VERS TUILE DAM NE MODIFIE PAS LE NUMERO D'EVENEMENT, ON DOIT DONC LE METTRE MANUELLEMENT AFIN DE GARANTIR LA R2USSITE DU TEST
            dam.numero = i;
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        QCOMPARE(freader.getNextEvent(new TuileDam_l()), FileReader::END_OF_FILE);
        freader.closeFiles();
    }

    void readCorruptedFile()
    {
        //Lecture d'un fichier corrompu (numéro d'événement NOK et taille OK), Comme la primitive de lecture pour tout les fichiers binaires
        //est la même, on peut ne faire qu'un seul test (cqua vers PaireDam
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/corrupted.cqua");
        QCOMPARE(freader.getNextEvent(new PaireDam), FileReader::FILE_CORRUPTED);
        freader.closeFiles();
    }

    void readCorruptedJson()
    {
        //Lecture d'un fichier corrompu (numéro d'événement NOK et taille OK), Comme la primitive de lecture pour tout les fichiers binaires
        //est la même, on peut ne faire qu'un seul test (cqua vers PaireDam
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/corrupted.json");
        QCOMPARE(freader.getNextEvent(new PaireDam_l), FileReader::FILE_CORRUPTED);
        freader.closeFiles();

        //Lecture d'un fichier corrompu (numéro d'événement NOK et taille OK), Comme la primitive de lecture pour tout les fichiers binaires
        //est la même, on peut ne faire qu'un seul test (cqua vers PaireDam
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/corrupted_complete.json");
        QCOMPARE(freader.getNextEvent(new PaireDam_l), FileReader::FILE_CORRUPTED);
        freader.closeFiles();

        //Lecture d'un fichier corrompu (numéro d'événement NOK et taille OK), Comme la primitive de lecture pour tout les fichiers binaires
        //est la même, on peut ne faire qu'un seul test (cqua vers PaireDam
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/corrupted_champs.json");
        QCOMPARE(freader.getNextEvent(new PaireDam_l), FileReader::FILE_CORRUPTED);
        freader.closeFiles();
    }

    void readCorruptedCsv()
    {
        //Lecture d'un fichier corrompu (numéro d'événement NOK et taille OK), Comme la primitive de lecture pour tout les fichiers binaires
        //est la même, on peut ne faire qu'un seul test (cqua vers PaireDam
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/corrupted.csv");
        QCOMPARE(freader.getNextEvent(new PaireDam_l), FileReader::FILE_CORRUPTED);
        freader.closeFiles();

        //Lecture d'un fichier corrompu (numéro d'événement NOK et taille OK), Comme la primitive de lecture pour tout les fichiers binaires
        //est la même, on peut ne faire qu'un seul test (cqua vers PaireDam
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/corrupted_complete.csv");
        QCOMPARE(freader.getNextEvent(new PaireDam_l), FileReader::FILE_CORRUPTED);
        freader.closeFiles();

        //Lecture d'un fichier corrompu (numéro d'événement NOK et taille OK), Comme la primitive de lecture pour tout les fichiers binaires
        //est la même, on peut ne faire qu'un seul test (cqua vers PaireDam
        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/corrupted_header_missing.csv");
        QCOMPARE(freader.getNextEvent(new PaireDam_l), FileReader::FILE_CORRUPTED);
        freader.closeFiles();
    }

    void writeCquaFile()
    {
        //Test écriture Cqua avec un PaireDam
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua");
        for(int i=0; i < 10; i++)
        {
            QCOMPARE(freader.writeNextEvent(&testPaireDam[i]) , FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua");
        for(int i=0; i < 10; i++)
        {
            PaireDam dam;
            freader.getNextEvent(&dam);
            QCOMPARE(dam, testPaireDam[i]);
        }
        freader.closeFiles();
        QFile file(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua");
        file.remove();

        // Test écriture Cqua avec un PaireDam_l
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua");
        QCOMPARE(freader.writeNextEvent(&testPaireDam_l[0]), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();
        file.remove();

        //Test écriture Cqua avec un TuileDam
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua");
        QCOMPARE(freader.writeNextEvent(&testTuileDam[0]), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();
        file.remove();

        // Test écriture Cqua avec un TuileDam_l
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cqua");
        QCOMPARE(freader.writeNextEvent(&testTuileDam_l[0]), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();
        file.remove();


    }

    void writeLquaFile()
    {
        // Test écriture Lqua avec un PaireDam
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lqua");
        for(int i=0; i < 10; i++)
        {
            QCOMPARE(freader.writeNextEvent(&testPaireDam[i]) , FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lqua");
        for(int i=0; i<20;i++)
        {
            TuileDam dam;
            freader.getNextEvent(&dam);
            QCOMPARE(dam, testTuileDam[i]);
            //QCOMPARE(freader.getNextEvent(new PaireDam()), FileReader::END_OF_FILE);
        }
        freader.closeFiles();
        QFile file(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lqua");
        file.remove();

        // Test écriture Lqua avec un PaireDam_l
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lqua");
        QCOMPARE(freader.writeNextEvent(&testPaireDam_l[0]), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();
        file.remove();

        // Test écriture Lqua avec un TuileDam
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lqua");
        for(int i=0; i< 20; i++){
            QCOMPARE(freader.writeNextEvent(&testTuileDam[i]), FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lqua");
        for(int i=0; i<20; i++){
            TuileDam dam;
            freader.getNextEvent(&dam);
            QCOMPARE(dam, testTuileDam[i]);
        }
        freader.closeFiles();
        file.remove();


        // Test écriture Lqua avec un TuileDam_l
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lqua");
        QCOMPARE(freader.writeNextEvent(&testTuileDam_l[0]), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();
        file.remove();
    }

    void writeCqlFile()
    {
        // Test écriture Cql avec un PaireDam
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cql");
        for(int i=0; i<10; i++)
        {
            QCOMPARE(freader.writeNextEvent(&testPaireDam[i]), FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cql");
        for(int i=0; i<10; i++)
        {
            PaireDam_l dam;
            freader.getNextEvent(&dam);
            QCOMPARE(dam, testPaireDam_l[i]);
        }
        freader.closeFiles();
        QFile file(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cql");
        file.remove();

        // Test écriture Cql avec un PaireDam_l
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cql");
        for(int i=0; i<10; i++)
        {
            QCOMPARE(freader.writeNextEvent(&testPaireDam_l[i]), FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cql");
        for(int i=0; i<10; i++)
        {
            PaireDam_l dam;
            freader.getNextEvent(&dam);
            QCOMPARE(dam, testPaireDam_l[i]);
        }
        freader.closeFiles();
        file.remove();

        // Test écriture Cql avec un TuileDam
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cql");
        QCOMPARE(freader.writeNextEvent(&testTuileDam[0]), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();
        file.remove();

        // Test écriture Cql avec un TuileDam_l
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.cql");
        QCOMPARE(freader.writeNextEvent(&testTuileDam_l[0]), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();
        file.remove();
    }

    void writeLqlTest()
    {
        // Test écriture Lql avec un PaireDam
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql");
        for(int i=0; i<10;i++)
        {
            QCOMPARE(freader.writeNextEvent(&testPaireDam[i]), FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql");
        for(int i=0; i<20; i++)
        {
            TuileDam_l dam;
            freader.getNextEvent(&dam);
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        freader.closeFiles();
        QFile file(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql");
        file.remove();


        // Test écriture Lql avec un PaireDam_l
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql");
        for(int i=0; i<10;i++)
        {
            QCOMPARE(freader.writeNextEvent(&testPaireDam_l[i]), FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql");
        for(int i=0; i<20; i++)
        {
            TuileDam_l dam;
            freader.getNextEvent(&dam);
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        freader.closeFiles();
        file.remove();


        // Test écriture Lql avec un TuileDam
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql");
        for(int i=0; i<20;i++)
        {
            QCOMPARE(freader.writeNextEvent(&testTuileDam[i]), FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql");
        for(int i=0; i<20; i++)
        {
            TuileDam_l dam;
            freader.getNextEvent(&dam);
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        freader.closeFiles();
        file.remove();


        // Test écriture Lql avec un TuileDam_l
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql");
        for(int i=0; i<20;i++)
        {
            QCOMPARE(freader.writeNextEvent(&testTuileDam_l[i]), FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.lql");
        for(int i=0; i<20; i++)
        {
            TuileDam_l dam;
            freader.getNextEvent(&dam);
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        freader.closeFiles();
        file.remove();

    }

    void writeJsonTest()
    {
        // Test écriture Json avec un PaireDam
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.json");
        for(int i=0; i<10; i++)
        {
            QCOMPARE(freader.writeNextEvent(&testPaireDam[i]), FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.json");
        for(int i=0; i<10; i++){
            TuileDam_l dam;
            freader.getNextEvent(&dam);
            dam.numero = i;
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        QFile file(QString(SOURCES_DIR) + "/testFiles/test_ecriture.json");
        file.remove();

        // Test écriture Json avec un PaireDam_l
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.json");
        for(int i=0; i<10; i++)
        {
            QCOMPARE(freader.writeNextEvent(&testPaireDam_l[i]), FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.json");
        for(int i=0; i<20; i++){
            TuileDam_l dam;
            freader.getNextEvent(&dam);
            dam.numero = i;
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        file.remove();

        // Test écriture Json avec un TuileDam
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.json");
        QCOMPARE(freader.writeNextEvent(&testTuileDam[0]), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();
        file.remove();

        // Test écriture Json avec un TuileDam_l
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.json");
        QCOMPARE(freader.writeNextEvent(&testTuileDam_l[0]), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();
        file.remove();
    }

    void writeCsvTest()
    {
        // Test écriture Json avec un PaireDam
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.csv");
        for(int i=0; i<10; i++)
        {
            QCOMPARE(freader.writeNextEvent(&testPaireDam[i]), FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.csv");
        for(int i=0; i<10; i++){
            TuileDam_l dam;
            freader.getNextEvent(&dam);
            dam.numero = i;
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        QFile file(QString(SOURCES_DIR) + "/testFiles/test_ecriture.csv");
        file.remove();

        // Test écriture Json avec un PaireDam_l
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.csv");
        for(int i=0; i<10; i++)
        {
            QCOMPARE(freader.writeNextEvent(&testPaireDam_l[i]), FileReader::OK);
        }
        freader.closeFiles();

        freader.openReadFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.csv");
        for(int i=0; i<20; i++){
            TuileDam_l dam;
            freader.getNextEvent(&dam);
            dam.numero = i;
            QCOMPARE(dam, testTuileDam_l[i]);
        }
        file.remove();

        // Test écriture Json avec un TuileDam
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.csv");
        QCOMPARE(freader.writeNextEvent(&testTuileDam[0]), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();
        file.remove();

        // Test écriture Json avec un TuileDam_l
        freader.openWriteFile(QString(SOURCES_DIR) + "/testFiles/test_ecriture.csv");
        QCOMPARE(freader.writeNextEvent(&testTuileDam_l[0]), FileReader::CONVERSION_NOT_ALLOWED);
        freader.closeFiles();
        file.remove();
    }

    void writeNoDocTest()
    {
        QCOMPARE(freader.writeNextEvent(&testPaireDam[0]), FileReader::FILE_NOT_OPENNED);

        QCOMPARE(freader.writeNextEvent(&testPaireDam_l[0]), FileReader::FILE_NOT_OPENNED);

        QCOMPARE(freader.writeNextEvent(&testTuileDam[0]), FileReader::FILE_NOT_OPENNED);

        QCOMPARE(freader.writeNextEvent(&testTuileDam_l[0]), FileReader::FILE_NOT_OPENNED);
    }

    void readNoDocTest()
    {
        QCOMPARE(freader.getNextEvent(&testPaireDam[0]), FileReader::FILE_NOT_OPENNED);
        QCOMPARE(freader.getNextEvent(&testPaireDam[0], true), FileReader::FILE_NOT_OPENNED);
        QCOMPARE(freader.getNextEvent(&testPaireDam_l[0]), FileReader::FILE_NOT_OPENNED);
        QCOMPARE(freader.getNextEvent(&testTuileDam[0]), FileReader::FILE_NOT_OPENNED);
        QCOMPARE(freader.getNextEvent(&testTuileDam_l[0]), FileReader::FILE_NOT_OPENNED);
    }
};


#endif //TST_FILEREADER
