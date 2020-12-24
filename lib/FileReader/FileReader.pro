#-------------------------------------------------
#
# Project created by QtCreator 2018-12-14T17:55:39
#
#-------------------------------------------------

QT       -= gui

TARGET = filereader
TEMPLATE = lib

CONFIG+=staticlib

DEFINES += FILEREADER_LIBRARY

SOURCES += filereader.cpp

HEADERS += filereader.h\
        filereader_global.h
INCLUDEPATH += lib/CommonQuadratorHeader/510/

QMAKE_CXXFLAGS += -fvisibility=hidden
QMAKE_CXXFLAGS += -std=c++11


contains(DEFINES, LOCAL_BUILD){
    unix {
        lib.path = FileReader/include
        lib.files = *.h
        target.path = FileReader/lib
#oui.path = Vektor/lib
#        oui.extra = ln -s libvektor.so.$$GIT_VERSION_FIN $$PWD/$$target.path/libvektor.so
        INSTALLS += target lib
    }
}else{
    unix {
        lib.path = /usr/include
        lib.files = *.h
        target.path = /usr/lib
#oui.path = Vektor/lib
#        oui.extra = ln -s libvektor.so.$$GIT_VERSION_FIN $$PWD/$$target.path/libvektor.so
        INSTALLS += target lib
    }
}
