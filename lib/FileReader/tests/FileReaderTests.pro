#-------------------------------------------------
#
# Project created by QtCreator 2019-09-16T14:50:24
#
#-------------------------------------------------

QT       += testlib

QT       -= gui

CONFIG += c++11
TARGET = tst_filereader
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += tst_filereader.cpp \
    ../filereader.cpp

INCLUDEPATH += ../../CommonQuadratorHeader/510 ../
DEFINES += SRCDIR=\\\"$$PWD/\\\"
