#ifndef MISC_H
#define MISC_H

#include <QCoreApplication>
#include <QString>
#include <stdio.h>
#include <iostream>

inline bool check_version(QString teststring, QString versionstring) 
{
  teststring=teststring.section('-',0,0);

  int premdigit=teststring.section('.',0,0).toInt();
  int secodigit=teststring.section('.',1,1).toInt();
  int troisdigit=teststring.section('.',2,2).toInt();

  int premdigiv=versionstring.section('.',0,0).toInt();
  int secodigiv=versionstring.section('.',1,1).toInt();
  int troisdigiv=versionstring.section('.',2,2).toInt();

  if (premdigit>premdigiv){
    return true;
  }else if (secodigit>secodigiv && premdigit==premdigiv)
  {
    return true;
  }else if(troisdigit>=troisdigiv && secodigit == secodigiv && premdigit==premdigiv){
    return true;
  }
  else {
    return false;
  }
}

#endif
