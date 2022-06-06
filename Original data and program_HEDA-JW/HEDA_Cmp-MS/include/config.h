
#ifndef CSTCHANGE_CONFIG_H
#include "common.h"
#define CSTCHANGE_CONFIG_H
void DeleteFirstLineInFile(string fileName);
int ReadID(string id);
void ReadFile(string XmlFile,string RscAlcFile);
void TaskCombine();
void ClearALL();
void ConfigParameter_CGA();
void ConfigParameter_HGA();
void ConfigParameter_MOELS();
void ConfigParameter_LWSGA();
void ConfigParameter_HPSO();
void ConfigParameter_TSEDA();
void ConfigParameter_HEDA();

void ConfigParameter_ADBRKGA();
#endif //CSTCHANGE_CONFIG_H
