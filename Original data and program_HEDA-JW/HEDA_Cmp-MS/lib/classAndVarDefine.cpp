#include "classAndVarDefine.h"

vector<Task> Tasks;
vector<Task> OriginalTasks;
vector<int> Oid;
vector<int> Nid;
vector<vector<int>> TskLstInLvl;
vector<int> LevelIdOfTask;
vector<vector<double>> ParChildTranFileSizeSum;
vector<Resource> Rscs;
double MinBW ;
vector<double> MaxLd;
vector<set<int> > Descendants;
vector<set<int> > Ancestors;
Paramet_HGA Parameter_HGA;
Paramet_MOELS Parameter_MOELS;
Paramet_LWSGA Parameter_LWSGA;
Paramet_CGA Parameter_CGA;
Paramet_HPSO Parameter_HPSO;
Paramet_TSEDA Parameter_TSEDA;
Paramet_HEDA Parameter_HEDA;
Paramet_ADBRKGA Parameter_ADBRKGA;
ComConst comConst;
double ModelScale;
