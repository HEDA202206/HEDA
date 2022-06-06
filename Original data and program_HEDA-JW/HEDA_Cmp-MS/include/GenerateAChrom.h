
#ifndef CSTCHANGE_GENERATEACHROM_H
#include "common.h"
#define CSTCHANGE_GENERATEACHROM_H

void W_Cal_Average(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& rankList, vector<vector<double>>& c, vector<double>& w);
void Calculate_Rank_t(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w);
void IntChr(chromosome& chrom);
vector<int> GnrSS_Lvl();
vector<double> GnrDecimalsByAscend();
double DcdEvl(chromosome& chrom, bool isForward);
double NrmDcd(chromosome& chrom, bool isForward);
double HrsDcd_CTP(chromosome& chrom);
double HrsDcd_EFT_ADBRKGA(chromosome& ch);
double GnrMS_Evl(chromosome& chrom);
double IFBDI(chromosome& ch);
void LBCAI(chromosome& ch);
chromosome GnrChr_HEFT_b(vector<double> Rank_b);
chromosome GnrPrtByRank_Rnd(vector<double>& Rank);
chromosome GnrPrtByRank_EFT(vector<double>& Rank);
void SeletRsc_EFT(chromosome& ch, vector<set<double>>& ITL, int& TaskId, int& RscId, double& FinalStartTime, double& FinalEndTime);
double FindIdleTimeSlot(set<double>& ITLofRscId, double& ExeTime, double& ReadyTime);
void UpdateITL(set<double>& ITLofRscId,double& StartTime,double& EndTime);
void RepairMapAndGnrRscAlcLst(chromosome& ch);
int FindNearestRscId(int TaskId, double value);
void RepairPriorityAndGnrSchOrd(chromosome& chrom);
void UpdateParticle(chromosome &ch,chromosome &Pbest, chromosome &Gbest, double &runtime, double &SchTime);
void InitProModelOfResAlc(vector<vector<double> >& PMR);
void InitProModelOfTskSch(vector<vector<double> >& PMS, vector<int>& NumOfAncestors,vector<int>& NumOfNonDescendants, vector<double>& Rank_b);
void GnrRscLstOfChr(chromosome& chrom, vector<vector<double> >& PMR);
chromosome GnrTskLstOfChr_prp(vector<vector<double> >& PMS, vector<double>& eta_TSO);
void UpdatePMR(vector<vector<double>>& PMR, chromosome& bstChrom);
void UpdatePMS(vector<vector<double>>& PMS, chromosome& bstChrom);

#endif //CSTCHANGE_GENERATEACHROM_H
