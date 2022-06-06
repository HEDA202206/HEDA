
#ifndef CSTCHANGE_CROSSOVER_H
#include "common.h"
#define CSTCHANGE_CROSSOVER_H

void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A,int& NumOfChormPerPop);
int SelectChrom(vector<double>& A);
void SelectionTournament(int& parent_1, int& parent_2 ,int& NumOfChormPerPop);
void MtnSS_TS(chromosome& a);
void MtnMS_SP(chromosome& a);
void Crossover_CGA(chromosome& pop1, chromosome& pop2);
void Mutation_CGA(chromosome& a);
void GnrTskSchLst_HGA(chromosome& chrom);
void Crossover_HGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_HGA(chromosome& chrom);
void RscLoadAdjust_HGA(vector<chromosome>& chromosomes);
void Crossover_LWSGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_LWSGA(chromosome& chrom);

#endif //CSTCHANGE_CROSSOVER_H
