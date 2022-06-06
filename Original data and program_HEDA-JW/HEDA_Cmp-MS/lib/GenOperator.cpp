#include <cstdlib>
#include <GenerateAChrom.h>
#include <unordered_set>
#include "tools.hpp"
#include "common.h"

using namespace std;

//{select two different chromosomes using the tournament method}
//it can only be used in the population where the chromosome have been sorted according fitness from good to bad
void SelectionTournament(int& parent_1, int& parent_2 , int& NumOfChromPerPop) {
    int P1 = rand() % NumOfChromPerPop;
    int P2 = rand() % NumOfChromPerPop;
    while (P1 == P2)
        P2 = rand() % NumOfChromPerPop;
    if (P1 < P2)
        parent_1 = P1;
    else
        parent_1 = P2;

    parent_2 = parent_1;
    while (parent_2 == parent_1) {
        P1 = rand() % NumOfChromPerPop;
        P2 = rand() % NumOfChromPerPop;
        while (P1 == P2)
            P2 = rand() % NumOfChromPerPop;
        if (P1 < P2)
            parent_2 = P1;
        else
            parent_2 = P2;
    }
}

void CrsMS_MP(chromosome& chrom1, chromosome& chrom2) {
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        int TaskId = TskLstInLvl[i][rand() % TskLstInLvl[i].size()]; //select a task randomly
        XY_SWAP(chrom1.RscAlcLst[TaskId], chrom2.RscAlcLst[TaskId], int);
    }
}

//{The MS mutation based multiple point}
void MtnMS_MP(chromosome& ch) {
    int gamma = 1 + rand() % (int(comConst.NumOfTsk / 4) + 1);
    while (gamma--) {
        int i = rand() % comConst.NumOfTsk;
        int j = rand() % Tasks[i].ElgRsc.size();
        ch.RscAlcLst[i] = Tasks[i].ElgRsc[j];
    }
}

//{HGA: single point crossover }
void CrsMS_SP(chromosome& ch1, chromosome& ch2) {
    int CrossPoint = rand() % (comConst.NumOfTsk - 1) + 1;
    for (int i = 0; i < CrossPoint; ++i) {
        XY_SWAP(ch1.RscAlcLst[i], ch2.RscAlcLst[i], int);
    }
}

//{HGA:two point crossover, HGA}
void CrsMS_DP(chromosome& ch1, chromosome& ch2) {
    int point1 = rand() % comConst.NumOfTsk;
    int point2 = rand() % comConst.NumOfTsk;
    while (point1 == point2) {
        point2 = rand() % comConst.NumOfTsk;
    }
    if (point1 > point2) {
        XY_SWAP(point1, point2, int);
    }
    for (int i = point1;  i <= point2; ++i ) {
        XY_SWAP(ch1.RscAlcLst[i], ch2.RscAlcLst[i], int);
    }
}

void GnrTskSchLst_HGA(chromosome& ch) {
    vector<double> w(comConst.NumOfTsk, 0);
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<int> ind(comConst.NumOfTsk);
    vector<vector<double>> TransferTime(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    //{calculate the transfer time between tasks when resource(Rsc) allocation has been determined}
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = ch.RscAlcLst[i];
        if(Tasks[i].parents.size() !=  0){
            for (int j = 0; j < Tasks[i].parents.size(); ++j) {
                int parent = Tasks[i].parents[j];
                int ParRsc = ch.RscAlcLst[parent];
                if(ParRsc != RscIndex){
                    TransferTime[parent][i] = ParChildTranFileSizeSum[parent][i] / VALUE * 8 / XY_MIN(Rscs[RscIndex].bw,Rscs[ParRsc].bw) ;
                }
            }
        }
        w[i] = Tasks[i].length / Rscs[RscIndex].pc;
    }
    Calculate_Rank_b(Rank_b,TransferTime, w);
    IndexSortByValueOnAscend(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ind[comConst.NumOfTsk - i - 1];
        ch.TskSchLst[i] = TaskIndex;
    }
}

//｛HGA crossover｝
void Crossover_HGA(chromosome& ch1, chromosome& ch2) {
    chromosome x1 = ch1;
    chromosome x2 = ch2;
    chromosome y1 = ch1;
    chromosome y2 = ch2;
    CrsMS_SP(x1, x2);
    GnrTskSchLst_HGA(x1);
    DcdEvl(x1, true);
    GnrTskSchLst_HGA(x2);
    DcdEvl(x2, true);

    CrsMS_DP(y1, y2);
    GnrTskSchLst_HGA(y1);
    DcdEvl(y1, true);
    GnrTskSchLst_HGA(y2);
    DcdEvl(y2, true);
    vector<chromosome> sub;
    sub.push_back(x1);
    sub.push_back(x2);
    sub.push_back(y1);
    sub.push_back(y2);
    sort(sub.begin(), sub.end(), SortPopOnFitValueByAscend);
    ch1 = sub[0];
    ch2 = sub[1];
}

void Mutation_HGA(chromosome& ch) {
    chromosome x = ch;
    chromosome y = ch;
    //{single point mutation on x}
    int point = rand() % comConst.NumOfTsk;
    x.RscAlcLst[point] = Tasks[point].ElgRsc[rand() % Tasks[point].ElgRsc.size()];
    //{double point mutation on y}
    int point1 = rand() % comConst.NumOfTsk;
    int point2 = rand() % comConst.NumOfTsk;
    while (point2 == point1) {
        point2 = rand() % comConst.NumOfTsk;
    }
    y.RscAlcLst[point1] = Tasks[point1].ElgRsc[rand() % Tasks[point1].ElgRsc.size()];
    y.RscAlcLst[point2] = Tasks[point2].ElgRsc[rand() % Tasks[point2].ElgRsc.size()];
    GnrTskSchLst_HGA(x);
    DcdEvl(x, true);
    GnrTskSchLst_HGA(y);
    DcdEvl(y, true);
    if ( y.FitnessValue + PrecisionValue < x.FitnessValue ) {
        ch = y;
    } else {
        ch = x;
    }
}

//load balance improvement for HGA
void RscLoadAdjust_HGA(vector<chromosome>& Pop) {
    vector<double> lb(Pop.size());
#pragma omp parallel for
    for(int n =0 ;n < Pop.size(); ++n){
        //calculate the finish times of all task for each resource and find out the maximum
        vector<double> FT(comConst.NumOfRsc,0);
        for(int j = 0 ;j < Pop[n].RscAlcLst.size(); ++j){
            int IndexRsc = Pop[n].RscAlcLst[j];
            if(FT[IndexRsc] < Pop[n].EndTime[j]){
                FT[IndexRsc] = Pop[n].EndTime[j];
            }
        }
        //{find the minimum in FT}
        double min = InfiniteValue;
        for(int j = 0; j < comConst.NumOfRsc; ++j){
            if(min>FT[j]){
                min = FT[j];
            }
        }
        lb[n] = Pop[n].FitnessValue - min;
    }
    //{sort lb}
    vector<int> IndexLb(Pop.size());
    IndexSortByValueOnAscend(IndexLb,lb);             //sorting chromosome by lb from small to large
    //{According to LB, select the last 50% from small to large to improve}
#pragma omp parallel for
    for (int n = Pop.size() / 2; n < Pop.size(); ++n) {
        chromosome chrom = Pop[IndexLb[n]];
        vector<double> ld (comConst.NumOfRsc,0);
        vector<vector<int>> TSK(comConst.NumOfRsc);
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            int RscIndex = chrom.RscAlcLst[i];
            ld[RscIndex] += (1.0 * Tasks[i].length) / Rscs[RscIndex].pc;
            TSK[RscIndex].push_back(i);
        }
        vector<int> Ind(comConst.NumOfRsc);
        IndexSortByValueOnAscend(Ind, ld);                                                           //load sort
        int BigRsc = Ind[Ind.size() - 1];                                                            //obtain the maximum
        int RandTask = TSK[BigRsc][rand() % TSK[BigRsc].size()];                                     //select a task from the Rsc with the largest load
        chrom.RscAlcLst[RandTask] = Tasks[RandTask].ElgRsc[rand() % Tasks[RandTask].ElgRsc.size()];  //reallocation
        GnrTskSchLst_HGA(chrom);
        DcdEvl(chrom, true);
        if ( chrom.FitnessValue + PrecisionValue < Pop[IndexLb[n]].FitnessValue ) {
            Pop[IndexLb[n]] = chrom;
        }
    }
}

//{LWSGA: level swapping (exchange all tasks in the level) }
void Crs_Lvl(chromosome& chrom1, chromosome& chrom2){
    int RandLevel = rand()%TskLstInLvl.size();
    //{Rsc swap}
    for (int i = 0; i < TskLstInLvl[RandLevel].size(); ++i) {
        int TaskIndex = TskLstInLvl[RandLevel][i];
        XY_SWAP(chrom1.RscAlcLst[TaskIndex], chrom2.RscAlcLst[TaskIndex], int);
    }
    //{find the start point of level }
    int pos = 0;
    for(int i = 0; i < RandLevel; ++i){
        pos += TskLstInLvl[i].size();
    }
    for (int i = 0; i < TskLstInLvl[RandLevel].size(); ++i) {
        XY_SWAP(chrom1.TskSchLst[pos], chrom2.TskSchLst[pos], int);
        ++pos;
    }
}

//{LWSGA: exchange two tasks in each level }
void CrsSS_ExcTskInLvl(chromosome& chrom1, chromosome& chrom2) {
    int p1 = 0, p2 = 0;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        int TaskId = TskLstInLvl[i][rand() % TskLstInLvl[i].size()];
        //{Since the task with small level must be in front of the task with large level, the search can be started from the last recorded position}
        for (int j = p1; j < comConst.NumOfTsk; ++j) {
            if (chrom1.TskSchLst[j] == TaskId) {
                p1 = j;
                break;
            }
        }
        for (int j = p2; j < comConst.NumOfTsk; ++j) {
            if (chrom2.TskSchLst[j] == TaskId) {
                p2 = j;
                break;
            }
        }
        XY_SWAP(chrom1.TskSchLst[p1], chrom1.TskSchLst[p2], int);
        XY_SWAP(chrom2.TskSchLst[p1], chrom2.TskSchLst[p2], int);
    }
}

void Crossover_LWSGA(chromosome& ch1, chromosome& ch2) {
    int method = rand() % 3;
    if (method == 0) {
        Crs_Lvl(ch1, ch2);
    } else if (method == 1) {
        CrsSS_ExcTskInLvl(ch1, ch2);
    } else {
        CrsMS_MP(ch1, ch2);
    }
}

//(mutation: exchange two tasks in level)
void MtnSS_ExcTskInLvl(chromosome& chrom) {
    int RandLevel = rand() % TskLstInLvl.size();
    int t1 = TskLstInLvl[RandLevel][rand() % TskLstInLvl[RandLevel].size()];
    int t2 = TskLstInLvl[RandLevel][rand() % TskLstInLvl[RandLevel].size()];
    int p1 = -1, p2 = -1;
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        if(chrom.TskSchLst[i] == t1) {
            p1 = i;
            break;
        }
    }
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        if(chrom.TskSchLst[i] == t2) {
            p2 = i;
            break;
        }
    }
    XY_SWAP(chrom.TskSchLst[p1], chrom.TskSchLst[p2], int);
}

//{select a level, rearrange these tasks in this level and reallocate the resources for these tasks in this level}
void Mtn_rebuild_level(chromosome& ch) {
    int SctLvl = rand() % TskLstInLvl.size();
    vector<int> TemTskLst = TskLstInLvl[SctLvl];
    random_shuffle(TemTskLst.begin(), TemTskLst.end()); // rearrange these tasks
    int StartIndex = 0;
    for(int i = 0; i < SctLvl; ++i) {
        StartIndex = StartIndex + TskLstInLvl[i].size();
    }
    for(int i = 0;i < TemTskLst.size(); ++i) {
        int index = StartIndex + i;
        int TaskId = TemTskLst[i];
        ch.TskSchLst[index] = TaskId;
        int RscIndex = Tasks[TaskId].ElgRsc[rand() % Tasks[TaskId].ElgRsc.size()];
        ch.RscAlcLst[TaskId] = RscIndex;
    }
}

void Mutation_LWSGA(chromosome& ch) {
    int method = rand() % 3;
    if (method == 0) {
        MtnSS_ExcTskInLvl(ch);
    } else if (method == 1) {
        MtnMS_MP(ch);
    } else {
        Mtn_rebuild_level(ch);
    }
}
