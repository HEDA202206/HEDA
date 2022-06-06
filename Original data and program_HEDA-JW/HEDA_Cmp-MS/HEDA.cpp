//
// Created by xieyi on 2022/5/29.
//

#include "HEDA.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"

double runHEDA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    ConfigParameter_HEDA();
    CalculateLevelList();
    CalculateDescendants();
    CalculateAncestors();

    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    W_Cal_Average(ww);
    C_Cal_Average(cc);
    Calculate_Rank_b(Rank_b, cc, ww);

    double MaxRank_b = Rank_b[0];
    for (int i = 1; i < comConst.NumOfTsk; ++i) {
        if (MaxRank_b + PrecisionValue < Rank_b[i]) {
            MaxRank_b = Rank_b[i];
        }
    }

    //初始化任务调度顺序概率模型
    vector<int> NumOfAncestors(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        NumOfAncestors[i] = Ancestors[i].size();
    }
    //递归计算，子孙任务的数量
    vector<int> NumOfDescendants(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        NumOfDescendants[i] = Descendants[i].size();
    }
    //总任务数-子孙任务数量
    vector<int> NumOfNonDescendants(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        NumOfNonDescendants[i] = comConst.NumOfTsk - NumOfDescendants[i];
    }

    vector<vector<double> > PMR(comConst.NumOfTsk, vector<double>(comConst.NumOfRsc, 0));
    InitProModelOfResAlc(PMR);
    vector<vector<double> > PMS(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    InitProModelOfTskSch(PMS, NumOfAncestors, NumOfNonDescendants, Rank_b); //PMS[i][k] represents the probability that the k-th scheduled task is task i

    vector<chromosome> Population(Parameter_HEDA.NumOfChromPerPop);
    chromosome Chrom_gb = GnrChr_HEFT_b(Rank_b);

    vector<double> eta_TSO(comConst.NumOfTsk);
    double RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;

    while (RunTime + PrecisionValue < Parameter_HEDA.RunTimeRatioOfStg1 * SchTime) {
        #pragma omp parallel for
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            eta_TSO[i] = pow(Rank_b[i]/MaxRank_b, (1-RunTime/SchTime) * Parameter_HEDA.eta);
        }
        #pragma omp parallel for
        for ( int n = 0; n < Parameter_HEDA.NumOfChromPerPop; ++n) {
            Population[n] = GnrTskLstOfChr_prp(PMS, eta_TSO);
            GnrMS_Evl(Population[n]);
        }

        for (int n = 0; n < Parameter_HEDA.NumOfChromPerPop; ++n) {
            if (Population[n].FitnessValue + PrecisionValue < Chrom_gb.FitnessValue)
                Chrom_gb = Population[n];
        }

        UpdatePMR(PMR, Chrom_gb); UpdatePMS(PMS, Chrom_gb);

        ++iteration;
        RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    }
    while (RunTime + PrecisionValue < SchTime){
        #pragma omp parallel for
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            eta_TSO[i] = pow(Rank_b[i]/MaxRank_b, (1-RunTime/SchTime) * Parameter_HEDA.eta);
        }
        #pragma omp parallel for
        for (int n = 0; n < Parameter_HEDA.NumOfChromPerPop; ++n) {
            Population[n] = GnrTskLstOfChr_prp(PMS, eta_TSO);
            GnrRscLstOfChr(Population[n],PMR);
            DcdEvl(Population[n],true);
        }

        chromosome Chrom_lb = Population[0];
        for (int n = 1; n < Parameter_HEDA.NumOfChromPerPop; ++n) {
            if (Population[n].FitnessValue + PrecisionValue < Chrom_lb.FitnessValue)
                Chrom_lb = Population[n];
        }

        IFBDI(Chrom_lb); LBCAI(Chrom_lb);

        if (Chrom_lb.FitnessValue + PrecisionValue < Chrom_gb.FitnessValue)
            Chrom_gb = Chrom_lb;

        UpdatePMR(PMR, Chrom_gb); UpdatePMS(PMS, Chrom_gb);

        ++iteration;
        RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    }
    SchTime = RunTime;
    return Chrom_gb.FitnessValue;
}
