#include "common.h"
#include "GenOperator.h"
#include "GenerateAChrom.h"

using namespace std;

//{calculate the level of tasks}
void CalculateLevelList() {
    LevelIdOfTask.resize(comConst.NumOfTsk);
    vector<int> InDegree;   //variables for recording the number of parent tasks whose level have not been calculated;
    vector<int> stk;        //a set for recording the index of tasks whose inDegree is equal to 0;
    InDegree.assign(comConst.NumOfTsk, 0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        InDegree[i] = Tasks[i].parents.size();
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if (InDegree[i] == 0) stk.push_back(i);
    }
    int MaxLevel = 0;
    while (!stk.empty()) {
        int v = stk[0];
        LevelIdOfTask[v] = 0;
        for (int i = 0; i < Tasks[v].parents.size(); ++i) {
            if (LevelIdOfTask[Tasks[v].parents[i]] >= LevelIdOfTask[v]) {
                LevelIdOfTask[v] = LevelIdOfTask[Tasks[v].parents[i]] + 1;
            }
        }
        if(LevelIdOfTask[v] + 1> MaxLevel) {
            MaxLevel = LevelIdOfTask[v] + 1;
            TskLstInLvl.resize(MaxLevel);
        }
        TskLstInLvl[LevelIdOfTask[v]].push_back(v);
        stk.erase(stk.begin());
        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            InDegree[Tasks[v].children[i]]--;
            if (InDegree[Tasks[v].children[i]] == 0) {
                stk.push_back(Tasks[v].children[i]);
            }
        }
    }
}

void CalculateDescendants() {
    Descendants.resize(comConst.NumOfTsk);
    for(int i = TskLstInLvl.size()-2; i >= 0; --i) {
        for(int taskId : TskLstInLvl[i]) {
            for(int childId : Tasks[taskId].children) {
                Descendants[taskId].insert(childId);
                Descendants[taskId].insert(Descendants[childId].begin(),Descendants[childId].end());
            }
        }
    }
}

void CalculateAncestors() {
    Ancestors.resize(comConst.NumOfTsk);
    for(int i = 1; i < TskLstInLvl.size(); ++i) {
        for(int taskId : TskLstInLvl[i]) {
            for(int parentId : Tasks[taskId].parents) {
                Ancestors[taskId].insert(parentId);
                Ancestors[taskId].insert(Ancestors[parentId].begin(),Ancestors[parentId].end());
            }
        }
    }
}

bool SortPopOnFitValueByAscend(chromosome& a, chromosome& b) {
    return a.FitnessValue + PrecisionValue < b.FitnessValue;
}

bool SortValueByDescend(pair<int,double>& a, pair<int,double>& b) {
    return a.second > b.second + PrecisionValue;
}

//{sorting the elements in "value" and the indexs are recorded in "ind"; for example, the index of the element with minimal value is recorded in "ind[0]"}
void IndexSortByValueOnAscend(vector<int>& ind, vector<int>& value) {
    for (int i = 0; i < ind.size(); ++i) {
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&value](int v1, int v2) { return value[v1] < value[v2]; });
}

void IndexSortByValueOnAscend(vector<int>& ind, vector<double>& fitness) {
    for (int i = 0; i < ind.size(); ++i) {
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&fitness](int v1, int v2) { return fitness[v1] < fitness[v2]; });
}

void IndexSortByValueOnDescend(vector<int>& ind, vector<int>& value) {
//    vector<double> result;
    for (int i = 0; i < ind.size(); ++i) {
//        result.push_back(value[i]);
        ind[i] = i;
    }
//    sort(ind.begin(), ind.end(), [&result](int v1, int v2) { return result[v1] > result[v2]; });
    sort(ind.begin(), ind.end(), [&value](int v1, int v2) { return value[v1] > value[v2]; });
}

void IndexSortByValueOnDescend(vector<int>& ind, vector<double>& fitness) {
    for (int i = 0; i < ind.size(); ++i) {
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&fitness](int v1, int v2) { return fitness[v1] > fitness[v2]; });
}

//{generate a random number in [Start,End) }
double RandomDouble(int start, int end) {
    double ret = start + rand() % ((end - start) * 1000) / 1000.0;
    return ret;
}

//{generate a random number in [Start,End] }
double RandomDouble2(int start, int end) {
    double ret = start + rand() % (end - start) +  rand() %  1001 / 1000.0;
    return ret;
}

int CalculateParNum(int taskIndex, vector<int>& markList) {
    if(Tasks[taskIndex].parents.empty()) {
        return 0;
    }
    int ret = 0;
    for(int i = 0; i < Tasks[taskIndex].parents.size(); ++i) {
        if(markList[Tasks[taskIndex].parents[i]])
            continue;
        markList[Tasks[taskIndex].parents[i]] = 1;
        ret += 1 + CalculateParNum(Tasks[taskIndex].parents[i], markList);
    }

    return ret;
}

int CalculateSonNum(int taskIndex, vector<int>& markList) {
    if(Tasks[taskIndex].children.empty()) {
        return 0;
    }
    int ret = 0;
    for(int i = 0; i < Tasks[taskIndex].children.size(); ++i) {
        if(markList[Tasks[taskIndex].children[i]])
            continue;
        markList[Tasks[taskIndex].children[i]] = 1;
        ret += 1 + CalculateSonNum(Tasks[taskIndex].children[i], markList);
    }

    return ret;
}
