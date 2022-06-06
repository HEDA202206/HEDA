#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "HEFT.h"
#include "HGA.h"
#include "LWSGA.h"
#include "HPSO.h"
#include "HEDA.h"

using namespace std;

int main() {
    srand((int) time(0));
    //{set the runtime (termination) time of the algorithm} -xy4
    map<string, double> SchTime;
    SchTime["Montage25_1.0"] = 3.545;
    SchTime["Montage50_1.0"] = 10.647;
    SchTime["Montage100_1.0"] = 19.263;

    SchTime["CyberShake30_1.0"] = 4.105;
    SchTime["CyberShake50_1.0"] = 7.637;
    SchTime["CyberShake100_1.0"] = 22.135;

    SchTime["Epigenomics24_1.0"] = 3.063;
    SchTime["Epigenomics47_1.0"] = 7.177;
    SchTime["Epigenomics100_1.0"] = 27.421;

    SchTime["Ligo30_1.0"] = 3.806;
    SchTime["Ligo50_1.0"] = 5.499;
    SchTime["Ligo100_1.0"] = 19.602;

    //{clear "result"}
    ofstream outfile("../result.txt", ios::out);
    outfile.close();

    string Model, NumOfTask, RscAvlRatio;
    do {
        string StrLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(0);
        }
        getline(iFile, StrLine);
        if (StrLine.size() < 1) {
            cout << "Empty input file" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile;
        string RscAlcFile;
        istringstream is(StrLine);
        is >> Model >> NumOfTask >> RscAvlRatio;  //NumOfTask, RscAvlRatio
        XmlFile = Model + "_" + NumOfTask + "_0.xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        cout <<endl<< Model << " " << NumOfTask << " " << RscAvlRatio << " ";

        double HEFT_SchTime  = 0 ;
        double HEFT_Result = runHEFT(XmlFile, RscAlcFile, HEFT_SchTime);
        ClearALL();

        double HGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int HGA_Iteration = 0;
        double HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration);
        ClearALL();

        double LWSGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int LWSGA_Iteration = 0;
        double LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration);
        ClearALL();

        double HPSO_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int HPSO_Iteration = 0;
        double HPSO_Result = runHPSO(XmlFile, RscAlcFile, HPSO_SchTime, HPSO_Iteration);
        ClearALL();

        double HEDA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int HEDA_Iteration = 0;
        double HEDA_Result = runHEDA(XmlFile, RscAlcFile, HEDA_SchTime, HEDA_Iteration);
        ClearALL();

        //results are written into the file
        outfile.open("../result.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(5);
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                << HEFT_Result << " " << HEFT_SchTime << " "
                << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " "
                << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " "
                << HPSO_Result << " " << HPSO_SchTime << " " << HPSO_Iteration << " "
                << HEDA_Result  << " " << HEDA_SchTime << " " << HEDA_Iteration  << " "
                << endl;
        outfile.close();
        DeleteFirstLineInFile("../fileList.txt"); //delete the first line in the file
    } while (1);
    //return 0;
}
