#include "FiniteElement.h"

void printTotalStiff(FiniteElement &);
void printLoadVector(FiniteElement &);
void printDisplacement(FiniteElement &);
void printUnitStiff(double *);

int main()
{
    clock_t start1 = 0, end1 = 0;
    DWORD start, end;
    start1 = clock();
    start = GetTickCount();
    double *us = new double[4]();
    FiniteElement fe;
    fe.feCircularStructure(100, 100);
    fe.feInput("source&result/fe_test.csv");
    fe.feCalculate();
    fe.feOutput();

    // fe.cstInitialize();
    // fe.feAllocate();
    // fe.cstBuildUnitStiff(1, 1, 1, us);
    // fe.feBuildTotalStiff();
    // fe.feBuildLoadVector();
    // fe.feConjugateGradientPar(fe.TotalStiffness, fe.LoadVector, fe.Displacement, fe.DOF);
    // fe.cstStrainStress();
    // fe.feOutput();

    // printUnitStiff(us);
    // printTotalStiff(fe);
    // printLoadVector(fe);
    // printDisplacement(fe);


    // for (int i = 0; i < 2; i++)
    //     for (int j = 0; j < 3; j++)
    //         cout << fe.CSTriangles[i].strain[j] << "\n";

    // for (int i = 0; i < 2; i++)
    //     for (int j = 0; j < 3; j++)
    //         cout << fe.CSTriangles[i].stress[j] << "\n";

    // ofstream fout("totalstiffness.csv", ios::out);
    // for (int i = 0; i < fe.DOF; i++)
    // {
    //     for (int j = 0; j < fe.DOF; j++)
    //     {
    //         fout << setw(5) << fe.TotalStiffness(i, j) << ",";
    //     }
    //     fout << "\n";
    // }

    end1 = clock();
    cout << "TIME:\n" << (double)(end1 - start1) / CLOCKS_PER_SEC << endl;
    end = GetTickCount();
    cout << (double)(end - start) / 1000 << endl;
    
    return 0;
}

void printTotalStiff(FiniteElement &fe)
{
    for (int i = 0; i < fe.DOF; i++)
    {
        for (int j = 0; j < fe.DOF; j++)
        {
            cout << setw(10) << fe.TotalStiffness(i, j) << "|";
        }
        cout << "\n";
    }
}

void printLoadVector(FiniteElement &fe)
{
    for (int i = 0; i < fe.DOF; i++)
        cout << fe.LoadVector[i] << "\n";
}

void printDisplacement(FiniteElement &fe)
{
    for (int i = 0; i < fe.DOF; i++)
        cout << fe.Displacement[i] << "\n";
}

void printUnitStiff(double *us)
{
    for (int k = 0; k < 4; k++)
        cout << setw(10) << us[k] << " | ";
    cout << "\n";
}