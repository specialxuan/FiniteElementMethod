#include "FiniteElement.h"

int main()
{
    clock_t start1 = 0, end1 = 0;
    DWORD start, end;
    start1 = clock();
    start = GetTickCount();
    // double *us = new double[4]();
    FiniteElement fe;
    fe.feCircularStructure(4, 5);
    // fe.feInput();
    fe.feInput("source&result/fe_test.csv");
    // fe.feCalculate();
    // fe.feOutput();
    
    // end1 = clock();
    // cout << (double)(end1 - start1) / CLOCKS_PER_SEC << endl;
    // end = GetTickCount();
    // cout << (double)(end - start) / 1000 << endl;
    fe.cstInitialize();
    // fe.cstBuildUnitStiff(1, 0, 0, us);
    fe.feAllocate();
    fe.feBuildTotalStiff();
    // fe.feBuildLoadVector();
    // fe.feConjugateGradientPar(fe.TotalStiffness, fe.LoadVector, fe.Displacement, fe.DOF);
    // fe.cstStrainStress();

    // for (int i = 0; i < 2; i++)
    //     for (int j = 0; j < 3; j++)
    //         cout << fe.CSTriangles[i].strain[j] << "\n";

    // for (int i = 0; i < 2; i++)
    //     for (int j = 0; j < 3; j++)
    //         cout << fe.CSTriangles[i].stress[j] << "\n";
    

    // for (int i = fe.NFIN * 2; i < fe.DOF; i++)
    // {
    //     cout << fe.Displacement[i] << "\n";
    // }

    // for (int i = fe.NFIN * 2; i < fe.DOF; i++)
    // {
    //     cout << fe.LoadVector[i] << "\n";
    // }
    

    // for (int i = fe.NFIN * 2; i < fe.DOF; i++)
    // {
    //     for (int j = fe.NFIN * 2; j < fe.DOF; j++)
    //     {
    //         cout << setw(5) << fe.TotalStiffness(i, j) << "|";
    //     }
    //     cout << "\n";
    // }

    // ofstream fout("totalstiffness.csv", ios::out);
    // for (int i = fe.NFIN * 2; i < fe.DOF; i++)
    // {
    //     for (int j = fe.NFIN * 2; j < fe.DOF; j++)
    //     {
    //         fout << setw(5) << fe.TotalStiffness(i, j) << ",";
    //     }
    //     fout << "\n";
    // }

    // for (int i = 0; i < 4; i++)
    // {
    //     cout << us[i] << "\n";
    // }
    
    // for (int i = 0; i < 3; i++)
    // {
    //     cout << fe.CSTriangles[0].a[i] << "\n";
    //     cout << fe.CSTriangles[0].b[i] << "\n";
    //     cout << fe.CSTriangles[0].c[i] << "\n";
    // }

    // cout << fe.CSTriangles[0].area << "\n";

    return 0;
}