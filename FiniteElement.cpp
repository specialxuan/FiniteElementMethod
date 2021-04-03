#include "FiniteElement.h"

int main()
{
    double *us = new double[4]();
    FiniteElement fe;
    fe.feInput();
    fe.cstInitialize();
    // fe.cstBuildUnitStiff(1, 0, 0, us);
    fe.feAllocate();
    fe.feBuildTotalStiff();

    for (int i = 0; i < fe.DOF; i++)
    {
        for (int j = 0; j < fe.DOF; j++)
        {
            cout << setw(10) << fe.TotalStiffness(i, j) << "|";
        }
        cout << "\n";
    }
    


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