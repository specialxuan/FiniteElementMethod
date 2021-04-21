#include "FiniteElement.h"

int main()
{
    clock_t start1 = 0, end1 = 0;
    DWORD start, end;
    start1 = clock();
    start = GetTickCount();
    FiniteElement fe;
    fe.feCircularStructure(121, 121);
    fe.feInput("source&result/fe_test.csv");
    fe.feCalculate();
    fe.feOutput();

    end1 = clock();
    cout << "TIME:\n" << (double)(end1 - start1) / CLOCKS_PER_SEC << endl;
    end = GetTickCount();
    cout << (double)(end - start) / 1000 << endl;
    
    return 0;
}