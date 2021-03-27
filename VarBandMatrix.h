#include <iostream>
#include <math.h>
using namespace std;

class VarBandMatrix
{
private:
    double *matrix;
    int *IV; // the location of diagonal element
    int DIM;
    int NSI; // upper limit
    double zero = 0;

public:
    VarBandMatrix();
    VarBandMatrix(VarBandMatrix &);
    ~VarBandMatrix();

    void initialize(int *, int);
    void initialize(VarBandMatrix &);
    double maximum();
    // VarBandMatrix &operator/(); TODO
    double &operator()(int, int);
};

VarBandMatrix::VarBandMatrix()
{
    matrix = NULL;
    IV = NULL;
    DIM = 0;
    NSI = 0;
}

VarBandMatrix::VarBandMatrix(VarBandMatrix &vbm)
{
    NSI = vbm.NSI;
    DIM = vbm.DIM;

    if (vbm.matrix != NULL)
    {
        matrix = new double[NSI]();
        memcpy(matrix, vbm.matrix, vbm.NSI * sizeof(double));
    }
    if (vbm.IV != NULL)
    {
        IV = new int[DIM]();
        memcpy(IV, vbm.IV, vbm.DIM * sizeof(int));
    }
}

VarBandMatrix::~VarBandMatrix()
{
    delete[] matrix;
    delete[] IV;
}

void VarBandMatrix::initialize(int *iv, int dim)
{
    DIM = dim;
    if (iv != NULL)
    {
        IV = new int[DIM]();
        memcpy(IV, iv, DIM * sizeof(int));
    }
    NSI = IV[DIM - 1];
    matrix = new double[NSI]();
}

void VarBandMatrix::initialize(VarBandMatrix &vbm)
{
    NSI = vbm.NSI;
    DIM = vbm.DIM;

    if (vbm.matrix != NULL)
    {
        matrix = new double[NSI]();
        memcpy(matrix, vbm.matrix, vbm.NSI * sizeof(double));
    }
    if (vbm.IV != NULL)
    {
        IV = new int[DIM]();
        memcpy(IV, vbm.IV, vbm.DIM * sizeof(int));
    }
}

double VarBandMatrix::maximum()
{
    double Max = 0;
    for (int i = 0; i < NSI; i++)
        if (fabs(matrix[i]) > Max)
            Max = matrix[i];
    return Max;
}

double &VarBandMatrix::operator()(int i, int j)
{
    if (zero != 0)
    {
        cout << "ERROR: The last assignment was out of bandwith!\n";
        zero = 0;
    }

    if (i == j)
    {
        return matrix[IV[i] - 1];
    }
    else if (j > i)
    {
        if ((IV[j] - j + i) > IV[j - 1])
            return matrix[IV[j] - j + i - 1];
        else
            return zero;
    }
    else if (i > j)
    {
        if ((IV[i] - i + j) > IV[i - 1])
            return matrix[IV[i] - i + j - 1];
        else
            return zero;
    }
} // TODO
