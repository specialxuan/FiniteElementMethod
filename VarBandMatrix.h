#include <iostream>
#include <math.h>
using namespace std;

class VarBandMatrix
{
private:
    double *matrix;
    int *diagelement; // the location of diagonal element
    int DIM;
    int lastelement; // upper limit
    double zero = 0;
    double Max;
    bool normalized = 0;

public:
    VarBandMatrix();
    VarBandMatrix(VarBandMatrix &);
    ~VarBandMatrix();

    void initialize(int *, int);
    void initialize(VarBandMatrix &);
    double normalize();
    double denormalize();
    double &operator()(int, int);

    friend class FiniteElement;
};

VarBandMatrix::VarBandMatrix()
{
    matrix = NULL;
    diagelement = NULL;
    DIM = 0;
    lastelement = 0;
}

VarBandMatrix::VarBandMatrix(VarBandMatrix &vbm)
{
    lastelement = vbm.lastelement;
    DIM = vbm.DIM;

    if (vbm.matrix != NULL)
    {
        matrix = new double[lastelement]();
        memcpy(matrix, vbm.matrix, vbm.lastelement * sizeof(double));
    }
    if (vbm.diagelement != NULL)
    {
        diagelement = new int[DIM]();
        memcpy(diagelement, vbm.diagelement, vbm.DIM * sizeof(int));
    }
}

VarBandMatrix::~VarBandMatrix()
{
    delete[] matrix;
    delete[] diagelement;
}

void VarBandMatrix::initialize(int *iv, int dim)
{
    DIM = dim;
    if (iv != NULL)
    {
        diagelement = new int[DIM]();
        memcpy(diagelement, iv, DIM * sizeof(int));
    }
    lastelement = diagelement[DIM - 1];
    matrix = new double[lastelement]();
}

void VarBandMatrix::initialize(VarBandMatrix &vbm)
{
    lastelement = vbm.lastelement;
    DIM = vbm.DIM;

    if (vbm.matrix != NULL)
    {
        matrix = new double[lastelement]();
        memcpy(matrix, vbm.matrix, vbm.lastelement * sizeof(double));
    }
    if (vbm.diagelement != NULL)
    {
        diagelement = new int[DIM]();
        memcpy(diagelement, vbm.diagelement, vbm.DIM * sizeof(int));
    }
}

double VarBandMatrix::normalize()
{
    if (normalized == 1)
        return Max;

    for (int i = 0; i < lastelement; i++)
        if (fabs(matrix[i]) > Max)
            Max = matrix[i];
    
    for (int i = 0; i < lastelement; i++)
        matrix[i] = matrix[i] / Max;

    normalized = 1;

    return Max;
}

double VarBandMatrix::denormalize()
{
    if (normalized == 0)
        return Max;

    for (int i = 0; i < lastelement; i++)
        matrix[i] = matrix[i] * Max;

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
        return matrix[diagelement[i] - 1];
    }
    else if (j > i)
    {
        if ((diagelement[j] - j + i) > diagelement[j - 1])
            return matrix[diagelement[j] - j + i - 1];
        else
            return zero;
    }
    else if (i > j)
    {
        if ((diagelement[i] - i + j) > diagelement[i - 1])
            return matrix[diagelement[i] - i + j - 1];
        else
            return zero;
    }
}
