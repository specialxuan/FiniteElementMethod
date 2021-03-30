#include "VarBandMatrix.h"
#include <Windows.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
using namespace std;

class FiniteElement
{
public:
    double EPS;   // epsilon, calculation accuracy
    double MaxTS; // maximum value in total stiffness matrix
    double MaxLV; // maximum value in load vector

    int TNN; // total number of nodes
    int DOF; // degree of freedom
    int NCST; // number of constant strain triangle unit
    int NOL; // number of loads

    struct Node // parameters of nodes
    {
        double xcn; // X coordinate of nodes
        double ycn; // Y coordinate of nodes
        double zcn; // Z coordinate of nodes
    } * Nodes;      // parameters of nodes

    struct ConstantStrainTriangle // constant strain triangle unit
    {
        int nodes[3];     // 3 nodes of the unit
        double a[3];      // paramaters of shape function
        double b[3];      // paramaters of shape function
        double c[3];      // paramaters of shape function
        double area;      // paramaters of shape function
        double mu;        // poisson ratio
        double elastic;   // elastic module
        double thickness; // thickness of the unit
        double strain[3]; // strain of the unit
        double stress[3]; // stress of the unit
        double rou;       // density of the unit
    } * CSTriangles;      // constant strain triangle unit

    struct Load // parameters of loads
    {
        int unit;      // the kind of unit with load
        int num;       // the number of unit with load
        int node[2];   // the node with load
        int pli;       // the plane of the load's in
        int kol;       // the kind of load
        double vol[3]; // the value of load
    } * Loads;         // parameters of loads

    VarBandMatrix TotalStiffness; // total stiffness matrix
    double *LoadVector;           // load vector
    double *Displacement;         // displacement of nodes

    bool ProgressBar; // open progress bar
    bool Parallel;    // open parallel
    int status;

    // allocate total stiffness matrix, load vector and displacement vector
    bool feAllocate()
    {
        return 0;
    }

    // initialize constant strain triangle unit
    bool cstInitialize()
    {
        for (int k = 0; k < NCST; k++)
        {
            double &a1 = CSTriangles[k].a[0]; // paramaters of shape function
            double &a2 = CSTriangles[k].a[1];
            double &a3 = CSTriangles[k].a[2];
            double &b1 = CSTriangles[k].b[0];
            double &b2 = CSTriangles[k].b[1];
            double &b3 = CSTriangles[k].b[2];
            double &c1 = CSTriangles[k].c[0];
            double &c2 = CSTriangles[k].c[1];
            double &c3 = CSTriangles[k].c[2];
            double x1 = Nodes[CSTriangles[k].nodes[0]].xcn; // cooridnates of nodes of the unit
            double x2 = Nodes[CSTriangles[k].nodes[1]].xcn;
            double x3 = Nodes[CSTriangles[k].nodes[2]].xcn;
            double y1 = Nodes[CSTriangles[k].nodes[0]].ycn;
            double y2 = Nodes[CSTriangles[k].nodes[1]].ycn;
            double y3 = Nodes[CSTriangles[k].nodes[2]].ycn;

            a1 = x2 * y3 - x3 * y2; // paramaters of shape function
            a2 = x3 * y1 - x1 * y3;
            a3 = x1 * y2 - x2 * y1;
            b1 = y2 - y3;
            b2 = y3 - y1;
            b3 - y1 - y2;
            c1 = x3 - x2;
            c2 = x1 - x3;
            c3 = x2 - x1;

            CSTriangles[k].area = (a1 + a2 + a3) / 2; // area of the unit
        }

        return 0;
    }

    // build unit stiffness matrix of constant strain triangle unit
    bool cstBuildUnitStiff(int k, int r, int s, double *us) // k is unit number, r and s is section number, us is a section in unit stiffness matrix
    {
        if (k < 0 || r < 0 || r > 2 || s < 0 || s > 2 || us == NULL)
            return sfPrintError(16);

        double *a = CSTriangles[k].a; // paramaters of shape function
        double *b = CSTriangles[k].b;
        double *c = CSTriangles[k].c;

        double Area = CSTriangles[k].area; // paramaters of the unit
        double Mu = CSTriangles[k].mu;
        double E = CSTriangles[k].elastic;
        double T = CSTriangles[k].thickness;

        double tmp = T * E * T / (4 * (1 - Mu * Mu)); // temperatary coefficient

        us[0 * 2 + 0] = (b[r] * b[s] + (1 - Mu) / 2 * c[r] * c[s]) * tmp; // unit stiffness matrix
        us[0 * 2 + 1] = (Mu * c[r] * b[s] + (1 - Mu) / 2 * b[r] * c[s]) * tmp;
        us[1 * 2 + 0] = ((1 - Mu) / 2 * c[r] * b[s] + Mu * b[r] * c[s]) * tmp;
        us[1 * 2 + 1] = ((1 - Mu) / 2 * b[r] * b[s] + c[r] * c[s]) * tmp;

        return 0;
    }

    // build unit load vector of constant strain triangle unit
    bool cstBuildUnitLoad(int i, double *llv) // k is load number, llv is local load vector
    {
        Load load = Loads[i]; // paramaters of load
        int n1 = load.node[0] - 1;
        int n2 = load.node[1] - 1;
        double A = CSTriangles[load.num].area;
        double t = CSTriangles[load.num].thickness;
        double x1 = Nodes[CSTriangles[load.num].nodes[n1]].xcn;
        double x2 = Nodes[CSTriangles[load.num].nodes[n2]].xcn;
        double y1 = Nodes[CSTriangles[load.num].nodes[n1]].ycn;
        double y2 = Nodes[CSTriangles[load.num].nodes[n2]].ycn;
        double l = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

        switch (load.kol)
        {
        case 0: // nodal force
            llv[n1 * 2 + 0] = load.vol[0];
            llv[n1 * 2 + 1] = load.vol[1];
            break;
        case 1: // body force
            llv[0] = llv[2] = llv[4] = -A * t * load.vol[0] / 3;
            llv[1] = llv[3] = llv[5] = -A * t * load.vol[1] / 3;
            break;
        case 2: // horizontal distribution force
            llv[n1 * 2 + 0] = llv[n2 * 2 + 0] = load.vol[0] * t * l / 2;
            llv[n1 * 2 + 1] = llv[n2 * 2 + 1] = load.vol[1] * t * l / 2;
            break;
        case 3: // triangle distribution force
            llv[n1 * 2 + 0] = load.vol[0] * t * l / 3;
            llv[n1 * 2 + 1] = load.vol[1] * t * l / 3;
            llv[n2 * 2 + 0] = load.vol[0] * t * l / 6;
            llv[n2 * 2 + 1] = load.vol[1] * t * l / 6;
            break;
        default:
            break;
        }
    }

    // calculate strain and stress of constant strain triangle unit
    bool cstStrainStress()
    {
        int p[3] = {0}; // p is a temperary vector for i0j0

        for (int k = 0; k < NCST; k++)
        {
            p[0] = 2 * (CSTriangles[k].nodes[0] - 1); // match the displacement with nods
            p[1] = 2 * (CSTriangles[k].nodes[1] - 1);
            p[2] = 2 * (CSTriangles[k].nodes[2] - 1);

            double u1 = Displacement[p[0] + 0]; // displacements
            double v1 = Displacement[p[0] + 1];
            double u2 = Displacement[p[1] + 0];
            double v2 = Displacement[p[1] + 1];
            double u3 = Displacement[p[2] + 0];
            double v3 = Displacement[p[2] + 1];

            double b1 = CSTriangles[k].b[0]; // paramaters of shape function
            double b2 = CSTriangles[k].b[1];
            double b3 = CSTriangles[k].b[2];
            double c1 = CSTriangles[k].c[0];
            double c2 = CSTriangles[k].c[1];
            double c3 = CSTriangles[k].c[2];

            CSTriangles[k].strain[0] = u1 * b1 + u2 * b2 + u3 * b3; // strain
            CSTriangles[k].strain[1] = v1 * c1 + v2 * c2 + v3 * c3;
            CSTriangles[k].strain[2] = u1 * c1 + u2 * c2 + u3 * c3 + v1 * b1 + v2 * b2 + v3 * b3;

            CSTriangles[k].stress[0] = CSTriangles[k].elastic * CSTriangles[k].strain[0]; //stress
            CSTriangles[k].stress[1] = CSTriangles[k].elastic * CSTriangles[k].strain[1];
            CSTriangles[k].stress[2] = CSTriangles[k].elastic * CSTriangles[k].strain[2];
        }
    }

    // build total stiffness matrix
    bool feBuildTotalStiff()
    {
        int dofNode = 2, numNode = 3;
        double us[dofNode * dofNode] = {0}; // unit stiffness matrix
        int p[3] = {0};                     // p is a temperary vector for i0j0

        for (int k = 0; k < NCST; k++) // for every nodes
        {
            p[0] = dofNode * (CSTriangles[k].nodes[0] - 1); // match the displacement with nods
            p[1] = dofNode * (CSTriangles[k].nodes[1] - 1);
            p[2] = dofNode * (CSTriangles[k].nodes[2] - 1);

            for (int i = 0; i < numNode; i++) // for every nodes of the unit
            {
                for (int j = i; j < numNode; j++) // for every nodes of the unit
                {
                    if (cstBuildUnitStiff(k, i, j, us)) // build unit stiffness matrix
                        return sfPrintError(7);
                    for (int m = 0; m < dofNode; m++)
                        for (int n = 0; n <= dofNode; n++)
                            TotalStiffness(p[i] + m, p[j] + n) += us[m * dofNode + n]; // superpose
                }
            }
        }
        return 0;
    }

    // build load vector
    bool sfBuildLoadVector()
    {
        for (int i = 0; i < NOL; i++) // for every loads
        {
            int dofNode = 2, p[3] = {0};
            double llv[6] = {0}; // local load vector

            if (cstBuildUnitLoad(i, llv)) // build unit stiffness matrix
                return sfPrintError(7);

            p[0] = dofNode * (CSTriangles[Loads[i].num].nodes[0] - 1); // match the displacement with nods
            p[1] = dofNode * (CSTriangles[Loads[i].num].nodes[1] - 1);
            p[2] = dofNode * (CSTriangles[Loads[i].num].nodes[2] - 1);

            for (int j = 0; j < 2; j++) // add local load vector to load vector
                for (int m = 0; m < dofNode; m++)
                    LoadVector[p[j] + m] += llv[m];
        }

        return 0;
    }

    // print error
    bool sfPrintError(int error)
    {
        cout << "ERROR:\t";
        switch (error)
        {
        case 1:
            cout << "\n";
            break;
        case 2:
            cout << "\n";
            break;
        case 3:
            cout << "\n";
            break;
        case 4:
            cout << "\n";
            break;
        case 5:
            cout << "\n";
            break;
        case 6:
            cout << "\n";
            break;
        case 7:
            cout << "\n";
            break;
        case 8:
            cout << "\n";
            break;
        case 9:
            cout << "\n";
            break;
        case 10:
            cout << "\n";
            break;
        case 11:
            cout << "\n";
            break;
        case 12:
            cout << "\n";
            break;
        case 13:
            cout << "\n";
            break;
        case 14:
            cout << "\n";
            break;
        case 15:
            cout << "\n";
            break;
        case 16:
            cout << "\n";
            break;
        case 17:
            cout << "\n";
            break;
        case 18:
            cout << "\n";
            break;
        case 19:
            cout << "\n";
            break;
        case 20:
            cout << "\n";
            break;
        case 21:
            cout << "\n";
            break;
        case 22:
            cout << "\n";
            break;
        case 23:
            cout << "\n";
            break;
        case 24:
            cout << "\n";
            break;
        case 25:
            cout << "\n";
            break;

        default:
            break;
        }
        status = 3; //error

        return 1;
    }

    FiniteElement();
    ~FiniteElement();
};

FiniteElement::FiniteElement()
{
    Nodes = NULL;
    CSTriangles = NULL;
    Loads = NULL;
    LoadVector = NULL;
    Displacement = NULL;
}

FiniteElement::~FiniteElement()
{
    delete[] Nodes;
    Nodes = NULL;
    delete[] CSTriangles;
    CSTriangles = NULL;
    delete[] Loads;
    Loads = NULL;
    delete[] LoadVector;
    LoadVector = NULL;
    delete[] Displacement;
    Displacement = NULL;

    status = 0; // initialization is completed
}
