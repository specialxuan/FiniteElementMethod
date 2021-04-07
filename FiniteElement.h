#include "VarBandMatrix.h"
#include <Windows.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
using namespace std;

class FiniteElement
{
public:
    double EPS = 1e-15; // epsilon, calculation accuracy
    double MaxTS;       // maximum value in total stiffness matrix
    double MaxLV;       // maximum value in load vector

    int TNN; // total number of nodes
    int NFIN;
    int NFRN;
    int DOF;  // degree of freedom
    int NCST; // number of constant strain triangle unit
    int NOL;  // number of loads

    struct Node // parameters of nodes
    {
        double xcn; // X coordinate of nodes
        double ycn; // Y coordinate of nodes
        bool fixed; // fixed or not
    } * Nodes;      // parameters of nodes

    struct ConstantStrainTriangle // constant strain triangle unit
    {
        int nodes[3];     // 3 nodes of the unit
        double a[3];      // paramaters of shape function
        double b[3];      // paramaters of shape function
        double c[3];      // paramaters of shape function
        double area;      // paramaters of shape function
        double elastic;   // elastic module
        double mu;        // poisson ratio
        double thick;     // thickness of the unit
        double rou;       // density of the unit
        double strain[3]; // strain of the unit
        double stress[3]; // stress of the unit
    } * CSTriangles;      // constant strain triangle unit

    struct Load // parameters of loads
    {
        int unit;      // the number of unit with load
        int node[2];   // the node with load
        int kol;       // the kind of load
        double vol[2]; // the value of load
    } * Loads;         // parameters of loads

    VarBandMatrix TotalStiffness; // total stiffness matrix
    double *LoadVector;           // load vector
    double *Displacement;         // displacement of nodes

    bool ProgressBar = 1; // open progress bar
    bool Parallel;    // open parallel
    int status;

    // allocate total stiffness matrix, load vector and displacement vector
    bool feAllocate()
    {
        int *IV = new int[DOF](); // the location of diagonal element
        int it = 0;
        int bandwidth1 = 0, bandwidth2 = 0, bandwidth3 = 0;
        int *perband = new int[TNN](); // bandwidth per line in total stiffness matrix
        // int perband[20] = {0};

        for (int i = 0; i < NCST; i++)
        {
            // int n1 = CSTriangles[i].nodes[0] - 1;
            // int n2 = CSTriangles[i].nodes[1] - 1;
            // int n3 = CSTriangles[i].nodes[2] - 1;
            // bandwidth1 = n3 - n1;
            // bandwidth2 = n2 - n1;
            // bandwidth3 = n3 - n2;

            int n[3] = {0};
            n[0] = CSTriangles[i].nodes[0] - 1;
            n[1] = CSTriangles[i].nodes[1] - 1;
            n[2] = CSTriangles[i].nodes[2] - 1;

            sort(n, n + 3);

            bandwidth1 = n[2] - n[0];
            bandwidth2 = n[1] - n[0];
            bandwidth3 = n[2] - n[1];

            if (bandwidth1 > perband[n[2]] && n[0] > NFIN - 1)
                perband[n[2]] = bandwidth1;
            if (bandwidth2 > perband[n[1]] && n[0] > NFIN - 1)
                perband[n[1]] = bandwidth2;
            if (bandwidth3 > perband[n[2]] && n[1] > NFIN - 1)
                perband[n[2]] = bandwidth3;
        }

        // for (int i = 0; i < TNN; i++)
        //     cout << perband[i] << "\n";

        for (int i = NFIN; i < TNN; i++)
            for (int j = 1; j <= 2; j++)
            {
                it++;
                if (it == 1)
                    IV[it - 1] = 2 * perband[i] + j;
                else
                    IV[it - 1] = IV[it - 2] + 2 * perband[i] + j;
            }

        TotalStiffness.initialize(IV, DOF); // allocate memory for total stiffness matrix
        LoadVector = new double[DOF]();     // allocate memory for load vector
        Displacement = new double[DOF]();   // allocate memory for displacement vector
        // for (int i = 0; i < DOF; i++)
        // {
        //     cout << IV[i] << "\n";
        // }

        delete[] perband;
        delete[] IV;

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
            double x1 = Nodes[CSTriangles[k].nodes[0] - 1].xcn; // cooridnates of nodes of the unit
            double x2 = Nodes[CSTriangles[k].nodes[1] - 1].xcn;
            double x3 = Nodes[CSTriangles[k].nodes[2] - 1].xcn;
            double y1 = Nodes[CSTriangles[k].nodes[0] - 1].ycn;
            double y2 = Nodes[CSTriangles[k].nodes[1] - 1].ycn;
            double y3 = Nodes[CSTriangles[k].nodes[2] - 1].ycn;

            // cout << x1 << "," << y1 << " " << x2 << "," << y2 << " " << x3 << "," << y3 << "\n";

            a1 = x2 * y3 - x3 * y2; // paramaters of shape function
            a2 = x3 * y1 - x1 * y3;
            a3 = x1 * y2 - x2 * y1;
            b1 = y2 - y3;
            b2 = y3 - y1;
            b3 = y1 - y2;
            c1 = x3 - x2;
            c2 = x1 - x3;
            c3 = x2 - x1;

            // cout << a1 << " " << a2 << " " << a3 << "\n";
            // cout << b1 << " " << b2 << " " << b3 << "\n";
            // cout << c1 << " " << c2 << " " << c3 << "\n";

            CSTriangles[k].area = (a1 + a2 + a3) / 2; // area of the unit

            // cout << CSTriangles[k].area << " " << CSTriangles[k].elastic << " " << CSTriangles[k].mu << " " << CSTriangles[k].thick << " " << CSTriangles[k].rou << "\n";
        }

        return 0;
    }

    // build unit stiffness matrix of constant strain triangle unit
    bool cstBuildUnitStiff(int k, int r, int s, double *us) // k is unit number, r and s is section number, us is a section in unit stiffness matrix
    {
        if (k < 0 || r < 0 || r > 2 || s < 0 || s > 2 || us == NULL)
            return fePrintError(7);

        // double *a = CSTriangles[k].a; // paramaters of shape function
        double *b = CSTriangles[k].b;
        double *c = CSTriangles[k].c;

        double Area = CSTriangles[k].area; // paramaters of the unit
        double Mu = CSTriangles[k].mu;
        double E = CSTriangles[k].elastic;
        double T = CSTriangles[k].thick;

        double tmp = E * T/ (4 * Area * (1 - Mu * Mu)); // temperatary coefficient
        // cout << tmp << "   tmp   \n";
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
        int num = load.unit - 1;
        int n1 = load.node[0] - 1;
        int n2 = load.node[1] - 1;
        double A = CSTriangles[num].area;
        double t = CSTriangles[num].thick;
        double x1 = Nodes[CSTriangles[num].nodes[n1] - 1].xcn;
        double x2 = Nodes[CSTriangles[num].nodes[n2] - 1].xcn;
        double y1 = Nodes[CSTriangles[num].nodes[n1] - 1].ycn;
        double y2 = Nodes[CSTriangles[num].nodes[n2] - 1].ycn;
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

        return 0;
    }

    // calculate strain and stress of constant strain triangle unit
    bool cstStrainStress()
    {
        int p[3] = {0}; // p is a temperary vector for i0j0

        for (int i = 0; i < NCST; i++)
        {
            p[0] = 2 * (CSTriangles[i].nodes[0] - NFIN - 1); // match the displacement with nods
            p[1] = 2 * (CSTriangles[i].nodes[1] - NFIN - 1);
            p[2] = 2 * (CSTriangles[i].nodes[2] - NFIN - 1);

            double u1 = 0, v1 = 0, u2 = 0, v2 = 0, u3 = 0, v3 = 0;

            if (p[0] < 0)
                u1 = 0, v1 = 0;
            else
                u1 = Displacement[p[0] + 0], v1 = Displacement[p[0] + 1];
            if (p[1] < 0)
                u2 = 0, v2 = 0;
            else
                u2 = Displacement[p[1] + 0], v2 = Displacement[p[1] + 1];
            if (p[2] < 0)
                u3 = 0, v3 = 0;
            else
                u3 = Displacement[p[2] + 0], v3 = Displacement[p[2] + 1];

            // double u1 = Displacement[p[0] + 0]; // displacements
            // double v1 = Displacement[p[0] + 1];
            // double u2 = Displacement[p[1] + 0];
            // double v2 = Displacement[p[1] + 1];
            // double u3 = Displacement[p[2] + 0];
            // double v3 = Displacement[p[2] + 1];

            double b1 = CSTriangles[i].b[0]; // paramaters of shape function
            double b2 = CSTriangles[i].b[1];
            double b3 = CSTriangles[i].b[2];
            double c1 = CSTriangles[i].c[0];
            double c2 = CSTriangles[i].c[1];
            double c3 = CSTriangles[i].c[2];

            double area2 = 2 * CSTriangles[i].area;
            double E = CSTriangles[i].elastic;
            double mu = CSTriangles[i].mu;
            double tmp = E / (1 - mu * mu);

            // CSTriangles[i].strain[0] = u1 * b1 + u2 * b2 + u3 * b3 / area2; // strain
            // CSTriangles[i].strain[1] = v1 * c1 + v2 * c2 + v3 * c3 / area2;
            // CSTriangles[i].strain[2] = u1 * c1 + u2 * c2 + u3 * c3 + v1 * b1 + v2 * b2 + v3 * b3 / area2;

            double *strain = CSTriangles[i].strain;
            double *stress = CSTriangles[i].stress;

            strain[0] = (u1 * b1 + u2 * b2 + u3 * b3) / area2; // strain
            strain[1] = (v1 * c1 + v2 * c2 + v3 * c3) / area2;
            strain[2] = (u1 * c1 + u2 * c2 + u3 * c3 + v1 * b1 + v2 * b2 + v3 * b3) / area2;

            stress[0] = tmp * (strain[0] + mu * strain[1]); //stress
            stress[1] = tmp * (mu * strain[0] + strain[1]);
            stress[2] = tmp * (1 - mu) / 2 * strain[2];
        }

        return 0;
    }

    // build total stiffness matrix
    bool feBuildTotalStiff()
    {
        int dofNode = 2, numNode = 3;
        double us[dofNode * dofNode] = {0}; // unit stiffness matrix
        int p[3] = {0};                     // p is a temperary vector for i0j0

        for (int k = 0; k < NCST; k++) // for every nodes
        {
            p[0] = dofNode * (CSTriangles[k].nodes[0] - NFIN - 1); // match the displacement with nods
            p[1] = dofNode * (CSTriangles[k].nodes[1] - NFIN - 1);
            p[2] = dofNode * (CSTriangles[k].nodes[2] - NFIN - 1);

            for (int i = 0; i < numNode; i++) // for every nodes of the unit
            {
                for (int j = i; j < numNode; j++) // for every nodes of the unit
                {
                    if (p[i] >= 0 && p[j] >= 0)
                    {
                        if (cstBuildUnitStiff(k, i, j, us)) // build unit stiffness matrix
                            return fePrintError(8);
                        
                        // for (int k = 0; k < 4; k++)
                        //     cout << setw(10) << us[k] << " | ";
                        // cout << "\n";

                        for (int m = 0; m < dofNode; m++)
                            for (int n = 0; n < dofNode; n++)
                                TotalStiffness(p[i] + m, p[j] + n) += us[m * dofNode + n]; // superpose
                    }
                }
            }
        }

        MaxTS = TotalStiffness.normalize();

        // for (int i = 0; i < TNN; i++)
        //     if (Nodes[i].fixed)
        //         TotalStiffness(i * 2, i * 2) = TotalStiffness(i * 2 + 1, i * 2 + 1) = 1;

        return 0;
    }

    // build load vector
    bool feBuildLoadVector()
    {
        for (int i = 0; i < NOL; i++) // for every loads
        {
            int dofNode = 2, p[3] = {0};
            double llv[6] = {0}; // local load vector

            if (cstBuildUnitLoad(i, llv)) // build unit stiffness matrix
                return fePrintError(9);

            p[0] = dofNode * (CSTriangles[Loads[i].unit - 1].nodes[0] - NFIN - 1); // match the displacement with nods
            p[1] = dofNode * (CSTriangles[Loads[i].unit - 1].nodes[1] - NFIN - 1);
            p[2] = dofNode * (CSTriangles[Loads[i].unit - 1].nodes[2] - NFIN - 1);

            for (int j = 0; j < 3; j++) // add local load vector to load vector
                for (int m = 0; m < dofNode; m++)
                    if (p[j] >= 0)
                        LoadVector[p[j] + m] += llv[j * 2 + m];
        }

        for (int i = 0; i < DOF; i++)
            if (fabs(LoadVector[i]) > MaxLV)
                MaxLV = fabs(LoadVector[i]);

        for (int i = 0; i < DOF; i++)
            LoadVector[i] = LoadVector[i] / MaxLV;

        return 0;
    }

    bool feConjugateGradient(VarBandMatrix &A, double *b, double *x, int N)
    {
        if (b == NULL || x == NULL || N == 0)
            return fePrintError(12);

        double *r = NULL, *p = NULL, *z = NULL;
        double gamma = 0, gamma_new = 0, gamma_new_sqrt = 0, alpha = 0, beta = 0;
        int percent = 0, percent_new = 0;

        if (ProgressBar)
            cout << "\rSolving equation      [ 0%% ][                                                 ]";

        r = (double *)malloc(N * sizeof(double));
        memset(r, 0, sizeof(double));
        p = (double *)malloc(N * sizeof(double));
        memset(p, 0, sizeof(double));
        z = (double *)malloc(N * sizeof(double));
        memset(z, 0, sizeof(double));

        // TODO: Max
        // for (int i = 0; i < NSI; i++)
        //     A[i] = A[i] / MaxTS;
        // for (int i = 0; i < N; i++)
        //     b[i] = b[i] / MaxLV;

        // x = [0 ... 0]
        // r = b - A * x
        // p = r
        // gamma = r' * r
        gamma = 0.0;
        for (int i = 0; i < N; ++i)
        {
            x[i] = 0.0;
            r[i] = b[i];
            p[i] = r[i];
            gamma += r[i] * r[i];
        }

        for (int n = 0; 1; ++n)
        {
            // z = A * p
            for (int i = 0; i < N; i++)
            {
                z[i] = 0.0;
                for (int j = 0; j < N; j++)
                    z[i] += A(i, j) * p[j];
            }

            // alpha = gamma / (p' * z)
            alpha = 0.0;
            for (int i = 0; i < N; ++i)
                alpha += p[i] * z[i];
            alpha = gamma / alpha;

            // x = x + alpha * p
            // r = r - alpha * z
            // gamma_new = r' * r
            gamma_new = 0.0;
            for (int i = 0; i < N; ++i)
            {
                x[i] += alpha * p[i];
                r[i] -= alpha * z[i];
                gamma_new += r[i] * r[i];
            }

            gamma_new_sqrt = sqrt(gamma_new);
            if (gamma_new_sqrt < EPS)
                break;

            if (ProgressBar)
            {
                percent_new = (int)((1 - log10(gamma_new_sqrt / EPS) / (-log10(EPS))) * 100);
                if (percent_new > percent)
                {
                    percent = percent_new;
                    cout << "\rSolving equation ";
                    for (int i = 0; i <= 4; i++)
                        if (i <= n % 4)
                            cout << ".";
                        else
                            cout << " ";
                    cout << "[ " << percent << "%% ]";
                    cout << "[";
                    for (int i = 0; i < 49; i++)
                        if (i < percent / 2)
                            cout << "=";
                        else
                            cout << " ";
                    cout << "]";
                }
                else
                {
                    cout << "\rSolving equation ";
                    for (int i = 0; i <= 4; i++)
                        if (i <= n % 4)
                            cout << ".";
                        else
                            cout << " ";
                }
            }

            beta = gamma_new / gamma;

            // p = r + (gamma_new / gamma) * p;
            for (int i = 0; i < N; ++i)
                p[i] = r[i] + beta * p[i];

            // gamma = gamma_new
            gamma = gamma_new;
        }

        // TODO: Max
        // for (int i = 0; i < NSI; i++)
        //     A[i] = A[i] * MaxTS;
        // for (int i = 0; i < N; i++)
        //     b[i] = b[i] * MaxLV;
        for (int i = 0; i < N; i++)
            x[i] = x[i] * MaxLV / MaxTS;

        if (ProgressBar)
            cout << "\rSolving equation done [ 100%% ][=================================================]\n";

        free(r);
        free(p);
        free(z);

        return 0;
    }

    // solve equation of matrix by conjugate gradient parallel
    bool feConjugateGradientPar(VarBandMatrix &A, double *b, double *x, int N)
    {
        if (b == NULL || x == NULL || N == 0)
            return fePrintError(10);

        double *r = NULL, *p = NULL, *z = NULL;
        double gamma = 0, gamma_new = 0, gamma_new_sqrt = 0, alpha = 0, beta = 0;
        int percent = 0, percent_new = 0;

        if (ProgressBar)
            cout << "\rSolving equation      [ 0% ][                                                 ]";

        r = (double *)malloc(N * sizeof(double));
        memset(r, 0, sizeof(double));
        p = (double *)malloc(N * sizeof(double));
        memset(p, 0, sizeof(double));
        z = (double *)malloc(N * sizeof(double));
        memset(z, 0, sizeof(double));

        // for (int i = 0; i < NSI; i++)
        //     A[i] = A[i] / MaxTS;
        // for (int i = 0; i < N; i++)
        //     b[i] = b[i] / MaxLV;

        //  x = [0 ... 0]
        //  r = b - A * x
        //  p = r
        //  gamma = r' * r
        gamma = 0.0;
#pragma omp parallel for reduction(+ \
                                   : gamma)
        for (int i = 0; i < N; ++i)
        {
            x[i] = 0.0;
            r[i] = b[i];
            p[i] = r[i];
            gamma += r[i] * r[i];
        }

        for (int n = 0; true; ++n)
        {
//  z = A * p
#pragma omp parallel for
            for (int i = 0; i < N; i++)
            {
                z[i] = 0.0;
                for (int j = 0; j < N; j++)
                {
                    if (i == j)
                    {
                        z[i] += A.matrix[A.diagelement[i] - 1] * p[j];
                    }
                    else if (j > i)
                    {
                        if ((A.diagelement[j] - j + i) > A.diagelement[j - 1])
                            z[i] += A.matrix[A.diagelement[j] - j + i - 1] * p[j];
                        else
                            z[i] += 0;
                    }
                    else if (i > j)
                    {
                        if ((A.diagelement[i] - i + j) > A.diagelement[i - 1])
                            z[i] += A.matrix[A.diagelement[i] - i + j - 1] * p[j];
                        else
                            z[i] += 0;
                    }
                }
            }

            //  alpha = gamma / (p' * z)
            alpha = 0.0;
#pragma omp parallel for reduction(+ \
                                   : alpha)
            for (int i = 0; i < N; ++i)
                alpha += p[i] * z[i];
            alpha = gamma / alpha;

            //  x = x + alpha * p
            //  r = r - alpha * z
            //  gamma_new = r' * r
            gamma_new = 0.0;
#pragma omp parallel for reduction(+ \
                                   : gamma_new)
            for (int i = 0; i < N; ++i)
            {
                x[i] += alpha * p[i];
                r[i] -= alpha * z[i];
                gamma_new += r[i] * r[i];
            }

            gamma_new_sqrt = sqrt(gamma_new);
            if (gamma_new_sqrt < EPS)
                break;

            if (ProgressBar)
            {
                percent_new = (int)((1 - log10(gamma_new_sqrt / EPS) / (-log10(EPS))) * 100);
                if (percent_new > percent)
                {
                    percent = percent_new;
                    cout << "\rSolving equation ";
                    for (int i = 0; i <= 4; i++)
                        if (i <= n % 4)
                            cout << ".";
                        else
                            cout << " ";
                    cout << "[ " << percent << "% ]";
                    cout << "[";
                    for (int i = 0; i < 49; i++)
                        if (i < percent / 2)
                            cout << "=";
                        else
                            cout << " ";
                    cout << "]";
                }
                else
                {
                    cout << "\rSolving equation ";
                    for (int i = 0; i <= 4; i++)
                        if (i <= n % 4)
                            cout << ".";
                        else
                            cout << " ";
                }
            }

            beta = gamma_new / gamma;

//  p = r + (gamma_new / gamma) * p;
#pragma omp parallel for
            for (int i = 0; i < N; ++i)
                p[i] = r[i] + beta * p[i];

            //  gamma = gamma_new
            gamma = gamma_new;
        }

        // for (int i = 0; i < NSI; i++)
        //     A[i] = A[i] * MaxTS;
        // for (int i = 0; i < N; i++)
        //     b[i] = b[i] * MaxLV;
        for (int i = 0; i < N; i++)
            x[i] = x[i] * MaxLV / MaxTS;

        if (ProgressBar)
            cout << "\rSolving equation done [ 100% ][=================================================]\n";

        free(r);
        free(p);
        free(z);

        return 0;
    }

    // print error
    bool fePrintError(int error)
    {
        cout << "ERROR:\t";
        switch (error)
        {
        case 1:
            cout << "Initializing const strain triangle unit failed!\n";
            break;
        case 2:
            cout << "Allocating memeories failed!\n";
            break;
        case 3:
            cout << "Building total stiffness matrix failed!\n";
            break;
        case 4:
            cout << "Building load vector failed!\n";
            break;
        case 5:
            cout << "Solving equation failed!\n";
            break;
        case 6:
            cout << "Calculating strain and stress failed!\n";
            break;
        case 7:
            cout << "There is something wrong in building unit stiffness matrix!\n";
            break;
        case 8:
            cout << "There is something wrong in building total stiffness matrix!\n";
            break;
        case 9:
            cout << "There is something wrong in building load vector!\n";
            break;
        case 10:
            cout << "There is something wrong in the equation!\n";
            break;
        case 11:
            cout << "There is no such file!\n";
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

    bool fePrintError(int row, int column)
    {
        if (column == 1)
            cout << "ERROR:\tRow: " << row << " Column: 1 : head is mismathced!\n";
        else
            cout << "ERROR:\tRow: " << row << " Column: " << column << " : data input failed!\n";
        status = 3; // error

        return 1;
    }

    FiniteElement();
    ~FiniteElement();

    bool feInput(const char *);
    bool feOutput(const char *);
    bool feCalculate();
    bool feCircularStructure(int, int);
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

bool FiniteElement::feInput(const char *inputFile = "source&result/fe.csv")
{
    if (status)
        this->~FiniteElement();
    const int one = 1;
    struct Row
    {
        char head[10];
        const int &cnt;
    } rows[20] = {
        {"TNN", one},
        {"NFIN", one},
        {"NCST", one},
        {"NOL", one},
        {"XCN", TNN},
        {"YCN", TNN},
        {"FIX", TNN},
        {"NODE1", NCST},
        {"NODE2", NCST},
        {"NODE3", NCST},
        {"ELASTIC", NCST},
        {"MU", NCST},
        {"THICK", NCST},
        {"ROU", NCST},
        {"UNIT", NOL},
        {"NODE1", NOL},
        {"NODE2", NOL},
        {"KOL", NOL},
        {"VOLX", NOL},
        {"VOLY", NOL},
    };

    int rowIndex = 0;   // Reset the number of rows to zero
    char buf[10] = {0}; // buffer and data string

    ifstream fin(inputFile, ios::in);
    if (!fin)
        return fePrintError(11);

    rowIndex = 1;
    fin.ignore(1000000, '\n'); // skip first line

    rowIndex = 2;
    fin.getline(buf, 10, ',');
    if (strcmp(rows[rowIndex - 2].head, buf))
        return fePrintError(rowIndex, 1);
    for (int i = 0; i < rows[rowIndex - 2].cnt; i++)
        if (!(fin >> TNN))
            return fePrintError(rowIndex, i + 1);
    fin.ignore(1000000, '\n');

    rowIndex = 3;
    fin.getline(buf, 10, ',');
    if (strcmp(rows[rowIndex - 2].head, buf))
        return fePrintError(rowIndex, 1);
    for (int i = 0; i < rows[rowIndex - 2].cnt; i++)
        if (!(fin >> NFIN))
            return fePrintError(rowIndex, i + 1);
    fin.ignore(1000000, '\n');

    rowIndex = 4;
    fin.getline(buf, 10, ',');
    if (strcmp(rows[rowIndex - 2].head, buf))
        return fePrintError(rowIndex, 1);
    for (int i = 0; i < rows[rowIndex - 2].cnt; i++)
        if (!(fin >> NCST))
            return fePrintError(rowIndex, i + 1);
    fin.ignore(1000000, '\n');

    rowIndex = 5;
    fin.getline(buf, 10, ',');
    if (strcmp(rows[rowIndex - 2].head, buf))
        return fePrintError(rowIndex, 1);
    for (int i = 0; i < rows[rowIndex - 2].cnt; i++)
        if (!(fin >> NOL))
            return fePrintError(rowIndex, i + 1);
    fin.ignore(1000000, '\n');

    NFRN = TNN - NFIN;
    DOF = 2 * NFRN;
    Nodes = new Node[TNN]();
    CSTriangles = new ConstantStrainTriangle[NCST]();
    Loads = new Load[NOL]();

    for (rowIndex = 6; rowIndex <= 21; rowIndex++)
    {
        fin.getline(buf, 10, ',');
        if (strcmp(rows[rowIndex - 2].head, buf))
            return fePrintError(rowIndex, 1);
        for (int i = 0; i < rows[rowIndex - 2].cnt; i++)
        {
            switch (rowIndex)
            {
            case 6:
                if (!(fin >> Nodes[i].xcn))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 7:
                if (!(fin >> Nodes[i].ycn))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 8:
                if (!(fin >> Nodes[i].fixed))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 9:
                if (!(fin >> CSTriangles[i].nodes[0]))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 10:
                if (!(fin >> CSTriangles[i].nodes[1]))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 11:
                if (!(fin >> CSTriangles[i].nodes[2]))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 12:
                if (!(fin >> CSTriangles[i].elastic))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 13:
                if (!(fin >> CSTriangles[i].mu))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 14:
                if (!(fin >> CSTriangles[i].thick))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 15:
                if (!(fin >> CSTriangles[i].rou))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 16:
                if (!(fin >> Loads[i].unit))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 17:
                if (!(fin >> Loads[i].node[0]))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 18:
                if (!(fin >> Loads[i].node[1]))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 19:
                if (!(fin >> Loads[i].kol))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 20:
                if (!(fin >> Loads[i].vol[0]))
                    return fePrintError(rowIndex, i + 1);
                break;
            case 21:
                if (!(fin >> Loads[i].vol[1]))
                    return fePrintError(rowIndex, i + 1);
                break;
            }
            fin.get();
        }
        fin.ignore(1000000, '\n');
    }
    fin.close();
    status = 1; // input procedure is completed

    return 0;

    // EPS = 1e-15;

    // TNN = 4;
    // DOF = 2 * TNN;
    // NCST = 2;
    // NOL = 1;

    // Nodes = new Node[TNN]();
    // CSTriangles = new ConstantStrainTriangle[NCST]();
    // Loads = new Load[NOL]();

    // Nodes[0].xcn = 0;
    // Nodes[0].ycn = 0;
    // Nodes[1].xcn = 2;
    // Nodes[1].ycn = 0;
    // Nodes[2].xcn = 2;
    // Nodes[2].ycn = 1;
    // Nodes[3].xcn = 0;
    // Nodes[3].ycn = 1;

    // Nodes[0].fixed = 1;
    // Nodes[1].fixed = 0;
    // Nodes[2].fixed = 0;
    // Nodes[3].fixed = 1;

    // CSTriangles[0].nodes[0] = 1;
    // CSTriangles[0].nodes[1] = 2;
    // CSTriangles[0].nodes[2] = 4;
    // CSTriangles[0].elastic = 1;
    // CSTriangles[0].mu = 1.0 / 3;
    // CSTriangles[0].thick = 1;

    // CSTriangles[1].nodes[0] = 2;
    // CSTriangles[1].nodes[1] = 3;
    // CSTriangles[1].nodes[2] = 4;
    // CSTriangles[1].elastic = 1;
    // CSTriangles[1].mu = 1.0 / 3;
    // CSTriangles[1].thick = 1;

    // Loads[0].kol = 0;
    // Loads[0].node[0] = 2;
    // Loads[0].node[1] = 1;
    // Loads[0].unit = 1;
    // Loads[0].vol[0] = 0;
    // Loads[0].vol[1] = 1;

    // return 0;
}

bool FiniteElement::feOutput(const char *outputFile = "source&result/feResult.csv")
{
    // if (status != 2)
    // {
    //     cout << "ERROR:\tCalculation is not completed!\n";
    //     return 0;
    // }

    ofstream fout(outputFile, ios::out);
    fout << setw(66) << "Calculation Of Finite Element Method,\n";
    fout << "Unit   , EPSILON_X       , EPSILON_Y       , GAMMA_XY        , SIGMA_X         , SIGMA_Y         , TAU_XY          ,\n";
    for (int i = 0; i < NCST; i++)
        fout << setw(6) << i + 1 << " , "
             << setw(15) << CSTriangles[i].strain[0] << " , "
             << setw(15) << CSTriangles[i].strain[1] << " , "
             << setw(15) << CSTriangles[i].strain[2] << " , "
             << setw(15) << CSTriangles[i].stress[0] << " , "
             << setw(15) << CSTriangles[i].stress[1] << " , "
             << setw(15) << CSTriangles[i].stress[2] << " ,\n";

    return 0;
}

bool FiniteElement::feCalculate()
{
    if (status == 0 || status == 3)
        return fePrintError(0);

    if (cstInitialize())
        return fePrintError(1);
    else
        cout << "Initializing const strain triangle unit succeed!\n";

    if (feAllocate())
        return fePrintError(2);
    else
        cout << "Allocating memories succeed!\n";

    if (feBuildTotalStiff())
        return fePrintError(3);
    else
        cout << "Building total stiffness matrix succeed!\n";

    if (feBuildLoadVector())
        return fePrintError(4);
    else
        cout << "Building load vector succeed!\n";

    if (feConjugateGradientPar(TotalStiffness, LoadVector, Displacement, DOF))
        return fePrintError(5);
    else
        cout << "Solving equation succeed!\n";

    if (cstStrainStress())
        return fePrintError(6);
    else
        cout << "Calculating strain and stress succeed!\n";

    status = 2;

    return 0;
}

bool FiniteElement::feCircularStructure(int m, int n)
{
    int NCST = (m - 1) * (n - 1) * 2;
    ofstream fout("source&result/fe_test.csv", ios::out);

    fout << "INPUT, Degree of freedom is " << 2 * m * n << ",\n";
    fout << "TNN," << m * n << ",\n";
    fout << "NFIN," << m << ",\n";
    fout << "NCST," << NCST << ",\n";
    fout << "NOL," << m - 1 << ",\n";

    double delta_x = 10.0 / (n - 1), delta_y = 1.0 / (m - 1);

    fout << "XCN,";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            fout << i * delta_x << ",";
    fout << "\n";

    fout << "YCN,";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            fout << j * delta_y << ",";
    fout << "\n";

    fout << "FIX,";
    for (int i = 0; i < m; i++)
        fout << "1,";
    for (int i = m; i < m * n; i++)
        fout << "0,";
    fout << "\n";

    fout << "NODE1,";
    for (int i = 0; i < n - 1; i++)
        for (int j = 0; j < m - 1; j++)
            fout << i * m + j + 1 << "," << i * m + j + 1 << ",";
    fout << "\n";

    fout << "NODE2,";
    for (int i = 0; i < n - 1; i++)
        for (int j = 0; j < m - 1; j++)
            fout << i * m + j + 1 + m + 1 << "," << i * m + j + 1 + m << ",";
    fout << "\n";

    fout << "NODE3,";
    for (int i = 0; i < n - 1; i++)
        for (int j = 0; j < m - 1; j++)
            fout << i * m + j + 1 + 1 << "," << i * m + j + 1 + m + 1<< ",";
    fout << "\n";

    fout << "ELASTIC,";
    for (int i = 0; i < NCST; i++)
        fout << "1,";
    fout << "\n";

    fout << "MU,";
    for (int i = 0; i < NCST; i++)
        fout << 1.0 / 3 << ",";
    fout << "\n";

    fout << "THICK,";
    for (int i = 0; i < NCST; i++)
        fout << "1,";
    fout << "\n";

    fout << "ROU,";
    for (int i = 0; i < NCST; i++)
        fout << "1,";
    fout << "\n";

    fout << "UNIT,";
    for (int j = 0; j < m - 1; j++)
        fout << ((n - 2) * (m - 1) + j + 1) * 2 << ",";
    fout << "\n";

    fout << "NODE1,";
    for (int i = 0; i < m - 1; i++)
        fout << "2,";
    fout << "\n";

    fout << "NODE2,";
    for (int i = 0; i < m - 1; i++)
        fout << "3,";
    fout << "\n";

    fout << "KOL,";
    for (int i = 0; i < m - 1; i++)
        fout << "2,";
    fout << "\n";

    fout << "VOLX,";
    for (int i = 0; i < m - 1; i++)
        fout << "0,";
    fout << "\n";

    fout << "VOLY,";
    for (int i = 0; i < m - 1; i++)
        fout << "1,";
    fout << "\n";

    fout << "END,\n";
    return 0;
}