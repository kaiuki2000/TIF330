#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

/* Auxiliary methods */
double Trapz_Scalar_Prod(vector<double> &f1, vector<double> &f2, double dx)
{
    double sum = 0;
    for (int n = 0; n < (size(f1) - 1); n++) sum += ((f1[n+1] * f2[n+1]) + (f1[n] * f2[n]))/2 * dx;
    return(sum);
}

void Write_Vec_to_file(vector<double> &Vec, string filename)
{
    cout << "Entering \"Write_Vec_to_file\"... " << endl;
    fstream w_File;
    w_File.open(filename, fstream::out);
    for (int i = 0; i < size(Vec); i++)
    {     
        if(i < (size(Vec) - 1)) w_File << Vec[i] << ",";
        else                    w_File << Vec[i] << endl;
    }
}

int main()
{
    /* Variables */
    const int    N    = 50; /* Number of basis functions */
    const int    SIZE = 2049;
    const double dx   = 2/double(SIZE - 1);
    const double a    = M_PI/2.;
    int          n    = 0;

    vector<double> x(SIZE, 0.);
    vector<double> y(SIZE, 0.);
    vector<double> RHS(SIZE, 0.);

    vector<double> f(N, 0.);
    vector<double> w(N, 0.);

    vector<double>         Tmp(SIZE, 0.);
    vector<vector<double>> Chi(N, Tmp);
    vector<double>         Ly(SIZE, 0.);
    
    /* Filling up x vector */
    for (auto &val : x)
    {
        val += -1. + double(n)*dx;
        n++;
    }

    /* Filling up RHS vector: Initial guess - parabola */
    for (size_t s = 0; s < SIZE; s++) RHS[s] = sqrt(1. + 4.*a*x[s]*x[s]);
    
    /* Filling up Chi matrix */
    for (size_t k = 0; k < N; k++)
    {
        for (size_t s = 0; s < SIZE; s++) /* 's' is a spatial index */
        {
            Chi[k][s] = cos((2.*k + 1.) * (M_PI/2.) * x[s]);
        }
    }

    /* Galerkin method loop */
    for (size_t it = 0; it < 20; it++) /* 20 iterations */
    {
        /* Resetting from previous iteration */
        fill(y.begin(), y.end(), 0.);
        fill(Ly.begin(), Ly.end(), 0.);

        for (size_t k = 0; k < N; k++)
        {
            f[k] = Trapz_Scalar_Prod(Chi[k], RHS, dx);
            w[k] = - (f[k])/(pow((2.*k + 1.) * M_PI/2., 2.));

            for (size_t s = 0; s < SIZE; s++) /* 's' is a spatial index */
            {
                /* Important: 'y' has to always start at ALL 0. */
                y[s] += w[k] * Chi[k][s];
            }
        }

        /* New RHS */
        for (size_t k = 0; k < N; k++)
        {
            for (size_t s = 0; s < SIZE; s++) /* 's' is a spatial index */
            {
                /* Important: 'Ly' has to always start at ALL 0. */
                Ly[s] += w[k] * (-1.) * ((2.*k + 1.) * M_PI/2.) * sin((2.*k + 1.) * (M_PI/2.) * x[s]);
            }
        }
        for (size_t s = 0; s < SIZE; s++)
        {
            Ly[s] = Ly[s] * Ly[s];
            RHS[s] = sqrt(1. + a*Ly[s]);
        }
    }
    Write_Vec_to_file(y, "yVec.txt");
    Write_Vec_to_file(x, "xVec.txt");
    cout << "Result [a = M_PI/2]: x = " << x[SIZE/2] << "; y = " << y[SIZE/2] << "." << endl;
    
    return(0);
}