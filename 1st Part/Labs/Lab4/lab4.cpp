#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "Prof_lab4.cpp"

#define SIZE 11

using namespace std;

/* Auxiliary methods */
double Trapz_Scalar_Prod(double f1[SIZE], double f2[SIZE], double dx)
{
    double sum = 0;
    for (int n = 0; n < (SIZE - 1); n++) sum += ((f1[n+1] * f2[n+1]) + (f1[n] * f2[n]))/2 * dx;
    return(sum);
}

void Write_to_file(vector<vector<double>> &Res, string filename)
{
    cout << "Entering \"Write_to_file\"... " << endl;
    cout << "Matrix dimensions: " << size(Res) << " " << size(Res[0]) << endl;
    fstream w_File;
    w_File.open(filename, fstream::out);
    for (int i = 0; i < size(Res); i++)
    {     
        for (int j = 0; j < size(Res[0]); j++)
        {
        if(j < (size(Res[0]) - 1)) w_File << Res[i][j] << ",";
        else                       w_File << Res[i][j] << endl;
        }
    }
    w_File.close();
}

void Write_to_file(vector<double[SIZE]> &Res, string filename)
{
    cout << "Entering \"Write_to_file\"... " << endl;
    cout << "Matrix dimensions: " << size(Res) << " " << SIZE << endl;
    fstream w_File;
    w_File.open(filename, fstream::out);
    for (int i = 0; i < size(Res); i++)
    {     
        for (int j = 0; j < SIZE; j++)
        {
        if(j < (SIZE - 1)) w_File << Res[i][j] << ",";
        else                       w_File << Res[i][j] << endl;
        }
    }
    w_File.close();
}

int main(int argc, const char** argv) {
    /* Variables: 1D case */
    const double dx   = 1./double(SIZE - 1.);
    const double T    = 1.;
    const double tau  = 0.01;
    /* Right-hand side */
    double d[SIZE] = {0.};
    double u[SIZE] = {0.};
    std::vector<std::vector<double>> All;

    for (size_t i = 0; i < SIZE; i++)
    {
        d[i] = 1.;
        cout << "d[" << i << "] = " << d[i] << endl;
    }
    d[0] = 0.; d[SIZE - 1] = 0.; /* Boundary conditions */
    double S0 = sqrt(Trapz_Scalar_Prod(d, d, dx));
    for (size_t i = 0; i < SIZE; i++) /* Normalization */
        d[i] /= S0;
    std::vector<double> Tmp1;
    Tmp1.assign(d, d + SIZE);
    All.push_back(Tmp1);
    Tmp1.clear();

    /* Tridiagonal matrix */
    double a[SIZE] = {0.};
    double b[SIZE] = {0.};
    double c[SIZE] = {0.};
 
    for (size_t i = 0; i < SIZE; i++)
    {
        a[i]     = -tau/(dx*dx);
        b[i]     = 1 + 2*tau/(dx*dx); 
        c[i]     = -tau/(dx*dx);
    }

    /* Initial step */
    thomasAlgorithm(a, b, c, d, u, SIZE);
    for (size_t i = 0; i < SIZE; i++)
    {
        a[i]     = -tau/(dx*dx);
        b[i]     = 1 + 2*tau/(dx*dx); 
        c[i]     = -tau/(dx*dx);
    }
    double S = sqrt(Trapz_Scalar_Prod(u, u, dx));
    for (size_t i = 0; i < SIZE; i++) /* Normalization */
        u[i] /= S;
    std::vector<double> Tmp;
    Tmp.assign(u, u + SIZE);
    All.push_back(Tmp);
    Tmp.clear();

    for (size_t i = 0; i < int(T/tau); i++)
    {
        memcpy(d, u, sizeof(d));
        thomasAlgorithm(a, b, c, d, u, SIZE);
        for (size_t i = 0; i < SIZE; i++)
        {
            a[i]     = -tau/(dx*dx);
            b[i]     = 1 + 2*tau/(dx*dx); 
            c[i]     = -tau/(dx*dx);
        }
        double S = sqrt(Trapz_Scalar_Prod(u, u, dx));
        for (size_t i = 0; i < SIZE; i++) /* Normalization */
            u[i] /= S;
    std::vector<double> Tmp;
    Tmp.assign(u, u + SIZE);
    All.push_back(Tmp);
    Tmp.clear();
    }
    
    cout << endl;
    for (size_t i = 0; i < SIZE; i++)
        cout << "u[" << i << "] = " << u[i] << endl;
    Write_to_file(All, "Lab4.txt");

    return 0;
}