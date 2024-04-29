#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <string>
#include <fstream>

/* Potential: Independent of time! */
double V(double x, double y){
    // return(-5. * pow((1. + pow(x/5., 2.) + pow(y/4., 2.)), -4.));
    return(0.); // Trying particle in a 2D box.
}
/* Primitive normalization (Riemann sum type). Ideally, should implement more sophisticated 2D integration techniques. */
template <size_t gridX, size_t gridY>
double Get_Norm(std::complex<double> (&grid)[gridX][gridY], double dx, double dy){
    double norm = 0;
    for (size_t m = 0; m < gridX - 0; m++)
        for (size_t n = 0; n < gridY - 0; n++)
            norm += (std::abs(grid[m][n])*std::abs(grid[m][n])) * (dx*dy);
    return(norm);
}

template <size_t gridX, size_t gridY>
void Update_psiGrid(std::complex<double> (&grid)[gridX][gridY], double dx, double dy, double dkx, double dky, double tMax, double dt, int ind, std::vector<std::complex<double>> &psiXY){
    fftw_plan p, pR;
    p   = fftw_plan_dft_2d(gridX, gridY, reinterpret_cast<fftw_complex*>(&grid[0]), reinterpret_cast<fftw_complex*>(&grid[0]),  FFTW_FORWARD, FFTW_ESTIMATE);
    pR  = fftw_plan_dft_2d(gridX, gridY, reinterpret_cast<fftw_complex*>(&grid[0]), reinterpret_cast<fftw_complex*>(&grid[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
    const std::complex<double> imU(0., 1.); /* Imaginary unit*/
    std::vector<double> kxVec(gridX, 0.), kyVec(gridY, 0.);

    /* Defining kVec's for later */
    for (size_t i = 0; i < gridX; i++)
    {
        kxVec[i] = (i < gridX/2 + 1) ? dkx * i : dkx * (i - gridX);
        kyVec[i] = (i < gridX/2 + 1) ? dky * i : dky * (i - gridX);
    }

    /* Time loop */
    for (size_t tIndex = 0; tIndex < int(tMax/dt); tIndex++)
    {
        /* Applying exp(-i*t/2*V(0)): V = 0, so we don't need to do anything. */

        /* Applying 2D FFT: Check the other script for help */
        fftw_execute(p);

        // Enforcing Psi = 0 at boundaries. For the particle in a 2D box case.
        for (size_t m = 0; m < gridX; m++)
        {
            for (size_t n = 0; n < gridY; n++)
            {
                if(m == 0)       {grid[m][n] = 0;};
                if(m == gridX-1) {grid[m][n] = 0;};
                if(n == 0)       {grid[m][n] = 0;};
                if(n == gridY-1) {grid[m][n] = 0;};
            }
        }

        /* There's some room for improvement here. No need to reassign everything to 'grid' right away. Can use 'out' */
        /* Applying exp(-i*k^2*t) */
        for (size_t i = 0; i < gridX; i++)
            for (size_t j = 0; j < gridY; j++)
                grid[i][j] *= std::exp(-1. * imU * (kxVec[i]*kxVec[i] + kyVec[j]*kyVec[j]) * dt);

        /* Applying Inverse 2D FFT: Check the other script for help */
        fftw_execute(pR);

        for (size_t i = 0; i < gridX; i++)
            for (size_t j = 0; j < gridY; j++)
                grid[i][j] /= double(gridX*gridY); /* Normalization */

        // Enforcing Psi = 0 at boundaries. For the particle in a 2D box case.
        for (size_t m = 0; m < gridX; m++)
        {
            for (size_t n = 0; n < gridY; n++)
            {
                if(m == 0)       {grid[m][n] = 0;};
                if(m == gridX-1) {grid[m][n] = 0;};
                if(n == 0)       {grid[m][n] = 0;};
                if(n == gridY-1) {grid[m][n] = 0;};
            }
        }

        /* Finally, applying exp(-i*t/2*V(t)); Same as before, since V != V(t): V = 0, so we don't need to do anything. */

        double norm = Get_Norm(grid, dx, dy); /* Normalizing */
        for (size_t m = 0; m < gridX; m++)
            for (size_t n = 0; n < gridY; n++)
                grid[m][n] /= sqrt(norm);

        std::cout << "Iteration " << tIndex << " of " << int(tMax/dt) - 1 << "... (Norm_Check) = " << Get_Norm(grid, dx, dy) << std::endl;
        psiXY.push_back(grid[ind][gridY/2]);
    }
    fftw_destroy_plan(p);
    fftw_destroy_plan(pR);
}

void Write_Vec_to_file(std::vector<std::complex<double>> &Vec, std::string filename)
{
    using namespace std;
    cout << "Entering \"Write_Vec_to_file\"... " << endl;
    fstream w_File;
    w_File.open(filename, fstream::out);
    for (int i = 0; i < size(Vec); i++)
    {     
        if(i < (size(Vec) - 1)) w_File << real(Vec[i]) << "," << imag(Vec[i]) << endl;
        else                    w_File << real(Vec[i]) << "," << imag(Vec[i]);
    }
}

void Write_k_Vec_to_file(std::vector<double> &Vec, std::string filename)
{
    using namespace std;
    cout << "Entering \"Write_T_Vec_to_file\"... " << endl;
    fstream w_File;
    w_File.open(filename, fstream::out);
    for (int i = 0; i < size(Vec); i++)
    {     
        if(i < (size(Vec) - 1)) w_File << Vec[i] << endl;
        else                    w_File << Vec[i];
    }
}

/* Reason why the norm check isn't = 0: Due to my deefinition of the norm! (Primitive integral...
Doesn't consider the extreme points. ) */
int main(int argc, const char** argv) {
    const int gridX    = 201; /* How should I choose these values? 2^i, 3^i, etc. ?*/
    const int gridY    = 201;
    const double dx    = 1./(gridX - 1);
    const double dy    = 1./(gridY - 1);
    const double dkx   = 2.*M_PI/1.;
    const double dky   = 2.*M_PI/1.;
    const double tMax  = 100.;
    const double dt    = 0.15;
    double       norm  = 0.;
    int          ind   = 0;
    bool         Flag  = true;
    static std::complex<double> psiGrid[gridX][gridY];
    std::vector<double> FreqVec;
    std::vector<std::complex<double>> psiXY;

    // Testing
    const std::complex<double> imU(0., 1.);

    /* Note: 'm' indexes x (rows); 'n' indexes y (columns) */
    for (size_t m = 0; m < gridX; m++)
    {
        if (fabs(-0. + m*dx - 0.1) < 0.01 && Flag == true)
        {
            ind = m;
            std::cout << "ind = " << m << std::endl;
            Flag = false;
        }
        for (size_t n = 0; n < gridY; n++)
        {
                psiGrid[m][n] = 1; // Trying particle in a 2D box.
                // psiGrid[m][n] = 1./sqrt(M_PI) * exp(-(-0. + m*dx-1)*(-0. + m*dx-1)-(-0. + n*dy-1)*(-0. + n*dy-1)); /* Initial condition */
                if(m == 0)       {psiGrid[m][n] = 0;};
                if(m == gridX-1) {psiGrid[m][n] = 0;};
                if(n == 0)       {psiGrid[m][n] = 0;};
                if(n == gridY-1) {psiGrid[m][n] = 0;};
        }
    }

    norm = Get_Norm(psiGrid, dx, dy);
    for (size_t m = 0; m < gridX; m++)
        for (size_t n = 0; n < gridY; n++)
            psiGrid[m][n] /= sqrt(norm); /* Primitive normalization */

    std::cout << "Before normalization: " << norm << "\nAfter normalization : " << Get_Norm(psiGrid, dx, dy) << std::endl;

    /* Specific point: x = 0.1, y = 0; */
    std::cout << "x = 0.1: " << -0. + ind*dx << std::endl;
    std::cout << "y = 0  : " << -0. + gridY/2*dy << std::endl;

    Update_psiGrid(psiGrid, dx, dy, dkx, dky, tMax, dt, ind, psiXY);
    Write_Vec_to_file(psiXY, "PsiXY.txt");

    /* Final step: FFT psiXY (in time domain) */
    const int gridT = int(tMax/dt);
    const double dw = 2 * M_PI / tMax; 
    std::vector<std::complex<double>> psiXY_FT;
    fftw_complex *outT;
    fftw_plan pT;
    outT = (fftw_complex *) fftw_malloc(gridT * sizeof(fftw_complex));
    pT   = fftw_plan_dft_1d(gridT, reinterpret_cast<fftw_complex*>(&psiXY[0]), outT, FFTW_FORWARD, FFTW_ESTIMATE);
    
    fftw_execute(pT);
    fftw_destroy_plan(pT);
    for (size_t i = 0; i < gridT; i++)
    {
            psiXY_FT.push_back(std::complex(outT[i][0], outT[i][1])/double(sqrt(gridT)));
            if(i < gridT/2)
                FreqVec.push_back(dw * i);
            else
                FreqVec.push_back(dw * (int(i) - gridT));
    }
    fftw_free(outT);
    Write_Vec_to_file(psiXY_FT, "PsiXY_FT.txt");
    Write_k_Vec_to_file(FreqVec, "wVec.txt");

    /* I think I'm just missing saving psiXY to a file, so I can then do the spectral analysis.
    Whatever that means... Should I compute this for |Psi|^2? Also, I should check for the normalization of Psi:
    Initial condition is enough, since the propagator mantains the magnitude! It should...
    Do the peaks correspond to eigen-frequencies (<=> eigen-energies)?
    Such that the lowest one would be the ground-state energy? Maybe... */

    /* Why do w just normalize after the inverse transform? */
    /* (Checklist of what is done) Done:
    - Normalization; */
    return 0;
}