#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <string>
#include <fstream>

/* Potential: Independent of time! */
double V(double x, double y){
    return(-5. * pow((1. + pow(x/5., 2.) + pow(y/4., 2.)), -4.));
}
/* Primitive normalization (Riemann sum type). Ideally, should implement more sophisticated 2D integration techniques. */
template <size_t gridX, size_t gridY>
double Get_Norm(std::complex<double> (&grid)[gridX][gridY], double dx, double dy){
    double norm = 0;
    for (size_t m = 0; m < gridX - 0; m++)
        for (size_t n = 0; n < gridY - 0; n++)
            norm += (abs(grid[m][n])*abs(grid[m][n])) * (dx*dy);
    return(norm);
}

template <size_t gridX, size_t gridY>
void Update_psiGrid(std::complex<double> (&grid)[gridX][gridY], double dx, double dy, double dkx, double dky, double tMax, double dt, int ind, std::vector<std::complex<double>> &psiXY){
    fftw_complex *out, *outR;
    fftw_plan p, pR;
    std::complex<double> *oneDarrayPointer = (std::complex<double> *)grid;
    out  = (fftw_complex *) fftw_malloc(gridX * gridY * sizeof(fftw_complex));
    outR = (fftw_complex *) fftw_malloc(gridX * gridY * sizeof(fftw_complex));
    p   = fftw_plan_dft_2d(gridX, gridY, (fftw_complex *)oneDarrayPointer, out,  FFTW_FORWARD,  FFTW_ESTIMATE);
    pR  = fftw_plan_dft_2d(gridX, gridY, (fftw_complex *)oneDarrayPointer, outR, FFTW_BACKWARD, FFTW_ESTIMATE);
    const std::complex<double> imU(0., 1.); /* Imaginary unit*/
    std::vector<double> kxVec, kyVec;

    /* Defining kVec's for later */
    for (size_t i = 0; i < gridX; i++)
    {
        if(i < gridX/2 + 1) /* Since I use an odd number of points, I add +1 */
            {
                kxVec.push_back(dkx * i);
                kyVec.push_back(dky * i);
            }
        else /* Due to "fftw3's" weirdness... */
            {
                kxVec.push_back(dkx * (i - gridX));
                kyVec.push_back(dky * (i - gridX));
            }
    }
    /* Time loop */
    for (size_t tIndex = 0; tIndex < int(tMax/dt); tIndex++)
    {
        /* Applying exp(-i*t/2*V(0)) */
        for (size_t m = 0; m < gridX; m++)
            for (size_t n = 0; n < gridY; n++)
                grid[m][n] *= std::exp(-imU * dt/2. * V(-10. + m*dx, -10. + n*dy));

        /* Applying 2D FFT: Check the other script for help */
        fftw_execute(p);
        for (size_t i = 0; i < gridX; i++)
            for (size_t j = 0; j < gridY; j++)
              grid[i][j] = std::complex(out[i*gridX + j][0], out[i*gridX + j][1]);

        /* There's some room for improvement here. No need to reassign everything to 'grid' right away. Can use 'out' */
        /* Applying exp(-i*k^2*t) */
        for (size_t i = 0; i < gridX; i++)
            for (size_t j = 0; j < gridY; j++)
                grid[i][j] *= std::exp(-imU * (kxVec[i]*kxVec[i] + kyVec[j]*kyVec[j]) * dt);

        /* Applying Inverse 2D FFT: Check the other script for help */
        fftw_execute(pR);
        for (size_t i = 0; i < gridX; i++)
            for (size_t j = 0; j < gridY; j++)
                grid[i][j] = std::complex(outR[i*gridX + j][0], outR[i*gridX + j][1])/double(gridX*gridY);

        /* Finally, applying exp(-i*t/2*V(t)); Same as before, since V != V(t) */
        for (size_t m = 0; m < gridX; m++)
            for (size_t n = 0; n < gridY; n++)
                grid[m][n] *= std::exp(-imU * dt/2. * V(-10. + m*dx, -10. + n*dy));

        std::cout << "Iteration " << tIndex << " of " << int(tMax/dt) - 1 << "... (Norm_Check) = " << Get_Norm(grid, dx, dy) << std::endl;
        psiXY.push_back(grid[ind][gridY/2]); /* Thing to FT in time, later. Choose, here, between Psi and |Psi|^2 */
        /* psiXY.push_back(abs(grid[ind][gridY/2]) * abs(grid[ind][gridY/2])); */
    }
    fftw_destroy_plan(p);
    fftw_destroy_plan(pR);
    fftw_free(out);
    fftw_free(outR);
    /* Periodic boundary conditions are automatically verified? How/Why? Ask professor... */
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

int main(int argc, const char** argv) {
    const int gridX    = 201; /* How should I choose these values? 2^i, 3^i, etc. ?*/
    const int gridY    = 201;
    const double dx    = 20./(gridX - 1);
    const double dy    = 20./(gridY - 1);
    const double dkx   = 2.*M_PI/20.;
    const double dky   = 2.*M_PI/20.;
    const double tMax  = 100.;
    const double dt    = 0.01;
    double       norm  = 0.;
    int          ind   = 0;
    bool         Flag  = true;
    std::complex<double> psiGrid[gridX][gridY];
    std::vector<double> FreqVec;
    std::vector<std::complex<double>> psiXY;

    /* Note: 'm' indexes x (rows); 'n' indexes y (columns) */
    for (size_t m = 0; m < gridX; m++)
    {
        if (fabs(-10. + m*dx - 0.1) < 0.1 && Flag == true)
        {
            ind = m;
            std::cout << "ind = " << m << std::endl;
            Flag = false;
            // break;
        }
        for (size_t n = 0; n < gridY; n++)
            psiGrid[m][n] = 1./sqrt(M_PI) * exp(-(-10. + m*dx-1)*(-10. + m*dx-1)-(-10. + n*dy-1)*(-10. + n*dy-1)); /* Initial condition */
    }

    norm = Get_Norm(psiGrid, dx, dy);
    for (size_t m = 0; m < gridX; m++)
        for (size_t n = 0; n < gridY; n++)
            psiGrid[m][n] /= sqrt(norm); /* Primitive normalization */

    std::cout << "Before normalization: " << norm << "\nAfter normalization : " << Get_Norm(psiGrid, dx, dy) << std::endl;

    /* Specific point: x = 0.1, y = 0; */
    std::cout << "x = 0.1: " << -10. + ind*dx << std::endl;
    std::cout << "y = 0  : " << -10. + gridY/2*dy << std::endl;

    Update_psiGrid(psiGrid, dx, dy, dkx, dky, tMax, dt, ind, psiXY);
    Write_Vec_to_file(psiXY, "PsiXY.txt");

    /* Final step: FFT psiXY (in time domain) */
    const int gridT = int(tMax/dt);
    const double dw = 2 * M_PI /tMax; 
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

/* I tried to solve the particle in a box case... But doesn't seem to give me correct results...
10.Smth != 9.8Smth ... */