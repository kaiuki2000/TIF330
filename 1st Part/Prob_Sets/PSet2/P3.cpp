#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <string>
#include <fstream>

/* I'm declaring this here, because I'm lazy... */
static std::vector<std::complex<double>> psiXY;

/* Potential: Independent of time! */
double V(double x, double y){
    return(-5. * pow((1. + pow(x/5., 2.) + pow(y/4., 2.)), -4.));
}

template <size_t gridX, size_t gridY>
void Update_psiGrid(std::complex<double> (&grid)[gridX][gridY], double dx, double dy, double dkx, double dky, double t, int ind){
    /* Applying exp(-i*t/2*V(0)) */
    /* Printing for debugging */
    // std::cout << std::endl;
    for (size_t m = 0; m < gridX; m++)
    {
        for (size_t n = 0; n < gridY; n++)
        {
            using namespace std::complex_literals;
            grid[m][n] *= std::exp(-1i * t/2. * V(-10 + m*dx, -10. + n*dy));
            // std::cout << grid[m][n] << ", ";
        }
    // std::cout << std::endl;
    }

    /* Applying 2D FFT: Check the other script for help */
    fftw_complex *out;
    fftw_plan p;
    std::complex<double> *oneDarrayPointer = (std::complex<double> *)grid;
    out = (fftw_complex *) fftw_malloc(gridX * gridY * sizeof(fftw_complex));
    p   = fftw_plan_dft_2d(gridX, gridY, (fftw_complex *)oneDarrayPointer, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    /* Printing for debugging */
    // std::cout << std::endl;
    for (size_t i = 0; i < gridX; i++)
    {
        for (size_t j = 0; j < gridY; j++)
        {
            grid[i][j] = std::complex(out[i*gridX + j][0], out[i*gridX + j][1]);
            // std::cout << grid[i][j] << ", ";
        }
        // std::cout << std::endl;
    }
    fftw_free(out);

    /* Applying exp(-i*k^2*t) */
    /* Printing for debugging */
    // std::cout << std::endl;
    for (size_t i = 0; i < gridX; i++)
    {
        for (size_t j = 0; j < gridY; j++)
        {
            double ky = 0.;
            double kx = 0.;
            
            if(i < gridX/2)
                kx = dkx * i;
            else /* fftw3 weirdness */
                kx = dkx * (i - gridX);

            if(j < gridY/2)
                ky = dky * j;
            else
                ky = dky * (j - gridY);

            using namespace std::complex_literals;
            grid[i][j] *= std::exp(-1i * (kx*kx + ky*ky) * t);
            // std::cout << grid[i][j] << ", ";
        }
        // std::cout << std::endl;
    }

    /* Applying Inverse 2D FFT: Check the other script for help */
    fftw_complex *outR;
    fftw_plan pR;
    std::complex<double> *oneDarrayPointerR = (std::complex<double> *)grid;
    outR = (fftw_complex *) fftw_malloc(gridX * gridY * sizeof(fftw_complex));
    pR   = fftw_plan_dft_2d(gridX, gridY, (fftw_complex *)oneDarrayPointerR, outR, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(pR);
    fftw_destroy_plan(pR);
    
    /* Printing for debugging */
    // std::cout << std::endl;
    for (size_t i = 0; i < gridX; i++)
    {
        for (size_t j = 0; j < gridY; j++)
        {
            grid[i][j] = std::complex(outR[i*gridX + j][0], outR[i*gridX + j][1])/double(gridX*gridY);
            // std::cout << grid[i][j] << ", ";
        }
        // std::cout << std::endl;
    }
    fftw_free(outR);

    /* Finally, applying exp(-i*t/2*V(t)); Same as before, since V != V(t) */
    /* Printing for debugging */
    // std::cout << std::endl;
    for (size_t m = 0; m < gridX; m++)
    {
        for (size_t n = 0; n < gridY; n++)
        {
            using namespace std::complex_literals;
            grid[m][n] *= std::exp(-1i * t/2. * V(-10. + m*dx, -10. + n*dy));
            // std::cout << grid[m][n] << ", ";
        }
        // std::cout << std::endl;
    }
    psiXY.push_back(grid[ind][gridY/2]);
}

template <size_t gridX, size_t gridY>
double Get_Norm(std::complex<double> (&grid)[gridX][gridY]){
    double norm = 0;
    for (size_t m = 0; m < gridX; m++)
        for (size_t n = 0; n < gridY; n++)
            norm += abs(grid[m][n]) * abs(grid[m][n]);
    return(norm);
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
    const double dy    = 20./(gridX - 1);
    const double dkx   = 2.*M_PI/20.;
    const double dky   = 2.*M_PI/20.;
    const double tMax  = 100.;
    const double dt    = 0.01;
    double       norm  = 0.;
    int          ind   = 0;
    std::complex<double> psiGrid[gridX][gridY];
    std::vector<double> FreqVec;

    /* Note: 'm' indexes x (rows); 'n' indexes y (columns) */
    for (size_t m = 0; m < gridX; m++)
    {
        if (fabs(-10. + m*dx - 0.1) < 0.1)
        {
            ind = m;
            std::cout << "ind = " << m << std::endl;
            break;
        }
        
        for (size_t n = 0; n < gridY; n++)
        {
            psiGrid[m][n] = 1./M_PI * exp(-(-10. + m*dx-1)*(-10. + m*dx-1)-(-10. + n*dy-1)*(-10. + n*dy-1)); /* Initial condition */
            norm += abs(psiGrid[m][n]) * abs(psiGrid[m][n]);
            // std::cout << psiGrid[m][n] << ", ";
        }
        // std::cout << std::endl;
    }

    for (size_t m = 0; m < gridX; m++)
        for (size_t n = 0; n < gridY; n++)
            psiGrid[m][n] /= sqrt(norm); /* Normalization */

    std::cout << "Before normalization: " << norm << "\nAfter normalization : " << Get_Norm(psiGrid) << std::endl;

    /* Specific point: x = 0.1, y = 0; */
    std::cout << "x = 0.1: " << -10. + ind*dx << std::endl;
    std::cout << "y = 0  : " << -10. + gridY/2*dy << std::endl;

    for (size_t i = 1; i < int(tMax/dt); i++){
        Update_psiGrid(psiGrid, dx, dy, dkx, dky, i*dt, ind);
        std::cout << "Iteration " << i << "..." << std::endl;
    }

    /* Final step: FFT psiXY (in time domain) */
    std::vector<std::complex<double>> psiXY_FT;
    Write_Vec_to_file(psiXY, "psiXY.txt");
    fftw_complex *outT;
    fftw_plan pT;
    const int gridT = int(tMax/dt) - 1;
    const double dk = 2 * M_PI /tMax; 
    outT = (fftw_complex *) fftw_malloc((int(tMax/dt) - 1) * sizeof(fftw_complex));
    pT   = fftw_plan_dft_1d((int(tMax/dt) - 1), reinterpret_cast<fftw_complex*>(&psiXY[0]), outT, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(pT);
    fftw_destroy_plan(pT);
     for (size_t i = 0; i < gridT; i++)
    {
            psiXY_FT.push_back(std::complex(outT[i][0], outT[i][1])/double(gridT));
            if(i < gridT/2 + 1)
                FreqVec.push_back(dk * i);
            else{
                FreqVec.push_back(dk * (int(i) - gridT));
                std::cout << dk * (int(i) - gridT) << std::endl;
            }
    }
    fftw_free(outT);
    Write_Vec_to_file(psiXY_FT, "psiXY_FT.txt");
    Write_k_Vec_to_file(FreqVec, "kVec.txt"); /* Change name to *kVec. I think it's better... /

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