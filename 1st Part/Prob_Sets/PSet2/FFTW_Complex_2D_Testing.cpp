#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <complex>

int main(int argc, const char** argv) {
    using namespace std;
    int n0 = 4;
    int n1 = 4;
    fftw_complex *out, *outR;
    fftw_plan p;

    std::complex<double> twoDarray[n0][n1];
    std::complex<double> twoDarrayOUT[n0][n1];

    /* Initial values */
    for (size_t i = 0; i < n0; i++)
    {
        for (size_t j = 0; j < n1; j++)
        {
            twoDarray[i][j].real(i);
            twoDarray[i][j].imag(j);
            std::cout << twoDarray[i][j] << ", ";
        }
        std::cout << std::endl;
    }

    std::complex<double> *oneDarrayPointer = (std::complex<double> *)twoDarray;
    out  = (fftw_complex *) fftw_malloc(n0 * n1 * sizeof(fftw_complex));
    outR = (fftw_complex *) fftw_malloc(n0 * n1 * sizeof(fftw_complex));

    /* Forward transform */
    p = fftw_plan_dft_2d(n0, n1, (fftw_complex *)oneDarrayPointer, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */
    fftw_destroy_plan(p);

    std::cout << std::endl;
    for (size_t i = 0; i < n0; i++)
    {
        for (size_t j = 0; j < n1; j++)
        {
            twoDarrayOUT[i][j] = std::complex(out[i*n0 + j][0], out[i*n0 + j][1]);
            std::cout << twoDarrayOUT[i][j] << ", ";
        }
        std::cout << std::endl;
    }

    /* Backward transform */
    fftw_plan pR;
    pR = fftw_plan_dft_2d(n0, n1, out, outR, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(pR); /* repeat as needed */
    fftw_destroy_plan(pR);

    std::cout << std::endl;
    for (size_t i = 0; i < n0; i++)
    {
        for (size_t j = 0; j < n1; j++)
        {
            twoDarrayOUT[i][j] = std::complex(outR[i*n0 + j][0], outR[i*n0 + j][1])/double(n0*n1);
            std::cout << twoDarrayOUT[i][j] << ", ";
        }
        std::cout << std::endl;
    }

    fftw_free(out);
    return 0;
};