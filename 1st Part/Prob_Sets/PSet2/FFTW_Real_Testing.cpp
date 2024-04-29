#include <iostream>
#include <cmath>
#include <fftw3.h>

/* Potential */
double V(double x, double y){
    return(-5. * pow((1. + pow(x/5., 2.) + pow(y/4., 2.)), -4.));
}

int main(int argc, const char** argv) {
    int          N = 11;
    double       *in, *outR;
    fftw_complex *out;
    fftw_plan    p, p_Inverse;

    in   = (double*)       malloc(sizeof(double) * N);
    out  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    outR = (double*)       malloc(sizeof(double) * N);

    std::cout << "Initial data: " << std::endl;
    for (size_t i = 0; i < N; i++){
        in[i] = sin(2.*M_PI/10. * i);
        std::cout << "in[" << i << "]: " << in[i]  << std::endl;
    }

    /* Forward transform*/
    p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    fftw_execute(p); /* Repeat as needed */
    fftw_destroy_plan(p);

    /* Padding? Idk what this is called... Re-introducing the "redundant" points */
    std::cout << "\nAfter FFT: " << std::endl;
    for (size_t i = 0; i < N; i++){
        if(i > N/2)
        {
            out[i][0] =   out[N - i][0];
            out[i][1] = - out[N - i][1];
        }
        std::cout << "out[" << i << "]: {" << out[i][0] << ", " << out[i][1] << "i}" << std::endl;
    }

    /* Inverse transform*/
    p_Inverse = fftw_plan_dft_c2r_1d(N, out, outR, FFTW_ESTIMATE);
    fftw_execute(p_Inverse); /* Repeat as needed */
    fftw_destroy_plan(p_Inverse);

    std::cout << "\nAfter FFT and Inverse FFT: " << std::endl;
    for (int i = 0; i < N; i++) {
    printf("outR[%d] = {%f, %fi}\n", i, outR[i]/N, 0.); /* Note the normalization here! */
    }
    std::cout << "Goes back to initial data, which is good!" << std::endl;

    free(in); free(outR); fftw_free(out);
    return 0;
}
