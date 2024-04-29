#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <complex>

int main(int argc, const char** argv) {
    using namespace std;
    int N = 10;
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (size_t i = 0; i < N; i++)
    {
        using namespace complex_literals;
        std::complex<double> temp(in[i][0], in[i][1]);
        temp += 1i * double(i);
        in[i][0] = real(temp);
        in[i][1] = imag(temp);
        std::cout << temp << std::endl;
    }
    

    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed */
    fftw_destroy_plan(p);

    std::cout << std::endl;
    for (size_t i = 0; i < N; i++)
    {
        std::cout << std::complex(out[i][0], out[i][1])/double(N) << std::endl;
    }

    fftw_free(in); fftw_free(out);
    return 0;
};
