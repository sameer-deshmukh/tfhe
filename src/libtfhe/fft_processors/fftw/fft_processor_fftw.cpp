#include <complex>
//#define complex _Complex
#include <fftw3.h>
#include "polynomials.h"
#include "lagrangehalfc_impl.h"
#include <cassert>
#include <cmath>
#include <mutex>

#include <fstream>

FFT_Processor_fftw::FFT_Processor_fftw(const int32_t N): _2N(2*N),N(N),Ns2(N/2) {
    rev_in = (double*) malloc(sizeof(double) * _2N);
    out = (double*) malloc(sizeof(double) * _2N);
    rev_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N+1));
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N+1));
    plan_fftw();
    omegaxminus1 = new cplx[_2N];
    for (int32_t x=0; x<_2N; x++) {
	omegaxminus1[x]=cplx(cos(x*M_PI/N)-1.,-sin(x*M_PI/N)); // instead of cos(x*M_PI/N)-1. + sin(x*M_PI/N) * I
	//exp(i.x.pi/N)-1
    }
}

void FFT_Processor_fftw::plan_fftw() {
    //ensure fftw plan thread safety
    static std::mutex mutex;
    std::lock_guard<std::mutex> lock(mutex);
    rev_p = fftw_plan_dft_r2c_1d(_2N, rev_in, rev_out, FFTW_ESTIMATE);
    p = fftw_plan_dft_c2r_1d(_2N, in, out, FFTW_ESTIMATE);
}

int execute_reverse_int_counter = 0;
void FFT_Processor_fftw::execute_reverse_int(cplx* res, const int* a) {
    cplx* rev_out_cplx = (cplx*) rev_out; //fftw_complex and cplx are layout-compatible
    for (int32_t i=0; i<N; i++) rev_in[i]=a[i]/2.;
    for (int32_t i=0; i<N; i++) rev_in[N+i]=-rev_in[i];

    std::ofstream file;
    file.open(std::string("execute_reverse_int_") + std::to_string(execute_reverse_int_counter) +
              std::string("_input.txt"), std::ios::app);
    file << std::setprecision(40) << N << std::endl;
    for (int i = 0; i < N; ++i) {
      file << rev_in[i] << std::endl;
    }
    file.close();

    fftw_execute(rev_p);
    for (int32_t i=0; i<Ns2; i++) res[i]=rev_out_cplx[2*i+1];
    for (int32_t i=0; i<=Ns2; i++) assert(abs(rev_out_cplx[2*i])<1e-20);

    file.open(std::string("execute_reverse_int_") + std::to_string(execute_reverse_int_counter) +
              std::string("_output.txt"), std::ios::app);
    file << std::setprecision(40) << N << std::endl;
    for (int i = 0; i < N; ++i) {
      file << res[i] << std::endl;
    }
    file.close();

    execute_reverse_int_counter++;
}
int reverse_torus32_counter = 0;
void FFT_Processor_fftw::execute_reverse_torus32(cplx* res, const Torus32* a) {
    static const double _2pm33 = 1./double(INT64_C(1)<<33);
    int32_t* aa = (int32_t*) a;
    cplx* rev_out_cplx = (cplx*) rev_out; //fftw_complex and cplx are layout-compatible
    for (int32_t i=0; i<N; i++) rev_in[i]=aa[i]*_2pm33;
    for (int32_t i=0; i<N; i++) rev_in[N+i]=-rev_in[i];

    std::ofstream file;
    file.open(std::string("execute_reverse_torus32_") + std::to_string(reverse_torus32_counter) +
              std::string("_input.txt"), std::ios::app);

    file << std::setprecision(40) << N << std::endl;
    for (int i = 0; i < N; ++i) {
      file << rev_in[i] << std::endl;
    }
    file.close();

    fftw_execute(rev_p);
    for (int32_t i=0; i<Ns2; i++) res[i]=rev_out_cplx[2*i+1];
    for (int32_t i=0; i<=Ns2; i++) assert(abs(rev_out_cplx[2*i])<1e-20);

    file.open(std::string("execute_reverse_torus32_") + std::to_string(reverse_torus32_counter) +
              std::string("_output.txt"), std::ios::app);

    file << std::setprecision(40) << N << std::endl;
    for (int i = 0; i < N; ++i) {
      file << res[i] << std::endl;
    }
    file.close();

    reverse_torus32_counter++;
}

int Torus32_counter = 0;
void FFT_Processor_fftw::execute_direct_Torus32(Torus32* res, const cplx* a) {
    static const double _2p32 = double(INT64_C(1)<<32);
    static const double _1sN = double(1)/double(N);
    cplx* in_cplx = (cplx*) in; //fftw_complex and cplx are layout-compatible
    for (int32_t i=0; i<=Ns2; i++) in_cplx[2*i]=0;
    for (int32_t i=0; i<Ns2; i++) in_cplx[2*i+1]=a[i];

    std::ofstream file;
    file.open(std::string("execute_direct_Torus32_") + std::to_string(Torus32_counter) +
              std::string("_input.txt"), std::ios::app);
    file << std::setprecision(40) << Ns2 * 2 << std::endl;

    for (int32_t i = 0; i < Ns2 * 2; ++i) {
      file << in_cplx[i] << std::endl;
    }
    file.close();

    fftw_execute(p);
    for (int32_t i=0; i<N; i++) res[i]=Torus32(int64_t(out[i]*_1sN*_2p32));
    //pas besoin du fmod... Torus32(int64_t(fmod(rev_out[i]*_1sN,1.)*_2p32));
    for (int32_t i=0; i<N; i++) assert(fabs(out[N+i]+out[i])<1e-20);

    file.open(std::string("execute_direct_Torus32_") + std::to_string(Torus32_counter) +
              std::string("_output.txt"), std::ios::app);
    file << std::setprecision(40) << Ns2 * 2 << std::endl;
    for (int i = 0; i < N; ++i) {
      file << res[i] << std::endl;
    }

    file.close();
    Torus32_counter++;
}

FFT_Processor_fftw::~FFT_Processor_fftw() {
    fftw_destroy_plan(p);
    fftw_destroy_plan(rev_p);
    fftw_free(in); fftw_free(rev_out);
    free(rev_in); free(out);
    delete[] omegaxminus1;
}


/**
 * FFT functions
 */
EXPORT void IntPolynomial_ifft(LagrangeHalfCPolynomial* result, const IntPolynomial* p) {
    fp1024_fftw.execute_reverse_int(((LagrangeHalfCPolynomial_IMPL*)result)->coefsC, p->coefs);
}
EXPORT void TorusPolynomial_ifft(LagrangeHalfCPolynomial* result, const TorusPolynomial* p) {
    fp1024_fftw.execute_reverse_torus32(((LagrangeHalfCPolynomial_IMPL*)result)->coefsC, p->coefsT);
}
EXPORT void TorusPolynomial_fft(TorusPolynomial* result, const LagrangeHalfCPolynomial* p) {
    fp1024_fftw.execute_direct_Torus32(result->coefsT, ((LagrangeHalfCPolynomial_IMPL*)p)->coefsC);
}
