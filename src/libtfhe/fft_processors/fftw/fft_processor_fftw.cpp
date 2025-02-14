#include <complex>
//#define complex _Complex
#include <fftw3.h>
#include "polynomials.h"
#include "lagrangehalfc_impl.h"
#include <cassert>
#include <cmath>
#include <mutex>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "next_posits/posit_32.h"

unsigned long long num_inputs = 0, num_outputs = 0, exact_inputs = 0, exact_outputs = 0;

#define TOLERANCE 1e-6

bool
mpfr_p32_is_exact(mpfr_srcptr exact_input, posit32 converted_posit) {
  mpfr_t conversion, rel_error, diff, fp32_rel_error, fp64_rel_error;
  mpfr_init2(conversion, precision);
  mpfr_init2(rel_error, precision);
  mpfr_init2(diff, precision);

  mpfr_init2(fp32_rel_error, precision);
  mpfr_init2(fp64_rel_error, precision);
  bool is_exact = false;

  posit2mpfr(conversion, converted_posit);

  mpfr_sub(diff, exact_input, conversion, MPFR_RNDN);
  mpfr_div(rel_error, diff, exact_input, MPFR_RNDN);
  mpfr_abs(rel_error, rel_error, MPFR_RNDN);

  if (mpfr_cmp_d(rel_error, TOLERANCE) <= 0) {
    is_exact = true;
  }
  else {
    mpfr_printf("input: %.32Rf converted: %.32Rf rel_error: %.32Rf\n",
                exact_input, conversion, rel_error);
  }

  float converted_float = mpfr_get_flt(exact_input, MPFR_RNDN);
  mpfr_set_flt(conversion, converted_float, MPFR_RNDN);
  mpfr_sub(diff, exact_input, conversion, MPFR_RNDN);
  mpfr_div(fp32_rel_error, diff, exact_input, MPFR_RNDN);
  mpfr_abs(fp32_rel_error, fp32_rel_error, MPFR_RNDN);

  double converted_double = mpfr_get_d(exact_input, MPFR_RNDN);
  mpfr_set_d(conversion, converted_double, MPFR_RNDN);
  mpfr_sub(diff, exact_input, conversion, MPFR_RNDN);
  mpfr_div(fp64_rel_error, diff, exact_input, MPFR_RNDN);
  mpfr_abs(fp64_rel_error, fp64_rel_error, MPFR_RNDN);

  FILE *file;
  double a[4];
  a[0] = mpfr_get_d(exact_input, MPFR_RNDN);
  a[1] = mpfr_get_d(rel_error, MPFR_RNDN);
  a[2] = mpfr_get_d(fp32_rel_error, MPFR_RNDN);
  a[3] = mpfr_get_d(fp64_rel_error, MPFR_RNDN);

  file = fopen("fft_accuracy_long.dat", "ab");
  fwrite(a, sizeof(double), 4, file);
  // fprintf(file, "%lf %lf %lf %lf\n", a, b, c, d);
  // mpfr_fprintf(file, "%.16Rf,%.16Rf,%.16Rf,%.16Rf\n",
  //              exact_input, rel_error, fp32_rel_error, fp64_rel_error);
  fclose(file);

  mpfr_clear(conversion);
  mpfr_clear(rel_error);
  mpfr_clear(diff);
  mpfr_clear(fp32_rel_error);
  mpfr_clear(fp64_rel_error);

  return is_exact;
}

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
    // IFFT
    rev_p = fftw_plan_dft_r2c_1d(_2N, rev_in, rev_out, FFTW_ESTIMATE);
    // FFT
    p = fftw_plan_dft_c2r_1d(_2N, in, out, FFTW_ESTIMATE);
}

// accepts an integer and outputs a complex.
int execute_reverse_int_counter = 0;
void FFT_Processor_fftw::execute_reverse_int(cplx* res, const int* a) {
    mpfr_t exact_input;
    mpfr_init2(exact_input, precision);
    cplx* rev_out_cplx = (cplx*) rev_out; //fftw_complex and cplx are layout-compatible
    for (int32_t i=0; i<N; i++) rev_in[i]=a[i]/2.;
    for (int32_t i=0; i<N; i++) rev_in[N+i]=-rev_in[i];

    // std::ofstream file;
    // file.open(std::string("execute_reverse_int_") + std::to_string(execute_reverse_int_counter) +
    //           std::string("_input.txt"), std::ios::app);
    // file << std::setprecision(40) << N << std::endl;
    // for (int i = 0; i < N; ++i) {
    //   file << rev_in[i] << std::endl;
    // }
    // file.close();

    for (int32_t i = 0; i < N; ++i) {
      // count number of inputs
      num_inputs++;
      mpfr_set_d(exact_input, rev_in[i], MPFR_RNDN);
      posit32 converted_posit = mpfr2posit(exact_input); // convert the input into a p32.
      // check if the posit exactly matches the int32 input with MPFR.
      if (mpfr_p32_is_exact(exact_input, converted_posit)) {
        exact_inputs++;
      }
    }

    fftw_execute(rev_p);
    for (int32_t i=0; i<Ns2; i++) res[i]=rev_out_cplx[2*i+1];
    for (int32_t i=0; i<=Ns2; i++) assert(abs(rev_out_cplx[2*i])<1e-20);

    // file.open(std::string("execute_reverse_int_") + std::to_string(execute_reverse_int_counter) +
    //           std::string("_output.txt"), std::ios::app);
    // file << std::setprecision(40) << N << std::endl;
    // for (int i = 0; i < N; ++i) {
    //   file << res[i] << std::endl;
    // }
    // file.close();

    // execute_reverse_int_counter++;
    mpfr_clear(exact_input);

    // std::cout << "inputs: " << num_inputs << " outputs: " << num_outputs
    //           << " exact ip: " << exact_inputs << " exact op: " << exact_outputs << std::endl;
}

// accepts an integer and outputs a complex.
int reverse_torus32_counter = 0;
void FFT_Processor_fftw::execute_reverse_torus32(cplx* res, const Torus32* a) {
    mpfr_t exact_input;
    mpfr_init2(exact_input, precision);
    static const double _2pm33 = 1./double(INT64_C(1)<<33);
    int32_t* aa = (int32_t*) a;
    cplx* rev_out_cplx = (cplx*) rev_out; //fftw_complex and cplx are layout-compatible
    for (int32_t i=0; i<N; i++) rev_in[i]=aa[i]*_2pm33;
    for (int32_t i=0; i<N; i++) rev_in[N+i]=-rev_in[i];

    // std::ofstream file;
    // file.open(std::string("execute_reverse_torus32_") + std::to_string(reverse_torus32_counter) +
    //           std::string("_input.txt"), std::ios::app);

    // file << std::setprecision(40) << N << std::endl;
    // for (int i = 0; i < N; ++i) {
    //   file << rev_in[i] << std::endl;
    // }
    // file.close();

    for (int32_t i = 0; i < N; ++i) {
      // count number of inputs
      num_inputs++;

      mpfr_set_d(exact_input, rev_in[i], MPFR_RNDN);
      posit32 converted_posit = mpfr2posit(exact_input);
      // check if the posit exactly matches the int32 input with MPFR.
      if (mpfr_p32_is_exact(exact_input, converted_posit)) {
        exact_inputs++;
      }
    }

    fftw_execute(rev_p);
    for (int32_t i=0; i<Ns2; i++) res[i]=rev_out_cplx[2*i+1];
    for (int32_t i=0; i<=Ns2; i++) assert(abs(rev_out_cplx[2*i])<1e-20);

    // file.open(std::string("execute_reverse_torus32_") + std::to_string(reverse_torus32_counter) +
    //           std::string("_output.txt"), std::ios::app);
    // file << std::setprecision(40) << N << std::endl;
    // for (int i = 0; i < N; ++i) {
    //   file << res[i] << std::endl;
    // }
    // file.close();
    // reverse_torus32_counter++;
    mpfr_clear(exact_input);
}

// accepts a complex and outputs an integer.
int Torus32_counter = 0;
void FFT_Processor_fftw::execute_direct_Torus32(Torus32* res, const cplx* a) {
    mpfr_t exact_output;
    mpfr_init2(exact_output, precision);
    static const double _2p32 = double(INT64_C(1)<<32);
    static const double _1sN = double(1)/double(N);
    cplx* in_cplx = (cplx*) in; //fftw_complex and cplx are layout-compatible
    for (int32_t i=0; i<=Ns2; i++) in_cplx[2*i]=0;
    for (int32_t i=0; i<Ns2; i++) in_cplx[2*i+1]=a[i];

    // std::ofstream file;
    // file.open(std::string("execute_direct_Torus32_") + std::to_string(Torus32_counter) +
    //           std::string("_input.txt"), std::ios::app);
    // file << std::setprecision(40) << Ns2 * 2 << std::endl;

    // for (int32_t i = 0; i < Ns2 * 2; ++i) {
    //   file << in_cplx[i] << std::endl;
    // }
    // file.close();

    fftw_execute(p);
    for (int32_t i=0; i<N; i++) res[i]=Torus32(int64_t(out[i]*_1sN*_2p32));
    //pas besoin du fmod... Torus32(int64_t(fmod(rev_out[i]*_1sN,1.)*_2p32));
    for (int32_t i=0; i<N; i++) assert(fabs(out[N+i]+out[i])<1e-20);

    // file.open(std::string("execute_direct_Torus32_") + std::to_string(Torus32_counter) +
    //           std::string("_output.txt"), std::ios::app);
    // file << std::setprecision(40) << Ns2 * 2 << std::endl;
    // for (int i = 0; i < N; ++i) {
    //   file << res[i] << std::endl;
    // }

    for (int32_t i = 0; i < N; ++i) {
      // count number of inputs
      num_outputs++;
      mpfr_set_d(exact_output, res[i], MPFR_RNDN);
      posit32 converted_posit = mpfr2posit(exact_output);       // check if the posit exactly matches the int32 input with MPFR.
      if (mpfr_p32_is_exact(exact_output, converted_posit)) {
        exact_outputs++;
      }
    }


    // file.close();
    // Torus32_counter++;

    mpfr_clear(exact_output);
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
