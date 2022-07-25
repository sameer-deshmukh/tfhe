cmake_minimum_required(VERSION 3.0)

# This is the nayuki fft processor for the tfhe library
#enable_language(ASM)

set(SRCS
    fft_processor_fftw.cpp
    lagrangehalfc_impl.cpp
    )

set(HEADERS
    lagrangehalfc_impl.h
    )

include_directories(${FFTW_INCLUDES})
add_library(tfhe-fft-fftw OBJECT ${SRCS} ${HEADERS})
set_property(TARGET tfhe-fft-fftw PROPERTY POSITION_INDEPENDENT_CODE ON)

pkg_check_modules(NEXT_POSITS REQUIRED IMPORTED_TARGET next_posits)
pkg_check_modules(GMP         REQUIRED IMPORTED_TARGET gmp)
pkg_check_modules(MPFR        REQUIRED IMPORTED_TARGET mpfr)

# target_link_libraries(tfhe-fft-fftw
#         PkgConfig::NEXT_POSITS)

# target_link_libraries(tfhe-fft-fftw ${NEXT_POSITS_LDFLAGS})
# target_link_libraries(tfhe-fft-fftw ${NEXT_POSITS_LIBRARIES})
# target_link_libraries(tfhe-fft-fftw "-rpath ${NEXT_POSITS_LIBRARY_DIRS}")
target_include_directories(tfhe-fft-fftw PUBLIC ${NEXT_POSITS_INCLUDE_DIRS})
target_include_directories(tfhe-fft-fftw PUBLIC ${GMP_INCLUDE_DIRS})
target_include_directories(tfhe-fft-fftw PUBLIC ${MPFR_INCLUDE_DIRS})
