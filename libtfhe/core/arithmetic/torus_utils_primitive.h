#ifndef TFHE_TORUS_UTILS_PRIMITIVE_H
#define TFHE_TORUS_UTILS_PRIMITIVE_H

#include "torus.h"
#include "zmodule_param.h"
#include "core/allocator/allocator.h"

#include <limits>
#include <type_traits>
#include <cstdint>
#include <cmath>

/**
 * Numeric utils for Torus data types
 * @tparam TORUS type of torus (must be an integral type)
 */
template<typename TORUS>
class TorusUtils {
    static_assert(std::is_integral<TORUS>::value and std::is_signed<TORUS>::value,
                  "TorusUtils<T> defined only for native signed integer types");

public:
    typedef typename std::make_signed<TORUS>::type INT_TYPE;
    static constexpr TORUS min_val = std::numeric_limits<TORUS>::min();
    static constexpr TORUS max_val = std::numeric_limits<TORUS>::max();

private:
    typedef typename std::make_unsigned<TORUS>::type UTORUS; ///< unsigned version of the torus
    typedef typename std::make_signed<TORUS>::type STORUS;   ///< signed version of the torus

    static constexpr uint32_t bit_cnt = sizeof(TORUS) * 8;       ///< bit size of the Torus

    /**
     * @brief 2^p as double where p = bit_size
     * if x is a double between -0.5 and 0.5, x*dtot_factor is the torus representation of x
     */
    static constexpr double dtot_factor = 4.0 * (TORUS(1) << (bit_cnt - 2));
    /**
     * @brief 2^{p-1} as UTORUS where p = bit_size
     */
    static constexpr UTORUS half_utorus = UTORUS(1) << (bit_cnt-1);

public:
    /**
     * @brief Converts real number to torus (mainly for printing)
     * @param d real to convert
     * @return torus value
     */
    static TORUS from_double(const double d, const ZModuleParams<TORUS> *params);

    /**
     * @brief Converts real number to torus (mainly for printing)
     * @param d real to convert
     * @param reps torus value
     */
    static void from_double(TORUS *reps, const double d, const ZModuleParams<TORUS> *params);

    /**
     * @brief Converts torus to real number (mainly for printing)
     * @param x torus element to convert
     * @return real number
    */
    static double to_double(const TORUS *x, const ZModuleParams<TORUS> *params);

    /**
     * @brief Rounds the torus value to the nearest multiple of 1/Msize
     * @param phase torus value
     * @param Msize discrete space size
     * @return approximated torus value
     */
    static void approxPhase(TORUS &res, const TORUS phase, const uint64_t Msize, const ZModuleParams<TORUS> *params, Allocator alloc);

    /**
     * @brief Mod-Rescale from the torus to Z/Msize.Z
     *  i.e. computes roundToNearestInteger(Msize*phase)
     * @param phase torus value to convert
     * @param Msize discrete space size
     * @return discrete space value in [0, MSize[
     */
    static uint64_t modSwitchFromTorus(const TORUS phase, const uint64_t Msize, const ZModuleParams<TORUS> *params, Allocator alloc);

    /**
     * @brief Converts discrete message space to torus
     *  i.e. value mu/Msize to a torus for mu
     * @param mu discrete space value (from [0,Msize[) to convert
     * @param Msize discrete space size
     * @return torus value
     */
    static void modSwitchToTorus(TORUS &res, const uint64_t mu, const uint64_t Msize, const ZModuleParams<TORUS> *params, Allocator alloc);

    /**
     * @brief Return absolute distance between 2 torus elements
     *
     * @param t1 first torus element
     * @param t2 second torus element
     * @return double value of the infinity norm
     */
    static double distance(const TORUS *t1, const TORUS *t2, const ZModuleParams<TORUS> *params, Allocator alloc);
};

#endif //TFHE_TORUS_UTILS_PRIMITIVE_H
