//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_CONCEPTS_H
#define DSP_SIMULATION_CONCEPTS_H

#include "frequency_t.hpp"

namespace mechdancer {
    constexpr static auto PI = 3.1415926535897932384626433832795;
    using floating_seconds = std::chrono::duration<float>;
    
    constexpr floating_seconds operator ""_sf(unsigned long long value) {
        return floating_seconds(value);
    }
    
    constexpr floating_seconds operator ""_sf(long double value) {
        return floating_seconds(value);
    }
    
    constexpr floating_seconds operator ""_msf(unsigned long long value) {
        return floating_seconds(value * 1e-3);
    }
    
    constexpr floating_seconds operator ""_msf(long double value) {
        return floating_seconds(value * 1e-3);
    }
    
    template<class t>
    concept Integer = std::is_integral_v<t>;
    
    template<class t>
    concept Floating = std::is_floating_point_v<t>;
    
    template<class t>
    concept Number = std::is_arithmetic_v<t>;
    
    template<class t>
    concept Frequency = requires(t f){ t::value; f.template cast_to<Hz_t>(); };
    
    template<class t>
    concept Time = requires(t time) { t::zero; std::chrono::duration_cast<floating_seconds>(time); };
}

#endif //DSP_SIMULATION_CONCEPTS_H
