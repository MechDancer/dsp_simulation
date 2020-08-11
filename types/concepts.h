//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_CONCEPTS_H
#define DSP_SIMULATION_CONCEPTS_H

#include <chrono>
#include <complex>

#include "frequency_t.hpp"

namespace mechdancer {
    template<class t>
    concept Integer = std::is_integral_v<t>;
    
    template<class t>
    concept Floating = std::is_floating_point_v<t>;
    
    template<class t>
    concept Number = std::is_arithmetic_v<t>;
    
    template<class t>
    concept Frequency = requires(t f){ t::value; f.template cast_to<Hz_t>(); };
    
    template<class t>
    concept Time = requires(t time) { t::zero; std::chrono::duration_cast<std::chrono::duration<double>>(time); };
    
    template<class t>
    concept Signal = requires(t signal){
        typename t::value_t;
        signal.values;
        signal.sampling_frequency;
        signal.begin_time;
    };
    
    template<class t>
    concept RealSignal = Signal<t> && Number<typename t::value_t>;
    
    template<class t, class u>
    concept ComplexSignal = Signal<t> && Number<u> && std::same_as<typename t::value_t, std::complex<u>>;
}

#endif //DSP_SIMULATION_CONCEPTS_H
