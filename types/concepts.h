//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_CONCEPTS_H
#define DSP_SIMULATION_CONCEPTS_H

#include "frequency_t.hpp"

namespace mechdancer {
    using floating_seconds = std::chrono::duration<float>;
    
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
