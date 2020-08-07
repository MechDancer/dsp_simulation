//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_CONCEPTS_H
#define DSP_SIMULATION_CONCEPTS_H

#include <chrono>
#include "frequency_t.hpp"

namespace mechdancer {
    template<class t>
    concept Number = std::is_arithmetic<t>::value;
    
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
}

#endif //DSP_SIMULATION_CONCEPTS_H
