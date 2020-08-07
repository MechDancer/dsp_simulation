//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_SIGNAL_T_HPP
#define DSP_SIMULATION_SIGNAL_T_HPP

#include <chrono>
#include <vector>

#include "concepts.h"

namespace mechdancer {
    enum class signal_domain { time, frequency };
    
    template<signal_domain domain, class _value_t, Frequency frequency_t, Time time_t>
    struct signal_t {
        using value_t = _value_t;
        
        std::vector<value_t> values;
        
        frequency_t sampling_frequency;
        time_t      begin_time;
    };
    
    template<signal_domain domain, class _value_t, Frequency frequency_t, Time time_t>
    auto signal_of(size_t size, frequency_t frequency, time_t time) {
        return signal_t<domain, _value_t, frequency_t, time_t>{
            .values = std::vector<_value_t>(size),
            .sampling_frequency = frequency,
            .begin_time = time,
        };
    }
}

#endif // DSP_SIMULATION_SIGNAL_T_HPP
