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
    
    template<signal_domain domain, Number _value_t, Frequency frequency_t, Time time_t>
    struct signal_t {
        using value_t = _value_t;
        
        std::vector<value_t> values;
        
        frequency_t sampling_frequency;
        time_t      begin_time;
    };
}

#endif // DSP_SIMULATION_SIGNAL_T_HPP
