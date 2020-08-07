//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_SIGNAL_T_HPP
#define DSP_SIMULATION_SIGNAL_T_HPP

#include <chrono>
#include <vector>
#include <complex>

#include "concepts.h"

namespace mechdancer {
    /// �ź�����
    /// \tparam _value_t ��������
    /// \tparam frequency_t Ƶ�����ͣ�frequency_t��
    /// \tparam time_t ʱ�����ͣ�std::chrono::duration��
    template<class _value_t, Frequency frequency_t, Time time_t>
    struct signal_t {
        using value_t = _value_t;
        
        std::vector<value_t> values;
        frequency_t          sampling_frequency;
        time_t               begin_time;
    };
    
    template<class _value_t = float, Frequency frequency_t, Time time_t>
    auto real_signal_of(size_t size, frequency_t frequency, time_t time) {
        return signal_t<_value_t, frequency_t, time_t>{
            .values = std::vector<_value_t>(size),
            .sampling_frequency = frequency,
            .begin_time = time,
        };
    }
    
    template<class _value_t = float, Frequency frequency_t, Time time_t>
    auto complex_signal_of(size_t size, frequency_t frequency, time_t time) {
        return real_signal_of<std::complex<_value_t>>(size, frequency, time);
    }
}

#endif // DSP_SIMULATION_SIGNAL_T_HPP
