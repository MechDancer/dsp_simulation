//
// Created by ydrml on 2020/8/13.
//

#ifndef DSP_SIMULATION_PROCESS_COMPLEX_H
#define DSP_SIMULATION_PROCESS_COMPLEX_H

#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "functions.h"

namespace mechdancer {
    /// 转换到特定类型计算实信号频谱
    /// \tparam target_t 复数值类型，标准库似乎不支持整数
    /// \tparam t 实信号类型
    /// \param signal 原信号
    /// \param size 最小计算谱长度
    /// \return 频谱
    template<Number target_t, RealSignal t>
    auto fft(t const &signal, size_t size = 0) {
        using frequency_t = typename t::frequency_t;
        using time_t = typename t::time_t;
        using result_t = signal_t <complex_t<target_t>, frequency_t, time_t>;
        
        size = enlarge_to_2_power(std::max(signal.values.size(), size));
        auto result = result_t{
            .values = std::vector<complex_t < target_t>>(size),
            .sampling_frequency = signal.sampling_frequency,
            .begin_time = signal.begin_time,
        };
        std::transform(signal.values.begin(), signal.values.end(), result.values.begin(),
                       [](auto x) { return complex_t < target_t > {static_cast<target_t>(x), 0}; });
        fft(result.values);
        return result;
    }
    
    template<class t, class u> requires ComplexSignal<t> && Frequency<u>
    void bandpass(t &signal, u min, u max) {
        using f_t = typename t::frequency_t;
        using complex_t = typename t::value_t;
        
        if (min >= max) throw std::invalid_argument("");
        
        auto n_min = std::round(signal.values.size() * min.template cast_to<f_t>().value / signal.sampling_frequency.value);
        auto n_max = std::round(signal.values.size() * max.template cast_to<f_t>().value / signal.sampling_frequency.value);
        
        signal.values[signal.values.size() / 2] = complex_t{};
        if (n_min != 0) signal.values[0] = complex_t{};
        
        if (n_min >= signal.values.size() / 2) return;
        std::fill(signal.values.begin(), signal.values.begin() + n_min, complex_t{});
        if (n_min > 0)
            std::fill(signal.values.end() - n_min + 1, signal.values.end(), complex_t{});
        if (n_max < signal.values.size() / 2)
            std::fill(signal.values.begin() + n_max, signal.values.end() - n_max + 1, complex_t{});
    }
}

#endif //DSP_SIMULATION_PROCESS_COMPLEX_H
