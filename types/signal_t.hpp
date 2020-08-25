//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_SIGNAL_T_HPP
#define DSP_SIMULATION_SIGNAL_T_HPP

#include <chrono>
#include <vector>

#include "concepts.h"
#include "complex_t.hpp"

namespace mechdancer {
    template<class t>
    concept Signal = requires(t signal){
        typename t::value_t;
        typename t::frequency_t;
        typename t::time_t;
        signal.values;
        signal.sampling_frequency;
        signal.begin_time;
    };
    
    template<class t>
    concept RealSignal = Signal<t> && Number<typename t::value_t>;
    
    template<class t>
    concept ComplexSignal = Signal<t> && std::same_as<typename t::value_t, complex_t<typename t::value_t::value_t>>;
    
    /// 信号类型
    /// \tparam _value_t 数据类型
    /// \tparam _frequency_t 频率类型（frequency_t）
    /// \tparam _time_t 时间类型（std::chrono::duration）
    template<class _value_t, Frequency _frequency_t, Time _time_t>
    struct signal_t {
        using value_t = _value_t;
        using frequency_t = _frequency_t;
        using time_t = _time_t;
        
        std::vector<value_t> values;
        _frequency_t sampling_frequency;
        _time_t begin_time;
    };
    
    /// 构造空的实信号
    /// \tparam _value_t 数据类型
    /// \tparam frequency_t 频率类型
    /// \tparam time_t 时间类型
    /// \param size 信号长度
    /// \param frequency 采样频率
    /// \param time 起始时间
    /// \return 实信号对象
    template<class _value_t = float, Frequency frequency_t, Time time_t>
    auto real_signal_of(size_t size, frequency_t frequency, time_t time) {
        return signal_t<_value_t, frequency_t, time_t>{
            .values = std::vector<_value_t>(size),
            .sampling_frequency = frequency,
            .begin_time = time,
        };
    }
    
    /// 构造空的复信号
    /// \tparam _value_t 复数数据类型
    /// \tparam frequency_t 频率类型
    /// \tparam time_t 时间类型
    /// \param size 信号长度
    /// \param frequency 采样频率
    /// \param time 起始时间
    /// \return 复信号对象
    template<class _value_t = float, Frequency frequency_t, Time time_t>
    auto complex_signal_of(size_t size, frequency_t frequency, time_t time) {
        return real_signal_of<complex_t<_value_t>>(size, frequency, time);
    }
    
    /// 实信号作为实部生成复信号
    /// \tparam t 信号类型
    /// \param signal 实信号
    /// \return 复信号
    template<RealSignal t, Number value_t = typename t::value_t>
    auto complex(t const &signal) {
        auto result = complex_signal_of<value_t>(signal.values.size(), signal.sampling_frequency, signal.begin_time);
        std::transform(signal.values.begin(), signal.values.end(), result.values.begin(), [](auto x) { return complex_t<value_t>::from_real(x); });
        return result;
    }
    
    /// 复信号以取实部的方式转换为实信号
    /// \tparam _value_t 复数值类型
    /// \tparam t 信号类型
    /// \param signal 复信号
    /// \return 复信号实部组成的实信号
    template<ComplexSignal t>
    auto real(t const &signal) {
        using value_t = typename t::value_t::value_t;
        auto result = real_signal_of<value_t>(signal.values.size(), signal.sampling_frequency, signal.begin_time);
        std::transform(signal.values.begin(), signal.values.end(), result.values.begin(), [](auto z) { return z.re; });
        return result;
    }
    
    /// 复信号以取模的方式转换为实信号
    /// \tparam _value_t 复数值类型
    /// \tparam t 信号类型
    /// \param signal 复信号
    /// \return 复信号模组成的实信号
    template<ComplexSignal t, class _value_t = typename t::value_t::value_t>
    auto abs(t const &signal) {
        auto result = real_signal_of<_value_t>(signal.values.size(), signal.sampling_frequency, signal.begin_time);
        std::transform(signal.values.begin(), signal.values.end(), result.values.begin(),
                       [](auto z) { return z.norm(); });
        return result;
    }
}

#endif // DSP_SIMULATION_SIGNAL_T_HPP
