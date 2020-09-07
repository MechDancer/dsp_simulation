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
        
        /// 转换信号类型
        /// \tparam __value_t 新的值类型
        /// \tparam __frequency_t 新的频率类型
        /// \tparam __time_t 新的时间类型
        /// \tparam converter_t 转换器类型
        /// \param new_size 新的离散信号长度
        /// \param converter 转换器
        /// \return 目标类型的信号
        template<class __value_t = _value_t, Frequency __frequency_t = _frequency_t, Time __time_t = _time_t, class converter_t>
        auto cast(size_t new_size, converter_t converter) const {
            auto result = signal_t<__value_t, __frequency_t, __time_t>{
                .values = std::vector<__value_t>(new_size > 0 ? new_size : values.size(), __value_t{}),
                .sampling_frequency = sampling_frequency.template cast_to<__frequency_t>(),
                .begin_time = std::chrono::duration_cast<__time_t>(begin_time),
            };
            if constexpr (std::is_same_v<converter_t, nullptr_t>)
                std::copy(values.begin(), values.end(), result.values.begin());
            else
                std::transform(values.begin(), values.end(), result.values.begin(), converter);
            return result;
        }
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
    auto signal_of(size_t size, frequency_t frequency, time_t time) {
        return signal_t<_value_t, frequency_t, time_t>{
            .values = std::vector<_value_t>(size),
            .sampling_frequency = frequency,
            .begin_time = time,
        };
    }
    
    /// 实信号作为实部生成复信号
    /// \tparam t 信号类型
    /// \param signal 实信号
    /// \return 复信号
    template<RealSignal t, Number value_t = typename t::value_t>
    auto complex(t const &signal) {
        return signal.template cast<complex_t<value_t>>(0, [](auto x) -> complex_t<value_t> { return x; });
    }
    
    /// 复信号以取实部的方式转换为实信号
    /// \tparam _value_t 复数值类型
    /// \tparam t 信号类型
    /// \param signal 复信号
    /// \return 复信号实部组成的实信号
    template<ComplexSignal t, class value_t = typename t::value_t::value_t>
    auto real(t const &signal) {
        return signal.template cast<value_t>(0, [](auto z) -> value_t { return z.re; });
    }
    
    /// 复信号以取模的方式转换为实信号
    /// \tparam _value_t 复数值类型
    /// \tparam t 信号类型
    /// \param signal 复信号
    /// \return 复信号模组成的实信号
    template<ComplexSignal t, class value_t = typename t::value_t::value_t>
    auto abs(t const &signal) {
        return signal.template cast<value_t>(0, [](auto z) -> value_t { return z.norm(); });
    }
}

#endif // DSP_SIMULATION_SIGNAL_T_HPP
