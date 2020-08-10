//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_SIGNAL_T_HPP
#define DSP_SIMULATION_SIGNAL_T_HPP

#include <chrono>
#include <vector>
#include <complex>

#include "concepts.h"
#include "../functions/fft.h"

namespace mechdancer {
    /// 信号类型
    /// \tparam _value_t 数据类型
    /// \tparam _frequency_t 频率类型（frequency_t）
    /// \tparam time_t 时间类型（std::chrono::duration）
    template<class _value_t, Frequency _frequency_t, Time time_t>
    struct signal_t {
        using value_t = _value_t;
        
        std::vector<value_t> values;
        _frequency_t         sampling_frequency;
        time_t               begin_time;
        
        /// 改变采样率重采样
        /// \tparam new_frequency_t 新采样频率类型
        /// \tparam new_signal_t 新信号类型
        /// \param new_fs 新采样率
        /// \param _times 处理倍率，倍率越高，精度越高
        /// \return 新信号
        template<Number times_t, Frequency new_frequency_t, Signal new_signal_t = signal_t<_value_t, new_frequency_t, time_t>>
        new_signal_t resample(new_frequency_t new_fs, times_t _times) const {
            using complex_t = std::complex<_value_t>;
            
            // 检查新频率，若新频率与旧频率相同，直接返回
            auto old_fs = sampling_frequency.template cast_to<new_frequency_t>();
            if (old_fs == new_fs)
                return new_signal_t{
                    .values = values,
                    .sampling_frequency = new_fs,
                    .begin_time = begin_time,
                };
            // 在基 2 条件下计算实际重采样倍率和降采样间隔，若间隔过小，要求提高倍率
            auto times    = enlarge_to_2_power(values.size() * _times) / enlarge_to_2_power(values.size());
            auto interval = static_cast<size_t>(static_cast<double>(old_fs.value * times) / new_fs.value + .5);
            if (interval == 0u)
                throw std::invalid_argument("processing times is too little");
            // 按是否进行升采样分类
            if (times > 1) {
                // 基于 FFT 升采样
                auto spectrum = std::vector<complex_t>(values.size());
                std::transform(values.begin(), values.end(), spectrum.begin(), [](value_t x) { return complex_t{x, 0}; });
                
                fft(spectrum);
                auto size = spectrum.size();
                spectrum.resize(size * times, complex_t{});
                std::copy_n(spectrum.begin() + size / 2, size / 2, spectrum.end() - size / 2);
                std::fill(spectrum.begin() + size / 2, spectrum.end() - size, complex_t{});
                ifft(spectrum);
                // 再降采样到目标采样率
                new_signal_t result{
                    .values = std::vector<value_t>((spectrum.size() + interval - 1) / interval),
                    .sampling_frequency = new_fs,
                    .begin_time = begin_time,
                };
                if (interval == 1u)
                    std::transform(spectrum.begin(), spectrum.end(), result.values.begin(), [](complex_t z) { return z.real(); });
                else {
                    auto p = spectrum.begin();
                    auto q = result.values.begin();
                    *q++ = p->real();
                    while (q != result.values.end()) {
                        p += interval;
                        *q++ = p->real();
                    }
                }
                return result;
            } else {
                // 从原信号直接降采样
                new_signal_t result{
                    .values = std::vector<value_t>((values.size() + interval - 1) / interval),
                    .sampling_frequency = new_fs,
                    .begin_time = begin_time,
                };
                if (interval == 1u)
                    std::copy(values.begin(), values.end(), result.values.begin());
                else {
                    auto p = values.begin();
                    auto q = result.values.begin();
                    *q++ = *p;
                    while (q != result.values.end()) {
                        p += interval;
                        *q++ = *p;
                    }
                }
                return result;
            }
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
        return real_signal_of<std::complex<_value_t>>(size, frequency, time);
    }
}

#endif // DSP_SIMULATION_SIGNAL_T_HPP
