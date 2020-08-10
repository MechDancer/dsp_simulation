//
// Created by ydrml on 2020/8/7.
//

#ifndef DSP_SIMULATION_PROCESSING_H
#define DSP_SIMULATION_PROCESSING_H

#include <vector>
#include <numeric>
#include <functional>

#include "fft.h"

namespace mechdancer {
    /// 求数值向量均值
    /// \tparam t 数据类型
    /// \param values 向量
    /// \return 均值
    template<class t>
    t mean(std::vector<t> const &values) {
        return std::accumulate(values.begin(), values.end(), t{}) / values.size();
    }
    
    /// 快速卷积
    /// \tparam t 实信号类型
    /// \param a 信号 1
    /// \param b 信号 2
    /// \param size 计算长度
    /// \return 卷积信号
    template<Signal _signal_t>
    _signal_t convolution(_signal_t const &a, _signal_t const &b, size_t size = 0) {
        using value_t = typename _signal_t::value_t;
        using complex_t = std::complex<value_t>;
        using spectrum_t = std::vector<complex_t>;
        
        if (a.sampling_frequency != b.sampling_frequency)
            throw std::invalid_argument("the two signals should be with same sampling_frequency");
        
        size = enlarge_to_2_power(std::max(a.values.size() + b.values.size() - 1, size));
        auto A = spectrum_t(size, complex_t{}),
             B = spectrum_t(size, complex_t{});
        
        std::transform(a.values.begin(), a.values.end(), A.begin(), [](value_t x) { return complex_t{x, 0}; });
        std::transform(b.values.begin(), b.values.end(), B.begin(), [](value_t x) { return complex_t{x, 0}; });
        
        fft(A);
        fft(B);
        for (auto p = A.begin(), q = B.begin(); p < A.end(); ++p, ++q) *p *= *q;
        ifft(A);
        
        _signal_t result{
            .values = std::vector<value_t>(size),
            .sampling_frequency = a.sampling_frequency,
            .begin_time = a.begin_time + b.begin_time,
        };
        std::transform(A.begin(), A.end(), result.values.begin(), [](complex_t z) { return z.real(); });
        return result;
    }
    
    /// 互相关模式
    enum class correlation_mode { basic, phat, noise_reduction };
    
    template<Number t, class complex_t = std::complex<t>>
    static complex_t correlation_basic(complex_t r, complex_t s) {
        return std::conj(r) * s;
    }
    
    template<Number t, class complex_t = std::complex<t>>
    static complex_t correlation_phat(complex_t r, complex_t s) {
        auto product = std::conj(r) * s;
        return product / std::abs(product);
    }
    
    template<Number t, class complex_t = std::complex<t>>
    static complex_t correlation_noise_reduction(complex_t r, complex_t s) {
        return std::conj(r) * s / std::abs(s);
    }
    
    /// 频域互相关
    /// \tparam _signal_t 信号类型
    /// \param ref 参考信号
    /// \param signal 目标信号
    /// \param size 计算长度
    /// \return 互相关谱
    template<correlation_mode mode = correlation_mode::basic, Signal _signal_t>
    _signal_t xcorr(_signal_t const &ref, _signal_t const &signal, size_t size = 0) {
        using value_t = typename _signal_t::value_t;
        using complex_t = std::complex<value_t>;
        using spectrum_t = std::vector<complex_t>;
        
        constexpr static auto
            fun = mode == correlation_mode::basic
                  ? correlation_basic<value_t>
                  : mode == correlation_mode::phat
                    ? correlation_phat<value_t>
                    : correlation_noise_reduction<value_t>;
        
        if (ref.sampling_frequency != signal.sampling_frequency)
            throw std::invalid_argument("the two signals should be with same sampling_frequency");
        
        size = enlarge_to_2_power(std::max(ref.values.size() + signal.values.size() - 1, size));
        auto R = spectrum_t(size, complex_t{}),
             S = spectrum_t(size, complex_t{});
        
        std::transform(ref.values.begin(), ref.values.end(), R.begin(), [](value_t x) { return complex_t{x, 0}; });
        std::transform(signal.values.begin(), signal.values.end(), S.begin(), [](value_t x) { return complex_t{x, 0}; });
        
        fft(R);
        fft(S);
        for (auto p = S.begin(), q = R.begin(); p < S.end(); ++p, ++q)
            if (*q == complex_t{})
                *p = complex_t{};
            else if (*p != complex_t{})
                *p = fun(*q, *p);
        ifft(S);
        
        using namespace std::chrono;
        using namespace std::chrono_literals;
        auto      lr = ref.values.size();
        auto      ls = signal.values.size();
        _signal_t result{
            .values = std::vector<value_t>(lr + ls - 1),
            .sampling_frequency = ref.sampling_frequency,
            .begin_time = duration_cast<typename _signal_t::time_t>(
                1s / ref.sampling_frequency.template cast_to<Hz_t>().value - ref.begin_time),
        };
        std::transform(S.end() - lr + 1, S.end(), result.values.begin(), [](complex_t z) { return z.real(); });
        std::transform(S.begin(), S.begin() + ls, result.values.begin() + lr - 1, [](complex_t z) { return z.real(); });
        return result;
    }
    
    /// 希尔伯特变换
    /// \tparam _signal_t 信号类型
    /// \tparam _value_t 信号值类型
    /// \tparam complex_t 复数值类型
    /// \tparam new_signal_t 新信号类型（复信号）
    /// \param signal 原信号
    /// \return 希尔伯特谱
    template<Signal _signal_t,
        class _value_t = typename _signal_t::value_t,
        class complex_t = std::complex<_value_t>,
        class new_signal_t = signal_t<
            complex_t,
            typename _signal_t::frequency_t,
            typename _signal_t::time_t>>
    new_signal_t hilbert(_signal_t const &signal) {
        auto size   = enlarge_to_2_power(signal.values.size());
        auto result = std::vector<complex_t>(size, complex_t{});
        std::transform(signal.values.begin(), signal.values.end(), result.begin(),
                       [](_value_t x) { return complex_t{x, 0}; });
        // 生成超前 90° 的信号（虚部）
        fft(result);
        {
            auto p = result.begin();
            ++p; // 避开 0 频率点，前一半，正频率部分，超前 90°
            while (p < result.begin() + size / 2)
                *p++ = {p->imag(), -p->real()};
            ++p; // 避开 0 频率点，后一半，负频率部分，滞后 90°
            while (p < result.end())
                *p++ = {-p->imag(), p->real()};
        }
        ifft(result);
        // 与原信号合并为复信号
        result.resize(signal.values.size());
        auto p = signal.values.begin();
        auto q = result.begin();
        while (p < signal.values.end()) *q++ = {*p++, q->real()};
        return new_signal_t{
            .values = result,
            .sampling_frequency = signal.sampling_frequency,
            .begin_time = signal.begin_time,
        };
    }
}

#endif // DSP_SIMULATION_PROCESSING_H
