//
// Created by ydrml on 2020/8/7.
//

#ifndef DSP_SIMULATION_PROCESSING_H
#define DSP_SIMULATION_PROCESSING_H

#include <vector>
#include <numeric>

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
        
        for (auto p = A.begin(), q = B.begin(); p < A.end(); ++p, ++q)
            *p *= *q;
        
        ifft(A);
        
        _signal_t result{
            .values = std::vector<value_t>(size),
            .sampling_frequency = a.sampling_frequency,
            .begin_time = a.begin_time + b.begin_time,
        };
        std::transform(A.begin(), A.end(), result.values.begin(), [](complex_t z) { return z.real(); });
        return result;
    }
}

#endif // DSP_SIMULATION_PROCESSING_H
