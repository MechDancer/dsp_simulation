//
// Created by ydrml on 2020/8/7.
//

#ifndef DSP_SIMULATION_FFT_H
#define DSP_SIMULATION_FFT_H

#include <cmath>
#include <vector>
#include <limits>

#include "../types/concepts.h"
#include "functions.h"

namespace mechdancer {
    template<Integer t>
    constexpr t omega_times = std::numeric_limits<t>::max() >> (sizeof(t) * 4);
    
    /// 用于正变换的 Ω_n^k
    /// \tparam t 复数值数据类型
    /// \param k k
    /// \param n n
    /// \return Ω_n^k
    template<Number t>
    static complex_t<t> omega(unsigned k, unsigned n) {
        if constexpr (std::is_integral_v<t>) {
            auto theta = 2 * PI * k / n;
            return {static_cast<t>(omega_times<t> * std::cos(theta)),
                    static_cast<t>(omega_times<t> * std::sin(theta))};
        } else {
            t theta = 2 * PI * k / n;
            return {static_cast<t>(std::cos(theta)),
                    static_cast<t>(std::sin(theta))};
        }
    }
    
    /// 用于反变换的 Ω_n^k
    /// \tparam t 复数值数据类型
    /// \param k k
    /// \param n n
    /// \return Ω_n^k
    template<Number t>
    static complex_t<t> i_omega(unsigned k, unsigned n) {
        return omega<t>(k, n).conjugate();
    }
    
    /// fft 操作
    enum class fft_operation { fft, ifft };
    
    /// 基 2 快速傅里叶变换
    /// \tparam operation fft 操作
    /// \tparam t 复数数据类型
    /// \param memory 信号数据空间
    template<fft_operation operation = fft_operation::fft, Number t = float>
    void fft(std::vector<complex_t<t>> &memory) {
        // 编译期断定使用哪种 Ω
        constexpr static auto Ω = operation == fft_operation::fft ? omega<t> : i_omega<t>;
        
        // 扩大尺寸到 2 的幂（以进行基 2 FFT）
        auto n = enlarge_to_2_power(memory.size());
        memory.resize(n, complex_t<t>::zero);
        // 错序
        for (size_t i = 0, j = 0; i < n; ++i) {
            if (i > j) std::swap(memory[i], memory[j]);
            for (size_t l = n >> 1u; (j ^= l) < l; l >>= 1u);
        }
        // 变换
        for (size_t m = 1; m < n; m <<= 1u) {
            auto a = memory.data(), b = a + m;
            const auto s = n / m / 2;
            for (size_t i = 0; i < s; ++i) {
                for (size_t j = 0; j < m; ++j, ++a, ++b)
                    if (b->re == 0 && b->im == 0)
                        *b = *a;
                    else {
                        complex_t<t> c;
                        if constexpr (std::is_integral_v<t>)
                            c = *b * Ω(s * j, n) / omega_times<t>;
                        else
                            c = *b * Ω(s * j, n);
                        *b = *a - c;
                        *a += c;
                    }
                a = b;
                b += m;
            }
        }
    }
    
    /// 反 fft
    /// \tparam t 复数数据类型
    /// \param memory 信号数据空间
    template<Number t = float>
    void ifft(std::vector<complex_t<t>> &memory) {
        fft<fft_operation::ifft>(memory);
        for (auto n = memory.size(); auto &p : memory) p /= n;
    }
}

#endif //DSP_SIMULATION_FFT_H
