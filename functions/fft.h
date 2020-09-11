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
    /// 使用整型进行傅里叶变换时的放大倍数
    /// \tparam t 整型
    template<Integer t>
    constexpr static t omega_times = std::numeric_limits<t>::max() >> (sizeof(t) * 4);
    
    /// 缓存 Ω 以加速 FFT 运算
    static std::vector<double> MEMORY;
    
    /// 用于正变换的 Ω_n^k
    /// \tparam t 复数值数据类型
    /// \param k k
    /// \param n n
    /// \return Ω_n^k
    template<Number t>
    static complex_t<t> omega(unsigned k, unsigned n) {
        auto quarter = n / 4;
        if (MEMORY.size() < quarter) {
            MEMORY.resize(quarter);
            for (unsigned i = 0; i < quarter; ++i)
                MEMORY[i] = std::cos(2 * PI * i / n);
        } else
            while (quarter < MEMORY.size()) {
                k <<= 1u;
                quarter <<= 1u;
            }
        auto i = k % quarter;
        auto result = complex_t<double>{MEMORY[i], i ? MEMORY[quarter - i] : 0};
        switch (k / quarter) {
            case 0:
                break;
            case 1:
                result = {-result.im, result.re};
                break;
            case 2:
                result = -result;
                break;
            case 3:
                result = {result.im, -result.re};
                break;
            default:
                throw std::runtime_error("");
        }
        
        if constexpr (std::is_integral_v<t>)
            return {omega_times<t> * result.re + .5f, omega_times<t> * result.im + .5f};
        else
            return {result.re, result.im};
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
        memory.resize(n, memory.back());
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
    
    /// 用于 FFT 的前后颠倒
    /// \tparam t 数据类型
    /// \param memory 数据
    template<class t>
    void fft_shift(std::vector<t> &memory) {
        auto n = enlarge_to_2_power(memory.size());
        memory.resize(n, memory.back());
        std::swap_ranges(memory.begin(), memory.begin() + n / 2, memory.begin() + n / 2);
    }
    
    template<Number t, unsigned ulp = 1>
    bool almost_equal(t x, t y) {
        constexpr static auto epsilon = std::numeric_limits<t>::epsilon() * ulp;
        auto diff = std::abs(x - y);
        return diff <= std::abs(x + y) * epsilon || diff < std::numeric_limits<t>::min();
    }
    
    template<unsigned order, Number t = float>
    void frft_spectial(std::vector<complex_t<t>> &signal, double sqrt_n) {
        static_assert(order == 1 || order == 2 || order == 3, "only for order 1, 2 or 3");
        if constexpr (order == 2) {
            std::reverse(signal.begin(), signal.end());
        } else if constexpr (order == 1) {
            fft_shift(signal);
            fft<fft_operation::fft>(signal);
            for (auto &x : signal) x /= sqrt_n;
            fft_shift(signal);
        } else {
            fft_shift(signal);
            fft<fft_operation::ifft>(signal); // 免除 ifft 幅度变换
            for (auto &x : signal) x /= sqrt_n;
            fft_shift(signal);
        }
    }
    
    template<Number t = float, Number order_t>
    void frft(std::vector<complex_t<t>> &signal, order_t order) {
        auto n = enlarge_to_2_power(signal.size());
        signal.resize(n, signal.back());
        auto sqrt_n = std::sqrt(n);
        
        // 解周期性
        while (order < 0) order += 4;
        while (order >= 4) order -= 4;
        // 处理特殊情况
        if (almost_equal<double>(order, 0))
            return;
        if (almost_equal<double>(order, 1)) {
            frft_spectial<1>(signal, sqrt_n);
            return;
        }
        if (almost_equal<double>(order, 2)) {
            frft_spectial<2>(signal, sqrt_n);
            return;
        }
        if (almost_equal<double>(order, 3)) {
            frft_spectial<3>(signal, sqrt_n);
            return;
        }
        // 归入 [.5, 1.5)
        if (order > 2) {
            order -= 2;
            frft_spectial<2>(signal, sqrt_n);
        }
        if (order >= 1.5) {
            order -= 1;
            frft_spectial<1>(signal, sqrt_n);
        } else if (order < .5) {
            order += 1;
            frft_spectial<3>(signal, sqrt_n);
            fft_shift(signal);
        }
        { // 时频域同时插值
            signal.resize(2 * n);
            auto p = signal.end(), q = signal.begin() + n;
            while (p != q) {
                *--p = {};
                *--p = *--q;
            }
            fft(signal);
            std::fill(signal.begin() + n / 2, signal.end() - n / 2 + 1, complex_t<t>{});
            fft<fft_operation::ifft>(signal); // 此处免除一次 ifft 的幅度变换，合并到下面
            signal.resize(8 * n, complex_t<t>{}); // 直接扩张到 8n，用于后续做卷积，插值后的信号只使用前面的 4n
            std::transform(signal.begin(), signal.begin() + 2 * n, signal.begin() + n, [n](auto z) { return z / n; });
            std::fill(signal.begin(), signal.begin() + n, complex_t<t>{});
        }
        { // 乘 + 卷 + 乘
            auto alpha = order * PI / 2;
            auto c1 = PI / 4 / n * -std::tan(alpha / 2);
            auto c2 = PI / 4 / n / std::sin(alpha);
            auto chirp = std::vector<complex_t<float>>(8 * n, complex_t<t>{});
            // 第一次乘基啁啾，同时构造另一个啁啾
            auto k = .5 - 2 * n;
            for (unsigned i = 0; i < 4 * n; ++i, ++k) {
                signal[i] *= complex_t<t>::exp(c1 * k * k);
                chirp[i] = complex_t<t>::exp(c2 * k * k);
            }
            // 卷一个啁啾
            fft(signal);
            fft(chirp);
            for (unsigned i = 0; i < 8 * n; ++i)
                signal[i] *= chirp[i];
            ifft(signal);
            { // 抽取，同时乘第三个啁啾和校正项
                // 计算固定的校正项
                const auto z = complex_t<t>::exp(alpha / 2 - PI / 4) / (2 * std::sqrtf(n * std::sinf(alpha)));
                // 初始化用到的循环变量
                auto p = signal.begin();
                auto q = signal.begin() + 3 * n;
                auto e = signal.begin() + n;
                k = .5 - n;
                // 这里重新算一次啁啾，因为第一次算的随着 FFT 毁灭了
                for (; p != e; ++p, q += 2, k += 2)
                    *p = *q * complex_t<t>::exp(c1 * k * k) * z;
                signal.resize(n);
            }
        }
    }
    
    template<class t0, class t1, class f0, class f1> requires Time<t0> && Time<t1> && Frequency<f0> && Frequency<f1>
    auto best_order(t0 t, f0 fs, t1 tl, f1 df) {
        auto x = std::sqrt(floating_seconds(t).count() / fs.template cast_to<Hz_t>().value);
        return std::atan2(floating_seconds(tl).count() / x, -df.template cast_to<Hz_t>().value * x) / (PI / 2);
    }
}

#endif //DSP_SIMULATION_FFT_H
