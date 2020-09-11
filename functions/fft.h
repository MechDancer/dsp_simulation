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
    /// ʹ�����ͽ��и���Ҷ�任ʱ�ķŴ���
    /// \tparam t ����
    template<Integer t>
    constexpr static t omega_times = std::numeric_limits<t>::max() >> (sizeof(t) * 4);
    
    /// ���� �� �Լ��� FFT ����
    static std::vector<double> MEMORY;
    
    /// �������任�� ��_n^k
    /// \tparam t ����ֵ��������
    /// \param k k
    /// \param n n
    /// \return ��_n^k
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
    
    /// ���ڷ��任�� ��_n^k
    /// \tparam t ����ֵ��������
    /// \param k k
    /// \param n n
    /// \return ��_n^k
    template<Number t>
    static complex_t<t> i_omega(unsigned k, unsigned n) {
        return omega<t>(k, n).conjugate();
    }
    
    /// fft ����
    enum class fft_operation { fft, ifft };
    
    /// �� 2 ���ٸ���Ҷ�任
    /// \tparam operation fft ����
    /// \tparam t ������������
    /// \param memory �ź����ݿռ�
    template<fft_operation operation = fft_operation::fft, Number t = float>
    void fft(std::vector<complex_t<t>> &memory) {
        // �����ڶ϶�ʹ������ ��
        constexpr static auto �� = operation == fft_operation::fft ? omega<t> : i_omega<t>;
        
        // ����ߴ絽 2 ���ݣ��Խ��л� 2 FFT��
        auto n = enlarge_to_2_power(memory.size());
        memory.resize(n, memory.back());
        // ����
        for (size_t i = 0, j = 0; i < n; ++i) {
            if (i > j) std::swap(memory[i], memory[j]);
            for (size_t l = n >> 1u; (j ^= l) < l; l >>= 1u);
        }
        // �任
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
                            c = *b * ��(s * j, n) / omega_times<t>;
                        else
                            c = *b * ��(s * j, n);
                        *b = *a - c;
                        *a += c;
                    }
                a = b;
                b += m;
            }
        }
    }
    
    /// �� fft
    /// \tparam t ������������
    /// \param memory �ź����ݿռ�
    template<Number t = float>
    void ifft(std::vector<complex_t<t>> &memory) {
        fft<fft_operation::ifft>(memory);
        for (auto n = memory.size(); auto &p : memory) p /= n;
    }
    
    /// ���� FFT ��ǰ��ߵ�
    /// \tparam t ��������
    /// \param memory ����
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
            fft<fft_operation::ifft>(signal); // ��� ifft ���ȱ任
            for (auto &x : signal) x /= sqrt_n;
            fft_shift(signal);
        }
    }
    
    template<Number t = float, Number order_t>
    void frft(std::vector<complex_t<t>> &signal, order_t order) {
        auto n = enlarge_to_2_power(signal.size());
        signal.resize(n, signal.back());
        auto sqrt_n = std::sqrt(n);
        
        // ��������
        while (order < 0) order += 4;
        while (order >= 4) order -= 4;
        // �����������
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
        // ���� [.5, 1.5)
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
        { // ʱƵ��ͬʱ��ֵ
            signal.resize(2 * n);
            auto p = signal.end(), q = signal.begin() + n;
            while (p != q) {
                *--p = {};
                *--p = *--q;
            }
            fft(signal);
            std::fill(signal.begin() + n / 2, signal.end() - n / 2 + 1, complex_t<t>{});
            fft<fft_operation::ifft>(signal); // �˴����һ�� ifft �ķ��ȱ任���ϲ�������
            signal.resize(8 * n, complex_t<t>{}); // ֱ�����ŵ� 8n�����ں������������ֵ����ź�ֻʹ��ǰ��� 4n
            std::transform(signal.begin(), signal.begin() + 2 * n, signal.begin() + n, [n](auto z) { return z / n; });
            std::fill(signal.begin(), signal.begin() + n, complex_t<t>{});
        }
        { // �� + �� + ��
            auto alpha = order * PI / 2;
            auto c1 = PI / 4 / n * -std::tan(alpha / 2);
            auto c2 = PI / 4 / n / std::sin(alpha);
            auto chirp = std::vector<complex_t<float>>(8 * n, complex_t<t>{});
            // ��һ�γ˻���ౣ�ͬʱ������һ�����
            auto k = .5 - 2 * n;
            for (unsigned i = 0; i < 4 * n; ++i, ++k) {
                signal[i] *= complex_t<t>::exp(c1 * k * k);
                chirp[i] = complex_t<t>::exp(c2 * k * k);
            }
            // ��һ�����
            fft(signal);
            fft(chirp);
            for (unsigned i = 0; i < 8 * n; ++i)
                signal[i] *= chirp[i];
            ifft(signal);
            { // ��ȡ��ͬʱ�˵�������౺�У����
                // ����̶���У����
                const auto z = complex_t<t>::exp(alpha / 2 - PI / 4) / (2 * std::sqrtf(n * std::sinf(alpha)));
                // ��ʼ���õ���ѭ������
                auto p = signal.begin();
                auto q = signal.begin() + 3 * n;
                auto e = signal.begin() + n;
                k = .5 - n;
                // ����������һ����ౣ���Ϊ��һ��������� FFT ������
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
