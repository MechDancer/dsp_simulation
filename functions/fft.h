//
// Created by ydrml on 2020/8/7.
//

#ifndef DSP_SIMULATION_FFT_H
#define DSP_SIMULATION_FFT_H

#include <cmath>
#include <complex>
#include <vector>
#include "../types/concepts.h"

namespace mechdancer {
    constexpr static auto PI = 3.1415926535897932384626433832795;
    
    template<Number t = float>
    std::complex<t> omega(unsigned k, unsigned n) {
        double theta = 2 * PI * k / n;
        return {static_cast<t>(std::cos(theta)),
                static_cast<t>(std::sin(theta))};
    }
    
    template<Number t = float>
    std::complex<t> i_omega(unsigned k, unsigned n) {
        double theta = 2 * PI * k / n;
        return {static_cast<t>(std::cos(theta)),
                -static_cast<t>(std::sin(theta))};
    }
    
    enum class fft_operation { fft, ifft };
    
    template<fft_operation operation = fft_operation::fft, Number t = float>
    void fft(std::vector<std::complex<t>> &memory) {
        constexpr static auto ¶ÿ = operation == fft_operation::fft ? omega<t> : i_omega<t>;
        
        const size_t n = memory.size();
        
        // ¥Ì–Ú
        for (size_t i = 0, j = 0; i < n; ++i) {
            if (i > j) std::swap(memory[i], memory[j]);
            for (size_t l = n >> 1u; (j ^= l) < l; l >>= 1u);
        }
        // ±‰ªª
        for (size_t m = 1; m < n; m <<= 1u) {
            auto        a = memory.data(), b = a + m;
            const auto  s = n / m / 2;
            for (size_t i = 0; i < s; ++i) {
                for (size_t j = 0; j < m; ++j, ++a, ++b)
                    if (b->real() == 0 && b->imag() == 0)
                        *b = *a;
                    else {
                        const auto c = *b * ¶ÿ(s * j, n);
                        *b = *a - c;
                        *a += c;
                    }
                a = b;
                b += m;
            }
        }
    }
}

#endif //DSP_SIMULATION_FFT_H
