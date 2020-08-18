//
// Created by ydrml on 2020/8/7.
//

#ifndef DSP_SIMULATION_FFT_H
#define DSP_SIMULATION_FFT_H

#include <cmath>
#include <vector>

#include "../types/concepts.h"
#include "functions.h"

namespace mechdancer {
    /// �������任�� ��_n^k
    /// \tparam t ����ֵ��������
    /// \param k k
    /// \param n n
    /// \return ��_n^k
    template<Number t = float>
    static complex_t<t> omega(unsigned k, unsigned n) {
        double theta = 2 * PI * k / n;
        return {static_cast<t>(std::cos(theta)),
                static_cast<t>(std::sin(theta))};
    }
    
    /// ���ڷ��任�� ��_n^k
    /// \tparam t ����ֵ��������
    /// \param k k
    /// \param n n
    /// \return ��_n^k
    template<Number t = float>
    static complex_t<t> i_omega(unsigned k, unsigned n) {
        double theta = 2 * PI * k / n;
        return {static_cast<t>(std::cos(theta)),
                -static_cast<t>(std::sin(theta))};
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
        memory.resize(n, complex_t<t>::zero);
        // ����
        for (size_t i = 0, j = 0; i < n; ++i) {
            if (i > j) std::swap(memory[i], memory[j]);
            for (size_t l = n >> 1u; (j ^= l) < l; l >>= 1u);
        }
        // �任
        for (size_t m = 1; m < n; m <<= 1u) {
            auto        a = memory.data(), b = a + m;
            const auto  s = n / m / 2;
            for (size_t i = 0; i < s; ++i) {
                for (size_t j = 0; j < m; ++j, ++a, ++b)
                    if (b->re == 0 && b->im == 0)
                        *b = *a;
                    else {
                        const auto c = *b * ��(s * j, n);
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
}

#endif //DSP_SIMULATION_FFT_H
