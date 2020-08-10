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
    /// �ź�����
    /// \tparam _value_t ��������
    /// \tparam _frequency_t Ƶ�����ͣ�frequency_t��
    /// \tparam time_t ʱ�����ͣ�std::chrono::duration��
    template<class _value_t, Frequency _frequency_t, Time time_t>
    struct signal_t {
        using value_t = _value_t;
        
        std::vector<value_t> values;
        _frequency_t         sampling_frequency;
        time_t               begin_time;
        
        /// �ı�������ز���
        /// \tparam new_frequency_t �²���Ƶ������
        /// \tparam new_signal_t ���ź�����
        /// \param new_fs �²�����
        /// \param _times �����ʣ�����Խ�ߣ�����Խ��
        /// \return ���ź�
        template<Number times_t, Frequency new_frequency_t, Signal new_signal_t = signal_t<_value_t, new_frequency_t, time_t>>
        new_signal_t resample(new_frequency_t new_fs, times_t _times) const {
            using complex_t = std::complex<_value_t>;
            
            // �����Ƶ�ʣ�����Ƶ�����Ƶ����ͬ��ֱ�ӷ���
            auto old_fs = sampling_frequency.template cast_to<new_frequency_t>();
            if (old_fs == new_fs)
                return new_signal_t{
                    .values = values,
                    .sampling_frequency = new_fs,
                    .begin_time = begin_time,
                };
            // �ڻ� 2 �����¼���ʵ���ز������ʺͽ�����������������С��Ҫ����߱���
            auto times    = enlarge_to_2_power(values.size() * _times) / enlarge_to_2_power(values.size());
            auto interval = static_cast<size_t>(static_cast<double>(old_fs.value * times) / new_fs.value + .5);
            if (interval == 0u)
                throw std::invalid_argument("processing times is too little");
            // ���Ƿ��������������
            if (times > 1) {
                // ���� FFT ������
                auto spectrum = std::vector<complex_t>(values.size());
                std::transform(values.begin(), values.end(), spectrum.begin(), [](value_t x) { return complex_t{x, 0}; });
                
                fft(spectrum);
                auto size = spectrum.size();
                spectrum.resize(size * times, complex_t{});
                std::copy_n(spectrum.begin() + size / 2, size / 2, spectrum.end() - size / 2);
                std::fill(spectrum.begin() + size / 2, spectrum.end() - size, complex_t{});
                ifft(spectrum);
                // �ٽ�������Ŀ�������
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
                // ��ԭ�ź�ֱ�ӽ�����
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
    
    /// ����յ�ʵ�ź�
    /// \tparam _value_t ��������
    /// \tparam frequency_t Ƶ������
    /// \tparam time_t ʱ������
    /// \param size �źų���
    /// \param frequency ����Ƶ��
    /// \param time ��ʼʱ��
    /// \return ʵ�źŶ���
    template<class _value_t = float, Frequency frequency_t, Time time_t>
    auto real_signal_of(size_t size, frequency_t frequency, time_t time) {
        return signal_t<_value_t, frequency_t, time_t>{
            .values = std::vector<_value_t>(size),
            .sampling_frequency = frequency,
            .begin_time = time,
        };
    }
    
    /// ����յĸ��ź�
    /// \tparam _value_t ������������
    /// \tparam frequency_t Ƶ������
    /// \tparam time_t ʱ������
    /// \param size �źų���
    /// \param frequency ����Ƶ��
    /// \param time ��ʼʱ��
    /// \return ���źŶ���
    template<class _value_t = float, Frequency frequency_t, Time time_t>
    auto complex_signal_of(size_t size, frequency_t frequency, time_t time) {
        return real_signal_of<std::complex<_value_t>>(size, frequency, time);
    }
}

#endif // DSP_SIMULATION_SIGNAL_T_HPP
