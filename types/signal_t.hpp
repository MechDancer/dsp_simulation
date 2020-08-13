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
    /// \tparam _time_t ʱ�����ͣ�std::chrono::duration��
    template<class _value_t, Frequency _frequency_t, Time _time_t>
    struct signal_t {
        using value_t = _value_t;
        using frequency_t = _frequency_t;
        using time_t = _time_t;
        
        std::vector<value_t> values;
        _frequency_t sampling_frequency;
        _time_t begin_time;
        
      
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
    
    /// ʵ�ź���Ϊʵ�����ɸ��ź�
    /// \tparam t �ź�����
    /// \param signal ʵ�ź�
    /// \return ���ź�
    template<RealSignal t>
    auto complex(t const &signal) {
        using value_t = typename t::value_t;
        
        auto result = complex_signal_of<value_t>(signal.values.size(), signal.sampling_frequency, signal.begin_time);
        std::transform(signal.values.begin(), signal.values.end(), result.values.begin(),
                       [](value_t x) { return std::complex<value_t>{x, 0}; });
        return result;
    }
    
    /// ���ź���ȡʵ���ķ�ʽת��Ϊʵ�ź�
    /// \tparam _value_t ����ֵ����
    /// \tparam t �ź�����
    /// \param signal ���ź�
    /// \return ���ź�ʵ����ɵ�ʵ�ź�
    template<class t, class _value_t = typename t::value_t::value_type> requires ComplexSignal<t, _value_t>
    auto real(t const &signal) {
        auto result = real_signal_of<_value_t>(signal.values.size(), signal.sampling_frequency, signal.begin_time);
        std::transform(signal.values.begin(), signal.values.end(), result.values.begin(),
                       [](std::complex<_value_t> z) { return z.real(); });
        return result;
    }
    
    /// ���ź���ȡģ�ķ�ʽת��Ϊʵ�ź�
    /// \tparam _value_t ����ֵ����
    /// \tparam t �ź�����
    /// \param signal ���ź�
    /// \return ���ź�ģ��ɵ�ʵ�ź�
    template<class t, class _value_t = typename t::value_t::value_type> requires ComplexSignal<t, _value_t>
    auto abs(t const &signal) {
        auto result = real_signal_of<_value_t>(signal.values.size(), signal.sampling_frequency, signal.begin_time);
        std::transform(signal.values.begin(), signal.values.end(), result.values.begin(),
                       [](std::complex<_value_t> z) { return std::abs(z); });
        return result;
    }
}

#endif // DSP_SIMULATION_SIGNAL_T_HPP
