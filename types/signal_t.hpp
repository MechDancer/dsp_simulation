//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_SIGNAL_T_HPP
#define DSP_SIMULATION_SIGNAL_T_HPP

#include <chrono>
#include <vector>

#include "concepts.h"
#include "complex_t.hpp"

namespace mechdancer {
    template<class t>
    concept Signal = requires(t signal){
        typename t::value_t;
        typename t::frequency_t;
        typename t::time_t;
        signal.values;
        signal.sampling_frequency;
        signal.begin_time;
    };
    
    template<class t>
    concept RealSignal = Signal<t> && Number<typename t::value_t>;
    
    template<class t>
    concept ComplexSignal = Signal<t> && std::same_as<typename t::value_t, complex_t<typename t::value_t::value_t>>;
    
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
        return real_signal_of<complex_t<_value_t>>(size, frequency, time);
    }
    
    /// ʵ�ź���Ϊʵ�����ɸ��ź�
    /// \tparam t �ź�����
    /// \param signal ʵ�ź�
    /// \return ���ź�
    template<RealSignal t, Number value_t = typename t::value_t>
    auto complex(t const &signal) {
        auto result = complex_signal_of<value_t>(signal.values.size(), signal.sampling_frequency, signal.begin_time);
        std::transform(signal.values.begin(), signal.values.end(), result.values.begin(), [](auto x) { return complex_t<value_t>::from_real(x); });
        return result;
    }
    
    /// ���ź���ȡʵ���ķ�ʽת��Ϊʵ�ź�
    /// \tparam _value_t ����ֵ����
    /// \tparam t �ź�����
    /// \param signal ���ź�
    /// \return ���ź�ʵ����ɵ�ʵ�ź�
    template<ComplexSignal t>
    auto real(t const &signal) {
        using value_t = typename t::value_t::value_t;
        auto result = real_signal_of<value_t>(signal.values.size(), signal.sampling_frequency, signal.begin_time);
        std::transform(signal.values.begin(), signal.values.end(), result.values.begin(), [](auto z) { return z.re; });
        return result;
    }
    
    /// ���ź���ȡģ�ķ�ʽת��Ϊʵ�ź�
    /// \tparam _value_t ����ֵ����
    /// \tparam t �ź�����
    /// \param signal ���ź�
    /// \return ���ź�ģ��ɵ�ʵ�ź�
    template<ComplexSignal t, class _value_t = typename t::value_t::value_t>
    auto abs(t const &signal) {
        auto result = real_signal_of<_value_t>(signal.values.size(), signal.sampling_frequency, signal.begin_time);
        std::transform(signal.values.begin(), signal.values.end(), result.values.begin(),
                       [](auto z) { return z.norm(); });
        return result;
    }
}

#endif // DSP_SIMULATION_SIGNAL_T_HPP
