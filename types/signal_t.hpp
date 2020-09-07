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
        
        /// ת���ź�����
        /// \tparam __value_t �µ�ֵ����
        /// \tparam __frequency_t �µ�Ƶ������
        /// \tparam __time_t �µ�ʱ������
        /// \tparam converter_t ת��������
        /// \param new_size �µ���ɢ�źų���
        /// \param converter ת����
        /// \return Ŀ�����͵��ź�
        template<class __value_t = _value_t, Frequency __frequency_t = _frequency_t, Time __time_t = _time_t, class converter_t>
        auto cast(size_t new_size, converter_t converter) const {
            auto result = signal_t<__value_t, __frequency_t, __time_t>{
                .values = std::vector<__value_t>(new_size > 0 ? new_size : values.size(), __value_t{}),
                .sampling_frequency = sampling_frequency.template cast_to<__frequency_t>(),
                .begin_time = std::chrono::duration_cast<__time_t>(begin_time),
            };
            if constexpr (std::is_same_v<converter_t, nullptr_t>)
                std::copy(values.begin(), values.end(), result.values.begin());
            else
                std::transform(values.begin(), values.end(), result.values.begin(), converter);
            return result;
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
    auto signal_of(size_t size, frequency_t frequency, time_t time) {
        return signal_t<_value_t, frequency_t, time_t>{
            .values = std::vector<_value_t>(size),
            .sampling_frequency = frequency,
            .begin_time = time,
        };
    }
    
    /// ʵ�ź���Ϊʵ�����ɸ��ź�
    /// \tparam t �ź�����
    /// \param signal ʵ�ź�
    /// \return ���ź�
    template<RealSignal t, Number value_t = typename t::value_t>
    auto complex(t const &signal) {
        return signal.template cast<complex_t<value_t>>(0, [](auto x) -> complex_t<value_t> { return x; });
    }
    
    /// ���ź���ȡʵ���ķ�ʽת��Ϊʵ�ź�
    /// \tparam _value_t ����ֵ����
    /// \tparam t �ź�����
    /// \param signal ���ź�
    /// \return ���ź�ʵ����ɵ�ʵ�ź�
    template<ComplexSignal t, class value_t = typename t::value_t::value_t>
    auto real(t const &signal) {
        return signal.template cast<value_t>(0, [](auto z) -> value_t { return z.re; });
    }
    
    /// ���ź���ȡģ�ķ�ʽת��Ϊʵ�ź�
    /// \tparam _value_t ����ֵ����
    /// \tparam t �ź�����
    /// \param signal ���ź�
    /// \return ���ź�ģ��ɵ�ʵ�ź�
    template<ComplexSignal t, class value_t = typename t::value_t::value_t>
    auto abs(t const &signal) {
        return signal.template cast<value_t>(0, [](auto z) -> value_t { return z.norm(); });
    }
}

#endif // DSP_SIMULATION_SIGNAL_T_HPP
