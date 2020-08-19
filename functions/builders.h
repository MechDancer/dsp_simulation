//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_BUILDERS_H
#define DSP_SIMULATION_BUILDERS_H

#include <string>
#include <fstream>
#include <functional>

#include "../types/signal_t.hpp"

namespace mechdancer {
    template<Number value_t>
    static value_t chirp_value(float f0, float k, floating_seconds t) {
        return static_cast<value_t>(std::sin(2 * PI * (f0 + k * t.count()) * t.count()));
    }
    
    template<Number value_t = float, Frequency f_t, Time t_t>
    auto chirp(f_t f0, f_t f1, t_t time) {
        using namespace std::placeholders;
        
        auto f0_Hz = f0.template cast_to<Hz_t>().value;
        auto k = (f1.template cast_to<Hz_t>().value - f0_Hz)
                 / floating_seconds(time).count()
                 / 2;
        
        return std::bind(chirp_value<value_t>, f0_Hz, k, _1);
    }
    
    /// �������źŲ���
    /// \tparam size ����
    /// \tparam value_t ��������
    /// \tparam frequency_t Ƶ������
    /// \tparam time_t ʱ������
    /// \tparam origin_t �����ź�����
    /// \param origin �����ź�
    /// \param fs ����Ƶ��
    /// \param t0 ��ʼʱ��
    /// \return ��ɢ�ź�
    template<class value_t = float, Frequency frequency_t, Time time_t, class origin_t>
    auto sample(size_t size, origin_t origin, frequency_t fs, time_t t0 = time_t::zero) {
        using float_s_t = std::chrono::duration<float>;
        
        signal_t<value_t, frequency_t, time_t> result{
            .values = std::vector<value_t>(size),
            .sampling_frequency = fs,
            .begin_time = t0,
        };
        
        auto dt = float_s_t(1) / fs.template cast_to<Hz_t>().value;
        auto t = float_s_t(t0);
        
        for (auto &x : result.values) {
            x = origin(t);
            t += dt;
        }
        
        return result;
    }
    
    /// ���ļ�����ʱ���ź�
    /// \tparam value_t ��������
    /// \tparam frequency_t Ƶ������
    /// \tparam time_t ʱ������
    /// \param file_name �ļ���
    /// \param fs ����Ƶ��
    /// \param t0 ��ʼʱ��
    /// \return ��ɢ�ź�
    template<class value_t = float, Frequency frequency_t, Time time_t>
    auto load(std::string const &file_name, frequency_t fs, time_t t0 = time_t::zero) {
        signal_t<value_t, frequency_t, time_t> result{
            .sampling_frequency = fs,
            .begin_time = t0,
        };
        
        std::ifstream file(file_name);
        value_t value;
        while (file >> value) result.values.push_back(value);
        return result;
    }
    
    /// �����źŵ��ļ�
    /// \tparam _signal_t �ź�����
    /// \param file_name �ļ���������·����
    /// \param signal �ź�
    /// \param formatter �����ʽ������
    template<Signal _signal_t, class sava_formatter_t>
    void save(std::string const &file_name, _signal_t const &signal, sava_formatter_t const &formatter) {
        for (std::ofstream file(file_name); auto x : signal.values) formatter(file, x);
    }
    
    #define SAVE_SIGNAL_FORMAT(PATH, S, TF) save(PATH, S, [](std::ofstream &file, typename decltype(S)::value_t x) { file << TF; })
    #define SAVE_SIGNAL_TF(PATH, S, TF) SAVE_SIGNAL_FORMAT(PATH, S, (TF) << std::endl)
    #define SAVE_SIGNAL(PATH, S) SAVE_SIGNAL_TF(PATH, S, x)
    #define SAVE_SIGNAL_AUTO(SB, S) SAVE_SIGNAL(SB.save(#S), S)
}

#endif // DSP_SIMULATION_BUILDERS_H
