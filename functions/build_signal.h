//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_BUILD_SIGNAL_H
#define DSP_SIMULATION_BUILD_SIGNAL_H

#include <string>
#include <fstream>
#include "../types/signal_t.hpp"

namespace mechdancer {
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
        using namespace std::chrono;
        
        signal_t<signal_domain::time, value_t, frequency_t, time_t> result{
            .values = std::vector<value_t>(size),
            .sampling_frequency = fs,
            .begin_time = t0,
        };
        
        auto dt = seconds(1) / fs.template cast_to<Hz_t>().value;
        auto t  = duration_cast<duration<double>>(t0);
        
        for (auto &x : result.values) {
            x = origin(t);
            t += dt;
        }
        
        return result;
    }
    
    /// ���ļ����أ�ʱ���ź�
    /// \tparam value_t ��������
    /// \tparam frequency_t Ƶ������
    /// \tparam time_t ʱ������
    /// \param file_name �ļ���
    /// \param fs ����Ƶ��
    /// \param t0 ��ʼʱ��
    /// \return ��ɢ�ź�
    template<class value_t = float, Frequency frequency_t, Time time_t>
    auto load(std::string const &file_name, frequency_t fs, time_t t0 = time_t::zero) {
        signal_t<signal_domain::time, value_t, frequency_t, time_t> result{
            .sampling_frequency = fs,
            .begin_time = t0,
        };
        
        std::ifstream file(file_name);
        value_t       value;
        
        while (file >> value)
            result.values.push_back(value);
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
    
    #define SAVE_SIGNAL(PATH, S) \
    save(PATH, S, [](std::ofstream &file, typename decltype(S)::value_t x) { file << x << std::endl; })
    
    #define SAVE_SIGNAL_TF(PATH, S, TF) \
    save(PATH, S, [](std::ofstream &file, typename decltype(S)::value_t x) { file << (TF) << std::endl; })
    
    #define SAVE_SIGNAL_FORMAT(PATH, S, TF) \
    save(PATH, S, [](std::ofstream &file, typename decltype(S)::value_t x) { file << TF; })
}

#endif // DSP_SIMULATION_BUILD_SIGNAL_H
