//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_BUILD_SIGNAL_H
#define DSP_SIMULATION_BUILD_SIGNAL_H

#include <string>
#include <fstream>
#include "../types/signal_t.hpp"

namespace mechdancer {
    /// 从连续信号采样
    /// \tparam size 长度
    /// \tparam value_t 数据类型
    /// \tparam frequency_t 频率类型
    /// \tparam time_t 时间类型
    /// \tparam origin_t 连续信号类型
    /// \param origin 连续信号
    /// \param fs 采样频率
    /// \param t0 起始时间
    /// \return 离散信号
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
    
    /// 从文件加载（时域）信号
    /// \tparam value_t 数据类型
    /// \tparam frequency_t 频率类型
    /// \tparam time_t 时间类型
    /// \param file_name 文件名
    /// \param fs 采样频率
    /// \param t0 起始时间
    /// \return 离散信号
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
    
    /// 保存信号到文件
    /// \tparam _signal_t 信号类型
    /// \param file_name 文件名（包括路径）
    /// \param signal 信号
    /// \param formatter 输出格式编码器
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
