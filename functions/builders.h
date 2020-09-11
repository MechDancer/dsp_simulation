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
    /// 从啁啾信号采样
    /// \tparam value_t 采样值类型
    /// \param f0 起始频率
    /// \param k 变频率
    /// \param t 时间
    /// \return 采样值
    template<Number value_t>
    static value_t sample_chirp(float f0, float k, floating_seconds t) {
        return static_cast<value_t>(std::sin(2 * PI * (f0 + k * t.count()) * t.count()));
    }
    
    /// 构造啁啾构造函数
    /// \tparam value_t 值类型
    /// \tparam f_t 频率类型
    /// \tparam t_t 时间类型
    /// \param f0 起始频率
    /// \param f1 终止频率
    /// \param time 时间
    /// \return 啁啾采样海曙
    template<Number value_t = float, Frequency f_t, Time t_t>
    auto chirp(f_t f0, f_t f1, t_t time) {
        using namespace std::placeholders;
        
        auto f0_Hz = f0.template cast_to<Hz_t>().value;
        auto k = (f1.template cast_to<Hz_t>().value - f0_Hz)
                 / floating_seconds(time).count()
                 / 2;
        
        return std::bind(sample_chirp<value_t>, f0_Hz, k, _1);
    }
    
    template<Number value_t = float, Frequency f_t>
    auto sin(f_t fs) {
        return [fs](auto t) { return static_cast<value_t>(std::sin(2 * PI * fs.template cast_to<Hz_t>().value * floating_seconds(t).count())); };
    }
    
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
        auto result = signal_of<value_t>(size, fs, t0);
        auto dt = floating_seconds(1) / fs.template cast_to<Hz_t>().value;
        auto t = floating_seconds(t0);
        for (auto &x : result.values) x = origin(std::exchange(t, t + dt));
        return result;
    }
    
    /// 按基本分隔符从 ASCII 文件加载时域信号
    /// \tparam value_t 数据类型
    /// \tparam frequency_t 频率类型
    /// \tparam time_t 时间类型
    /// \param file_name 文件名
    /// \param fs 采样频率
    /// \param t0 起始时间
    /// \return 离散信号
    template<class value_t = float, Frequency frequency_t, Time time_t>
    auto load(std::string const &file_name, frequency_t fs, time_t t0 = time_t::zero) {
        auto file = std::ifstream(file_name);
        auto result = signal_of(0, fs, t0);
        value_t value;
        while (file >> value) result.values.push_back(value);
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
    
    #define SAVE_SIGNAL_FORMAT(PATH, S, TF) save(PATH, S, [](std::ofstream &file, typename decltype(S)::value_t x) { file << TF; })
    #define SAVE_SIGNAL_TF(PATH, S, TF) SAVE_SIGNAL_FORMAT(PATH, S, (TF) << std::endl)
    #define SAVE_SIGNAL(PATH, S) SAVE_SIGNAL_TF(PATH, S, x)
    #define SAVE_SIGNAL_AUTO(SB, S) SAVE_SIGNAL(SB.save(#S), S)
}

#endif // DSP_SIMULATION_BUILDERS_H
