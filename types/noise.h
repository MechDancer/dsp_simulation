//
// Created by ydrml on 2020/7/23.
//

#ifndef SIMULATION_NOISE_H
#define SIMULATION_NOISE_H

#include <vector>
#include <cmath>
#include <numeric>
#include <random>

#include "signal_t.hpp"

namespace mechdancer {
    template<class t = float>
    struct db_t {
        t value;
        
        [[nodiscard]] t to_ratio() const {
            return std::pow(10, value / 10);
        }
        
        [[nodiscard]] db_t operator-() const {
            return {-value};
        }
        
        auto operator<=>(db_t const &others) const = default;
    };
    
    /// 构造信噪比值的分贝字面值
    constexpr db_t<> operator ""_db(long double db) {
        return {static_cast<float>(db)};
    }
    
    /// 构造信噪比值的分贝字面值
    constexpr db_t<> operator ""_db(unsigned long long db) {
        return {static_cast<float>(db)};
    }
    
    /// 计算实信号能量
    /// \tparam sample_t 信号类型
    /// \param signal 信号
    /// \return 能量值
    template<RealSignal t, class value_t = typename t::value_t>
    value_t energy(t const &signal) {
        return std::accumulate(signal.values.begin(), signal.values.end(), value_t{},
                               [](value_t sum, value_t x) { return sum + x * x; });
    }
    
    template<RealSignal t, Number snr_t>
    typename t::value_t sigma_noise(t const &signal, snr_t snr) {
        return std::sqrt(energy(signal) / snr / signal.values.size());
    }
    
    template<RealSignal t, Number snr_t>
    typename t::value_t sigma_noise(t const &signal, db_t<snr_t> snr) {
        return sigma_noise(signal, snr.to_ratio());
    }
    
    /// 给信号加上高斯白噪声
    /// \tparam t 实信号类型
    /// \tparam sigma_t 标准差值类型
    /// \param signal 实信号
    /// \param sigma 噪声标准差
    template<RealSignal t, Number sigma_t>
    void add_noise(t &signal, sigma_t sigma) {
        if (sigma != 0) {
            using value_t = typename t::value_t;
            std::random_device                rd{};
            std::mt19937                      gen{rd()};
            std::normal_distribution<value_t> d{value_t{}, sigma};
            for (auto &x : signal.values) x += d(gen);
        }
    }
    
    /// 给信号加上高斯白噪声
    /// \tparam snr_t 信噪比值类型
    /// \tparam t 实信号类型
    /// \tparam _value_t 信号数据类型
    /// \param signal 实信号
    /// \param snr 信噪比数值
    template<RealSignal t, Number snr_t>
    void add_noise_measured(t &signal, snr_t snr) {
        add_noise(signal, sigma_noise(signal, snr));
    }
    
    /// 给信号加上高斯白噪声
    /// \tparam snr_t 信噪比值类型
    /// \tparam t 实信号类型
    /// \tparam _value_t 信号数据类型
    /// \param signal 实信号
    /// \param snr 信噪比数值
    template<RealSignal t, Number snr_t>
    void add_noise_measured(t &signal, db_t<snr_t> snr) {
        add_noise(signal, snr.to_ratio());
    }
}

#endif // SIMULATION_NOISE_H
