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
    
    /// ���������ֵ�ķֱ�����ֵ
    constexpr db_t<> operator ""_db(long double db) {
        return {static_cast<float>(db)};
    }
    
    /// ���������ֵ�ķֱ�����ֵ
    constexpr db_t<> operator ""_db(unsigned long long db) {
        return {static_cast<float>(db)};
    }
    
    /// ����ʵ�ź�����
    /// \tparam sample_t �ź�����
    /// \param signal �ź�
    /// \return ����ֵ
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
    
    /// ���źż��ϸ�˹������
    /// \tparam t ʵ�ź�����
    /// \tparam sigma_t ��׼��ֵ����
    /// \param signal ʵ�ź�
    /// \param sigma ������׼��
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
    
    /// ���źż��ϸ�˹������
    /// \tparam snr_t �����ֵ����
    /// \tparam t ʵ�ź�����
    /// \tparam _value_t �ź���������
    /// \param signal ʵ�ź�
    /// \param snr �������ֵ
    template<RealSignal t, Number snr_t>
    void add_noise_measured(t &signal, snr_t snr) {
        add_noise(signal, sigma_noise(signal, snr));
    }
    
    /// ���źż��ϸ�˹������
    /// \tparam snr_t �����ֵ����
    /// \tparam t ʵ�ź�����
    /// \tparam _value_t �ź���������
    /// \param signal ʵ�ź�
    /// \param snr �������ֵ
    template<RealSignal t, Number snr_t>
    void add_noise_measured(t &signal, db_t<snr_t> snr) {
        add_noise(signal, snr.to_ratio());
    }
}

#endif // SIMULATION_NOISE_H
