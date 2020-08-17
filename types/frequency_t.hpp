//
// Created by ydrml on 2020/8/6.
//

#ifndef DSP_SIMULATION_FREQUENCY_T_HPP
#define DSP_SIMULATION_FREQUENCY_T_HPP

#include <cmath>
#include <chrono>

namespace mechdancer {
    /// 频率类型
    /// \tparam _value_t 内部字段类型
    /// \tparam _ratio 单位相对于赫兹的比例
    template<class _value_t, class _ratio = std::ratio<1>>
    struct frequency_t {
        using value_t = _value_t;
        using ratio = _ratio;
        
        value_t value;
        
        template<class target_t>
        target_t cast_to() const {
            constexpr static double k = static_cast<double>(ratio::num * target_t::ratio::den) / (ratio::den * target_t::ratio::num);
            return {static_cast<typename target_t::value_t>(value * k)};
        }
        
        template<class t>
        size_t index_of(t time) const {
            constexpr static double k = static_cast<double>(ratio::num) / (ratio::den);
            return std::round(floating_seconds(time).count() * k * value);
        }
        
        template<class t>
        t duration_of(size_t n) const {
            constexpr static double k = static_cast<double>(ratio::num) / (ratio::den);
            return std::chrono::duration_cast<t>(floating_seconds(n / k / value));
        }
        
        auto operator<=>(frequency_t const &others) const = default;
    };
    
    using Hz_t = frequency_t<float>;
    using kHz_t = frequency_t<float, std::kilo>;
    using MHz_t = frequency_t<float, std::mega>;
    using GHz_t = frequency_t<float, std::giga>;
    
    constexpr Hz_t operator ""_Hz(unsigned long long value) {
        return {static_cast<float>(value)};
    }
    
    constexpr Hz_t operator ""_Hz(long double value) {
        return {static_cast<float>(value)};
    }
    
    constexpr kHz_t operator ""_kHz(unsigned long long value) {
        return {static_cast<float>(value)};
    }
    
    constexpr kHz_t operator ""_kHz(long double value) {
        return {static_cast<float>(value)};
    }
    
    constexpr MHz_t operator ""_MHz(unsigned long long value) {
        return {static_cast<float>(value)};
    }
    
    constexpr MHz_t operator ""_MHz(long double value) {
        return {static_cast<float>(value)};
    }
    
    constexpr GHz_t operator ""_GHz(unsigned long long value) {
        return {static_cast<float>(value)};
    }
    
    constexpr GHz_t operator ""_GHz(long double value) {
        return {static_cast<float>(value)};
    }
}

#endif // DSP_SIMULATION_FREQUENCY_T_HPP
