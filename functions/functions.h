//
// Created by ydrml on 2020/8/10.
//

#ifndef DSP_SIMULATION_FUNCTIONS_H
#define DSP_SIMULATION_FUNCTIONS_H

#include "../types/concepts.h"

namespace mechdancer {
    /// 若两个类型相同，取此类型，否则取第三个类型
    template<class a, class b, class else_then>
    using same_or = std::conditional_t<std::is_same_v<a, b>, a, else_then>;
    
    /// 两个实信号类型的和类型
    template<class t, class u> requires RealSignal<t> && RealSignal<u>
    using common_type = signal_t<
        decltype(typename t::value_t{} + typename u::value_t{}),
        same_or<typename t::frequency_t, typename u::frequency_t, Hz_t>,
        same_or<typename t::time_t, typename u::time_t, floating_seconds>>;
    
    template<class t> requires Integer<t>
    t enlarge_to_2_power(t value) {
        if ((value - 1) & value) {
            t n = 4;
            while (n < value) n <<= 1u;
            return n;
        }
        return value;
    }
    
    template<class t> requires Floating<t>
    size_t enlarge_to_2_power(t value) {
        return enlarge_to_2_power(static_cast<size_t>(value + .5));
    }
    
    template<class t, class u, class v, class f>
    void zip(t p, t q, u e, v r, f fun) {
        while (p != e) *r++ = fun(*p++, *q++);
    }
}

#endif // DSP_SIMULATION_FUNCTIONS_H
