//
// Created by ydrml on 2020/5/21.
//

#ifndef FFT_COMPLEX_T_HPP
#define FFT_COMPLEX_T_HPP

#include <cmath>

#include "concepts.h"

namespace mechdancer {
    template<Number t = float>
    struct complex_t {
        using value_t = t;
        
        t re, im;
        
        constexpr complex_t() : re(0), im(0) {}
        
        template<Number u>
        constexpr complex_t(u re) : re(static_cast<t>(re)), im(0) {}
        
        template<class u, class v>
        constexpr complex_t(u re, v im) : re(static_cast<t>(re)), im(static_cast<t>(im)) {}
        
        template<Number u>
        constexpr static complex_t exp(u theta) {
            return {std::cos(theta), std::sin(theta)};
        }
        
        [[nodiscard]]
        t norm() const {
            return std::hypot(re, im);
        }
        
        [[nodiscard]]
        t arg() const {
            return std::atan2(im, re);
        }
        
        [[nodiscard]]
        complex_t conjugate() const {
            return {re, -im};
        }
        
        [[nodiscard]]
        complex_t normalize() const {
            auto l = norm();
            return l == 0 ? complex_t{} : complex_t{re / l, im / l};
        }
        
        [[nodiscard]]
        bool is_zero() const {
            return re == 0 && im == 0;
        }
        
        [[nodiscard]]
        complex_t operator+() const {
            return {re, im};
        }
        
        [[nodiscard]]
        complex_t operator-() const {
            return {-re, -im};
        }
        
        [[nodiscard]]
        complex_t operator+(const complex_t &others) const {
            return {re + others.re, im + others.im};
        }
        
        [[nodiscard]]
        complex_t operator-(const complex_t &others) const {
            return {re - others.re, im - others.im};
        }
        
        [[nodiscard]]
        complex_t operator*(const complex_t &others) const {
            return {re * others.re - im * others.im, re * others.im + im * others.re};
        }
        
        [[nodiscard]]
        complex_t operator/(const complex_t &others) const {
            auto k = others.re * others.re + others.im * others.im;
            return {(re * others.re + im * others.im) / k, (im * others.re - re * others.im) / k};
        }
        
        complex_t operator+=(const complex_t &others) {
            return *this = {re + others.re, im + others.im};
        }
        
        complex_t operator-=(const complex_t &others) {
            return *this = {re - others.re, im - others.im};
        }
        
        complex_t operator*=(const complex_t &others) {
            return *this = {re * others.re - im * others.im, re * others.im + im * others.re};
        }
        
        complex_t operator/=(const complex_t &others) {
            auto k = others.re * others.re + others.im * others.im;
            return *this = {(re * others.re + im * others.im) / k, (im * others.re - re * others.im) / k};
        }
        
        #define OPERATOR(WHAT, RE, IM) \
        template<Number num_t>         \
        [[nodiscard]]  auto operator WHAT(const num_t &others) const {  \
            using result_t = decltype(t{} WHAT num_t{}); \
            return complex_t<result_t>{                  \
                static_cast<result_t>(RE),               \
                static_cast<result_t>(IM),               \
            };                                           \
        }
        
        OPERATOR(+, re + others, im)
        
        OPERATOR(-, re - others, im)
        
        OPERATOR(*, re * others, im * others)
        
        OPERATOR(/, re / others, im / others)
        
        #undef OPERATOR
        
        template<Number num_t>
        complex_t &operator+=(const num_t &others) {
            re += others;
            return *this;
        }
        
        template<Number num_t>
        complex_t &operator-=(const num_t &others) {
            re -= others;
            return *this;
        }
        
        template<Number num_t>
        complex_t &operator*=(const num_t &others) {
            re *= others;
            im *= others;
            return *this;
        }
        
        template<Number num_t>
        complex_t &operator/=(const num_t &others) {
            re /= others;
            im /= others;
            return *this;
        }
    };
}

#endif //FFT_COMPLEX_T_HPP
