//
// Created by ydrml on 2020/8/7.
//

#ifndef DSP_SIMULATION_PROCESSING_H
#define DSP_SIMULATION_PROCESSING_H

#include <vector>
#include <numeric>

template<class t>
t mean(std::vector<t> const &values) {
    return std::accumulate(values.begin(), values.end(), t{}) / values.size();
}

#endif // DSP_SIMULATION_PROCESSING_H
