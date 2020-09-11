#include <iostream>

#include "../functions/builders.h"
#include "../functions/process_real.h"
#include "../functions/script_builder.hh"

using namespace mechdancer;
using namespace std::chrono;

/// 为多探头调制
/// \tparam _signal_t 信号类型
/// \param base 待调信号
/// \return 已调信号
template<class _signal_t> requires RealSignal<_signal_t>
auto modulate(_signal_t base) {
    auto excitation31 = base * sample(base.values.size(), sin(29_kHz), 1_MHz, 0_sf);
    auto excitation40 = base * sample(base.values.size(), sin(34_kHz), 1_MHz, 0_sf);
    bandpass(excitation31, 29_kHz, 33_kHz);
    bandpass(excitation40, 38_kHz, 42_kHz);
    return excitation31 + excitation40;
}

/// 为多探头解调
/// \tparam _signal_t 信号类型
/// \param received 已调信号
/// \return 解调信号
template<class _signal_t> requires RealSignal<_signal_t>
auto demodulate(_signal_t received) {
    auto received31 = received;
    auto received40 = received;
    bandpass(received31, 29_kHz, 33_kHz);
    bandpass(received40, 38_kHz, 42_kHz);
    auto temp1 = received31 * sample(received31.values.size(), sin(29_kHz), 1_MHz, 0_sf);
    auto temp2 = received40 * sample(received40.values.size(), sin(34_kHz), 1_MHz, 0_sf);
    bandpass(temp1, 0_kHz, 14.5_kHz);
    bandpass(temp2, 0_kHz, 17.0_kHz);
    return temp1 + temp2;
}

int main() {
    auto script_builder = script_builder_t("data");
    
    // 构造信道
    auto transceiver = load("../31+40_2048_1M.txt", 1_MHz, floating_seconds(0));
    auto multi_path = signal_of(800, 1_MHz, 0_sf);
    multi_path.values.front() = 1;
    multi_path.values[400] = 2;
    multi_path.values.back() = 3;
    
    // 单边带调制
    auto base = sample(4096, chirp(7_kHz, -7_kHz, 4096us), 1_MHz, 0_sf);
    auto excitation = modulate(base);
    
    // 解调
    auto reference = demodulate(convolution(excitation, transceiver));
    auto received = convolution(excitation, convolution(transceiver, multi_path));
    auto recovered = demodulate(received);
    
    auto size = 16384;
    base.values.resize(size, 0);
    recovered.values.resize(size, 0);
    reference.values.resize(size, 0);
    
    SAVE_SIGNAL_AUTO(script_builder, base);
    SAVE_SIGNAL_AUTO(script_builder, excitation);
    SAVE_SIGNAL_AUTO(script_builder, received);
    SAVE_SIGNAL_AUTO(script_builder, reference);
    SAVE_SIGNAL_AUTO(script_builder, recovered);
    
    auto extracting = 64;
    auto order = best_order(16384us, Hz_t{1e6f / extracting}, 4096us, -14_kHz);
    std::cout << "the best order = " << order << std::endl;
    std::cout << "transformation ratio = " << std::abs(std::sin(PI / 2 * (order - 1))) << std::endl;
    
    size /= extracting;
    auto reference16kHz = signal_of(size, Hz_t{1e6f / extracting}, 0_sf);
    auto recovered16kHz = signal_of(size, Hz_t{1e6f / extracting}, 0_sf);
    for (auto i = 0; i < size; ++i) {
        reference16kHz.values[i] = reference.values[i * extracting];
        recovered16kHz.values[i] = recovered.values[i * extracting];
    }
    SAVE_SIGNAL_AUTO(script_builder, reference16kHz);
    SAVE_SIGNAL_AUTO(script_builder, recovered16kHz);
    
    return 0;
}
