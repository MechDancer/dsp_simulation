#include <iostream>

#include "../functions/builders.h"
#include "../functions/process_real.h"
#include "../functions/script_builder.hh"

using namespace mechdancer;
using namespace std::chrono;

template<class _signal_t> requires RealSignal<_signal_t>
auto modulate(_signal_t base) {
    auto loader29 = sample(base.values.size(), [](auto t) { return std::sin(2 * PI * 29e3f * floating_seconds(t).count()); }, 1_MHz, 0_sf);
    auto loader34 = sample(base.values.size(), [](auto t) { return std::sin(2 * PI * 34e3f * floating_seconds(t).count()); }, 1_MHz, 0_sf);
    auto excitation31 = base * loader29;
    auto excitation40 = base * loader34;
    bandpass(excitation31, 29_kHz, 33_kHz);
    bandpass(excitation40, 38_kHz, 42_kHz);
    return excitation31 + excitation40;
}

template<class _signal_t> requires RealSignal<_signal_t>
auto demodulate(_signal_t received) {
    auto received31 = received;
    auto received40 = received;
    bandpass(received31, 29_kHz, 33_kHz);
    bandpass(received40, 38_kHz, 42_kHz);
    auto loader29 = sample(received31.values.size(), [](auto t) { return std::sin(2 * PI * 29e3f * floating_seconds(t).count()); }, 1_MHz, 0_sf);
    auto loader34 = sample(received40.values.size(), [](auto t) { return std::sin(2 * PI * 34e3f * floating_seconds(t).count()); }, 1_MHz, 0_sf);
    auto temp1 = received31 * loader29;
    auto temp2 = received40 * loader34;
    bandpass(temp1, 0_kHz, 14_kHz);
    bandpass(temp2, 0_kHz, 17_kHz);
    return temp1 + temp2;
}

int main() {
    auto script_builder = script_builder_t("data");
    
    auto extracting = 64;
    auto order = best_order(16384us, Hz_t{1e6f / extracting}, 8192us, -14_kHz);
    std::cout << "the best order = " << order << std::endl;
    std::cout << "transformation ratio = " << std::abs(std::sin(PI / 2 * (order - 1))) << std::endl;
    
    // 构造信道
    auto transceiver = load("../31+40_2048_1M.txt", 1_MHz, floating_seconds(0));
    auto multi_path = signal_of(800, 1_MHz, 0_sf);
    multi_path.values.front() = 1;
    multi_path.values[400] = 2;
    multi_path.values.back() = 3;
    
    // 单边带调制
    auto base = sample(8192, chirp(7_kHz, -7_kHz, 8192us), 1_MHz, 0_sf);
    auto excitation = modulate(base);
    
    // 解调
    auto reference = demodulate(convolution(excitation, transceiver));
    auto received = convolution(excitation, convolution(transceiver, multi_path));
    auto recovered = demodulate(received);
    
    base.values.resize(16384, 0);
    recovered.values.resize(16384, 0);
    reference.values.resize(16384, 0);
    
    SAVE_SIGNAL_AUTO(script_builder, base);
    SAVE_SIGNAL_AUTO(script_builder, excitation);
    SAVE_SIGNAL_AUTO(script_builder, received);
    SAVE_SIGNAL_AUTO(script_builder, reference);
    SAVE_SIGNAL_AUTO(script_builder, recovered);
    
    auto reference8kHz = signal_of(16384 / extracting, kHz_t{1000.0f / extracting}, 0_sf);
    auto recovered8kHz = signal_of(16384 / extracting, kHz_t{1000.0f / extracting}, 0_sf);
    for (auto i = 0; i < 16384 / extracting; ++i) {
        reference8kHz.values[i] = reference.values[i * extracting];
        recovered8kHz.values[i] = recovered.values[i * extracting];
    }
    SAVE_SIGNAL_AUTO(script_builder, reference8kHz);
    SAVE_SIGNAL_AUTO(script_builder, recovered8kHz);
    
    return 0;
}
