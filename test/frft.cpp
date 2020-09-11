#include <iostream>
#include <filesystem>

#include "../functions/builders.h"
#include "../functions/process_real.h"
#include "../functions/script_builder.hh"

using namespace mechdancer;

int main() {
    using namespace std::chrono_literals;
    constexpr static auto MAIN_FS = 500_kHz; // ²ÉÑùÂÊ
    auto order = best_order(16.384ms, MAIN_FS, 1ms, -1_kHz);
    std::cout << "the best order = " << order << std::endl;
    std::cout << "transformation ratio = " << std::abs(std::sin(PI / 2 * (order - 1))) << std::endl;
    
    script_builder_t script_builder("data");
    auto transceiver = resample(load("../31+40_2048_1M.txt", 1_MHz, floating_seconds(0)), MAIN_FS, 2);
    auto excitation0 = sample(2000, chirp(42_kHz, 38_kHz, 4ms), MAIN_FS, floating_seconds(0));
    auto excitation1 = sample(2000, chirp(42_kHz, 38_kHz, 4ms), MAIN_FS, floating_seconds(10e-3));
    
    excitation1 = signal_of(1, MAIN_FS, floating_seconds(0)) + excitation1;
    auto excitation = excitation0 + excitation1;
    SAVE_SIGNAL_AUTO(script_builder, excitation);
    SAVE_SIGNAL_AUTO(script_builder, transceiver);
    {
        auto r0 = excitation0; // convolution(excitation0, transceiver);
        auto r1 = excitation1; // convolution(excitation1, transceiver);
        r0.values.resize(8192, 0);
        r1.values.resize(8192, 0);
        SAVE_SIGNAL_AUTO(script_builder, r0);
        SAVE_SIGNAL_AUTO(script_builder, r1);
    }
    
    return 0;
}
