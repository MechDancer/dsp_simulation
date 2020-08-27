#include "../functions/builders.h"
#include "../functions/process_real.h"
#include "../functions/script_builder.hh"

using namespace mechdancer;
using namespace std::chrono_literals;

int main() {
    auto script_builder = script_builder_t("data");
    auto transceiver = load("../2048_1M_0.txt", 1_MHz, 0s);
    transceiver.values.erase(transceiver.values.begin() + 1600, transceiver.values.end());
    auto excitation0 = sample(1'000, chirp(38_kHz, 42_kHz, 1ms), 1_MHz, 0s);
    //    auto excitation1 = sample(1000, [](auto time) {
    //        using namespace std::chrono;
    //        auto t = duration_cast<floating_seconds>(time);
    //        auto x = chirp(38_kHz, 42_kHz, .8ms);
    //        return t < .8ms ? x(t) : t < .91ms ? -x(t - .8ms) : 0;
    //    }, 1_MHz, 0s);
    auto excitation1 = sample(1000, [](auto time) {
        using namespace std::chrono;
        auto t = duration_cast<floating_seconds>(time).count();
        return t < 5e-4 ? std::sin(2 * PI * 40e3 * t) : std::sin(2 * PI * 60e3 * t);
    }, 1_MHz, 0s);
    auto reference0 = convolution(transceiver, excitation0);
    auto reference1 = convolution(transceiver, excitation1);
    SAVE_SIGNAL_AUTO(script_builder, reference0);
    SAVE_SIGNAL_AUTO(script_builder, reference1);
    auto corr0 = correlation<correlation_mode::basic>(reference0, reference0);
    auto corr1 = correlation<correlation_mode::basic>(reference1, reference1);
    auto corr2 = correlation<correlation_mode::noise_reduction>(reference0, reference0);
    auto corr3 = correlation<correlation_mode::noise_reduction>(reference1, reference1);
    SAVE_SIGNAL_AUTO(script_builder, corr0);
    SAVE_SIGNAL_AUTO(script_builder, corr1);
    SAVE_SIGNAL_AUTO(script_builder, corr2);
    SAVE_SIGNAL_AUTO(script_builder, corr3);
    return 0;
}
