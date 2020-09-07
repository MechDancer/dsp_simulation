#include <iostream>

#include "../functions/builders.h"
#include "../functions/process_real.h"
#include "../functions/script_builder.hh"

using namespace mechdancer;

template<class t0, class t1, class f0, class f1> requires Time<t0> && Time<t1> && Frequency<f0> && Frequency<f1>
auto best_alpha(t0 t, f0 fs, t1 tl, f1 df) {
    auto x = std::sqrt(floating_seconds(t).count() / fs.template cast_to<Hz_t>().value);
    return std::atan2(floating_seconds(tl).count() / x, -df.template cast_to<Hz_t>().value * x) / (PI / 2);
}

int main() {
    using namespace std::chrono_literals;
    std::cout << best_alpha(4ms, 1_MHz, 1ms, 10_kHz) << std::endl;
    
    script_builder_t script_builder("data");
    constexpr static auto MAIN_FS = 1_MHz; // ²ÉÑùÂÊ
    auto reference = sample(1000, chirp(5_kHz, 15_kHz, 1ms), 1_MHz, floating_seconds(0));
    SAVE_SIGNAL_AUTO(script_builder, reference);
    return 0;
}
