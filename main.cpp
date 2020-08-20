#include <iostream>
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <thread>
#include <mutex>
#include <atomic>

#include "functions/builders.h"
#include "functions/process_real.h"
#include "functions/process_complex.h"
#include "types/noise.h"
#include "functions/script_builder.hh"

using namespace mechdancer;

int main() {
    // region 准备环境
    using namespace std::chrono;
    using namespace std::chrono_literals;
    script_builder_t script_builder("data");
    // endregion
    std::vector<std::vector<unsigned short>> slices;
    {
        std::ofstream read_me("../data/README.md");
        read_me << "# 说明\n"
                   "- 信号来源：016.BIN\n"
                   "- 算法：噪声抑制白化互相关";
    }
    { // 加载接收信号
        constexpr static auto path = "../016.BIN";
        auto size = std::filesystem::file_size(path);
        auto signal = std::vector<unsigned short>(size / 2);
        std::ifstream(path, std::ios_base::binary)
            .read(reinterpret_cast<char *>(signal.data()), size);
        for (auto p = signal.begin(); p < signal.end(); ++p)
            if (*p > 4096) {
                slices.emplace_back();
                slices.back().push_back(*p - 4096);
            } else
                slices.back().push_back(*p);
        std::cout << "parsed " << slices.size() << " groups of signal" << std::endl;
    }
    constexpr static auto MAIN_FS = 1_MHz; // 采样率
    // region 信源信道仿真
    auto transceiver = load("../2048_1M_0.txt", MAIN_FS, 0s);
    auto excitation = sample(1'000, chirp(28_kHz, 70_kHz, 1ms), MAIN_FS, 0s);
    auto reference = convolution(transceiver, excitation); // 构造接收信号
    {  // 构造发送信号
        auto sending = real_signal_of<unsigned short>(excitation.values.size(), excitation.sampling_frequency, excitation.begin_time);
        std::transform(excitation.values.begin(), excitation.values.end(), sending.values.begin(),
                       [](auto x) { return static_cast<unsigned short>(std::round(x * 2000 + 2048)); });
        sending.values.push_back(2048);
        SAVE_SIGNAL_AUTO(script_builder, sending);
    }
    
    std::vector<std::thread> tasks;
    std::atomic_int i{0};
    std::mutex mutex;
    for (auto &slice : slices) {
        tasks.emplace_back([&reference, &mutex, &i, &slice, &script_builder] {
            signal_t<unsigned short, decltype(MAIN_FS), floating_seconds>
                received{
                .values = std::move(slice),
                .sampling_frequency = MAIN_FS,
                .begin_time = floating_seconds(0),
            };
            auto spectrum = correlation<correlation_mode::noise_reduction>(reference, received);
            // auto spectrum = rceps(received + reference) - rceps(received - reference);
            spectrum.values.erase(spectrum.values.begin() + spectrum.values.size() / 2,
                                  spectrum.values.end());
            std::stringstream string_builder;
            string_builder << "spectrum" << i++;
            std::string name = string_builder.str();
            {
                std::lock_guard<decltype(mutex)> _(mutex);
                name = script_builder.save(name);
            }
            SAVE_SIGNAL(name, spectrum);
        });
    }
    auto j = 0;
    for (auto &task : tasks) {
        task.join();
        std::cout << "group" << j++ << " done" << std::endl;
    }
    return 0;
}


