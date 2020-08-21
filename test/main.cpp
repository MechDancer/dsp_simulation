#include <iostream>
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <thread>
#include <mutex>

#include "../functions/builders.h"
#include "../functions/process_real.h"
#include "../functions/script_builder.hh"

using namespace mechdancer;

#define 互相关
#ifndef 互相关
#define 倒谱
#endif

int main() {
    // region 准备环境
    using namespace std::chrono;
    using namespace std::chrono_literals;
    script_builder_t script_builder("data");
    // endregion
    std::vector<std::vector<unsigned short>> slices;
    std::ofstream("../data/README.md")
        << "# 说明\n"
           "- 信号来源：017.BIN\n"
           "- 信号说明：30 米"
           "- 算法：噪声抑制白化互相关";
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
            } else if (slices.back().size() < 100000)
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
    std::vector<long> result(slices.size());
    
    std::mutex mutex;
    for (auto i = 0; i < slices.size(); ++i)
        tasks.emplace_back([&, i] {
            signal_t<unsigned short, decltype(MAIN_FS), floating_seconds>
                received{
                .values = std::move(slices[i]),
                .sampling_frequency = MAIN_FS,
                .begin_time = floating_seconds(0),
            };
            #if defined(倒谱)
            auto spectrum = rceps(received + reference) - rceps(received - reference);
            auto &values = spectrum.values;
            std::fill(values.begin(), values.begin() + reference.values.size(), 0);
            #elif defined(互相关)
            auto spectrum = correlation<correlation_mode::noise_reduction>(reference, received);
            auto max = 1.2e4;
            auto p_count = 0;
            auto n_count = 0;
            auto j = 0;
            for (auto &x : spectrum.values) {
                if (x <= max * (1 + p_count / 500.0)) {
                    x = 0;
                    if (result[i] > 0) ++p_count;
                } else {
                    max = x;
                    p_count = 0;
                    ++n_count;
                    result[i] = j;
                }
                ++j;
            }
            if (n_count < 3) result[i] = 0;
            #endif
            std::stringstream string_builder;
            string_builder << "group" << i;
            std::string name = string_builder.str();
            {
                std::lock_guard<decltype(mutex)> _(mutex);
                std::cout << name << "(" << spectrum.values.size() << ") saving" << std::endl;
                name = script_builder.save(name);
            }
            SAVE_SIGNAL(name, spectrum);
        });
    for (auto &task : tasks) task.join();
    auto file = std::ofstream(script_builder.save("result"));
    for (auto j : result)
        file << j << '\t' << (j - 3046) * 343e-6 << std::endl;
    return 0;
}


