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

#define �����
#ifndef �����
#define ����
#endif

int main() {
    // region ׼������
    using namespace std::chrono;
    using namespace std::chrono_literals;
    script_builder_t script_builder("data");
    // endregion
    std::vector<std::vector<unsigned short>> slices;
    {
        std::ofstream read_me("../data/README.md");
        read_me << "# ˵��\n"
                   "- �ź���Դ��015.BIN\n"
                   "- �ź�˵����20 ��"
                   "- �㷨���������ư׻������";
    }
    { // ���ؽ����ź�
        constexpr static auto path = "../015.BIN";
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
    constexpr static auto MAIN_FS = 1_MHz; // ������
    // region ��Դ�ŵ�����
    auto transceiver = load("../2048_1M_0.txt", MAIN_FS, 0s);
    auto excitation = sample(1'000, chirp(28_kHz, 70_kHz, 1ms), MAIN_FS, 0s);
    auto reference = convolution(transceiver, excitation); // ��������ź�
    {  // ���췢���ź�
        auto sending = real_signal_of<unsigned short>(excitation.values.size(), excitation.sampling_frequency, excitation.begin_time);
        std::transform(excitation.values.begin(), excitation.values.end(), sending.values.begin(),
                       [](auto x) { return static_cast<unsigned short>(std::round(x * 2000 + 2048)); });
        sending.values.push_back(2048);
        SAVE_SIGNAL_AUTO(script_builder, sending);
    }
    
    std::vector<std::thread> tasks;
    std::vector<size_t> result(slices.size());
    
    std::mutex mutex;
    for (auto i = 0; i < 20; ++i)
        tasks.emplace_back([&, i] {
            signal_t<unsigned short, decltype(MAIN_FS), floating_seconds>
                received{
                .values = std::move(slices[i]),
                .sampling_frequency = MAIN_FS,
                .begin_time = floating_seconds(0),
            };
            #if defined(����)
            auto spectrum = rceps(received + reference) - rceps(received - reference);
            spectrum.values.erase(spectrum.values.begin() + spectrum.values.size() / 2,
                                  spectrum.values.end());
            #elif defined(�����)
            auto spectrum = correlation<correlation_mode::noise_reduction>(reference, received);
            auto max = 1.2e4;
            for (auto &x : spectrum.values)
                if (x <= max) x = 0;
                else max = x;
            max /= 3;
            auto max_i = 0;
            for (auto x : spectrum.values) {
                ++max_i;
                if (x > max) break;
            }
            #endif
            std::stringstream string_builder;
            string_builder << "group" << i;
            std::string name = string_builder.str();
            {
                std::lock_guard<decltype(mutex)> _(mutex);
                std::cout << name << "(" << spectrum.values.size() << ") saving" << std::endl;
                #ifdef �����
                result[i] = max_i - 1;
                #endif
                name = script_builder.save(name);
            }
            SAVE_SIGNAL(name, spectrum);
        });
    for (auto &task : tasks) task.join();
    auto file = std::ofstream(script_builder.save("result"));
    for (auto j = 0; j < result.size(); ++j) {
        std::cout << "group " << j << " = " << result[j] * 343e-6 << std::endl;
        file << result[j] * 343e-6 << std::endl;
    }
    return 0;
}


