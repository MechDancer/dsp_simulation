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

struct peak_t { size_t index; float value; };

int main() {
    // region 准备环境
    using namespace std::chrono;
    using namespace std::chrono_literals;
    script_builder_t script_builder("data");
    // endregion
    std::vector<std::vector<unsigned short>> slices;
    std::ofstream("../data/README.md")
        << "# 说明\n"
           "- 信号来源：107.BIN\n"
           "- 信号说明：13 米"
           "- 算法：噪声抑制白化互相关";
    { // 加载接收信号
        constexpr static auto path = "../107.BIN";
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
    transceiver.values.erase(transceiver.values.begin() + 1600, transceiver.values.end());
    // auto excitation = sample(1'000, chirp(38_kHz, 42_kHz, 1ms), MAIN_FS, 0s);
    auto excitation = sample(1000, [](auto time) {
        using namespace std::chrono;
        auto t = duration_cast<floating_seconds>(time);
        auto x = chirp(38_kHz, 42_kHz, .8ms);
        return t < .8ms ? x(t) : t < .91ms ? -x(t - .8ms) : 0;
    }, MAIN_FS, 0s);
    auto reference = convolution(transceiver, excitation); // 构造接收信号
    SAVE_SIGNAL_AUTO(script_builder, reference);
    auto self_corr = correlation<correlation_mode::basic>(reference, reference);
    SAVE_SIGNAL_AUTO(script_builder, self_corr);
    { // 构造发送信号
        auto sending = real_signal_of<unsigned short>(excitation.values.size(), excitation.sampling_frequency, excitation.begin_time);
        std::transform(excitation.values.begin(), excitation.values.end(), sending.values.begin(),
                       [](auto x) { return static_cast<unsigned short>(std::round(x * 2000 + 2048)); });
        sending.values.push_back(2048);
        SAVE_SIGNAL_AUTO(script_builder, sending);
    }
    {
        std::ofstream file(script_builder.save("ref_signal.c"));
        auto values = std::vector<float>();
        for (auto i = 0; i < reference.values.size(); i += 6)
            values.push_back(reference.values[i]);
        file << "#include \"correlation.h\"" << std::endl
             << std::endl
             << "const complex_value_t ref_signal[" << values.size() << "] = {" << std::endl;
        for (auto x : values)
            file << '\t' << x << ',' << std::endl;
        file << "};" << std::endl;
    }
    
    std::vector<std::thread> tasks;
    std::vector<peak_t> result(slices.size());
    
    std::mutex mutex;
    for (auto i = 0; i < slices.size(); ++i)
        tasks.emplace_back([&, i] {
            signal_t<unsigned short, decltype(MAIN_FS), floating_seconds>
                received{
                .values = std::move(slices[i]),
                .sampling_frequency = MAIN_FS,
                .begin_time = floating_seconds(0),
            };
            auto size = enlarge_to_2_power(std::max(reference.values.size(), received.values.size()));
            auto R = complex(reference);
            auto S = complex<decltype(received), float>(received);
            R.values.resize(size, R.values.back());
            S.values.resize(size, S.values.back());
            fft(R.values);
            fft(S.values);
            {
                auto p = R.values.begin() + size * (36e3f / 1e6f);
                auto q = S.values.begin() + size * (36e3f / 1e6f);
                auto e = S.values.begin() + size * (44e3f / 1e6f);
                std::fill(S.values.begin(), q, complex_t<float>::zero);
                std::fill(e, S.values.end(), complex_t<float>::zero);
                do {
                    if (p->is_zero())
                        *q = *p;
                    else if (!q->is_zero())
                        *q *= p->conjugate() / std::sqrt(p->norm()) / q->norm();
                    ++p;
                } while (++q < e);
            }
            ifft(S.values);
            S.values.erase(S.values.begin() + received.values.size(), S.values.end());
            auto spectrum = mechdancer::abs(S);
            {
                result[i] = {0, 0};
                // 不搜索超前部分
                // 搜索滞后部分中的完整信号
                auto end = spectrum.values.end() - reference.values.size() + 1;
                // 找到最大值
                auto max = spectrum.values.begin();
                for (auto p = max + 1; p < end; ++p)
                    if (*p > *max) max = p;
                // 找到首次达到最大值 40% 的位置
                auto p = max > spectrum.values.begin() + 6000 ? max - 6000 : spectrum.values.begin();
                result[i].value = *max * .26f;
                while (*p++ < result[i].value);
                // 找到这个位置后的第一个极大值
                while (*p > result[i].value)
                    result[i].value = *p++;
                *p = 0;
                result[i].index = p - spectrum.values.begin();
            }
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
    result.erase(result.end() - 1);
    for (auto j : result)
        file << j.index << '\t' << j.value << '\t' << j.index * 343e-6 << std::endl;
    return 0;
}


