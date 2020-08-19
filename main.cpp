#include <iostream>
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <cmath>

#include "functions/builders.h"
#include "functions/process_real.h"
#include "functions/process_complex.h"
#include "types/noise.h"

using namespace mechdancer;

constexpr static auto PATH = "../data"; // 数据保存路径
template<class t> requires Number<t>
constexpr auto sound_speed(t temperature) {
    return 20.048 * std::sqrt(temperature + 273.15);
}

/// 基于互相关的时延估计
/// \tparam t 参考信号类型
/// \tparam u 接收信号类型
template<class t, class u> requires RealSignal<t> && RealSignal<u>
[[maybe_unused]] void simulation(t, u);

int main() {
    // region 准备环境
    using namespace std::chrono;
    using namespace std::chrono_literals;
    
    std::filesystem::remove_all("../data");
    std::filesystem::create_directory("../data");
    // endregion
    // region 参数
    constexpr static auto MAIN_FS = 1_MHz;  // 仿真采样率（取决于测量脉冲响应的采样率）
    constexpr static auto DISTANCE = 3;     // 实际距离（米）
    constexpr static auto TEMPERATURE = 20; // 气温（℃）
    // endregion
    // region 信源信道仿真
    auto transceiver = load("../2048_1M_0.txt", MAIN_FS, 0s);
    auto excitation = sample(1'000, chirp(39_kHz, 61_kHz, 1ms), MAIN_FS, 0s);
    auto excitation_noise = excitation;
    add_noise_measured(excitation_noise, -3_db);
    auto reference = convolution(transceiver, excitation_noise); // 构造接收信号
    // 加噪
    auto DELAY = static_cast<size_t>(std::lround(DISTANCE * MAIN_FS.cast_to<Hz_t>().value / sound_speed(TEMPERATURE)));
    auto received = real_signal_of(DELAY + 3 * reference.values.size(), reference.sampling_frequency, 0s);
    std::copy(reference.values.begin(), reference.values.end(), received.values.begin() + DELAY);
    add_noise(received, sigma_noise(received, -12_db));
    // endregion
    // region 接收机仿真
    // 降低采样率重采样，模拟低采样率的嵌入式处理器
    auto sampling_float = resample(received, 600_kHz, 64);
    // 降低精度到 12 位，均值 1600
    using sample_t = unsigned short;
    auto sampling = real_signal_of<sample_t>(sampling_float.values.size(), sampling_float.sampling_frequency, sampling_float.begin_time);
    std::transform(sampling_float.values.begin(), sampling_float.values.end(), sampling.values.begin(),
                   [](auto x) { return static_cast<sample_t>(x / 15 + 1600); });
    SAVE_SIGNAL_AUTO({ PATH }, sampling)
    std::cout << "延迟点数 = " << DELAY * .6 << std::endl;
    { // 测试直接使用各种方法
        auto temp = resample(reference, sampling.sampling_frequency, 64);
        // 互相关
        auto test1 = correlation(temp, sampling_float);
        // 倒谱
        auto test2 = rceps(sampling + temp) - rceps(sampling - temp);
        SAVE_SIGNAL_AUTO({ PATH }, test1)
        SAVE_SIGNAL_AUTO({ PATH }, test2)
    }
    // endregion
    return 0;
}

template<class t, class u> requires RealSignal<t> && RealSignal<u>
[[maybe_unused]] void simulation(t reference, u sampling) {
    constexpr static auto FRAME_SIZE = 1024; // 帧长度
    // region 重叠分帧，模拟内存不足的嵌入式系统
    auto frames = std::vector<u>();
    using Tt = typename u::time_t;
    auto fs = sampling.sampling_frequency;
    for (auto i = 0; i < sampling.values.size(); i += FRAME_SIZE) {
        frames.push_back({decltype(sampling.values)(FRAME_SIZE), fs, fs.template duration_of<Tt>(i)});
        if (i + FRAME_SIZE < sampling.values.size())
            std::copy_n(sampling.values.begin() + i, FRAME_SIZE, frames.back().values.begin());
        else
            std::copy(sampling.values.begin() + i, sampling.values.end(), frames.back().values.begin());
    }
    // 分帧保存到文件
    std::ofstream file("../data/frames.txt");
    for (size_t i = 0; i < FRAME_SIZE; ++i) {
        for (auto const &frame : frames)
            file << frame.values[i] << '\t';
        file << std::endl;
    }
    // endregion
    // region 算法仿真
    auto buffer = real_signal_of<int>(768, 150_kHz, floating_seconds(0));
    auto reference150kHz = resample(reference, 150_kHz, 4);
    reference150kHz.values.erase(reference150kHz.values.begin() + 257, reference150kHz.values.end());
    auto result = std::vector<common_type<decltype(buffer), decltype(reference150kHz)>>();
    for (auto i = 0; i < frames.size(); ++i) {
        const auto &frame = frames[i];
        // 降采样
        auto spectrum = fft<typename decltype(buffer)::value_t>(frame);
        auto middle = spectrum.values[spectrum.values.size() / 2];
        spectrum.values.erase(spectrum.values.begin() + FRAME_SIZE / 8, spectrum.values.end() - FRAME_SIZE / 8);
        spectrum.values[spectrum.values.size() / 2] = middle;
        spectrum.sampling_frequency = 150_kHz;
        // 带通滤波 根本没必要带通滤波，通带之外本来就没有响应，滤狠了反为不美
        // bandpass(spectrum, 30_kHz, 70_kHz);
        ifft(spectrum.values);
        auto part = real(spectrum);
        // 保存为叠帧
        std::move(buffer.values.begin() + 256, buffer.values.end(), buffer.values.begin());
        std::copy(part.values.begin(), part.values.end(), buffer.values.begin() + 512);
        if (i % 2 || i < 2) continue;
        // 计算互相关
        auto corr = correlation<correlation_mode::noise_reduction>(reference150kHz, buffer);
        result.push_back(corr);
    }
    // endregion
    // region 生成结论
    auto sum = real_signal_of(512 * (result.size() + 1), 150_kHz, floating_seconds(0));
    for (auto i = 0; i < result.size(); ++i)
        for (auto j = 0; j < result[i].values.size(); ++j)
            sum.values[i * 512 + j] += result[i].values[j];
    SAVE_SIGNAL_AUTO({ PATH }, sum)
    auto max_i = 0;
    for (auto i = 1; i < sum.values.size(); ++i)
        if (sum.values[i] > sum.values[max_i])
            max_i = i;
    std::cout << "估计点数 = " << (max_i - 256) * 4 << std::endl;
    // endregion
}
