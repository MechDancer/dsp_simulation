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

constexpr static auto PATH = "../data"; // ���ݱ���·��
template<class t> requires Number<t>
constexpr auto sound_speed(t temperature) {
    return 20.048 * std::sqrt(temperature + 273.15);
}

/// ���ڻ���ص�ʱ�ӹ���
/// \tparam t �ο��ź�����
/// \tparam u �����ź�����
template<class t, class u> requires RealSignal<t> && RealSignal<u>
[[maybe_unused]] void simulation(t, u);

int main() {
    // region ׼������
    using namespace std::chrono;
    using namespace std::chrono_literals;
    
    std::filesystem::remove_all("../data");
    std::filesystem::create_directory("../data");
    // endregion
    // region ����
    constexpr static auto MAIN_FS = 1_MHz;  // ��������ʣ�ȡ���ڲ���������Ӧ�Ĳ����ʣ�
    constexpr static auto DISTANCE = 3;     // ʵ�ʾ��루�ף�
    constexpr static auto TEMPERATURE = 20; // ���£��棩
    // endregion
    // region ��Դ�ŵ�����
    auto transceiver = load("../2048_1M_0.txt", MAIN_FS, 0s);
    auto excitation = sample(1'000, chirp(39_kHz, 61_kHz, 1ms), MAIN_FS, 0s);
    auto excitation_noise = excitation;
    add_noise_measured(excitation_noise, -3_db);
    auto reference = convolution(transceiver, excitation_noise); // ��������ź�
    // ����
    auto DELAY = static_cast<size_t>(std::lround(DISTANCE * MAIN_FS.cast_to<Hz_t>().value / sound_speed(TEMPERATURE)));
    auto received = real_signal_of(DELAY + 3 * reference.values.size(), reference.sampling_frequency, 0s);
    std::copy(reference.values.begin(), reference.values.end(), received.values.begin() + DELAY);
    add_noise(received, sigma_noise(received, -12_db));
    // endregion
    // region ���ջ�����
    // ���Ͳ������ز�����ģ��Ͳ����ʵ�Ƕ��ʽ������
    auto sampling_float = resample(received, 600_kHz, 64);
    // ���;��ȵ� 12 λ����ֵ 1600
    using sample_t = unsigned short;
    auto sampling = real_signal_of<sample_t>(sampling_float.values.size(), sampling_float.sampling_frequency, sampling_float.begin_time);
    std::transform(sampling_float.values.begin(), sampling_float.values.end(), sampling.values.begin(),
                   [](auto x) { return static_cast<sample_t>(x / 15 + 1600); });
    SAVE_SIGNAL_AUTO({ PATH }, sampling)
    std::cout << "�ӳٵ��� = " << DELAY * .6 << std::endl;
    { // ����ֱ��ʹ�ø��ַ���
        auto temp = resample(reference, sampling.sampling_frequency, 64);
        // �����
        auto test1 = correlation(temp, sampling_float);
        // ����
        auto test2 = rceps(sampling + temp) - rceps(sampling - temp);
        SAVE_SIGNAL_AUTO({ PATH }, test1)
        SAVE_SIGNAL_AUTO({ PATH }, test2)
    }
    // endregion
    return 0;
}

template<class t, class u> requires RealSignal<t> && RealSignal<u>
[[maybe_unused]] void simulation(t reference, u sampling) {
    constexpr static auto FRAME_SIZE = 1024; // ֡����
    // region �ص���֡��ģ���ڴ治���Ƕ��ʽϵͳ
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
    // ��֡���浽�ļ�
    std::ofstream file("../data/frames.txt");
    for (size_t i = 0; i < FRAME_SIZE; ++i) {
        for (auto const &frame : frames)
            file << frame.values[i] << '\t';
        file << std::endl;
    }
    // endregion
    // region �㷨����
    auto buffer = real_signal_of<int>(768, 150_kHz, floating_seconds(0));
    auto reference150kHz = resample(reference, 150_kHz, 4);
    reference150kHz.values.erase(reference150kHz.values.begin() + 257, reference150kHz.values.end());
    auto result = std::vector<common_type<decltype(buffer), decltype(reference150kHz)>>();
    for (auto i = 0; i < frames.size(); ++i) {
        const auto &frame = frames[i];
        // ������
        auto spectrum = fft<typename decltype(buffer)::value_t>(frame);
        auto middle = spectrum.values[spectrum.values.size() / 2];
        spectrum.values.erase(spectrum.values.begin() + FRAME_SIZE / 8, spectrum.values.end() - FRAME_SIZE / 8);
        spectrum.values[spectrum.values.size() / 2] = middle;
        spectrum.sampling_frequency = 150_kHz;
        // ��ͨ�˲� ����û��Ҫ��ͨ�˲���ͨ��֮�Ȿ����û����Ӧ���˺��˷�Ϊ����
        // bandpass(spectrum, 30_kHz, 70_kHz);
        ifft(spectrum.values);
        auto part = real(spectrum);
        // ����Ϊ��֡
        std::move(buffer.values.begin() + 256, buffer.values.end(), buffer.values.begin());
        std::copy(part.values.begin(), part.values.end(), buffer.values.begin() + 512);
        if (i % 2 || i < 2) continue;
        // ���㻥���
        auto corr = correlation<correlation_mode::noise_reduction>(reference150kHz, buffer);
        result.push_back(corr);
    }
    // endregion
    // region ���ɽ���
    auto sum = real_signal_of(512 * (result.size() + 1), 150_kHz, floating_seconds(0));
    for (auto i = 0; i < result.size(); ++i)
        for (auto j = 0; j < result[i].values.size(); ++j)
            sum.values[i * 512 + j] += result[i].values[j];
    SAVE_SIGNAL_AUTO({ PATH }, sum)
    auto max_i = 0;
    for (auto i = 1; i < sum.values.size(); ++i)
        if (sum.values[i] > sum.values[max_i])
            max_i = i;
    std::cout << "���Ƶ��� = " << (max_i - 256) * 4 << std::endl;
    // endregion
}
