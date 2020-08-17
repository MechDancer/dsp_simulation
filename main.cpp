#include <iostream>
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <cmath>

#include "functions/builders.h"
#include "functions/process_real.h"
#include "functions/process_complex.h"
#include "types/noise.h"

int main() {
    // region ׼������
    using namespace mechdancer;
    using namespace std::chrono;
    using namespace std::chrono_literals;
    
    std::filesystem::remove_all("../data");
    std::filesystem::create_directory("../data");
    // endregion
    constexpr static auto PATH = "../data";  // ���ݱ���·��
    constexpr static auto MAIN_FS = 1_MHz;   // ��������ʣ�ȡ���ڲ���������Ӧ�Ĳ����ʣ�
    constexpr static auto DISTANCE = 3;      // ʵ�ʾ��루�ף�
    constexpr static auto TEMPERATURE = 20;  // ���£��棩
    constexpr static auto FRAME_SIZE = 1024; // ֡����
    // region ��Դ�ŵ�����
    // �շ�ϵͳ
    auto transceiver0 = load("../2048_1M_0.txt", MAIN_FS, 0s);
    auto transceiver1 = load("../2048_1M_1.txt", MAIN_FS, 0s);
    // �����źţ����룩
    auto excitation_pure = sample(1'000, chirp(39_kHz, 61_kHz, 1ms), MAIN_FS, 0s);
    auto excitation = excitation_pure;
    //    add_noise_measured(excitation, 0);
    // �ο�����
    auto reference0 = convolution(transceiver0, excitation);
    auto reference1 = convolution(transceiver1, excitation);
    // ����
    auto DELAY = static_cast<size_t>(std::lround(DISTANCE * MAIN_FS.cast_to<Hz_t>().value / (20.048 * std::sqrt(TEMPERATURE + 273.15))));
    auto received = real_signal_of(DELAY + 3 * reference0.values.size(), reference0.sampling_frequency, 0s);
    std::copy(reference0.values.begin(), reference0.values.end(), received.values.begin() + DELAY);
    //    add_noise(delay_received, sigma_noise(received, 0));
    // endregion
    // region ���ջ�����
    // ���Ͳ������ز�����ģ��Ͳ����ʵ�Ƕ��ʽ������
    auto sampling_float = resample(received, 600_kHz, 64);
    // ���;��ȵ� 12 λ����ֵ 1600
    using sample_t = unsigned short;
    auto sampling = real_signal_of<sample_t>(sampling_float.values.size(), sampling_float.sampling_frequency, sampling_float.begin_time);
    std::transform(sampling_float.values.begin(), sampling_float.values.end(), sampling.values.begin(),
                   [](auto x) { return static_cast<sample_t>(x / 15 + 1600); });
    { // ���Ի�����㷨
        auto reff = resample(reference0, 600_kHz, 64);
        auto corr = correlation(reff, sampling_float);
        std::cout << "����ط�λ�� = " << DELAY * .6 + reff.values.size() - 1 << std::endl;
        SAVE_SIGNAL_AUTO({ PATH }, reff)
        SAVE_SIGNAL_AUTO({ PATH }, sampling)
        SAVE_SIGNAL_AUTO({ PATH }, corr)
    }
    // �ص���֡��ģ���ڴ治���Ƕ��ʽϵͳ
    auto frames = std::vector<decltype(sampling)>();
    {
        using _time_t = typename decltype(sampling)::time_t;
        auto fs = sampling.sampling_frequency;
        for (auto i = 0; i < sampling.values.size(); i += FRAME_SIZE) {
            frames.push_back({std::vector<sample_t>(FRAME_SIZE), fs, fs.duration_of<_time_t>(i)});
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
    }
    // endregion
    // region �㷨����
    std::ofstream file("../data/spectrum.txt");
    auto buffer = real_signal_of(768, 150_kHz, 0s);
    auto reference150kHz = resample(reference0, 150_kHz, 4);
    reference150kHz.values.erase(reference150kHz.values.begin() + 256, reference150kHz.values.end());
    for (auto const &frame : frames) {
        constexpr static auto threshold = 1.6;
        // ������
        auto spectrum = fft<float>(frame);
        spectrum.values.erase(spectrum.values.begin() + FRAME_SIZE / 8, spectrum.values.begin() + FRAME_SIZE / 2);
        spectrum.values.erase(spectrum.values.begin() + FRAME_SIZE / 8 + 1, spectrum.values.end() - FRAME_SIZE / 8 + 1);
        spectrum.sampling_frequency = 150_kHz;
        // ��ͨ�˲�
        bandpass(spectrum, 30_kHz, 70_kHz);
        ifft(spectrum.values);
        auto part = real(spectrum);
        // ����Ϊ��֡
        std::move(buffer.values.begin() + 256, buffer.values.end(), buffer.values.begin());
        std::copy(part.values.begin(), part.values.end(), buffer.values.begin() + 512);
        // ���㻥���
        auto temp = correlation<correlation_mode::basic>(reference150kHz, buffer);
        for (auto value : temp.values) file << value << '\t';
        file << std::endl;
    }
    // endregion
    return 0;
}
