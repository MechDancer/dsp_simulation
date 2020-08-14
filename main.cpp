#include <iostream>
#include <algorithm>
#include <filesystem>

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
    constexpr static auto DISTANCE = 5;      // ʵ�ʾ��루�ף�
    constexpr static auto TEMPERATURE = 20;  // ���£��棩
    constexpr static auto FRAME_SIZE = 1024; // ֡����
    // region ��Դ�ŵ�����
    // �շ�ϵͳ
    auto transceiver = load("../2048_1M.txt", MAIN_FS, 0s);
    std::for_each(transceiver.values.begin(), transceiver.values.end(), [_mean = mean(transceiver.values)](auto &x) { x -= _mean; });
    // �����ź�
    auto excitation = sample(1'000, chirp(39_kHz, 61_kHz, 1ms), MAIN_FS, 0s);
    SAVE_SIGNAL_AUTO({ PATH }, excitation)
    // �ο�����
    auto reference = convolution(transceiver, excitation);
    SAVE_SIGNAL_AUTO({ PATH }, reference)
    // ����
    auto DELAY = floating_seconds(DISTANCE / (20.048 * std::sqrt(TEMPERATURE + 273.15)));
    auto received = real_signal_of(reference.values.size(), reference.sampling_frequency, DELAY);
    std::copy(reference.values.begin(), reference.values.end(), received.values.begin());
    auto delay_received = sum(real_signal_of(30000, MAIN_FS, 0s), received);
    add_noise(delay_received, sigma_noise(received, -5_db));
    SAVE_SIGNAL_AUTO({ PATH }, delay_received)
    // endregion
    // region ���ջ�����
    // ���Ͳ������ز�����ģ��Ͳ����ʵ�Ƕ��ʽ������
    auto sampling_float = resample(delay_received, 600_kHz, 6);
    // ���;��ȵ� 12 λ����ֵ 1600
    using sample_t = unsigned short;
    auto sampling = real_signal_of<sample_t>(sampling_float.values.size(), sampling_float.sampling_frequency, sampling_float.begin_time);
    std::transform(sampling_float.values.begin(), sampling_float.values.end(), sampling.values.begin(),
                   [](auto x) { return static_cast<sample_t>(x / 15 + 1600); });
    SAVE_SIGNAL_AUTO({ PATH }, sampling)
    auto last = 0;
    auto i = 0;
    auto ii = 0;
    std::ofstream counter_save("../data/counter.txt");
    auto counter = 0;
    for (short x : sampling.values) {
        x -= 1600;
        if ((x ^ std::exchange(last, x)) < 0) {
            if (unsigned(i - ii - 4) <= 5u)
                ++counter;
            else if (counter > 0)
                --counter;
            ii = i;
        }
        counter_save << counter << std::endl;
        ++i;
    }
    
    //    // �ص���֡��ģ���ڴ治���Ƕ��ʽϵͳ
    //    auto frames = std::vector<decltype(sampling)>();
    //    {
    //        using _time_t = typename decltype(sampling)::time_t;
    //        auto fs = sampling.sampling_frequency;
    //        for (auto i = 0; i < sampling.values.size(); i += FRAME_SIZE) {
    //            frames.push_back({std::vector<sample_t>(FRAME_SIZE), fs, fs.duration_of<_time_t>(i)});
    //            if (i + FRAME_SIZE < sampling.values.size())
    //                std::copy_n(sampling.values.begin() + i, FRAME_SIZE, frames.back().values.begin());
    //            else
    //                std::copy(sampling.values.begin() + i, sampling.values.end(), frames.back().values.begin());
    //        }
    //        // ��֡���浽�ļ�
    //        std::ofstream file("../data/frames.txt");
    //        for (size_t i = 0; i < FRAME_SIZE; ++i) {
    //            for (auto const &frame : frames)
    //                file << frame.values[i] << '\t';
    //            file << std::endl;
    //        }
    //    }
    //    // endregion
    //    // region �㷨���棨��ˮ�߲�����
    //    std::ofstream file("../data/spectrum.txt");
    //    auto buffer = real_signal_of(768, 150_kHz, 0s);
    //    auto last = 2e7;
    //    auto state = false;
    //    for (auto const &frame : frames) {
    //        // ������
    //        auto spectrum = fft<float>(frame);
    //        spectrum.values.erase(spectrum.values.begin() + FRAME_SIZE / 8, spectrum.values.end() - FRAME_SIZE / 8);
    //        spectrum.sampling_frequency = 150_kHz;
    //        // ��ͨ�˲�
    //        bandpass(spectrum, 39_kHz, 61_kHz);
    //        ifft(spectrum.values);
    //        auto part = real(spectrum);
    //        // ����Ϊ��֡
    //        std::move(buffer.values.begin() + 256, buffer.values.end(), buffer.values.begin());
    //        std::copy(part.values.begin(), part.values.end(), buffer.values.begin() + 512);
    //        // �����ź������仯
    //        auto e = energy(buffer);
    //        auto k = e / std::exchange(last, e);
    //        state = k > (state ? .6 : 3);
    //        std::cout << state << '\t' << k << '\t' << e << std::endl;
    //        // ��֡���浽�ļ�
    //        for (auto &value : buffer.values) file << value << '\t';
    //        file << std::endl;
    //    }
    //    // endregion
    return 0;
}
