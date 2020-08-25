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

struct peak_t { size_t index; float value; };

int main() {
    // region ׼������
    using namespace std::chrono;
    using namespace std::chrono_literals;
    script_builder_t script_builder("data");
    // endregion
    std::vector<std::vector<unsigned short>> slices;
    std::ofstream("../data/README.md")
        << "# ˵��\n"
           "- �ź���Դ��016.BIN\n"
           "- �ź�˵����10 ��"
           "- �㷨���������ư׻������";
    { // ���ؽ����ź�
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
    constexpr static auto MAIN_FS = 1_MHz; // ������
    // region ��Դ�ŵ�����
    auto transceiver = load("../2048_1M_0.txt", MAIN_FS, 0s);
    auto excitation = sample(1'000, chirp(28_kHz, 70_kHz, 1ms), MAIN_FS, 0s);
    auto reference = convolution(transceiver, excitation); // ��������ź�
    SAVE_SIGNAL_AUTO(script_builder, reference);
    { // ���췢���ź�
        auto sending = real_signal_of<unsigned short>(excitation.values.size(), excitation.sampling_frequency, excitation.begin_time);
        std::transform(excitation.values.begin(), excitation.values.end(), sending.values.begin(),
                       [](auto x) { return static_cast<unsigned short>(std::round(x * 2000 + 2048)); });
        sending.values.push_back(2048);
        SAVE_SIGNAL_AUTO(script_builder, sending);
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
            #if defined(����)
            auto spectrum = rceps(received + reference) - rceps(received - reference);
            auto &values = spectrum.values;
            std::fill(values.begin(), values.begin() + reference.values.size(), 0);
            #elif defined(�����)
            auto size = enlarge_to_2_power(std::max(reference.values.size(), received.values.size()));
            auto R = complex(reference);
            auto S = complex<decltype(received), float>(received);
            R.values.resize(size, R.values.back());
            S.values.resize(size, S.values.back());
            fft(R.values);
            fft(S.values);
            {
                auto p = R.values.begin() + size * (20e3f / 1e6f);
                auto q = S.values.begin() + size * (20e3f / 1e6f);
                auto e = S.values.begin() + size * (80e3f / 1e6f);
                std::fill(S.values.begin(), q, complex_t<float>::zero);
                std::fill(e, S.values.end(), complex_t<float>::zero);
                do {
                    if (p->is_zero())
                        *q = *p;
                    else if (!q->is_zero())
                        *q *= p->conjugate();
                    ++p;
                } while (++q < e);
            }
            ifft(S.values);
            auto spectrum = mechdancer::abs(S);
            {
                auto begin = spectrum.values.begin() + reference.values.size() - 1;
                result[i] = {0, 0};
                // �ҵ����ֵ
                for (auto x : spectrum.values)
                    if (x > result[i].value) result[i].value = x;
                result[i].value /= 2;
                // �ҵ��׸��壨���ֵ��һ�룩
                auto p = begin;
                while (*p++ < result[i].value);
                // �ҵ������ļ���ֵ��ֱ��һ���µ� 5% �Ĺ�
                for (auto q = p; *q > result[i].value * .95f; ++q)
                    if (*q > result[i].value)
                        result[i].value = *q;
                // �ҵ��״δﵽ����ֵ 98% ��λ��
                while (*p++ < result[i].value * .98f);
                // �ҵ����λ�ú�ĵ�һ������ֵ
                while (*p > result[i].value)
                    result[i].value = *p++;
                result[i].index = p - begin;
            }
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
        file << j.index + reference.values.size() - 1 << '\t' << j.value << '\t' << j.index * 343e-6 << std::endl;
    return 0;
}


