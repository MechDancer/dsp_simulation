#include <iostream>
#include <complex>
#include <algorithm>
#include <filesystem>

#include "functions/builders.h"
#include "functions/processing.h"
#include "types/db_t.hpp"

using namespace mechdancer;
using namespace std::chrono;
using namespace std::chrono_literals;

int main() {
    std::filesystem::remove_all("../data");
    std::filesystem::create_directory("../data");
    
    // �շ�ϵͳ
    auto transceiver = load("C:\\Users\\ydrml\\Desktop\\����\\2048_1M.txt", 1_MHz, 0s);
    std::for_each(transceiver.values.begin(), transceiver.values.end(), [_mean = mean(transceiver.values)](auto &x) { x -= _mean; });
    // �����ź�
    auto excitation = sample(1000, chirp(39_kHz, 61_kHz, 1ms), 1_MHz, 0s);
    SAVE_SIGNAL("../data/signal.txt", excitation);
    // �ο�����
    auto reference = convolution(transceiver, excitation);
    SAVE_SIGNAL("../data/reference.txt", reference);
    // ����
    auto received = reference;
    add_noise(received, -2_db);
    SAVE_SIGNAL("../data/received.txt", received);
    
    SAVE_SIGNAL("../data/resampled.txt", reference.resample(Hz_t{1e8f / 168 / 4}, 4));
    SAVE_SIGNAL_TF("../data/hilbert.txt", hilbert(reference), std::abs(x));
    SAVE_SIGNAL("../data/xcorr.txt", xcorr(reference, reference));
    
    auto f = complex(reference);
    fft(f.values);
    SAVE_SIGNAL_TF("../data/f.txt", f, std::abs(x));
    ifft(f.values);
    SAVE_SIGNAL_TF("../data/t.txt", f, x.real());
    
    return 0;
}
