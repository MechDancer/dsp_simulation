#include <iostream>
#include <complex>
#include <algorithm>
#include <filesystem>

#include "functions/builders.h"
#include "functions/processing.h"
#include "types/noise.h"

int main() {
    // region ׼������
    using namespace mechdancer;
    using namespace std::chrono;
    using namespace std::chrono_literals;
    
    std::filesystem::remove_all("../data");
    std::filesystem::create_directory("../data");
    // endregion
    // region ����
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
    add_noise_measured(received, -2_db);
    SAVE_SIGNAL("../data/received.txt", received);
    // endregion
    // region �����㷨
    // endregion
    return 0;
}
