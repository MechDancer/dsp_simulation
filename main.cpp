#include <iostream>
#include <complex>
#include <algorithm>
#include <filesystem>

#include "functions/builders.h"
#include "functions/processing.h"
#include "types/noise.h"

int main() {
    // region 准备环境
    using namespace mechdancer;
    using namespace std::chrono;
    using namespace std::chrono_literals;
    
    std::filesystem::remove_all("../data");
    std::filesystem::create_directory("../data");
    // endregion
    constexpr static auto PATH = "../data";
    constexpr static auto MAIN_FS = 10_MHz;
    // region 仿真
    // 收发系统
    auto transceiver = load("C:\\Users\\ydrml\\Desktop\\数据\\2048_1M.txt", MAIN_FS, 0s);
    std::for_each(transceiver.values.begin(), transceiver.values.end(), [_mean = mean(transceiver.values)](auto &x) { x -= _mean; });
    // 激励信号
    auto excitation = sample(100'000, chirp(39_kHz, 61_kHz, 10ms), MAIN_FS, 0s);
    SAVE_SIGNAL_AUTO({ PATH }, excitation)
    // 参考接收
    auto reference = convolution(transceiver, excitation);
    SAVE_SIGNAL_AUTO({ PATH }, reference)
    // 加噪
    auto received = reference;
    add_noise_measured(received, -2_db);
    SAVE_SIGNAL_AUTO({ PATH }, received)
    // endregion
    // region 测试算法
    // endregion
    return 0;
}
