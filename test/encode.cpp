#include "../functions/builders.h"
#include "../functions/process_real.h"
#include "../functions/script_builder.hh"

using namespace mechdancer;
using namespace std::chrono_literals;

// 这个测试尝试通过反相激励构造拖尾较小的参考信号

int main() {
    auto script_builder = script_builder_t("data");
    // 加载传输系统
    auto transceiver = load("../2048_1M_0.txt", 1_MHz, 0s);
    transceiver.values.erase(transceiver.values.begin() + 1600, transceiver.values.end());
    // 构造激励信号
    auto excitation = sample(535, [](auto time) {
        using namespace std::chrono;
        auto t = duration_cast<floating_seconds>(time).count();
        auto x = std::sin(2 * PI * 40e3 * t);
        return t < .4e-3 ? x : -x;
    }, 1_MHz, 0s);
    // 保存参考信号
    SAVE_SIGNAL(script_builder.save("signal"), convolution(excitation, transceiver));
    return 0;
}
