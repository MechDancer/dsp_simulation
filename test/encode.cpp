#include "../functions/builders.h"
#include "../functions/process_real.h"
#include "../functions/script_builder.hh"

using namespace mechdancer;
using namespace std::chrono_literals;

// ������Գ���ͨ�����༤��������β��С�Ĳο��ź�

int main() {
    auto script_builder = script_builder_t("data");
    // ���ش���ϵͳ
    auto transceiver = load("../2048_1M_0.txt", 1_MHz, 0s);
    transceiver.values.erase(transceiver.values.begin() + 1600, transceiver.values.end());
    // ���켤���ź�
    auto excitation = sample(535, [](auto time) {
        using namespace std::chrono;
        auto t = duration_cast<floating_seconds>(time).count();
        auto x = std::sin(2 * PI * 40e3 * t);
        return t < .4e-3 ? x : -x;
    }, 1_MHz, 0s);
    // ����ο��ź�
    SAVE_SIGNAL(script_builder.save("signal"), convolution(excitation, transceiver));
    return 0;
}
