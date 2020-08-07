#include <iostream>
#include <functional>

#include "functions/build_signal.h"

using namespace mechdancer;
using namespace std::chrono;
using namespace std::chrono_literals;

int main() {
    //    auto sin_10kHz = [](duration<double> t) { return std::sin(6.28 * 10e3 * t.count()); };
    //    auto signal    = sample(1024, sin_10kHz, 1_MHz, 5s);
    
    auto t0 = load<>("C:\\Users\\ydrml\\Desktop\\Êý¾Ý\\2048_1M.txt", 1_MHz, 0s);
    
    for (auto x : t0.values) std::cout << x << std::endl;
    SAVE_SIGNAL("../data/t0.txt", t0);
    
    return 0;
}
