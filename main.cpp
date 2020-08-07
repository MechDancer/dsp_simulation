#include <iostream>
#include <complex>
#include <algorithm>

#include "functions/builders.h"
#include "functions/processing.h"
#include "functions/fft.h"

using namespace mechdancer;
using namespace std::chrono;
using namespace std::chrono_literals;

int main() {
    auto t0 = load("C:\\Users\\ydrml\\Desktop\\Êý¾Ý\\2048_1M.txt", 1_MHz, 0s);
    auto f  = frequency_signal_of(2048, 1_MHz, 0s);
    
    auto _mean = mean(t0.values);
    std::transform(t0.values.begin(), t0.values.end(), f.values.begin(),
                   [=](float x) { return std::complex<float>{x - _mean, 0}; });
    
    fft(f.values);
    SAVE_SIGNAL_TF("../data/f.txt", f, std::abs(x));
    ifft(f.values);
    SAVE_SIGNAL_TF("../data/t0.txt", f, x.real());
    
    return 0;
}
