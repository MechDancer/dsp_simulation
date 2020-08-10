#include <iostream>
#include <complex>
#include <algorithm>

#include "functions/builders.h"
#include "functions/processing.h"

using namespace mechdancer;
using namespace std::chrono;
using namespace std::chrono_literals;

int main() {
    auto t0 = load("C:\\Users\\ydrml\\Desktop\\Êý¾Ý\\2048_1M.txt", 1_MHz, 0s);
    auto f  = complex_signal_of(3000, 1_MHz, 0s);
    
    for (auto _mean = mean(t0.values); auto &x : t0.values)
        x -= _mean;
    std::transform(t0.values.begin(), t0.values.end(), f.values.begin(),
                   [=](float x) { return std::complex<float>{x, 0}; });
    
    SAVE_SIGNAL("../data/up_sampling.txt", t0.resample(1.5_MHz, 99));
    
    SAVE_SIGNAL("../data/conv.txt", convolution(t0, t0));
    fft(f.values);
    SAVE_SIGNAL_TF("../data/f.txt", f, std::abs(x));
    ifft(f.values);
    SAVE_SIGNAL_TF("../data/t0.txt", f, x.real());
    
    return 0;
}
