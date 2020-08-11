#include <iostream>
#include <complex>
#include <algorithm>
#include <filesystem>

#include "functions/builders.h"
#include "functions/processing.h"

using namespace mechdancer;
using namespace std::chrono;
using namespace std::chrono_literals;

int main() {
    std::filesystem::remove_all("../data");
    std::filesystem::create_directory("../data");
    
    auto 收发系统 = load("C:\\Users\\ydrml\\Desktop\\数据\\2048_1M.txt", 1_MHz, 0s);
    
    std::for_each(收发系统.values.begin(), 收发系统.values.end(), [_mean = mean(收发系统.values)](auto &x) { x -= _mean; });
    
    SAVE_SIGNAL("../data/resampled.txt", 收发系统.resample(Hz_t{1e8f / 168 / 4}, 4));
    SAVE_SIGNAL("../data/conv.txt", convolution(收发系统, 收发系统));
    SAVE_SIGNAL_TF("../data/hilbert.txt", hilbert(收发系统), std::abs(x));
    SAVE_SIGNAL("../data/xcorr.txt", xcorr(收发系统, 收发系统));
    
    auto f = complex(收发系统);
    fft(f.values);
    SAVE_SIGNAL_TF("../data/f.txt", f, std::abs(x));
    ifft(f.values);
    SAVE_SIGNAL_TF("../data/t.txt", f, x.real());
    
    return 0;
}
