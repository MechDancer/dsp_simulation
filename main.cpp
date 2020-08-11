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
    
    auto �շ�ϵͳ = load("C:\\Users\\ydrml\\Desktop\\����\\2048_1M.txt", 1_MHz, 0s);
    
    std::for_each(�շ�ϵͳ.values.begin(), �շ�ϵͳ.values.end(), [_mean = mean(�շ�ϵͳ.values)](auto &x) { x -= _mean; });
    
    SAVE_SIGNAL("../data/resampled.txt", �շ�ϵͳ.resample(Hz_t{1e8f / 168 / 4}, 4));
    SAVE_SIGNAL("../data/conv.txt", convolution(�շ�ϵͳ, �շ�ϵͳ));
    SAVE_SIGNAL_TF("../data/hilbert.txt", hilbert(�շ�ϵͳ), std::abs(x));
    SAVE_SIGNAL("../data/xcorr.txt", xcorr(�շ�ϵͳ, �շ�ϵͳ));
    
    auto f = complex(�շ�ϵͳ);
    fft(f.values);
    SAVE_SIGNAL_TF("../data/f.txt", f, std::abs(x));
    ifft(f.values);
    SAVE_SIGNAL_TF("../data/t.txt", f, x.real());
    
    return 0;
}
