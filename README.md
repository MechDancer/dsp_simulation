# dsp_simulation
下一代信号处理仿真实验库，基于 C++20

> [上一版](https://github.com/YdrMaster/signal_processing_simulation)

- 当前支持功能
  - 频率类型 `frequency_t`
  - 信号类型 `signal_t`
  - 基 2 快速傅里叶变换 `fft`/`ifft`
  - 基于 fft 的快速卷积
  - 基于 fft 的快速互相关，和两种白化滤波模式
  - 希尔伯特变换

- 这一版目标：

  1. 用 concept 强化类型系统，增强类型安全性（类似 Kotlin inline class 的功能）
  2. 实现完整、与单片机上逻辑完全一致的超声波测距仿真
