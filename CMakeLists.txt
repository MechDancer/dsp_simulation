cmake_minimum_required(VERSION 3.17)
project(dsp_simulation)

set(CMAKE_CXX_STANDARD 20)

add_executable(dsp_simulation main.cpp
        types/concepts.h
        types/db_t.hpp
        types/frequency_t.hpp
        types/signal_t.hpp

        functions/builders.h
        functions/functions.h
        functions/fft.h
        functions/processing.h)
