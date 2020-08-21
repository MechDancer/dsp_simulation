cmake_minimum_required(VERSION 3.17)
project(dsp_simulation)

set(CMAKE_CXX_STANDARD 20)

add_executable(dsp_simulation test/simulation.cpp
        types/concepts.h
        types/complex_t.hpp
        types/frequency_t.hpp
        types/noise.h
        types/signal_t.hpp

        functions/builders.h
        functions/functions.h
        functions/fft.h
        functions/process_real.h
        functions/process_complex.h

        functions/script_builder.cc
        functions/script_builder.hh)
