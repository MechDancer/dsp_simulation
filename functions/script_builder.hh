//
// Created by ydrml on 2020/8/19.
//

#ifndef DSP_SIMULATION_SCRIPT_BUILDER_HH
#define DSP_SIMULATION_SCRIPT_BUILDER_HH

#include <string>
#include <vector>

namespace mechdancer {
    class script_builder_t {
        std::string data_file_path;
        std::vector<std::string> files;
    
    public:
        explicit script_builder_t(std::string const &);
        
        ~script_builder_t();
        
        std::string save(std::string const &);
    };
}

#endif //DSP_SIMULATION_SCRIPT_BUILDER_HH
