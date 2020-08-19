//
// Created by ydrml on 2020/8/19.
//

#include "script_builder.hh"
#include <sstream>
#include <fstream>
#include <iostream>
#include <filesystem>

mechdancer::script_builder_t::script_builder_t(std::string const &path) {
    data_file_path = __FILE__;
    data_file_path.erase(data_file_path.find_last_of('\\'));
    data_file_path.erase(data_file_path.find_last_of('\\') + 1);
    std::stringstream builder;
    builder << data_file_path << path << '\\';
    data_file_path = builder.str();
    std::filesystem::remove_all(data_file_path);
    std::filesystem::create_directory(data_file_path);
}

mechdancer::script_builder_t::~script_builder_t() {
    std::stringstream builder;
    builder << data_file_path << "script.txt";
    
    std::ofstream script(builder.str());
    builder.str("");
    
    for (auto const &file : files)
        builder << "load('" << data_file_path << file << ".txt');" << std::endl;
    
    auto text = builder.str();
    std::cout << std::endl << text;
    script << text;
}

std::string mechdancer::script_builder_t::save(std::string const &file) {
    files.push_back(file);
    std::stringstream builder;
    builder << data_file_path << file << ".txt";
    return builder.str();
}