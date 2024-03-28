#include "RGFM_parser.hpp"
#include <iostream>
#include <string>
int main(int argc, char** argv){
    int grid_num = 0;
    if (argc - 1 > 0)
    {
        std::string s = argv[1];
        int grid_num = std::stoi(s);
    }
    RGFM_parser parser(grid_num);
    return 0;
};