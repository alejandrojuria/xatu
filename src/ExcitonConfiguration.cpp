#include "ExcitonConfiguration.hpp"

ExcitonConfiguration::ExcitonConfiguration(){
    throw std::invalid_argument("Error: ExcitonConfiguration must be invoked with one argument (filename)");
};

ExcitonConfiguration::ExcitonConfiguration(std::string filename){
    parse
}