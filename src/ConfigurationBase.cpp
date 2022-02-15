#include <iostream>
#include <stdexcept>
#include <sstream>
#include <algorithm>

#include "ConfigurationBase.hpp"

ConfigurationBase::ConfigurationBase(){
    throw std::invalid_argument("ConfigurationBase must be called with one argument (filename)");
};

ConfigurationBase::ConfigurationBase(std::string file) : filename(file){
    if(file.empty()){
        throw std::invalid_argument("ConfigurationBase: Filename must not be empty");
    }
    
    m_file.open(file.c_str());
    std::cout << file.c_str() << std::endl;
    if(!m_file.is_open()){
        throw std::invalid_argument("ConfigurationBase: File does not exist");
    }
};

std::string ConfigurationBase::parseArgument(std::string line){
    std::size_t pos = line.find("#");
    std::string str = line.substr(pos + 1);
    str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
    return str;
};

void ConfigurationBase::extractArguments(){
    std::vector<std::string> arguments;
    std::string line;
    while (std::getline(m_file, line)){
        if (line.find("#") != std::string::npos && m_file.peek() != EOF){
            std::string arg = parseArgument(line);
            arguments.push_back(arg);
        }
    }
    restartFileStream();
    this->foundArguments = arguments;
};

void ConfigurationBase::checkArguments(){
    if(expectedArguments.empty()){
        throw std::logic_error("Expected arguments must be defined first");
    };
    for (auto arg = expectedArguments.begin(); arg!= expectedArguments.end(); arg++){
            if(!(std::find(foundArguments.begin(), foundArguments.end(), *arg) != foundArguments.end())){
                throw std::logic_error("Missing arguments in config. file");
            }
        }
};

void ConfigurationBase::extractRawContent(){
    std::string line, arg;
    std::vector<std::string> content;
    while (std::getline(m_file, line)){
        // Argument detection and parsing
        if (line.find("#") != std::string::npos){
            if (!content.empty()) {
                contents[arg] = content;
                content.clear();
            }
            arg = parseArgument(line);
        }
        // Empty line or commentary detection
        else if (!line.size() || (line.find("!") != std::string::npos)){
            continue;
        }
        // Store argument contents
        else{
            content.push_back(standarizeLine(line));
        }
        if (m_file.eof()) {
            contents[arg] = content;
        }
    }
    
    restartFileStream();
}

void ConfigurationBase::restartFileStream(){
    m_file.clear();
    m_file.seekg(0);
}

std::string ConfigurationBase::standarizeLine(std::string& line) {
    if (line.find(',') != std::string::npos){
        std::replace(line.begin(), line.end(), ',', ' ');
    }
    if (line.find(';') != std::string::npos) {
        std::replace(line.begin(), line.end(), ';', ' ');
    }
    return line;
}

void ConfigurationBase::printContent() {
    for (int i = 0; i != foundArguments.size(); i++) {
        std::string arg = foundArguments[i];
        std::cout << arg << std::endl;
        std::vector<std::string> section = contents[arg];
        printVector(section);
    }
}

double ConfigurationBase::parseFraction(std::string& content){
    std::istringstream iss(content);
    double numerator, denominator, fraction;
    std::string numeratorStr, denominatorStr;
    std::getline(iss, numeratorStr, '/');
    std::getline(iss, denominatorStr);
    numerator = std::stod(numeratorStr);
    denominator = std::stod(denominatorStr);
    try{
        fraction = numerator / denominator;
    }
    catch(std::exception& e){
        fraction = numerator;
    }
    return fraction;
}