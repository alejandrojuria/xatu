#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>


class ConfigurationBase {

    public:
        std::string filename;
        std::string text;
        std::vector<std::string> foundArguments;
        std::vector<std::string> expectedArguments;
        struct configuration;
        std::map<std::string, std::vector<std::string>> contents;

    protected:
        std::ifstream m_file;

    public:
        ConfigurationBase();
        ConfigurationBase(std::string);

        std::string parseArgument(std::string);
        void extractArguments();
        void extractRawContent();
        virtual void parseContent() = 0;
        void checkArguments();

        template<typename T>
        std::vector<T> parseLine(const std::string& line){
            std::vector<T> values;
            std::istringstream iss(line);
            T value;
            while (iss >> value){
                values.push_back(value);
            }
            return values;
        };

        std::string standarizeLine(std::string&);

        template<typename T>
        T parseScalar(std::string& line){
            T value;
            std::istringstream iss(line);
            iss >> value;
            return value;
        }

        template<typename T>
        void printVector(std::vector<T>& v){
            for (auto i = v.begin(); i != v.end(); i++){
                std::cout << *i << std::endl;
            }
        };

        void printContent();
    
    protected:
        void restartFileStream();
                
};