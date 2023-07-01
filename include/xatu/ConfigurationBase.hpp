#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>

namespace xatu {

/**
 * The ConfigurationBase class provides functionality to parse general configuration files. 
 * Its purpose is not be instantiated by itself, but to serve as a base class for more specialized
 * configuration classes.
 */
class ConfigurationBase {

    public:
        // Name of configuration file to be parsed.
        std::string filename;
        // Raw text read from the file without any parsing.
        std::string text;
        // List of arguments found in configuration file.
        std::vector<std::string> foundArguments;
        // List of arguments expected to be present in configuration file.
        std::vector<std::string> expectedArguments;
        // Parsed information of file already fully structured.
        struct configuration;
        // Semi-structured contents of the file. Each key of the dictionary corresponds
        // to one argument, while its value is the raw text corresponding to that argument,
        // which is to be processed further.
        std::map<std::string, std::vector<std::string>> contents;

    protected:
        std::ifstream m_file;

    protected:
        ConfigurationBase();
        ConfigurationBase(std::string);

        std::string parseArgument(std::string);
        void extractArguments();
        void extractRawContent();
        virtual void parseContent() = 0;
        void checkArguments();

        /**
         * Template method to parse a line made of some values of type T into a vector
         * with those values.
         * @param line String to be parsed into a vector.
         * @returns Vector of values T.
         */
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

        /**
         * Template method to parse one line containing only one value of a type T.
         * @param line String to parse.
         * @returns Value T.
         */
        template<typename T>
        T parseScalar(std::string& line){
            T value;
            std::istringstream iss(line);
            iss >> value;
            return value;
        }

        /**
         * Template method to print a vector of type T to screen.
         * @param v Vector to print.
         */
        template<typename T>
        void printVector(std::vector<T>& v){
            for (auto i = v.begin(); i != v.end(); i++){
                std::cout << *i << "\t"; 
            }
            std::cout << std::endl;
        };

        void printContent();
        double parseFraction(std::string&);
    
        void restartFileStream(); 
};

}