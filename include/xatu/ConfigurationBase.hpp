#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <map>

namespace xatu {

/**
 * The ConfigurationBase abstract class provides functionality to parse general configuration files. 
 * Its purpose is not be instantiated by itself, but to serve as a base class for more specialized
 * configuration classes
 */
class ConfigurationBase {

    public:
        // List of arguments found in configuration file
        std::vector<std::string> foundArguments;
        // List of arguments expected to be present in configuration file
        std::vector<std::string> expectedArguments;
        // Semi-structured contents of the file. Each key of the dictionary corresponds
        // to one argument, while its value is the raw text corresponding to that argument,
        // which is to be processed further
        std::map<std::string, std::vector<std::string>> contents;

    protected:
        // Configuration file in ifstream format, ready to be parsed
        std::ifstream m_file;

    protected:
        ConfigurationBase();
        ConfigurationBase(std::string);
        virtual ~ConfigurationBase(){};

        // Methods used in the constructor
        void extractArguments();
        void extractRawContent();
        void checkArguments();

        // Parsing methods
        std::string  parseArgument(std::string);
        virtual void parseContent() = 0;
        double       parseFraction(std::string&);
        std::string  parseWord(std::string&);
        std::string  standarizeLine(std::string&);

        // Utilities
        void printContent();
        void restartFileStream();

        // Templated methods
        /**
         * Template method to parse one line containing only one value of a type T
         * @param line String to parse
         * @return Value T
         */
        template<typename T>
        T parseScalar(std::string& line){
            T value;
            std::istringstream iss(line);
            iss >> value;
            return value;
        }

        /**
         * Template method to parse a line made of some values of type T into a vector
         * with those values
         * @param line String to be parsed into a vector
         * @return Vector of values T
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

        /**
         * Template method to print a vector of type T to screen
         * @param v Vector to print
         * @return void
         */
        template<typename T>
        void printVector(std::vector<T>& v){
            for (auto i = v.begin(); i != v.end(); i++){
                std::cout << *i << "\t"; 
            }
            std::cout << std::endl;
        };
 
};

}
