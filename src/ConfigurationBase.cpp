#include <stdexcept>
#include <sstream>
#include <algorithm>
#include "xatu/ConfigurationBase.hpp"

namespace xatu {

/**
 * Default constructor.
 * @details The default constructor throws an error as the class must always be
 * initialized providing a filename with the configuration to be parsed.
 */
ConfigurationBase::ConfigurationBase(){
    throw std::invalid_argument("ConfigurationBase must be called with one argument (filename)");
};

/**
 * File constructor.
 * @details This constructor takes a filename as input, which is passed to the class with
 * member initialization to avoid a call to the default constructor. Then, it checks if the string
 * is actually non-empty and if it specifies an existing file. Otherwise, it throws an error.
 * @param file Name of configuration file to be parsed.
 */
ConfigurationBase::ConfigurationBase(std::string file) : filename(file){
    if(file.empty()){
        throw std::invalid_argument("ConfigurationBase: file must not be empty");
    }
    
    m_file.open(file.c_str());
    if(!m_file.is_open()){
        throw std::invalid_argument("ConfigurationBase: file '" + file + "' does not exist");
    }
};

/**
 * Method to parse arguments.
 * @details Given some string, it looks for the # delimiter which signals that there is 
 * an argument. Then, it strips everything except for the argument, which is systematically
 * lowercased.
 * @param line Line containing the argument.
 * @returns Argument.
 */
std::string ConfigurationBase::parseArgument(std::string line){
    std::size_t pos = line.find("#");
    std::string str = line.substr(pos + 1);
    str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
};

/**
 * Method to extract all arguments from the raw text.
 * @details This method finds all lines containing an argument delimiter,
 * which are then parsed to extract the argument itself, and then are stored in the
 * attribute foundArguments.
 */
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

/**
 * Method to check if the arguments found contain at least the expected arguments.
 * @details Checks if all expected argument are present within the found ones. If not,
 * it throws an error.
 */
void ConfigurationBase::checkArguments(){
    if(expectedArguments.empty()){
        throw std::logic_error("Expected arguments must be defined first");
    };
    for (auto arg = expectedArguments.begin(); arg!= expectedArguments.end(); arg++){
            if(!(std::find(foundArguments.begin(), foundArguments.end(), *arg) != foundArguments.end())){
                throw std::logic_error("Missing arguments in config. file: missing " + *arg);
            }
        }
};

/**
 * Method to extract the information from the configuration file into
 * semi-structured form.
 * @details The file contents are parsed into a dictionary, whose
 * keys are the arguments and the values are the raw text between delimiters.
 */
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
            /* Skip empty or lines with comments */;
        }
        // Store argument contents
        else{
            content.push_back(standarizeLine(line));
        }
        // If EOF append lastly parsed content        
        if (m_file.eof()) {
            contents[arg] = content;
        }
    }
    
    restartFileStream();
}

/**
 * Method to go back to the beginning of the configuration file.
 * @details Useful to perform several reads of the file, to extract
 * different bits of information (e.g. to extract all arguments independently from
 * parsing the contents).
 */
void ConfigurationBase::restartFileStream(){
    m_file.clear();
    m_file.seekg(0);
}

/**
 * Method to standarize the value delimiters present in each row.
 * @details The value delimiters can be ",", ";" or simply " ".
 * All these possible delimiters are systematically converted to " " (blank) for simplicity.
 * @param line String to be standarized.
 * @returns Standarized line.
 */
std::string ConfigurationBase::standarizeLine(std::string& line) {
    if (line.find(',') != std::string::npos){
        std::replace(line.begin(), line.end(), ',', ' ');
    }
    if (line.find(';') != std::string::npos) {
        std::replace(line.begin(), line.end(), ';', ' ');
    }
    return line;
}

/**
 * Auxiliary method to print the contents of the dictionary with the semi-structured information.
 * @details Useful for debugging.
 */
void ConfigurationBase::printContent() {
    for (int i = 0; i != foundArguments.size(); i++) {
        std::string arg = foundArguments[i];
        std::cout << arg << std::endl;
        std::vector<std::string> section = contents[arg];
        printVector(section);
    }
}

/**
 * Method to parse a line containing a fraction written symbolically.
 * @details This method is intented to be used with a line where a fraction
 * is written symbolically, i.e. "a/b". Currently unused.
 * @param content String containing the symbolic fraction.
 * @returns Numeric value of fraction.
 */
double ConfigurationBase::parseFraction(std::string& content){
    std::istringstream iss(content);
    double numerator, denominator, fraction;
    std::string numeratorStr, denominatorStr;
    std::getline(iss, numeratorStr, '/');
    std::getline(iss, denominatorStr);
    numerator = std::stod(numeratorStr);
    try{
        denominator = std::stod(denominatorStr);
        fraction = numerator / denominator;
    }
    catch(std::exception& e){
        fraction = numerator;
    }
    return fraction;
}

/**
 * Method to parse a string from the configuration files.
 * @details This method removes all spaces and change the string to lowercase.
 * @param content String to be parsed.
 * @returns Parsed string.
 */
std::string ConfigurationBase::parseWord(std::string& content){
    std::string str = content;
    str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;

};

}