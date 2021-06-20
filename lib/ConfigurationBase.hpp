#include <string>
#include <vector>
#include <fstream>
#include <map>


class ConfigurationBase {
    public:
        ConfigurationBase();
        ConfigurationBase(std::string);

        std::string parseArgument(std::string);
        void extractArguments();
        void extractRawContent();
        virtual void parseContent() = 0;
        void checkArguments();

        template<typename T>
        std::vector<T> parseLine(const std::string&);
        std::string standarizeLine(std::string&);

        template<typename T>
        T parseScalar(std::string&);

        template<typename T>
        void printVector(std::vector<T>&);

        void printContent();
    
    protected:
        void restartFileStream();


    public:
        std::string filename;
        std::string text;
        std::vector<std::string> foundArguments;
        std::vector<std::string> expectedArguments;
        struct configuration;
        std::map<std::string, std::vector<std::string>> contents;

    protected:
        std::ifstream m_file;
        
        
};