#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <typeinfo>

using namespace std;

int main(){
    string line = "1 2 3";
    istringstream iss(line);
    int x, y, z;
    string species;

    if (!(iss >> x >> y >> z >> species).fail()){
        cout << "Works" << endl;
        cout << x << endl;
        cout << species << endl;
    }    
    else{
        cout << "Fail" << endl;
    }

    return 0;
}