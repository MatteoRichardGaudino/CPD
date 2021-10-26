#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char const *argv[]){

    if(argc == 1) return 1;

    ifstream fs(argv[1], ifstream::in);

    int nil;
    int count = 0;
    string s;

    double time = 0;
    double sum = 0;
    while (fs.good()) {
     fs >> nil;
     fs >> s;

     fs >> time;
     cout << time << endl;

     sum += time;
     count++;
    }

    fs.close();

    time /= count;
    cout << time << endl;


    return 0;
}