#include <iostream>
#include "FileUtils.h"

using namespace std;


int main() {
    CVRPInstance instance1 = FileUtils::readInstanceFile("./ins/Christofields/CMT1.in");
    CVRPInstance instance2 = FileUtils::readInstanceFile("./ins/Golden/Golden1.in");
    CVRPInstance instance3 = FileUtils::readInstanceFile("./ins/testeinstance.in");

    cout << "Instance 1: " << endl;
    instance1.printInstance();
    cout << "Instance 2: " << endl;
    instance2.printInstance();

    instance3.printDistanceMatrix();

    return 0;
}