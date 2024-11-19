#include <iostream>
#include "FileUtils.h"

using namespace std;

int main() {
    FileUtils fileUtils;
    InstanceObj instance1;
    InstanceObj instance2;

    instance1 = fileUtils.readInstanceFile("./ins/CMT1.in");

    cout << "Cris: " << endl; 
    cout << "Node coord: " << instance1.nodeCoords[50].node << " " << instance1.nodeCoords[50].x << " " << instance1.nodeCoords[50].y << endl;
    cout << "Node demand: " << instance1.demands[50].node << " " << instance1.demands[50].value << endl;

    instance2 = fileUtils.readInstanceFile("./ins/Golden1.in");
    cout << "Golden: " << endl; 

    cout << "Node coord: " << instance2.nodeCoords[50].node << " " << instance2.nodeCoords[50].x << " " << instance2.nodeCoords[50].y << endl;
    cout << "Node demand: " << instance2.demands[50].node << " " << instance2.demands[50].value << endl;

    return 0;
}