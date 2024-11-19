#include <iostream>
#include "FileUtils.h"

using namespace std;

int main() {
    FileUtils fileUtils;
    InstanceObj instance;

    instance = fileUtils.readInstanceFile("./ins/CMT1.in");

    cout << "Node coord: " << instance.nodeCoords[50].node << " " << instance.nodeCoords[50].x << " " << instance.nodeCoords[50].y << endl;
    cout << "Node demand: " << instance.demands[50].node << " " << instance.demands[50].value << endl;

    return 0;
}