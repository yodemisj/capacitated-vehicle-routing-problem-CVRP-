#include <iostream>
#include "FileUtils.h"
#include "CVRPSolver.h"
#include <string>

using namespace std;

int main() {
    // CVRPInstance instance1 = FileUtils::readInstanceFile("./ins/Christofields/CMT14.in");
    CVRPSolver crvpSolver = CVRPSolver();
    // CVRPInstance instance2 = FileUtils::readInstanceFile("./ins/Golden/Golden1.in");
    // CVRPInstance instance3 = FileUtils::readInstanceFile("./ins/testeinstance.in");

    // cout << "Instance 1: " << endl;
    // instance1.printInstance();
    // cout << "Instance 2: " << endl;
    // instance2.printInstance();

    // instance3.printDistanceMatrix();

    // for(int i = 1; i <= 14; i++) {
    //     string instancePath = "./ins/Christofields/CMT";
    //     instancePath.append(to_string(i)).append(".in");
    //     CVRPInstance instance1 = FileUtils::readInstanceFile(instancePath);
    //     cout << "\nINSTANCE : " << instance1.getName() << "\n"; 
    //     crvpSolver.solve(instance1);
    //     cout << "\n";
    //     crvpSolver.solveRCL(instance1, 0.02);
    // }

    for(int i = 1; i <= 20; i++) {
        string instancePath = "./ins/Golden/Golden";
        instancePath.append(to_string(i)).append(".in");
        CVRPInstance instance1 = FileUtils::readInstanceFile(instancePath);
        cout << "\nINSTANCE : " << instance1.getName() << "\n"; 
        crvpSolver.solve(instance1);
        cout << "\n";
        crvpSolver.solveRCL(instance1, 0.02);
    }

    return 0;
}