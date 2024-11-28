#include <iostream>
#include "FileUtils.h"
#include "CVRPSolver.h"
#include <string>
#include <chrono>

using namespace std;

void solveInstanceWithMetrics(const CVRPInstance &instance, double alpha, vector<CVRPResult> &results)
{
    CVRPSolver cvrpSolver = CVRPSolver();
    auto start = chrono::high_resolution_clock::now();
    double CK_WT = cvrpSolver.solve(instance);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> CK_WT_ExecTime = end - start;

    start = chrono::high_resolution_clock::now();
    double RCL = cvrpSolver.solveRCL(instance, alpha);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> RCL_ExecTime = end - start;

    // Salvar resultados no vetor
    results.push_back({instance.getName(),
                       instance.getDimension() - 1,
                       instance.getVehicleCapacity(),
                       instance.getDistance(),
                       instance.getServiceTime(),
                       CK_WT,
                       CK_WT_ExecTime.count(),
                       RCL,
                       RCL_ExecTime.count()});
}

int main()
{
    vector<CVRPResult> results;

    for (int i = 1; i <= 14; i++)
    {
        string instancePath = "./ins/Christofields/CMT";
        instancePath.append(to_string(i)).append(".in");
        CVRPInstance instance = FileUtils::readInstanceFile(instancePath);
        cout << "\nINSTANCE: " << instance.getName() << "\n";

        solveInstanceWithMetrics(instance, 0.02, results);
    }

    // for (int i = 1; i <= 20; i++) {
    //     string instancePath = "./ins/Golden/Golden";
    //     instancePath.append(to_string(i)).append(".in");
    //     CVRPInstance instance = FileUtils::readInstanceFile(instancePath);
    //     cout << "\nINSTANCE: " << instance.getName() << "\n";

    //     solveInstanceWithMetrics(instance, 0.02, results);
    // }

    // Salvar resultados em um arquivo CSV
    FileUtils::saveCVRPResultsToCSV("cvrp_results.csv", results);

    return 0;
}