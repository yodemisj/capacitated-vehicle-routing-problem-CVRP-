#include <iostream>
#include "FileUtils.h"
#include "CVRPSolver.h"
#include <string>
#include <chrono>

#define MAX_ITERATOR 1000
#define TABU_SIZE 100

using namespace std;

void solveInstanceWithMetrics(const CVRPInstance &instance, double alpha, vector<CVRPSHeuristicResult> &results)
{
    CVRPSolver cvrpSolver = CVRPSolver();

    auto start = chrono::high_resolution_clock::now();
    cvrpSolver.solveRCL(instance, alpha);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> RCL_ExecTime = end - start;

    vector<vector<int>> originalRoutes = cvrpSolver.getRoutes();

    // 2OPT FIRST IMPROVEMENT
    start = chrono::high_resolution_clock::now();
    double TS_2OPT_FI = cvrpSolver.runTabuSearch(instance, MAX_ITERATOR, TABU_SIZE);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> TS_2OPT_FI_ExecTime = end - start + RCL_ExecTime;
    cout << ">> TS_2OPT_FI: " << TS_2OPT_FI << endl;

    cvrpSolver.setRoutes(originalRoutes);

    // 2OPT BEST IMPROVEMENT
    start = chrono::high_resolution_clock::now();
    double TS_2OPT_BI = cvrpSolver.runTabuSearch(instance, MAX_ITERATOR, TABU_SIZE, true, false);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> TS_2OPT_BI_ExecTime = end - start + RCL_ExecTime;
    cout << ">> TS_2OPT_BI: " << TS_2OPT_BI << endl;
    
    cvrpSolver.setRoutes(originalRoutes);

    // 3OPT FIRST IMPROVEMENT
    start = chrono::high_resolution_clock::now();
    double TS_3OPT_FI = cvrpSolver.runTabuSearch(instance, MAX_ITERATOR, TABU_SIZE, false);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> TS_3OPT_FI_ExecTime = end - start + RCL_ExecTime;
    cout << ">> TS_3OPT_FI: " << TS_3OPT_FI << endl;

    cvrpSolver.setRoutes(originalRoutes);

    // 3OPT BEST IMPROVEMENT
    start = chrono::high_resolution_clock::now();
    double TS_3OPT_BI = cvrpSolver.runTabuSearch(instance, MAX_ITERATOR, TABU_SIZE, false, false);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> TS_3OPT_BI_ExecTime = end - start + RCL_ExecTime;
    cout << ">> TS_3OPT_BI: " << TS_3OPT_BI << endl;

    // Salvar resultados no vetor
    results.push_back({instance.getName(),
                       instance.getDimension() - 1,
                       instance.getVehicleCapacity(),
                       instance.getDistance(),
                       instance.getServiceTime(),
                       TS_2OPT_FI_ExecTime.count(),
                       TS_2OPT_BI_ExecTime.count(),
                       TS_2OPT_FI,
                       TS_2OPT_BI,
                       TS_3OPT_FI_ExecTime.count(),
                       TS_3OPT_BI_ExecTime.count(),
                       TS_3OPT_FI,
                       TS_3OPT_BI});
}

int main()
{
    // vector<CVRPSHeuristicResult> results;

    CVRPSolver cvrpSolver = CVRPSolver();
    // string instancePath = "./ins/Christofields/CMT";
    string instancePath = "./ins/Golden/Golden";
    instancePath.append(to_string(1)).append(".in");
    CVRPInstance instance = FileUtils::readInstanceFile(instancePath);
    cout << "\nINSTANCE: " << instance.getName() << "\n";

    cvrpSolver.runGeneticAlgorithm(instance, 0.02);

    // for (int i = 1; i <= 14; i++)
    // {
    //     string instancePath = "./ins/Christofields/CMT";
    //     instancePath.append(to_string(i)).append(".in");
    //     CVRPInstance instance = FileUtils::readInstanceFile(instancePath);
    //     cout << "\nINSTANCE: " << instance.getName() << "\n";

    //     solveInstanceWithMetrics(instance, 0.02, results);
    // }

    // for (int i = 1; i <= 20; i++)
    // {
    //     string instancePath = "./ins/Golden/Golden";
    //     instancePath.append(to_string(i)).append(".in");
    //     CVRPInstance instance = FileUtils::readInstanceFile(instancePath);
    //     cout << "\nINSTANCE: " << instance.getName() << "\n";

    //     solveInstanceWithMetrics(instance, 0.02, results);
    // }

    // Salvar resultados em um arquivo CSV
    // FileUtils::saveCVRPSHeuristicResultsToCSV("cvrp_s_heuristic_results.csv", results);

    return 0;
}