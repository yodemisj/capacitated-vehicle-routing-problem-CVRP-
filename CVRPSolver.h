#pragma once
#include "CVRPInstance.h"
#include <algorithm>
#include <tuple>
#include <random>
#include <thread>
#include <mutex>
using namespace std;

class CVRPSolver
{
private:
    vector<vector<int>> routes;
    vector<tuple<int, int, double>> savings;

    void initializeRoutes(int dimension, int depotIndex)
    {
        for (int i = 1; i <= dimension; i++)
        {
            if (i != depotIndex)
            {
                vector<int> route(1);
                route[0] = i;
                routes.push_back(route);
            }
        }
    }

    static void calculateSavingsForRange(int start, int end, int dimension, int depotIndex, const vector<vector<double>> &distanceMatrix, vector<tuple<int, int, double>> &localSavings, mutex &mtx)
    {
        for (int i = start; i <= end; i++)
        {
            for (int j = i + 1; j <= dimension; j++)
            {
                if (i != depotIndex && j != depotIndex)
                {
                    double saving = distanceMatrix[depotIndex][i] +
                                    distanceMatrix[depotIndex][j] -
                                    distanceMatrix[i][j];

                    lock_guard<mutex> lock(mtx);
                    localSavings.emplace_back(i, j, saving);
                }
            }
        }
    }

    void calculateSavings(const CVRPInstance &instance)
    {
        int dimension = instance.getDimension();
        int depotIndex = instance.getDepotIndex();
        vector<vector<double>> distanceMatrix = instance.getDistanceMatrix();

        vector<tuple<int, int, double>> localSavings;
        mutex mtx;

        int qtdThreads = thread::hardware_concurrency();
        int blockSize = dimension / qtdThreads;

        vector<thread> threads;

        for (int i = 0; i < qtdThreads; i++)
        {
            int start = i * blockSize + 1;
            int end = (i == qtdThreads - 1) ? dimension : (i + 1) * blockSize;

            // Lançar a thread para calcular savings no intervalo [start, end]
            threads.emplace_back(calculateSavingsForRange, start, end, dimension, depotIndex, ref(distanceMatrix), ref(savings), ref(mtx));
        }

        for (auto &t : threads)
        {
            t.join();
        }

        sort(savings.begin(), savings.end(), [](const auto &a, const auto &b)
             { return get<2>(a) > get<2>(b); });
    }

    int calculateRouteDemand(vector<int> route, vector<int> demands)
    {
        int totalDemand = 0;

        for (int node : route)
        {
            totalDemand += demands[node];
        }
        return totalDemand;
    }

    void mergeRoutes(int customerI, int customerJ, int vehicleCapacity, double distance, int depotIndex, double serviceTime, const vector<vector<double>> distanceMatrix, vector<int> demands)
    {
        int routeIndexI = -1, routeIndexJ = -1;

        for (size_t i = 0; i < routes.size(); i++)
        {
            if (find(routes[i].begin(), routes[i].end(), customerI) != routes[i].end())
            {
                routeIndexI = i;
            }
            if (find(routes[i].begin(), routes[i].end(), customerJ) != routes[i].end())
            {
                routeIndexJ = i;
            }

            if (routeIndexI != -1 && routeIndexJ != -1)
                break;
        }

        if ((routeIndexI == routeIndexJ) || (routeIndexI == -1 && routeIndexJ == -1))
            return;

        int totalDemand = 0;

        for (int customer : routes[routeIndexI])
        {
            totalDemand += demands[customer];
        }
        for (int customer : routes[routeIndexJ])
        {
            totalDemand += demands[customer];
        }

        if (distance > 0)
        {
            double totalDistance = 0;

            totalDistance += distanceMatrix[depotIndex][customerI];
            totalDistance += distanceMatrix[customerJ][depotIndex];
            for (size_t i = 0; i < routes[routeIndexI].size() - 1; i++)
            {
                totalDistance += distanceMatrix[routes[routeIndexI][i]][routes[routeIndexI][i + 1]];
            }
            for (size_t i = 0; i < routes[routeIndexJ].size() - 1; i++)
            {
                totalDistance += distanceMatrix[routes[routeIndexJ][i]][routes[routeIndexJ][i + 1]];
            }
            totalDistance += distanceMatrix[customerI][customerJ];

            if (serviceTime > 0)
            {
                int totalServiceTime = (routes[routeIndexI].size() * serviceTime) + (routes[routeIndexJ].size() * serviceTime);
                totalDistance += totalServiceTime;
            }
            // Se a distancia excede o limite
            if (totalDistance > distance)
                return;
        }

        if (totalDemand <= vehicleCapacity)
        {
            if (routes[routeIndexI].back() == customerI && routes[routeIndexJ].front() == customerJ)
            {
                // Concatenar rota J à rota I
                routes[routeIndexI].insert(routes[routeIndexI].end(), routes[routeIndexJ].begin(), routes[routeIndexJ].end());
                routes.erase(routes.begin() + routeIndexJ); // Remover rota J
            }
            else if (routes[routeIndexI].front() == customerI && routes[routeIndexJ].back() == customerJ)
            {
                // Concatenar rota I à rota J
                routes[routeIndexJ].insert(routes[routeIndexJ].end(), routes[routeIndexI].begin(), routes[routeIndexI].end());
                routes.erase(routes.begin() + routeIndexI); // Remover rota I
            }
        }
    }

    double calculateCost(const vector<vector<double>> distanceMatrix, int depotIndex, double serviceTime)
    {
        double cost = 0;
        for (vector<int> route : routes)
        {
            cost += distanceMatrix[depotIndex][route.front()];
            cost += distanceMatrix[route.back()][depotIndex];
            for (size_t i = 1; i < route.size(); i++)
            {
                cost += distanceMatrix[route[i - 1]][route[i]];
            }

            // if (serviceTime > 0)
            // {
            //     cost += route.size() * serviceTime;
            // }
        }

        return cost;
    }

    double calculateRouteCost(const vector<int> route, const vector<vector<double>> distanceMatrix, int depotIndex)
    {
        double cost = 0;
        cost += distanceMatrix[depotIndex][route.front()];
        cost += distanceMatrix[route.back()][depotIndex];

        for (size_t i = 1; i < route.size(); i++)
        {
            cost += distanceMatrix[route[i - 1]][route[i]];
        }

        return cost;
    }

    vector<tuple<int, int, double>> createRCL(double alpha)
    {
        vector<tuple<int, int, double>> rcl;
        if (savings.empty())
            return rcl;

        double maxSaving = get<2>(savings.front());
        double minSaving = get<2>(savings.back());
        double threshold = maxSaving - alpha * (maxSaving - minSaving);

        for (const auto &saving : savings)
        {
            if (get<2>(saving) >= threshold)
            {
                rcl.push_back(saving);
            }
            else
            {
                break;
            }
        }

        return rcl;
    }

    void clear()
    {
        routes.clear();
        savings.clear();
    }

    void printRoute(const vector<int> route)
    {
        for (int node : route)
        {
            cout << node << " ";
        }
        cout << endl;
    }

    vector<int> swap2Opt(vector<int> route, int i, int j)
    {
        vector<int> subRouteMid(route.begin() + i + 1, route.begin() + j + 1);
        reverse(subRouteMid.begin(), subRouteMid.end());

        vector<int> vectorResult;
        vectorResult.insert(vectorResult.end(), route.begin(), route.begin() + i + 1);
        vectorResult.insert(vectorResult.end(), subRouteMid.begin(), subRouteMid.end());
        vectorResult.insert(vectorResult.end(), route.begin() + j + 1, route.end());

        return vectorResult;
    }

    vector<int> generateNeighborhood_2Opt(vector<int> route, CVRPInstance instance)
    {
        vector<int> bestRoute = route;
        int totalDeamnd = 0;
        double newCost = 0;
        double bestCost = calculateRouteCost(route, instance.getDistanceMatrix(), instance.getDepotIndex());
        
        if(route.size() <= 3) {
            return route;
        }

        for (size_t i = 0; i < route.size() - 2; i++)
        {
            for (size_t j = i + 2; j < route.size() - 1; j++)
            {
                vector<int> newRoute = swap2Opt(route, i, j);
                newCost = calculateRouteCost(newRoute, instance.getDistanceMatrix(), instance.getDepotIndex());

                if (newCost < bestCost)
                {
                    bestRoute = newRoute;
                    bestCost = newCost;

                    cout << "Permutacao : ";
                    printRoute(newRoute);
                    cout << "COST : " << newCost << endl;
                }
            }
        }

        return bestRoute;
    }

public:
    CVRPSolver() : routes() {}

    double solve(CVRPInstance instance)
    {
        clear();
        initializeRoutes(instance.getDimension(), instance.getDepotIndex());
        calculateSavings(instance);

        for (const tuple<int, int, double> &saving : savings)
        {
            int customerI = get<0>(saving);
            int customerJ = get<1>(saving);

            mergeRoutes(customerI, customerJ, instance.getVehicleCapacity(), instance.getDistance(), instance.getDepotIndex(), instance.getServiceTime(), instance.getDistanceMatrix(), instance.getDemands());
        }

        // for (size_t i = 0; i < routes.size(); i++)
        // {
        //     cout << "Route #" << i + 1 << ": ";
        //     for (size_t j = 0; j < routes[i].size(); j++)
        //     {
        //         cout << routes[i][j] << " ";
        //     }
        //     cout << "\n";
        // }

        // cout << "Cost " << calculateCost(instance.getDistanceMatrix(), instance.getDepotIndex()) << "\n";

        // cout << "Savings:" << endl;
        // for (const auto& saving : savings) {
        //     cout << "Savin   g(" << get<0>(saving) << ", " << get<1>(saving)
        //         << ") = " << get<2>(saving) << endl;
        // }
        return calculateCost(instance.getDistanceMatrix(), instance.getDepotIndex(), instance.getServiceTime());
    }

    double solveRCL(CVRPInstance instance, double alpha)
    {
        clear();

        random_device rd;
        mt19937 rng(rd());

        initializeRoutes(instance.getDimension(), instance.getDepotIndex());
        calculateSavings(instance);

        vector<tuple<int, int, double>> rcl;

        while (!savings.empty())
        {
            rcl = createRCL(alpha);
            if (rcl.empty())
                break;

            uniform_int_distribution<int> dist(0, rcl.size() - 1);
            int selectedIndex = dist(rng);

            int customerI = get<0>(rcl[selectedIndex]);
            int customerJ = get<1>(rcl[selectedIndex]);

            mergeRoutes(customerI, customerJ, instance.getVehicleCapacity(), instance.getDistance(), instance.getDepotIndex(), instance.getServiceTime(), instance.getDistanceMatrix(), instance.getDemands());

            savings.erase(remove(savings.begin(), savings.end(), rcl[selectedIndex]), savings.end());
        }

        cout << "Solved by RCL \n";

        for (size_t i = 0; i < routes.size(); i++)
        {
            cout << "Route #" << i + 1 << ": ";
            for (size_t j = 0; j < routes[i].size(); j++)
            {
                cout << routes[i][j] << " ";
            }
            cout << "\n";
        }

        // cout << "Cost " << calculateCost(instance.getDistanceMatrix(), instance.getDepotIndex()) << "\n";

        // cout << "Savings:" << endl;
        // for (const auto& saving : savings) {
        //     cout << "Savin   g(" << get<0>(saving) << ", " << get<1>(saving)
        //         << ") = " << get<2>(saving) << endl;
        // }
        run2opt(instance);
        return calculateCost(instance.getDistanceMatrix(), instance.getDepotIndex(), instance.getServiceTime());
    }

    void run2opt(CVRPInstance instance)
    {
        int i = 1;
        double cost = 0;
        for (auto route : routes)
        {
            cout << "-------------------------------------------" << endl;
            cout << "Route #" << i++ << ": ";
            printRoute(route);
            cost = calculateRouteCost(route, instance.getDistanceMatrix(), instance.getDepotIndex());
            cout << "COST : " << cost << endl;

            route = generateNeighborhood_2Opt(route, instance);
        }
    }
};