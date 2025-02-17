#pragma once
#include "CVRPInstance.h"
#include <algorithm>
#include <tuple>
#include <random>
#include <thread>
#include <mutex>
#include <deque>
using namespace std;

typedef struct Chromosome
{
    vector<int> nodes;
    double fitness;

    Chromosome(vector<int> nodes, double fitness)
    {
        this->nodes = nodes;
        this->fitness = fitness;
    };
} Chromosome;

typedef struct Route
{
    vector<int> path;
    double cost;

    Route(vector<int> path, double cost)
    {
        this->path = path;
        this->cost = cost;
    };
} Route;

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

    double calculateCost(const vector<vector<double>> distanceMatrix, int depotIndex)
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

    double calculateRouteCostWithDepot(const vector<int> route, const vector<vector<double>> distanceMatrix)
    {
        double cost = 0;

        for (size_t i = 1; i < route.size(); i++)
        {
            cost += distanceMatrix[route[i - 1]][route[i]];
        }

        cost += distanceMatrix[route[0]][route.back()];

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

public:
    void printRoute(const vector<int> route)
    {
        for (int node : route)
        {
            cout << node << " ";
        }
        cout << endl;
    }

    int bestMove(int a, int b, int c, int d, int e, int f, CVRPInstance instance)
    {
        vector<vector<double>> distances = instance.getDistanceMatrix();
        double gains[8];
        double maxGain;
        int bestCase = 0;
        int i;

        gains[0] = 0;
        gains[1] = distances[a][e] + distances[b][f] - distances[a][b] - distances[e][f];
        gains[2] = distances[c][e] + distances[d][f] - distances[c][d] - distances[e][f];
        gains[3] = distances[a][c] + distances[b][d] - distances[a][b] - distances[c][d];

        int deletedEdges = distances[a][b] + distances[c][d] + distances[e][f];
        gains[4] = distances[a][c] + distances[b][e] + distances[d][f] - deletedEdges;
        gains[5] = distances[a][e] + distances[d][b] + distances[c][f] - deletedEdges;
        gains[6] = distances[a][d] + distances[e][c] + distances[b][f] - deletedEdges;
        gains[7] = distances[a][d] + distances[e][b] + distances[c][f] - deletedEdges;

        for (i = 1; i < 8; i++)
        {
            if (gains[i] < 0 && gains[i] < maxGain)
            {
                bestCase = i;
                maxGain = gains[i];
            }
        }

        return bestCase;
    }

    void reverseSubRoute(vector<int> &route, int start, int end)
    {
        size_t routeSize = route.size();
        size_t left = start;
        size_t right = end;
        int temp;

        size_t half = ((routeSize + end - start + 1) % routeSize) / 2;

        for (size_t i = 1; i < half + 1; i++)
        {
            temp = route[right];
            route[right] = route[left];
            route[left] = temp;

            left = (left + 1) % routeSize;
            right = (routeSize + right - 1) % routeSize;
        }
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

    vector<int> swap3Opt(vector<int> route, int bestCase, int i, int j, int k)
    {

        int routeSize = route.size();
        auto routeBegin = route.begin();
        auto routeEnd = route.end();

        switch (bestCase)
        {
        case 1:
            reverseSubRoute(route, ((k + 1) % routeSize), i);
            break;

        case 2:
            reverseSubRoute(route, ((j + 1) % routeSize), k);
            break;

        case 3:
            reverseSubRoute(route, ((i + 1) % routeSize), j);
            break;

        case 4:
            reverseSubRoute(route, ((j + 1) % routeSize), k);
            reverseSubRoute(route, ((i + 1) % routeSize), j);
            break;

        case 5:
            reverseSubRoute(route, ((k + 1) % routeSize), i);
            reverseSubRoute(route, ((i + 1) % routeSize), j);
            break;

        case 6:
            reverseSubRoute(route, ((k + 1) % routeSize), i);
            reverseSubRoute(route, ((j + 1) % routeSize), k);
            break;

        case 7:
            reverseSubRoute(route, ((k + 1) % routeSize), i);
            reverseSubRoute(route, ((i + 1) % routeSize), j);
            reverseSubRoute(route, ((j + 1) % routeSize), k);
            break;
        }
        return route;
    }

    vector<int> removeDepotFromRoute(vector<int> route, int depotIndex)
    {
        vector<int> routeResult;
        int index;
        for (size_t i = 0; i < route.size(); i++)
        {
            if (route[i] == depotIndex)
            {
                index = i;
                break;
            }
        }

        for (size_t i = index + 1; i < route.size() + index; i++)
        {
            routeResult.push_back(route[i % route.size()]);
        }

        return routeResult;
    }

    Route generateBestNeighborhood_2Opt(vector<int> route, CVRPInstance instance, deque<vector<int>> &tabuList)
    {
        double newNeighborCost = 0;
        double routeCost = calculateRouteCost(route, instance.getDistanceMatrix(), instance.getDepotIndex());

        if (route.size() <= 3)
        {
            return Route(route, routeCost);
        }

        double bestNeighborCost = numeric_limits<double>::infinity();
        vector<int> bestNeighbor;

        for (size_t i = 0; i < route.size() - 2; i++)
        {
            for (size_t j = i + 2; j < route.size() - 1; j++)
            {
                vector<int> newNeighbor = swap2Opt(route, i, j);

                if (find(tabuList.begin(), tabuList.end(), newNeighbor) == tabuList.end())
                {
                    newNeighborCost = calculateRouteCost(newNeighbor, instance.getDistanceMatrix(), instance.getDepotIndex());

                    if (newNeighborCost < bestNeighborCost)
                    {
                        bestNeighbor = newNeighbor;
                        bestNeighborCost = newNeighborCost;
                    }
                }
            }
        }

        return Route(bestNeighbor, bestNeighborCost);
    }

    Route generateNeighborhood_3Opt(vector<int> route, CVRPInstance instance, deque<vector<int>> &tabuList)
    {
        route.insert(route.begin(), instance.getDepotIndex());
        int i, j, k;
        int a, b, c, d, e, f;
        int bestCase;
        int routeSize = route.size();
        double newNeighborCost = 0;

        if (routeSize < 4)
        {
            return Route(route, calculateRouteCostWithDepot(route, instance.getDistanceMatrix()));
        }

        double bestNeighborCost = numeric_limits<double>::infinity();
        vector<int> bestNeighbor;

        size_t c1, c2, c3;

        for (c1 = 0; c1 < routeSize; c1++)
        {
            i = c1;
            a = route[i];
            b = route[(i + 1) % routeSize];

            for (c2 = 1; c2 < routeSize - 2; c2++)
            {
                j = (i + c2) % routeSize;
                c = route[j];
                d = route[(j + 1) % routeSize];

                for (c3 = c2 + 1; c3 < routeSize; c3++)
                {
                    k = (i + c3) % routeSize;
                    e = route[k];
                    f = route[(k + 1) % routeSize];

                    bestCase = bestMove(a, b, c, d, e, f, instance);
                    if (bestCase != 0)
                    {
                        vector<int> newNeighbor = swap3Opt(route, bestCase, i, j, k);

                        if (find(tabuList.begin(), tabuList.end(), newNeighbor) == tabuList.end())
                        {
                            newNeighborCost = calculateRouteCost(newNeighbor, instance.getDistanceMatrix(), instance.getDepotIndex());

                            if (newNeighborCost < bestNeighborCost)
                            {
                                bestNeighbor = newNeighbor;
                                bestNeighborCost = newNeighborCost;

                                cout << "BEST CASE : " << bestCase << endl;
                                cout << "Permutacao : ";
                                printRoute(bestNeighbor);
                                cout << "COST : " << newNeighborCost << endl;
                            }
                        }
                    }
                }
            }
        }

        return Route(removeDepotFromRoute(bestNeighbor, instance.getDepotIndex()), bestNeighborCost);
    }

public:
    CVRPSolver() : routes() {}

    vector<int> tabuSearch(vector<int> initialSolution, CVRPInstance instance, int maxIterations, int tabuListSize, bool is2Opt = true, bool fisrtImprovement = true)
    {
        vector<int> currentSolution = initialSolution;
        vector<int> bestSolution = initialSolution;
        const vector<vector<double>> distances = instance.getDistanceMatrix();
        int depotIndex = instance.getDepotIndex();

        double initialSolutionCost = calculateRouteCost(initialSolution, instance.getDistanceMatrix(), instance.getDepotIndex());

        double bestCost = calculateRouteCost(currentSolution, distances, depotIndex);

        deque<vector<int>> tabuList;
        int iterations = 0;

        while (iterations < maxIterations)
        {
            ++iterations;

            Route bestNeighbor = is2Opt ? generateBestNeighborhood_2Opt(currentSolution, instance, tabuList) : generateNeighborhood_3Opt(currentSolution, instance, tabuList);

            if (bestNeighbor.cost == numeric_limits<double>::infinity())
            {
                break;
            }

            currentSolution = bestNeighbor.path;
            if (bestNeighbor.cost < bestCost)
            {
                bestSolution = bestNeighbor.path;
                bestCost = bestNeighbor.cost;

                if (fisrtImprovement)
                {
                    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Primeira melhora!" << endl;
                    cout << "Iteração " << iterations << ": Custo Atual = " << bestCost << endl;
                    cout << "Initial cost: " << initialSolutionCost << endl;
                    break;
                }
            }
            // cout << "bestNeighbor: ";
            // printRoute(bestNeighbor.path);
            // cout << "COST: " << bestNeighbor.cost << endl;

            tabuList.push_back(currentSolution);
            if (tabuList.size() > tabuListSize)
            {
                tabuList.pop_front();
            }

            // cout << "Iteração " << iterations << ": Custo Atual = " << bestCost << endl;
        }

        return bestSolution;
    }

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
        return calculateCost(instance.getDistanceMatrix(), instance.getDepotIndex());
    }

    vector<vector<int>> getRoutes()
    {
        return routes;
    }

    void setRoutes(vector<vector<int>> newRoutes)
    {
        routes = newRoutes;
    }

    void solveRCL(CVRPInstance instance, double alpha)
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
    }

    double runTabuSearch(CVRPInstance instance, int maxIterator, int tabuListSize, bool is2Opt = true, bool isFirstImprovement = true)
    {
        for (size_t i = 0; i < routes.size(); i++)
        {

            routes[i] = tabuSearch(routes[i], instance, maxIterator, tabuListSize, is2Opt, isFirstImprovement);
        }

        return calculateCost(instance.getDistanceMatrix(), instance.getDepotIndex());
    }

    int splitProcedure(CVRPInstance instance, vector<int> &P, int capacity, int distance)
    {
        int n = instance.getDimension();
        vector<double> V(n + 1, numeric_limits<double>::infinity());
        V[0] = 0;

        for (int i = 1; i <= n; i++)
        {
            int load = 0;
            int cost = 0;
            int j = i;

            do
            {
                load += loadCost[j];
                if (i == j)
                {
                    cost += weights[0][j] + weights[j][0];
                }
                else
                {
                    cost = cost - weights[j - 1][0] + weights[j - 1][j] + weights[j][0];
                }

                if (load <= capacity && cost <= distance)
                {
                    if (V[i - 1] + cost < V[j])
                    { // Relaxa
                        V[j] = V[i - 1] + cost;
                        P[j] = i - 1;
                    }
                    j++;
                }
                else
                {
                    break;
                }

            } while (j <= n);
        }

        return V[n];
    }

    Chromosome parseToChromosome(vector<vector<int>> routes)
    {
        // PS: não consideramos o depósito
        vector<int> result;
        for (vector<int> subroute : routes)
        {
            result.insert(result.end(), subroute.begin(), subroute.end());
        }

        splitProcedure(result);
    }

    void runGeneticAlgorithm(CVRPInstance instance, double alpha)
    {
        int k, tryCount;

        const int maxIterations = 30000;
        const int maxNoImprovementIterations = 10000;
        int noImprovementCount = 0, improvementStreak = 0;
        int populationSize;

        solveRCL(instance, alpha);

        Chromosome solution_1 = Chromosome(); // Solução da heurística CW

        population[0] = solution_1;
        population[1] = solution_2;

        k = !isSpaced() ? 0 : 1;

        population[k + 1] = solution_3;

        if (isSpaced())
            k = k + 1;

        while (k < POPULATION_SIZE && tryCount <= TRY_NUMBER)
        {
            k = k + 1;
            tryCount = 0;

            while (isSpaced() && tryCount <= TRY_NUMBER)
            {
                tryCount += 1;
                Chromosome randomSolution = generateRandomSolution();
                population[k] = randomSolution;
            }
        }

        if (tryCount > TRY_NUMBER)
            populationSize = k - 1;

        sortSolutions();

        int iterations = 0;
        while (iterations < maxIterations && noImprovementCount < maxNoImprovementIterations && population[0].fitness != lowerBound())
        {
            Chromosome P1 = selectParent();
            Chromosome P2 = selectParent();

            vector<int> childSolution;
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> dis(0, 1);

            pair<vector<int>, vector<int>> children = crossoverOX(P1.nodes, P2.nodes);
            childSolution = dis(gen) == 0 ? children.first : children.second;

            k = getRandomNumberInRange(populationSize / 2, populationSize);
            double random = getRandomDouble(0.0, 1.0);
            double pm = 0.05;
            if (random < pm)
            {
                vector<int> mutation = localSearch(childSolution);
                vector<int> p(adj.size());
                int fit = splitProcedure(p, 10, INF);
                Chromosome mutatedChromosome = Chromosome(mutation, p, fit);

                Chromosome aux = population[k];
                population[k] = mutatedChromosome;

                if (isSpaced())
                {
                    childSolution = mutatedChromosome.nodes;
                }

                p.clear();
                fit = splitProcedure(p, 10, INF);
                Chromosome childChromosome = Chromosome(childSolution, p, fit);

                population[k] = childChromosome;

                if (isSpaced())
                {
                    iterations++;
                    if (childChromosome.fitness < population[0].fitness)
                    {
                        noImprovementCount = 0;
                    }
                    else
                    {
                        noImprovementCount += 1;
                    }
                    sortSolutions(); // Era melhor usar um shift para mover apenas novo elemento para a posição correta
                }
                else
                {
                    population[k] = aux;
                }
            }
        }
    }
};