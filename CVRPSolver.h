#pragma once
#include "CVRPInstance.h"
#include <algorithm>
#include <tuple>
#include <random>
#include <deque>
#include <chrono>
using namespace std;

#define POPULATION_SIZE 30
#define TRY_NUMBER 50

static mt19937 gen(chrono::steady_clock::now().time_since_epoch().count() + random_device{}());

typedef struct Chromosome
{
    vector<int> nodes;
    vector<int> P;
    double fitness;

    Chromosome(vector<int> nodes, vector<int> p, double fitness)
    {
        this->nodes = nodes;
        this->fitness = fitness;
        this->P = p;
    };

    bool operator<(const Chromosome &other) const
    {
        return fitness < other.fitness;
    }
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

    void calculateSavings(const CVRPInstance &instance)
    {
        int dimension = instance.getDimension();
        int depotIndex = instance.getDepotIndex();
        vector<vector<double>> distanceMatrix = instance.getDistanceMatrix();

        for (int i = 1; i <= dimension; i++)
        {
            for (int j = i + 1; j <= dimension; j++)
            {
                if (i != depotIndex && j != depotIndex)
                {
                    double saving = distanceMatrix[depotIndex][i] +
                                    distanceMatrix[depotIndex][j] -
                                    distanceMatrix[i][j];
                    savings.emplace_back(i, j, saving);
                }
            }
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

    bool isSpaced(vector<Chromosome> population, size_t idx, double delta = 0.5)
    {
        for (size_t i = 0; i < population.size(); i++)
        {
            if (i == idx)
                continue;
            else if (abs(population[i].fitness - population[idx].fitness) < delta)
            {
                return false;
            }
        }
        return true;
    }

    bool isSpaced(vector<Chromosome> population, Chromosome newChromosome, double delta = 0.5)
    {
        for (size_t i = 0; i < population.size(); i++)
        {
            if (abs(population[i].fitness - newChromosome.fitness) < delta)
            {
                return false;
            }
        }
        return true;
    }

    Chromosome generateRandomSolution(CVRPInstance instance, mt19937 &gen)
    {
        vector<int> solution;
        vector<int> P(instance.getDimension());

        for (int i = 1; i <= instance.getDimension(); ++i)
        {
            if (i == instance.getDepotIndex())
                continue;
            solution.push_back(i);
        }

        shuffle(solution.begin(), solution.end(), gen);
        solution.emplace(solution.begin(), instance.getDepotIndex());
        // printRoute(solution);
        int fitness = splitProcedure(instance, solution, P);
        // cout << "Fitness: " << fitness << endl;
        return Chromosome(solution, P, fitness);
    }

    void sortSolutions(vector<Chromosome> population)
    {
        sort(population.begin(), population.end());
    }

    Chromosome selectParent(const vector<Chromosome> &population)
    {
        int populationSize = population.size();

        int idx1 = getRandomInt(0, populationSize - 1, gen);
        int idx2;
        do
        {
            idx2 = getRandomInt(0, populationSize - 1, gen);
        } while (idx1 == idx2);

        if (population[idx1].fitness < population[idx2].fitness)
        {
            return population[idx1];
        }
        else
        {
            return population[idx2];
        }
    }

    double getRandomDouble(double min, double max, mt19937 &gen)
    {
        if (min > max)
        {
            swap(min, max);
        }

        // mt19937 rng(static_cast<unsigned int>(time(nullptr)));
        uniform_real_distribution<double> dist(min, max);
        return dist(gen);
    }

    int getRandomInt(int min, int max, mt19937 &gen)
    {
        // static thread_local mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> dist(min, max);
        return dist(gen);
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

    vector<int> swap2Opt(const vector<int> route, int i, int j)
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

    pair<vector<int>, vector<int>> crossoverOX(const vector<int> p1, const vector<int> p2, mt19937 &gen)
    {
        int p1Size = p1.size();
        int p2Size = p2.size();

        if (p1Size != p2Size)
        {
            return {};
        }

        vector<int> c1(p1Size), c2(p1Size);

        // random_device rd;
        // mt19937 gen(rd());
        uniform_int_distribution<> dis(1, p1Size - 2);

        int i = dis(gen);
        int j = dis(gen);

        while (i == j)
        {
            j = dis(gen);
        }

        if (j < i)
        {
            swap(i, j);
        }

        for (int k = i; k <= j; k++)
        {
            c1[k] = p1[k];
            c2[k] = p2[k];
        }

        int indexC1 = (j + 1) % p1Size;
        int indexC2 = (j + 1) % p1Size;

        for (int k = j + 1; k < p1Size + j + 1; k++)
        {
            int index = k % p1Size;

            auto itP1 = find(c2.begin(), c2.end(), p1[index]);
            auto itP2 = find(c1.begin(), c1.end(), p2[index]);

            if (itP1 == c2.end())
            {
                c2[indexC2] = p1[index];
                indexC2++;
                indexC2 %= p1Size;
            }

            if (itP2 == c1.end())
            {
                c1[indexC1] = p2[index];
                indexC1++;
                indexC1 %= p1Size;
            }
        }
        return make_pair(c1, c2);
    }

    void shiftSolution(vector<Chromosome> &population, size_t k)
    {
        size_t n = population.size();

        while (k > 0 && population[k].fitness < population[k - 1].fitness)
        {
            swap(population[k], population[k - 1]);
            k--;
        }

        while (k < n - 1 && population[k].fitness > population[k + 1].fitness)
        {
            swap(population[k], population[k + 1]);
            k++;
        }
    }

    Chromosome generateBestNeighborhood_2Opt_GA(Chromosome route, CVRPInstance instance)
    {
        double newNeighborCost = 0;
        vector<int> P(route.nodes.size());

        if (route.nodes.size() <= 3)
        {
            return route;
        }

        Chromosome bestNeighbor = route;

        for (size_t i = 0; i < route.nodes.size() - 2; i++)
        {
            for (size_t j = i + 2; j < route.nodes.size() - 1; j++)
            {
                vector<int> newNeighborRoute = swap2Opt(route.nodes, i, j);
                P.clear();
                newNeighborCost = splitProcedure(instance, route.nodes, P);
                Chromosome newNeighbor = Chromosome(newNeighborRoute, P, newNeighborCost);

                if (newNeighborCost < bestNeighbor.fitness)
                {
                    bestNeighbor = newNeighbor;
                    // if(bestNeighborCost < routeCost) return bestNeighbor; // Primeira melhora
                }
            }
        }

        return bestNeighbor;
    }

    Chromosome localSearch(CVRPInstance instance, Chromosome route)
    {
        Chromosome newRoute = generateBestNeighborhood_2Opt_GA(route, instance);
        return newRoute;
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

        // random_device rd;
        // mt19937 rng(rd());

        initializeRoutes(instance.getDimension(), instance.getDepotIndex());
        calculateSavings(instance);

        vector<tuple<int, int, double>> rcl;

        while (!savings.empty())
        {
            rcl = createRCL(alpha);
            if (rcl.empty())
                break;

            uniform_int_distribution<int> dist(0, rcl.size() - 1);
            // int selectedIndex = dist(rng);
            int selectedIndex = dist(gen);

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

    double splitProcedure(const CVRPInstance instance, const vector<int> route, vector<int> &P)
    {
        int capacity = instance.getVehicleCapacity();
        double distance = (instance.getDistance() == -1) ? numeric_limits<double>::infinity() : instance.getDistance();
        int n = instance.getDimension();

        vector<double> V(n + 1, numeric_limits<double>::infinity());
        vector<int> demands = instance.getDemands();
        vector<vector<double>> distanceMatrix = instance.getDistanceMatrix();
        V[route[0]] = 0;

        for (int i = 1; i < n; i++)
        {
            int load = 0;
            double cost = 0;
            int j = i;

            do
            {
                load += demands[route[j]];
                if (i == j)
                {
                    cost += distanceMatrix[route[0]][route[j]] + distanceMatrix[route[j]][route[0]];
                }
                else
                {
                    cost = cost - distanceMatrix[route[j - 1]][route[0]] + distanceMatrix[route[j - 1]][route[j]] + distanceMatrix[route[j]][route[0]];
                }
                if (load <= capacity && cost <= distance)
                {
                    if (V[route[i - 1]] + cost < V[route[j]])
                    { // Relaxa
                        V[route[j]] = V[route[i - 1]] + cost;
                        P[route[j]] = route[i - 1]; // to do: verificar isso aqui
                    }
                    j++;
                }
                else
                {
                    break;
                }

            } while (j < route.size());
        }
        return V[route[n - 1]];
    }

    Chromosome parseToChromosome(CVRPInstance instance, vector<vector<int>> routes)
    {
        vector<int> result;
        vector<int> P(instance.getDimension() + 1);
        result.push_back(instance.getDepotIndex());
        for (vector<int> subroute : routes)
        {
            result.insert(result.end(), subroute.begin(), subroute.end());
        }
        double fitness = splitProcedure(instance, result, P);
        return Chromosome(result, P, fitness);
    }

    void partialReplacement(vector<Chromosome> &population, CVRPInstance instance, int rho = 8, int maxRestartAttempts = 5)
    {
        int replacements = 0;
        int attempts = 0;
        int tryCount = 0;
        vector<Chromosome> newPopulation;
        newPopulation.reserve(8);

        while (replacements < rho && attempts < maxRestartAttempts)
        {
            attempts++;
            newPopulation.clear();
            int tryCount = 0;
            while (newPopulation.size() < 8 && tryCount < TRY_NUMBER)
            {
                Chromosome newChrom = generateRandomSolution(instance, gen);

                if (isSpaced(population, newChrom) && isSpaced(newPopulation, newChrom))
                {
                    newPopulation.push_back(newChrom);
                    shiftSolution(newPopulation, newPopulation.size() - 1);
                }
            }
            for (int k = 0; k < newPopulation.size(); k++)
            {
                Chromosome bestChild = newPopulation[k];

                for (const Chromosome &parent : population)
                {
                    if (bestChild.fitness < population.back().fitness)
                        break;
                    vector<int> P(instance.getDimension());
                    pair<vector<int>, vector<int>> children = crossoverOX(newPopulation[k].nodes, parent.nodes, gen);
                    int fit1 = splitProcedure(instance, children.first, P);
                    int fit2 = splitProcedure(instance, children.second, P);
                    vector<int> child = (fit1 <= fit2) ? children.first : children.second;
                    int fit = (fit1 < fit2) ? fit1 : fit2;
                    Chromosome childChromosome = Chromosome(child, P, fit);
                    if (isSpaced(population, childChromosome) && childChromosome.fitness < bestChild.fitness)
                    {
                        bestChild = childChromosome;
                    }
                }

                if (bestChild.fitness < population.back().fitness)
                {
                    population.back() = bestChild;
                    replacements++;
                    shiftSolution(population, population.size() - 1);
                }
            }
        }
    }

    void genCWSolutions(CVRPInstance instance, vector<Chromosome> &population, int n, double alpha)
    {
        int k = 0;
        int tryCount = 0;
        while (k < n && tryCount <= TRY_NUMBER)
        {
            tryCount = 0;
            while (tryCount <= TRY_NUMBER)
            {
                tryCount += 1;
                solveRCL(instance, alpha);
                Chromosome solution_1 = parseToChromosome(instance, this->routes);
                cout << solution_1.fitness << endl;
                if (tryCount == 1)
                    population.push_back(solution_1);
                else
                    population.back() = solution_1;
                if (isSpaced(population, k))
                    break;
            }
            shiftSolution(population, k);
            k++;
        }
        if (tryCount > TRY_NUMBER)
            population.pop_back();
    }

    void runPhase(CVRPInstance instance, vector<Chromosome> &population, int maxIterations, int maxNoImprovementIterations, int k, double pm)
    {
        int iterations = 0;
        int noImprovementCount = 0;
        int populationSize = population.size();
        while (iterations < maxIterations && noImprovementCount < maxNoImprovementIterations)
        {
            Chromosome P1 = selectParent(population);
            Chromosome P2 = selectParent(population);

            pair<vector<int>, vector<int>> children = crossoverOX(P1.nodes, P2.nodes, gen);
            vector<int> childSolution = (getRandomInt(0, 1, gen) == 0) ? children.first : children.second;

            vector<int> P(instance.getDimension());
            int fitness = splitProcedure(instance, childSolution, P);
            Chromosome childChromosome = Chromosome(childSolution, P, fitness);

            k = getRandomInt((populationSize - 1) / 2, populationSize - 1, gen);
            double random = getRandomDouble(0.0, 1.0, gen);
            Chromosome aux = population[k];
            if (random < pm)
            {
                Chromosome mutatedChromosome = localSearch(instance, childChromosome);

                population[k] = mutatedChromosome;

                if (isSpaced(population, k))
                {
                    childChromosome = mutatedChromosome;
                }
            }

            population[k] = childChromosome;

            if (isSpaced(population, k))
            {
                iterations++;
                if (childChromosome.fitness < population[0].fitness)
                {
                    noImprovementCount = 0;
                    cout << "Improvement: " << childChromosome.fitness << endl;
                }
                else
                {
                    noImprovementCount += 1;
                }
                shiftSolution(population, k);
            }
            else
            {
                population[k] = aux;
            }
        }
    }

    double runGeneticAlgorithm(CVRPInstance instance, double alpha)
    {
        int k, tryCount = 0;
        vector<Chromosome> population;
        int maxIterations = 30000;
        int maxNoImprovementIterations = 10000;
        int noImprovementCount = 0, improvementStreak = 0;
        int populationSize;
        population.reserve(POPULATION_SIZE);

        genCWSolutions(instance, population, 3, alpha);
        cout << "P[0]: " << population[0].fitness << endl;
        cout << "P[1]: " << population[1].fitness << endl;
        cout << "P[2]: " << population[2].fitness << endl;

        k = population.size() - 1;

        while (k < POPULATION_SIZE && tryCount <= TRY_NUMBER)
        {
            k++;
            tryCount = 0;

            while (tryCount <= TRY_NUMBER)
            {
                tryCount += 1;
                Chromosome randomSolution = generateRandomSolution(instance, gen);
                if (tryCount == 1)
                    population.push_back(randomSolution);
                else
                    population.back() = randomSolution;
                if (isSpaced(population, k))
                    break;
            }
            shiftSolution(population, k);
        }
        if (tryCount > TRY_NUMBER)
            population.pop_back();

        int maxPhases;
        int pm = 0.05;
        bool isCMT = instance.getName().find("CMT") != string::npos;
        if (isCMT)
        {
            runPhase(instance, population, maxIterations, maxNoImprovementIterations, k, pm);
            maxPhases = 10;
            pm = 0.1;
            partialReplacement(population, instance);
        }
        else
        {
            maxPhases = 8;
            pm = 0.2;
        }

        maxIterations = 2000;
        maxNoImprovementIterations = 2000;

        for (int phase = 0; phase < maxPhases; phase++)
        {
            runPhase(instance, population, maxIterations, maxNoImprovementIterations, k, pm);
            partialReplacement(population, instance);
        }

        return population[0].fitness;
    }
};