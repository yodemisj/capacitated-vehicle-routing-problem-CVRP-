#pragma once
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <limits>
#include <deque>
#include <ctime>

#define POPULATION_SIZE 30
#define TRY_NUMBER 50

using namespace std;

// (σ = 30, Δ = 0,5, ρ = 8)

struct Chromosome {
    vector<int> nodes;
    vector<int> P;
    int fitness;

    Chromosome(const vector<int>& n, const vector<int>& p, int f) : nodes(n), P(p), fitness(f) {}

    bool operator<(const Chromosome& other) const {
        return fitness > other.fitness;
    }
};

const int INF = numeric_limits<int>::max();

class Graph {
private: 
    vector<vector<int>> adj;
    vector<vector<int>> weights;
    vector<int> loadCost;
    vector<Chromosome> population;

public:
    Graph(int n) {
        adj.resize(n + 1, vector<int>());
        loadCost.resize(n + 1, 0);
        weights.resize(n + 1, vector<int>(n + 1, -1));
        population = vector<Chromosome>();
        population.reserve(POPULATION_SIZE);
    }

    ~Graph() {
        adj.clear();
        weights.clear();
    }

    void addEdge(const int u,const int v, const int weight) {
            adj[u].push_back(v);
            setWeigth(u, v, weight);
    }

    void setWeigth(const int u, const int v, const int weight) {
        weights[u][v] = weight;
    }

    void setLoadCost(const int u, const int cost) {
        loadCost[u] = cost;
    }

    int splitProcedure(vector<int>& P, int W, int L) {
        int n = adj.size() - 1;
        vector<int> V(n + 1);
        V[0] = 0;
        for(int i = 1; i <= n; i++) {
            V[i] = INF;
        }

        for(int i = 1; i <= n; i++ ) {
            int load = 0;
            int cost = 0;
            int j = i;

            do {
                load += loadCost[j];
                if (i == j) {
                    cost += weights[0][j] + weights[j][0];
                } else {
                    cost = cost - weights[j-1][0] + weights[j-1][j] + weights[j][0];
                }

                if(load <= W && cost <= L) {
                    if (V[i-1] + cost < V[j]) { // Relaxa
                        V[j] = V[i-1] + cost;
                        P[j] = i - 1;
                    }
                    j++;
                } else {
                    break;
                }

            } while (j <= n);
        }

        return V[n];
    }

    void extractSolution(const vector<int>& P, vector<vector<int>>& trip) {
        const int n = adj.size() - 1;
        trip.clear(); 
        trip.resize(n + 1);

        int t = 0;
        int j = n; 
        int i;

        do {
            t += 1;
            i = P[j];
            for (int k = i + 1; k <= j; k++) {
                trip[t].push_back(k);
            }
            j = i;
        } while (i > 0);
    }

    bool isSpaced(double delta = 0.5) {
        int populationSize = population.size();
        
        for (int i = 0; i < populationSize; ++i) {
            for (int j = i + 1; j < populationSize; ++j) {
                if (population[i].nodes.empty() || population[j].nodes.empty()) continue;

                if (abs(population[i].fitness - population[j].fitness) < delta) {
                    return false;
                }
            }
        }
        
        return true;
    }

    int getRandomNumberInRange(int min, int max) {
        if (min > max) {
            swap(min, max);
        }

        return min + rand() % (max - min + 1);
    }

    double getRandomDouble(double min, double max) {
        if (min > max) {
            swap(min, max);
        }

        mt19937 rng(static_cast<unsigned int>(time(nullptr))); 
        uniform_real_distribution<double> dist(min, max); 
        return dist(rng); 
    }

    Chromosome generateRandomSolution() {
        vector<int> solution;
        vector<int> P(adj.size());

        for (int i = 1; i < adj.size(); ++i) {
            solution.push_back(i);
        }

        random_device rd;
        mt19937 gen(rd());

        shuffle(solution.begin(), solution.end(), gen);
        int cost = splitProcedure(P, 10, INF);
        return Chromosome(solution, P, cost);
    }

    void sortSolutions() {
        sort(population.begin(), population.end());
    }

    int lowerBound() {
        return 0;
    }

    Chromosome selectParent() {
        int populationSize = population.size();

        int idx1 = rand() % populationSize;
        int idx2 = rand() % populationSize;

        if (population[idx1].fitness < population[idx2].fitness) {
            return population[idx1];
        } else {
            return population[idx2];
        }
    }

    pair<vector<int>, vector<int>> crossoverOX(const vector<int> p1, const vector<int> p2) {
        int p1Size = p1.size();
        int p2Size = p2.size();

        if (p1Size != p2Size) {
            return {};
        }

        vector<int> c1(p1Size), c2(p1Size);

        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dis(1, p1Size - 2);

        int i = dis(gen);
        int j = dis(gen);

        if (j <= i) {
            swap(i, j);
        }

        for(int k = i; k <= j; k++) {
            c1[k] = p1[k];
            c2[k] = p2[k];
        }
        
        int indexC1 = (j + 1) % p1Size; 
        int indexC2 = (j + 1) % p1Size; 

        for (int k = j + 1; k < p1Size + j + 1; k++) {
            int index = k % p1Size;

            auto itP1 = find(c2.begin(), c2.end(), p1[index]);
            auto itP2 = find(c1.begin(), c1.end(), p2[index]);

            if(itP1 == c2.end()) {
                c2[indexC2] = p1[index];
                indexC2++;
                indexC2 %= p1Size;
            }

            if(itP2 == c1.end()) {
                c1[indexC1] = p2[index];
                indexC1++;
                indexC1 %= p1Size;
            }
        }
    }

    vector<int> moveM1(vector<int> route, int u, int v) {
        int clientU = route[u];
        route.erase(route.begin() + u);
        route.insert(route.begin() + v + 1, clientU);
        return route;
    }

    vector<int> moveM2(vector<int> route, int u, int x, int v) {
        int clientU = route[u];
        int clientX = route[x];
        route.erase(route.begin() + x);
        route.erase(route.begin() + u);
        route.insert(route.begin() + v + 1, clientX);
        route.insert(route.begin() + v + 1, clientU);
        return route;
    }

    vector<int> moveM3(vector<int> route, int u, int x, int v) {
        int clientU = route[u];
        int clientX = route[x];
        route.erase(route.begin() + x); 
        route.erase(route.begin() + u);
        route.insert(route.begin() + v + 1, clientU);
        route.insert(route.begin() + v + 1, clientX);
        return route;
    }

    vector<int> moveM4(vector<int> route, int u, int v) {
        swap(route[u], route[v]);
        return route;
    }

    vector<int> moveM5(vector<int> route, int u, int x, int v) {
        if (u < x) {
            int clientV = route[v];
            int clientU = route[u];
            int clientX = route[x];
            
            route.erase(route.begin() + v);
            route.erase(route.begin() + x);
            route.erase(route.begin() + u);

            route.insert(route.begin() + u, clientV);
            route.insert(route.begin() + u, clientX);
            route.insert(route.begin() + u, clientU);
        } else {
            int clientV = route[v];
            int clientU = route[u];
            int clientX = route[x];
            
            route.erase(route.begin() + u);
            route.erase(route.begin() + x);
            route.erase(route.begin() + v);

            route.insert(route.begin() + v, clientU);
            route.insert(route.begin() + v, clientX);
            route.insert(route.begin() + v, clientV);
        }

        return route;
    }

    vector<int> moveM6(vector<int> route, int u, int x, int v, int y) {
        int clientU = route[u];
        int clientX = route[x];
        int clientV = route[v];
        int clientY = route[y];

        if (u < x && v < y) {
            if (x < v) {
                route.erase(route.begin() + y);
                route.erase(route.begin() + v);
                route.erase(route.begin() + x);
                route.erase(route.begin() + u);

                route.insert(route.begin() + u, clientY);
                route.insert(route.begin() + u, clientV);
                route.insert(route.begin() + v, clientX);
                route.insert(route.begin() + v, clientU);
            } else {
                route.erase(route.begin() + x);
                route.erase(route.begin() + u);
                route.erase(route.begin() + y);
                route.erase(route.begin() + v);

                route.insert(route.begin() + v, clientX);
                route.insert(route.begin() + v, clientU);
                route.insert(route.begin() + u, clientY);
                route.insert(route.begin() + u, clientV);
            }
        } else {
            route.erase(route.begin() + max(x, y));
            route.erase(route.begin() + min(x, y));
            route.erase(route.begin() + max(u, v));
            route.erase(route.begin() + min(u, v));

            route.insert(route.begin() + min(u, v), clientY);
            route.insert(route.begin() + min(u, v), clientV);
            route.insert(route.begin() + min(x, y), clientX);
            route.insert(route.begin() + min(x, y), clientU);
        }

        return route;
    }

    vector<int> moveM7(vector<int> route, int u, int x, int v, int y) {
        int clientU = route[u];
        int clientX = route[x];
        int clientV = route[v];
        int clientY = route[y];

        if (u < v) {
            route.erase(route.begin() + y);
            route.erase(route.begin() + v);
            route.erase(route.begin() + x);
            route.erase(route.begin() + u);
        } else {
            route.erase(route.begin() + x); 
            route.erase(route.begin() + u); 
            route.erase(route.begin() + y); 
            route.erase(route.begin() + v); 
        }

        route.insert(route.begin() + min(u, v), clientV);
        route.insert(route.begin() + min(u, v), clientU);

        route.insert(route.begin() + max(x, y), clientY);
        route.insert(route.begin() + max(x, y), clientX);

        return route;
    }

    vector<int> moveM8(vector<int> route, int u, int x, int v, int y) {
        int clientU = route[u];
        int clientX = route[x];
        int clientV = route[v];
        int clientY = route[y];

        if (u < v) {
            route.erase(route.begin() + y); 
            route.erase(route.begin() + v); 
            route.erase(route.begin() + x); 
            route.erase(route.begin() + u); 
        } else {
            route.erase(route.begin() + x); 
            route.erase(route.begin() + u); 
            route.erase(route.begin() + y);
            route.erase(route.begin() + v); 
        }

        route.insert(route.begin() + min(u, v), clientV);
        route.insert(route.begin() + min(u, v), clientU);

        route.insert(route.begin() + max(x, y), clientY);
        route.insert(route.begin() + max(x, y), clientX);

        return route;
    }

    vector<int> moveM9(vector<int> route, int u, int x, int v, int y) {
        int clientU = route[u];
        int clientX = route[x];
        int clientV = route[v];
        int clientY = route[y];

        if (u < v) {
            route.erase(route.begin() + y); 
            route.erase(route.begin() + v); 
            route.erase(route.begin() + x); 
            route.erase(route.begin() + u);
        } else {
            route.erase(route.begin() + x); 
            route.erase(route.begin() + u); 
            route.erase(route.begin() + y); 
            route.erase(route.begin() + v); 
        }

        route.insert(route.begin() + min(u, v), clientY);
        route.insert(route.begin() + min(u, v), clientU);

        route.insert(route.begin() + max(x, y), clientV);
        route.insert(route.begin() + max(x, y), clientX);

        return route;
    }

    vector<int> extractCromossome(vector<vector<int>> solution) {
        vector<int> chromosome;

        for (const auto& route : solution) {
            chromosome.insert(chromosome.end(), route.begin(), route.end());
        }

        return chromosome;
    }

    vector<int> localSearch(const vector<int>& childSolution) {
        vector<int> p(adj.size());
        vector<vector<int>> solution;
        int fit = splitProcedure(p, 10, INF);
        Chromosome childChromosome = Chromosome(childSolution, p, fit);
        extractSolution(p, solution);

        bool improved = true;
        while (improved) {
            improved = false;

            for (size_t i = 0; i < solution.size(); ++i) {
                for (size_t u = 0; u < solution[i].size(); ++u) {
                    for (size_t j = 0; j < solution.size(); ++j) {
                        for (size_t v = 0; v < solution[j].size(); ++v) {
                            if (i == j && u == v) continue;

                            size_t x = (u + 1 < solution[i].size()) ? u + 1 : 0;  
                            size_t y = (v + 1 < solution[j].size()) ? v + 1 : 0; 

                            for (int move = 1; move <= 9; ++move) {
                                vector<vector<int>> newSolution = solution;

                                switch (move) {
                                    case 1:
                                        if (solution[i][u] != 0) {
                                            newSolution[i] = moveM1(newSolution[i], u, v);
                                        }
                                        break;
                                    case 2:
                                        if (solution[i][u] != 0 && x < solution[i].size() && solution[i][x] != 0) {
                                            newSolution[i] = moveM2(newSolution[i], u, x, v);
                                        }
                                        break;
                                    case 3:
                                        if (solution[i][u] != 0 && x < solution[i].size() && solution[i][x] != 0) {
                                            newSolution[i] = moveM3(newSolution[i], u, x, v);
                                        }
                                        break;
                                    case 4:
                                        if (solution[i][u] != 0 && solution[j][v] != 0) {
                                            newSolution[i] = moveM4(newSolution[i], u, v);
                                        }
                                        break;
                                    case 5:
                                        if (solution[i][u] != 0 && x < solution[i].size() && solution[i][x] != 0 && solution[j][v] != 0) {
                                            newSolution[i] = moveM5(newSolution[i], u, x, v);
                                        }
                                        break;
                                    case 6:
                                        if (solution[i][u] != 0 && solution[i][x] != 0 && solution[j][v] != 0 && solution[j][y] != 0) {
                                            // Imcompleto
                                        }
                                        break;
                                    case 7:
                                            // Imcompleto
                                        break;
                                    case 8:
                                            // Imcompleto
                                        break;
                                    case 9:
                                            // Imcompleto
                                        break;
                                    default:
                                        break;
                                }

                                vector<int> p(adj.size());
                                int newCost = splitProcedure(p, 10, INF);
                                p.clear();
                                int currentCost = splitProcedure(p, 10, INF);

                                if (newCost < currentCost) {
                                    solution = newSolution;
                                    improved = true;
                                    break;
                                }
                            }
                            if (improved) break;
                        }
                        if (improved) break;
                    }
                    if (improved) break;
                }
            }
        }

        return extractCromossome(solution);
    }

    void geneticAlgorithm() {
        int k, tryCount;

        const int maxIterations = 30000; 
        const int maxNoImprovementIterations = 10000; 
        int noImprovementCount = 0, improvementStreak = 0;
        int populationSize;

        Chromosome solution_1 = Chromosome({}, {}, INF);// Solução da heurística CW
        Chromosome solution_2 = Chromosome({}, {}, INF); // Solução da heurística MJ
        Chromosome solution_3 = Chromosome({}, {},INF); // Solução da heurística GM
        
        population[0] =  solution_1;
        population[1] =  solution_2;

        k = !isSpaced() ? 0 : 1;

        population[k + 1] =  solution_3;

        if(isSpaced()) k = k + 1;

        while (k < POPULATION_SIZE && tryCount <= TRY_NUMBER) {
            k = k + 1;
            tryCount = 0;

            while (isSpaced() && tryCount <= TRY_NUMBER) {
                tryCount += 1;
                Chromosome randomSolution = generateRandomSolution();
                population[k] = randomSolution;
            }
        }

        if(tryCount > TRY_NUMBER) populationSize = k - 1;

        sortSolutions();

        int iterations = 0;
        while (iterations < maxIterations && noImprovementCount < maxNoImprovementIterations && population[0].fitness != lowerBound()) {
            Chromosome P1 = selectParent();
            Chromosome P2 = selectParent();

            vector<int> childSolution;
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> dis(0, 1);

            pair<vector<int>, vector<int>> children = crossoverOX(P1.nodes, P2.nodes);
            childSolution = dis(gen) == 0 ? children.first : children.second;

            k = getRandomNumberInRange(populationSize/2, populationSize);
            double random = getRandomDouble(0.0, 1.0);
            double pm = 0.05;
            if(random < pm) { 
                vector<int> mutation = localSearch(childSolution);
                vector<int> p(adj.size());
                int fit = splitProcedure(p, 10, INF);
                Chromosome mutatedChromosome = Chromosome(mutation, p, fit); 

                Chromosome aux = population[k];
                population[k] = mutatedChromosome;

                if (isSpaced()) {
                    childSolution = mutatedChromosome.nodes;
                }

                p.clear();
                fit = splitProcedure(p, 10, INF);
                Chromosome childChromosome = Chromosome(childSolution, p, fit);

                population[k] = childChromosome;

                if(isSpaced()) {
                    iterations++;
                    if(childChromosome.fitness < population[0].fitness) {
                        noImprovementCount = 0;
                    } else {
                        noImprovementCount += 1;
                    }
                    sortSolutions(); // Era melhor usar um shift para mover apenas novo elemento para a posição correta
                } else {
                    population[k] = aux;
                }

            }
        }
    }


};
