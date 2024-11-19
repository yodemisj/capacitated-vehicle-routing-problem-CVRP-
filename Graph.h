#pragma once
#include <vector>

using namespace std;

class Graph {
private: 
    vector<vector<int>> adj;
    vector<vector<int>> weights;
    vector<int> loadCost;
public: 
    Graph(int n) {
        adj.resize(n + 1, vector<int>());
        loadCost.resize(n + 1, 0);
        weights.resize(n + 1, vector<int>(n + 1, -1));
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
};