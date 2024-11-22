#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

typedef struct NodeCoord {
    double x;
    double y;
} NodeCoord;

class CVRPInstance {
private: 
    string name;
    int depotIndex;
    double comment;
    double distance;
    int dimension;                           
    int vehicleCapacity;
    vector<int> demands;                 
    vector<NodeCoord> coordinates;
    vector<vector<double>> distanceMatrix; 
public: 
    CVRPInstance(string name, double comment, int dimension, int vehicleCapacity, const vector<int> demands, const vector<NodeCoord> coordinates, double distance = -1, int depotIndex = 1)
        : name(name), depotIndex(depotIndex), comment(comment), distance(distance), dimension(dimension), vehicleCapacity(vehicleCapacity), demands(demands), coordinates(coordinates)
    {
        calculateDistanceMatrix();
    }

    void calculateDistanceMatrix(){ 
        int totalNodes = dimension;
        distanceMatrix.resize(totalNodes + 1, vector<double>(totalNodes + 1 , 0.0));

        for (int i = 1; i <= totalNodes; ++i) {
            for (int j = 1; j <= totalNodes; ++j) {
                if (i != j) {
                    distanceMatrix[i][j] = euclideanDistance(coordinates[i], coordinates[j]);
                }
            }
        }
    }

    double euclideanDistance(const NodeCoord& p1, const NodeCoord& p2) const {
        return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
    }

    string getName() const { return name; }
    int getDimension() const { return dimension; }
    int getVehicleCapacity() const { return vehicleCapacity; }
    const vector<int>& getDemands() const { return demands; }
    const vector<NodeCoord>& getCoordinates() const { return coordinates; }
    const vector<vector<double>>& getDistanceMatrix() const { return distanceMatrix; }

    void printInstance() const {
        cout << "Número de clientes: " << dimension << "\n";
        cout << "Capacidade do veículo: " << vehicleCapacity << "\n";
        cout << "Demandas: ";
        for (const int& d : demands) {
            cout << d << " ";
        }
        cout << "\nCoordenadas:\n";
        for (size_t i = 0; i < coordinates.size(); ++i) {
            cout << "  Nó " << i << ": (" << coordinates[i].x << ", " << coordinates[i].y << ")\n";
        }
    }
    
    void printDistanceMatrix() const {
        cout << "\nDistances:\n";
        for (size_t i = 1; i < distanceMatrix.size(); ++i) {
            for(size_t j = i; j < distanceMatrix[i].size(); j++){
                cout << distanceMatrix[i][j] << " ";
            }   
            cout << endl;
        }
    }
    
};