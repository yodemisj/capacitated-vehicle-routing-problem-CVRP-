#pragma once
#include "CVRPInstance.h"
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

class FileUtils{
    public: 

    static CVRPInstance readInstanceFile(const string filePath) {
        ifstream inputFile(filePath);
        string aux, name;
        double comment, distance = -1;
        int dimension, vehicleCapacity, depotIndex;
        vector<int> demands;                 
        vector<NodeCoord> coordinates;

        if(!inputFile.is_open()) throw runtime_error("Erro ao abrir o arquivo!");

        inputFile >> aux;
        inputFile >> aux;
        inputFile >> name;
        // cout << "instance name: " << name << "\n";

        inputFile >> aux;
        inputFile >> aux;
        inputFile >> comment;
        // cout << "instance comment: " << comment << "\n";
        
        inputFile >> aux;
        inputFile >> aux;
        inputFile >> aux; // Ignore type

        inputFile >> aux;
        inputFile >> aux;
        inputFile >> dimension;
        // cout << "instance dimension: " << dimension << "\n";

        inputFile >> aux;
        inputFile >> aux;
        inputFile >> aux; // Ignore edge_weight_type

        inputFile >> aux;
        inputFile >> aux;
        inputFile >> vehicleCapacity;
        // cout << "instance capacity: " << vehicleCapacity << "\n";

        inputFile >> aux; // NODE_COORD_SECTION or DISTANCE

        if (aux == "DISTANCE") {
            inputFile >> aux;
            inputFile >> distance;
            inputFile >> aux; 
        }

        cout << distance << "\n";


        coordinates = vector<NodeCoord>(dimension + 1);

        for(int i = 0; i < dimension; i++) {
            int nodeIndex;
            inputFile >> nodeIndex;
            inputFile >> coordinates[nodeIndex].x;
            inputFile >> coordinates[nodeIndex].y;
        }

        inputFile >> aux; // DEMAND_SECTION

        demands = vector<int>(dimension + 1);

        for(int i = 0; i < dimension; i++) {
            int nodeIndex;
            inputFile >> nodeIndex;
            inputFile >> demands[nodeIndex];
        }

        inputFile >> aux; // DEPOT_SECTION

        inputFile >> depotIndex;
        // cout << "instance depot: " << depotIndex << "\n";

        CVRPInstance instance = CVRPInstance(name, comment, dimension, vehicleCapacity, demands, coordinates, distance, depotIndex );

        return instance;
        
    }
};