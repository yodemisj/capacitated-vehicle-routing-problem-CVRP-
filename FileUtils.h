#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

typedef struct InstanceInfo{
    string name;
    double comment;
    string type;
    int dimension;
    string edgeWeightType;
    double distance;
    int capacity;
    int depot;
} InstanceInfo;

typedef struct NodeCoord {
    int node;
    double x;
    double y;
} NodeCoord;

typedef struct Demand {
    int node;
    int value;
} Demand;

typedef struct InstanceObj {
    InstanceInfo info;
    vector<NodeCoord> nodeCoords;
    vector<Demand> demands; 
} InstanceObj;

class FileUtils{
    public: 

    static InstanceObj readInstanceFile(const string filePath) {
        ifstream inputFile(filePath);
        string aux;
        InstanceObj instance;

        if(!inputFile.is_open()) throw runtime_error("Erro ao abrir o arquivo!");

        instance.info.distance = 0;

        inputFile >> aux;
        inputFile >> aux;
        inputFile >> instance.info.name;
        cout << "instance name: " << instance.info.name << "\n";

        inputFile >> aux;
        inputFile >> aux;
        inputFile >> instance.info.comment;
        cout << "instance comment: " << instance.info.comment << "\n";
        
        inputFile >> aux;
        inputFile >> aux;
        inputFile >> instance.info.type;
        cout << "instance type: " << instance.info.type << "\n";

        inputFile >> aux;
        inputFile >> aux;
        inputFile >> instance.info.dimension;
        cout << "instance dimension: " << instance.info.dimension << "\n";

        inputFile >> aux;
        inputFile >> aux;
        inputFile >> instance.info.edgeWeightType;
        cout << "instance edgeWeightType: " << instance.info.edgeWeightType << "\n";

        inputFile >> aux;
        inputFile >> aux;
        inputFile >> instance.info.capacity;
        cout << "instance capacity: " << instance.info.capacity << "\n";




        inputFile >> aux; // NODE_COORD_SECTION or DISTANCE

        if (aux == "DISTANCE") {
            inputFile >> aux;
            inputFile >> instance.info.distance;
            inputFile >> aux; 
        }

        cout << instance.info.distance << "\n";


        instance.nodeCoords = vector<NodeCoord>(instance.info.dimension);

        for(int i = 0; i < instance.info.dimension; i++) {
            inputFile >> instance.nodeCoords[i].node;
            inputFile >> instance.nodeCoords[i].x;
            inputFile >> instance.nodeCoords[i].y;
        }

        inputFile >> aux; // DEMAND_SECTION

        instance.demands = vector<Demand>(instance.info.dimension);

        for(int i = 0; i < instance.info.dimension; i++) {
            inputFile >> instance.demands[i].node;
            inputFile >> instance.demands[i].value;
        }

        inputFile >> aux; // DEPOT_SECTION

        inputFile >> instance.info.depot;
        cout << "instance depot: " << instance.info.depot << "\n";

        return instance;
        
    }
};