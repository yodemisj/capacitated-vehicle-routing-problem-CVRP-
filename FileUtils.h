#pragma once
#include "CVRPInstance.h"
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;

struct CVRPResult {
    string instanceName;
    int numClients;
    int vehicleCapacity;
    double distance;
    double serviceTime;
    double CK_WT;
    double CK_WT_ExecTime;
    double RCL;
    double RCL_ExecTime;
};

class FileUtils{
    public: 

    static CVRPInstance readInstanceFile(const string filePath) {
        ifstream inputFile(filePath);
        string aux, name;
        double comment, distance = -1, serviceTime = -1;
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

        if (aux == "SERVICE_TIME") { // SERVICE_TIME or NODE_COORD_SECTION
            inputFile >> aux;
            inputFile >> serviceTime;
            inputFile >> aux; 
        } 

        // cout << distance << "\n";

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

        CVRPInstance instance = CVRPInstance(name, comment, dimension, vehicleCapacity, demands, coordinates, distance, serviceTime, depotIndex );

        inputFile.close();

        return instance;
        
    }

    static void saveCVRPResultsToCSV(const string& filePath, const vector<CVRPResult>& results) {
        ofstream file(filePath);

        file << "Nome da Instância;Número de Clientes;Capacidade;Distância;Tempo de Serviço;CK WT;Tempo de Execucao CK_WT (s);RCL;Tempo de Execucao RCL (s)\n";

        for (const auto& result : results) {
            stringstream execTimeStream;
            execTimeStream << fixed << setprecision(7) << result.CK_WT_ExecTime;  // Formata o tempo de execução com 7 casas decimais
            string CK_WT_ExecTimeStr = execTimeStream.str();
            replace(CK_WT_ExecTimeStr.begin(), CK_WT_ExecTimeStr.end(), '.', ',');  // Substitui o ponto decimal por vírgula (formato europeu)

            execTimeStream << fixed << setprecision(7) << result.RCL_ExecTime;  // Formata o tempo de execução com 7 casas decimais
            string RCL_ExecTimeStr = execTimeStream.str();
            replace(RCL_ExecTimeStr.begin(), RCL_ExecTimeStr.end(), '.', ',');

            file << result.instanceName << ";"
                << result.numClients << ";"
                << result.vehicleCapacity << ";"
                << fixed << setprecision(2) << result.distance << ";"
                << fixed << setprecision(2) << result.serviceTime << ";"
                << fixed << setprecision(2) << result.CK_WT << ";"
                << CK_WT_ExecTimeStr << ";"
                << fixed << setprecision(2) << result.RCL << ";"
                << RCL_ExecTimeStr << "\n";
        }

        file.close();
    }
};