#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

/**
 * @brief Function to create initial 2d vector for adjacency matrix
 * @param A adjacency matrix
 */
void initial_agreement(vector<vector<int>>& A) {
    for (size_t i = 0; i < A.size(); i++) {
        A[i][i] = 0;
    }
}

/**
 * @brief Function to flip the node or link at a given position if the energy change is favorable
 * 
 * @param A Adjacency matrix
 * @param nodes Vector of nodes
 * @param T Temprature
 * @param gen Random number generator
 * @param ran_u Uniform real distribution
 * @param ran_pos Uniform int distribution
 */
void check_flip(vector<vector<int>>& A, vector<int>& nodes, double T, mt19937& gen, uniform_real_distribution<double>& ran_u, uniform_int_distribution<int>& ran_pos) {
    int N = nodes.size();
    double P = static_cast<double>(N) / (N + N * (N - 1) / 2);

    // Link choose
    if ( ran_u(gen) > P) {

        int row = ran_pos(gen);
        int col = ran_pos(gen);

        int dE = 2 * nodes[row] * nodes[col] * A[row][col];
        if (dE <= 0 || ran_u(gen) < exp(-dE / T)) {
            A[row][col] *= -1;
            A[col][row] *= -1;
        }
    }
    // Node choose
    else {
        int index = ran_pos(gen);
        int sum = 0;
        for (int i = 0; i < N; i++) {
            sum += A[index][i] * nodes[i];
        }
        int dE = 2 *nodes[index] * sum;
        
        if (dE <= 0 || (ran_u(gen)) < exp(-dE / T)) {
            nodes[index] *= -1;
        }
    }
}

/**
 * @brief Function to calculate node link corrolation
 * 
 * @param A Adjacency matrix
 * @param nodes Vector of nodes
 * @return double 
 */
double node_link_corr(vector<vector<int>> A, vector<int> nodes) {
    int N = nodes.size();
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            sum += nodes[i] * A[i][j];
        }
    }
    return (2 * sum) / (N * N);
}

/**
 * @brief Function to calculate the energy of network per number of triplets
 * 
 * @param A Adjacency matrix
 * @param nodes Vector of nodes
 * @return double 
 */
double node_link_node_corr(vector<vector<int>> A, vector<int> nodes) {
    int N = nodes.size();
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            sum += nodes[i] * A[i][j] * nodes[j];
        }
    }
    return (-2 * sum) / (N * N);
}


/**
 * @brief Function to write vectors to a CSV file
 * 
 * @param temp temprature vector
 * @param Q node_link_corr vector
 * @param Energy energy vector
 * @param filename file name
 */
void writeVectorsToCSV(const vector<double>& temp, const vector<double>& nlc, const vector<double>& energy, const string& filename) {
    ofstream file(filename);

    if (file.is_open()) {
        // Write column headers
        file << "Temperature,node link corr,Energy\n";

        // Write data to the file
        size_t size = temp.size();
        for (size_t i = 0; i < size; ++i) {
            file << temp[i] << "," << nlc[i] << "," << energy[i] << "\n";
        }

        // Close the file
        file.close();
        cout << "File saved successfully!" << endl;
    }
    else {
        cout << "Error opening the file!" << endl;
    }
}


int main() {
    const int N = 16;
    const int agreement = -1; 
    const int nens = 50; 
    const double minTemp = 0.001;
    const double maxTemp = 0.2 * sqrt(N-1) + sqrt(N-1);
    const int numTemps = 40;
    const int steps = pow(N, 3);
    double T_increment = (maxTemp - minTemp) / (numTemps - 1);


    mt19937 gen(random_device{}());
	uniform_int_distribution<int> ran_pos(0, N-1); //Get any random integer
	uniform_real_distribution<double> ran_u(0.0, 1.0); //Our uniform variable generator


    vector<double> temp;
    vector<double> Energy;
    vector<double> node_link;

    auto start = high_resolution_clock::now();

    for (double T = minTemp; T <= maxTemp; T+= T_increment){
        temp.push_back(T);

        double avgE = 0.0;
        double avgQ = 0.0;

        for (int ens = 0; ens < nens; ens++) {
            vector<int> nodes(N, agreement);
            vector<vector<int>> A(N, vector<int>(N, 1));
            initial_agreement(A);

            for (int i = 0; i < steps; i++) {
                check_flip(A, nodes, T, gen, ran_u, ran_pos);
            }

            double E = node_link_node_corr(A, nodes);
            double Q = node_link_corr(A, nodes);

            avgE += E;
            avgQ += Q;
        }
        avgE /= nens;
        avgQ /= nens;

        Energy.push_back(avgE);
        node_link.push_back(avgQ);
    }

    string filePrefix = "C:\\Users\\Asus\\Documents\\Programming\\Coev-balance\\Coevolutionary\\Coev_";
    string fileExtension = ".csv";
    string fileName = filePrefix+ to_string(N) + "_agreement_" + to_string(agreement)+ fileExtension;

    writeVectorsToCSV(temp, node_link, Energy, fileName);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);

    cout << "Simulation completed in " << duration.count() << " seconds." << endl;

    return 0;
}