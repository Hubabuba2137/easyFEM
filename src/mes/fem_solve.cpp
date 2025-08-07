#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <chrono>

#include "structs.h"
#include "matrix.h"

#include "gauss.cpp"
#include "fem_funcs.cpp"
#include "write_to_vtu.cpp"

void solve(Fem::Solution solution, bool wrtie_to_vtu=false, bool write_temperature=false){
    Fem::GlobalData conf = solution.conf;
    int node_number = conf.node_number;

    std::vector<double> t0(node_number, conf.init_temperature);
    std::vector<double> t1 = t0;

    Fem::Matrix H_C(node_number, node_number);
    for(int i=0; i<node_number; i++){
        for(int j=0; j<node_number;j++){
            H_C[i][j] = solution.Global_H[i][j]+solution.Global_C[i][j]/conf.time_step;
        }
    }

    std::chrono::steady_clock::time_point for_beg = std::chrono::steady_clock::now();
    //główna pętla
    std::cout<<"Beginning time integration...\n";
    for(double time = conf.time_step; time <= conf.total_time; time+=conf.time_step){
        for(int row=0; row<node_number; row++){
            t1[row] = 0;
            for(int col=0; col<node_number; col++){
                t1[row] += (solution.Global_C[row][col] / conf.time_step) * t0[col];
            }
            t1[row] += solution.Global_P[row][0];
        }

        t1 = Gauss(H_C, t1);

        if(write_temperature){
            std::cout << "Temperature at time " << time << "s:\n";
            std::cout << "\tMIN: " << *std::min_element(t1.begin(), t1.end()) << " MAX: " << *std::max_element(t1.begin(), t1.end()) << "\n";
        }

        if(wrtie_to_vtu){
            write_to_vtu_file((int)time, solution.nodes, t1, solution.elements);
        }
        t0=t1;
    }

    std::chrono::steady_clock::time_point for_end = std::chrono::steady_clock::now();
    std::cout << "Integration time = " << std::chrono::duration_cast<std::chrono::milliseconds>(for_end - for_beg).count() << "[ms]" << std::endl;
}
