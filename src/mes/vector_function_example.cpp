// Example program of how to use vector of funtion pointers
#include <iostream>
#include <vector>
#include <functional>

std::function<float(float)> N1 = [](float x){return 0.25*(1-x);};
std::function<float(float)> N2 = [](float x){return -0.25*(1-x);};
std::function<float(float)> N3 = [](float x){return 0.25*(1+x);};
std::function<float(float)> N4 = [](float x){return -0.25*(1+x);};


int main()
{
    std::vector<float> vals = {1,2,3,4};
    std::vector<std::function<float(float)>> funcs;
    
    funcs.push_back(N1);
    funcs.push_back(N2);
    funcs.push_back(N3);
    funcs.push_back(N4);
    
    for(int i=0; i<funcs.size();i++){
        std::cout<<funcs[i](vals[i])<<"\n";
    }
}