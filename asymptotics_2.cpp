#include <fstream>
#include <iostream>
#include <iomanip>
#include "QuantumGraph.hpp"

int main(int argc, char const* argv[])
{
    constexpr double upperBound = 25000.;
    constexpr double lowerBound = 24000.;
    constexpr double step = 50;
    constexpr size_t depth = 30;
    constexpr size_t numberOfNobes = 1000;

    QuantumGraph graph1{
        [](double) {return 0.;},
        [](double) {return 0.;},
        [](double) {return 0.;}
    };
    graph1.setNumberOfNodes(numberOfNobes);
    
    std::clog << "Starting eigenvalue calculation for graph1\n";
    std::vector<double> eigenvalues1 = graph1.calcEigenvalues(lowerBound, upperBound, step, depth);
    std::clog << "Finished eigenvalue calculation for graph1\n";

    QuantumGraph graph2{
        [](double x) {return x;},
        [](double x) {return -x;},
        [](double x) {return -sqrt(x);}
    };
    graph2.setNumberOfNodes(numberOfNobes);
    
    std::clog << "Starting eigenvalue calculation for graph2\n";
    std::vector<double> eigenvalues2 = graph2.calcEigenvalues(lowerBound, upperBound, step, depth);
    std::clog << "Finished eigenvalue calculation for graph2\n";

    QuantumGraph graph3{
        [](double x) {return sin(x);},
        [](double x) {return exp(-x);},
        [](double x) {return exp(-x * x);}
    };
    graph3.setNumberOfNodes(numberOfNobes);
    
    std::clog << "Starting eigenvalue calculation for graph3\n";
    std::vector<double> eigenvalues3 = graph3.calcEigenvalues(lowerBound, upperBound, step, depth);
    std::clog << "Finished eigenvalue calculation for graph3\n";

    std::ofstream ofs{ "ev-24000-25000.csv" };
	ofs << std::ios::fixed << std::setprecision(10);
    ofs << "graph#1";
    for (auto i : eigenvalues1)
    {
        ofs << ", " << i;
    }
    ofs << "\ngraph#2";
    for (auto i : eigenvalues2)
    {
        ofs << ", " << i;
    }
    ofs << "\ngraph#3";
    for (auto i : eigenvalues3)
    {
        ofs << ", " << i;
    }

    std::clog << "Program finished";
    return 0;
}
