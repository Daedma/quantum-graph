#include "QuantumGraph.hpp"
#include <chrono>
#include <fstream>
#include <thread>
#include <vector>
#include <functional>
#include <mutex>

void runExperiment(const QuantumGraph& graph, const std::string& graphName, double lowerBound, double higherBound, double step, double error, size_t windowSize, const std::string& algorithm, std::mutex& fileMutex, std::ofstream& outputFile)
{
	constexpr size_t numRuns = 20;
	size_t eigenvalueCount = 0;
	double totalTime = 0.0;

	for (size_t i = 0; i < numRuns; ++i)
	{
		auto start = std::chrono::high_resolution_clock::now();
		std::vector<double> eigenvalues;

		if (algorithm == "fixed_step")
		{
			eigenvalues = graph.calcEigenvalues(lowerBound, higherBound, step, error);
		}
		else if (algorithm == "adaptive_step")
		{
			eigenvalues = graph.calcEigenvalues(lowerBound, higherBound, windowSize, step, error);
		}

		auto end = std::chrono::high_resolution_clock::now();
		totalTime += std::chrono::duration<double>(end - start).count();
		eigenvalueCount = eigenvalues.size();
	}

	double avgTime = totalTime / numRuns;

	// Write results to CSV
	std::lock_guard<std::mutex> lock(fileMutex);
	outputFile << graphName << "," << algorithm << "," << lowerBound << "," << higherBound << "," << step << "," << eigenvalueCount << "," << avgTime << "\n";
}

int main()
{
// Define three graphs with different potentials
	QuantumGraph graph1([](double x) { return x * x; }, [](double x) { return x + 1; }, [](double x) { return std::sin(x); });
	QuantumGraph graph2([](double x) { return std::cos(x); }, [](double x) { return std::exp(-x); }, [](double x) { return x; });
	QuantumGraph graph3([](double x) { return 1.0; }, [](double x) { return 2.0; }, [](double x) { return 3.0; });

	// Define parameter sets (lowerBound, higherBound, step)
	std::vector<std::tuple<double, double, double>> parameterSets = {
		{ 0.0, 100.0, 0.1 },
		{ 100.0, 1000.0, 5.0 },
		{ 1000.0, 5000.0, 15.0 }
	};

	constexpr double error = 1e-6;
	constexpr size_t windowSize = 2;

	// Open CSV file for writing results
	std::ofstream outputFile("experiments/adaptive-step.csv");
	outputFile << "Graph,Algorithm,LowerBound,HigherBound,Step,EigenvalueCount,AvgTime\n";

	// Mutex for thread-safe file writing
	std::mutex fileMutex;

	// Run experiments in parallel
	std::vector<std::thread> threads;
	for (const auto& [lowerBound, higherBound, step] : parameterSets)
	{
		threads.emplace_back(runExperiment, graph1, "Graph1", lowerBound, higherBound, step, error, windowSize, "fixed_step", std::ref(fileMutex), std::ref(outputFile));
		threads.emplace_back(runExperiment, graph1, "Graph1", lowerBound, higherBound, step, error, windowSize, "adaptive_step", std::ref(fileMutex), std::ref(outputFile));

		threads.emplace_back(runExperiment, graph2, "Graph2", lowerBound, higherBound, step, error, windowSize, "fixed_step", std::ref(fileMutex), std::ref(outputFile));
		threads.emplace_back(runExperiment, graph2, "Graph2", lowerBound, higherBound, step, error, windowSize, "adaptive_step", std::ref(fileMutex), std::ref(outputFile));

		threads.emplace_back(runExperiment, graph3, "Graph3", lowerBound, higherBound, step, error, windowSize, "fixed_step", std::ref(fileMutex), std::ref(outputFile));
		threads.emplace_back(runExperiment, graph3, "Graph3", lowerBound, higherBound, step, error, windowSize, "adaptive_step", std::ref(fileMutex), std::ref(outputFile));

		for (auto& thread : threads)
		{
			thread.join();
		}

		threads.clear();
	}

	outputFile.close();
	return 0;
}