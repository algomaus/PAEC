#pragma once
#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include "external/gnuplot-iostream.h"

inline void plot(std::unordered_map<double, double> &data, const std::string &xLabel, const std::string &yLabel,
		const std::string &filename) {
	Gnuplot gp("tee " + filename + ".gnu | gnuplot -persist");
	//Gnuplot gp;
	double minX = std::numeric_limits<double>::infinity();
	double maxX = -std::numeric_limits<double>::infinity();
	double minY = std::numeric_limits<double>::infinity();
	double maxY = -std::numeric_limits<double>::infinity();
	std::vector<std::pair<double, double> > dataPoints;
	for (auto kv : data) {
		dataPoints.push_back(std::make_pair(kv.first, kv.second));
		if (kv.first < minX)
			minX = kv.first;
		if (kv.first > maxX)
			maxX = kv.first;
		if (kv.second < minY)
			minY = kv.second;
		if (kv.second > maxY)
			maxY = kv.second;
	}

	std::sort(dataPoints.begin(), dataPoints.end(),
			[](const std::pair<double,double> &left, const std::pair<double,double> &right) {
				return left.first < right.first;
			});

	if (minX == maxX) {
		minX--;
		maxX++;
	}
	if (minY == maxY) {
		minY--;
		maxY++;
	}

	// write data points to file
	std::ofstream fileData;
	fileData.open(filename + "_data.csv");
	for (auto kv : dataPoints) {
		fileData << kv.first << "," << kv.second << std::endl;
	}
	fileData.close();

	gp << "set terminal png size 800,600 enhanced font \"Helvetica,10\"" << std::endl;
	gp << "set output '" << filename << ".png'" << std::endl;
	gp << "set xrange [" << minX << ":" << maxX << "]\nset yrange [" << minY << ":" << maxY << "]\n";
	gp << "set xlabel \"" << xLabel << "\"" << std::endl;
	gp << "set ylabel \"" << yLabel << "\"" << std::endl;

	gp << "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5" << std::endl;
	gp << "plot '-' with linespoints ls 1 notitle" << std::endl;
	gp.send1d(dataPoints);
}

inline std::unordered_map<double, double> extractDatasetFromPlotscript(const std::string &filename) {
	std::unordered_map<double, double> data;
	std::ifstream file;
	file.open(filename);
	bool started = false;
	std::string line;
	while (std::getline(file, line)) {
		if (!started) {
			if (line.find("plot") == 0) {
				started = true;
			}
		} else {
			if (line.size() > 2) { // ignore the ending "e" line
				std::stringstream stream(line);
				double key;
				double value;
				stream >> key >> value;
				data[key] = value;
			}
		}
	}
	file.close();
	return data;
}
