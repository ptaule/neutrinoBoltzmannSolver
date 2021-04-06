/*
   io.cpp

   Created by Petter Taule on 29.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>

#include "../include/io.hpp"
#include "../include/quantity.hpp"


void read_file_and_interpolate(
        const std::string& filename,
        int column_a,
        int column_b,
        std::vector<double>& x,
        std::vector<double>& y
        )
{
    // Remove potential content of x and y vectors
    x.clear();
    y.clear();

    std::ifstream input(filename);
    std::string line;

    if (input.fail()) {
        std::cerr << "Could not read file: " << filename << ".\n";
        input.close();
        return;
    }

    while(getline(input, line)) {
        // Ignore lines beginning with #
        if (line[0] == '#') continue;

        std::stringstream ss(line);
        double value;

        // Counter
        int i = 0;
        while(ss >> value) {
            // If column is equal to column_a, store value in x
            if (i == column_a) {
                x.push_back(value);
            }
            if (i == column_b) {
                y.push_back(value);
            }
            ++i;
        }
    }

    // Sizes should be equal
    if (x.size() != y.size()) {
        std::cerr << "Dimension mismatch: x-grid and data do not have equal \
            sizes!. Exiting." << std::endl;
        input.close();
        return;
    }

    input.close();
}



void read_file_and_interpolate2d(
        const std::string& x_grid_file,
        const std::string& y_grid_file,
        const std::string& data_file,
        std::vector<double>& x,
        std::vector<double>& y,
        std::vector<std::vector<double>>& z
        )
{
    // Remove potential content of x,y,z vectors
    x.clear();
    y.clear();
    z.clear();

    // Read x-grid
    std::ifstream input(x_grid_file);
    std::string line;

    if (input.fail()) {
        std::cerr << "Could not read file: " << x_grid_file << ".\n";
        input.close();
        return;
    }

    while(getline(input, line)) {
        // Ignore lines beginning with #
        if (line[0] == '#') continue;

        std::stringstream ss(line);
        double value;
        // Counter
        int i = 0;

        while(ss >> value) ++i;

        if (i != 1) {
            std::cerr << "More than one column in " << x_grid_file <<
                ". Exiting." << std::endl;
            input.close();
            return;
        }
        x.push_back(value);
    }
    size_t x_size = x.size();
    input.close();

    // Read y-grid
    input.open(y_grid_file);
    if (input.fail()) {
        std::cerr << "Could not read file: " << y_grid_file << ".\n";
        input.close();
        return;
    }

    while(getline(input, line)) {
        // Ignore lines beginning with #
        if (line[0] == '#') continue;

        std::stringstream ss(line);
        double value;
        // Counter
        int i = 0;

        while(ss >> value) ++i;

        if (i != 1) {
            std::cerr << "More than one column in " << y_grid_file <<
                ". Exiting." << std::endl;
            input.close();
            return;
        }
        y.push_back(value);
    }
    size_t y_size = y.size();
    input.close();

    // Read z-grid (data-grid)
    input.open(data_file);
    if (input.fail()) {
        std::cerr << "Could not read file: " << data_file << ".\n";
        input.close();
        return;
    }

    z.resize(x_size);
    for (auto& el : z) {
        el.resize(y_size);
    }

    size_t x_idx = 0;
    while(getline(input, line)) {
        // Ignore lines beginning with #
        if (line[0] == '#') continue;

        if (x_idx >= x_size) {
            std::cerr << "Dimension mismatch between x-grid and first dimension \
                of z-grid in " << data_file << ". Exiting." << std::endl;
            input.close();
            return;
        }

        std::stringstream ss(line);
        double value;

        // Counter
        size_t y_idx = 0;

        while(ss >> value) {
            if (y_idx >= y_size) {
                std::cerr << "Dimension mismatch between x-grid and second \
                    dimension of z-grid in " << data_file << ". Exiting." <<
                    std::endl;
                input.close();
                return;
            }
            z[x_idx][y_idx] = value;
            ++y_idx;
        }
        ++x_idx;
    }
    input.close();
}



void write_results(
        const std::string& filename,
        const std::vector<double>& z_vals,
        const std::vector<Quantity>& data,
        double k
        )
{
    std::ofstream output(filename);

    if (output.fail()) {
        std::cerr << "Could not write to file: " << filename << ".\n";
        output.close();
        return;
    }

    output << "# k = " << k << " h/Mpc" << std::endl;

    size_t size = z_vals.size();
    if (size != data.size()) {
        std::cerr << "z_vals and data have unequal dimensions." << std::endl;
    }

    for (size_t i = 0; i < size; ++i) {
        output << z_vals[i] << "\t\t";
        output << data[i] << std::endl;
    }

    output.close();
}
