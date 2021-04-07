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
#include "../include/measurement.hpp"

#define RESERVE_SIZE 500

template <class T>
using Vec1D = std::vector<T>;
template <class T>
using Vec2D = std::vector<std::vector<T>>;


/* Read n_columns number of columns from file , if there are lines with less
 * columns, throw error. Skips empty lines and lines beginning with '#' */
void read_columns_from_file(
        const std::string& filename,  /* in, name of file to be read */
        size_t n_columns,             /* in, number of columns       */
        Vec2D<double>& columns        /* out, columns                */
        )
{
    std::ifstream input(filename);
    std::string line;

    if (input.fail()){
        throw(std::invalid_argument("Could not find " + filename + "."));
    }

    columns.clear();
    columns.resize(n_columns);
    for (size_t i = 0; i < n_columns; ++i) {
        columns.at(i).reserve(RESERVE_SIZE);
    }

    while (getline(input, line)) {
        /* Ignore empty lines or lines beginning with # */
        if (line.empty() || line.at(0) == '#') {
            continue;
        }

        std::stringstream ss(line);
        double value;

        // Column counter
        size_t i = 0;
        while (ss >> value && i < n_columns) {
            columns.at(i++).push_back(value);
        }
        if (i != n_columns) {
            throw(std::runtime_error("Number of columns less than n_columns in " +
                                     filename));
        }
    }
    input.close();
}



/* Reads data from filename, using indexing scheme:
 * index = row + column * n_rows
 * If the number of rows/columns found does not equal n_rows/n_columns, an
 * error is thrown */
void read_data_grid_from_file(
        const std::string& filename,
        Vec1D<double>& data,
        size_t n_rows,
        size_t n_columns
        )
{
    std::ifstream input;
    std::string line;

    input.open(filename);
    if (input.fail()){
        throw(std::invalid_argument("Could not find " + filename + "."));
    }

    size_t row_idx = 0;
    size_t column_idx = 0;
    size_t idx = 0;

    data.resize(n_rows * n_columns);

    while (getline(input, line)) {
        /* Ignore empty lines or lines beginning with # */
        if (line.empty() || line.at(0) == '#') {
            continue;
        }

        if (row_idx >= n_rows) {
            throw(std::runtime_error("Number of rows exceeds n_rows."));
        }

        double value;
        std::stringstream ss(line);

        column_idx = 0;
        idx = row_idx;

        while (ss >> value) {
            if (column_idx >= n_columns) {
                throw(
                    std::runtime_error("Number of columns exceeds n_columns."));
            }

            data.at(idx) = value;
            column_idx++;
            idx = row_idx + column_idx * n_rows;
        }
        if (column_idx != n_columns) {
            throw(std::runtime_error(
                "Number of columns does not equal n_columns."));
        }
        row_idx++;
    }
    if (row_idx != n_rows) {
        throw(std::runtime_error("Number of columns does not equal n_rows."));
    }
    input.close();
}



void write_results(
        const std::string& filename,
        const std::vector<double>& z_vals,
        const std::vector<Measurement>& data,
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
