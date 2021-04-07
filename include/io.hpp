/*
   io.hpp

   Created by Petter Taule on 29.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef IO_HPP
#define IO_HPP

#include <string>
#include <vector>


void read_columns_from_file(
        const std::string& filename,
        size_t n_columns,
        std::vector<std::vector<double>>& columns
        );


void read_data_grid_from_file(
        const std::string& filename,
        std::vector<double>& data,
        size_t n_rows,
        size_t n_columns
        );

class Measurement;

void write_results(
        const std::string& filename,
        const std::vector<double>& z_vals,
        const std::vector<Measurement>& data,
        double k
        );

#endif /* ifndef IO_HPP */
