/*
   io.hpp

   Created by Petter Taule on 29.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef IO_HPP
#define IO_HPP

#include <string>
#include <vector>

void read_file_and_interpolate(
        const std::string& filename,
        int column_a,
        int column_b,
        std::vector<double>& x,
        std::vector<double>& y
        );

void read_file_and_interpolate2d(
        const std::string& x_grid_file,
        const std::string& y_grid_file,
        const std::string& data_file,
        std::vector<double>& x,
        std::vector<double>& y,
        std::vector<std::vector<double>>& z
        );

class Measurement;

void write_results(
        const std::string& filename,
        const std::vector<double>& z_vals,
        const std::vector<Measurement>& data,
        double k
        );

#endif /* ifndef IO_HPP */
