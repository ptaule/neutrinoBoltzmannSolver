/*
   io.cpp

   Created by Petter Taule on 29.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "io.hpp"
#include "quantity.hpp"

#define INTERPOL_TYPE

void read_file_and_interpolate(
        const std::string& filename,
        int column_a,
        int column_b,
        gsl_interp_accel** acc,
        gsl_spline** spline,
        double (*x_op)(double), /* Optional function applied to x */
        double (*y_op)(double)  /* Optional function applied to y */
        )
{
    std::ifstream input(filename);
    std::string line;

    std::vector<double> x;
    std::vector<double> y;

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
                if (x_op == NULL) {
                    x.push_back(value);
                } else {
                    x.push_back(x_op(value));
                }
            }
            if (i == column_b) {
                if (y_op == NULL) {
                    y.push_back(value);
                } else {
                    y.push_back(y_op(value));
                }
            }
            ++i;
        }
    }

    // Sizes should be equal
    if (x.size() != y.size()) {
        std::cerr << "x and y do not have equal sizes!. Exiting.\n";
        input.close();
        return;
    }
    int size = x.size();

    // If x is descending (e.g. redshift, reverse vector)
    // We assume x strictly ascending/descending

    if (x[0] > x[1]) {
        std::reverse(x.begin(), x.end());
        std::reverse(y.begin(), y.end());
    }

    // Interpolation

    *acc = gsl_interp_accel_alloc();
    *spline = gsl_spline_alloc(gsl_interp_cspline, size);

    gsl_spline_init(*spline, x.data(), y.data(), size);

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


double redshift_to_scale_factor(double redshift) {
    return 1/(redshift + 1);
}
