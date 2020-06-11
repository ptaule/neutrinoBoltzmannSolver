/*
   io.hpp

   Created by Petter Taule on 29.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef IO_HPP
#define IO_HPP

#include <string>
#include <vector>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

void read_file_and_interpolate(
        const std::string& filename,
        int column_a,
        int column_b,
        gsl_interp_accel** acc,
        gsl_spline** spline,
        double (*x_op)(double)=NULL, /* Optional function applied to x */
        double (*y_op)(double)=NULL  /* Optional function applied to y */
        );

class Quantity;

void write_results(
        const std::string& filename,
        const std::vector<double>& z_vals,
        const std::vector<Quantity>& data,
        double k
        );

double redshift_to_scale_factor(double redshift);

#endif /* ifndef IO_HPP */
