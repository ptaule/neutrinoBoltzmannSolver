/*
   main.cpp

   Created by Petter Taule on 29.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <cmath>
#include <vector>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

#include "io.hpp"
#include "perturbations.hpp"


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Number of arguments should be 2! Exiting." << std::endl;
        return 1;
    }

    int k_index = atoi(argv[1]);
    int z_idx_shift = atoi(argv[2]);

    std::vector<double> z_vals = { 0.0, 0.54, 1.08, 1.62, 2.16, 2.7, 3.24,
        3.78, 4.32, 4.86, 5.4, 5.94, 6.48, 7.02, 7.56, 8.1, 8.64, 9.18, 9.72,
        10.3, 10.8, 11.3, 11.9, 12.4, 13.0, 13.5, 14.0, 14.6, 15.1, 15.7, 16.2,
        16.7, 17.3, 17.8, 18.4, 18.9, 19.4, 20.0, 20.5, 21.1, 21.6, 22.1, 22.7,
        23.2, 23.8, 24.3, 24.8, 25.4, 25.9, 26.5
    };

    std::vector<double>::const_iterator first = z_vals.begin() + z_idx_shift;
    std::vector<double>::const_iterator last = first + 10;
    std::vector<double> z_vals_subset(first, last);

    double k_vals[100] = {
        1.0e-5, 1.18e-5, 1.39e-5, 1.64e-5, 1.93e-5, 2.28e-5,
        2.69e-5, 3.17e-5, 3.73e-5, 4.4e-5, 5.19e-5, 6.12e-5,
        7.21e-5, 8.5e-5, 0.0001, 0.000118, 0.000139, 0.000164,
        0.000194, 0.000228, 0.000269, 0.000317, 0.000374, 0.000441,
        0.00052, 0.000613, 0.000723, 0.000852, 0.00101, 0.00118,
        0.0014, 0.00165, 0.00194, 0.00229, 0.0027, 0.00318,
        0.00375, 0.00442, 0.00522, 0.00615, 0.00725, 0.00855,
        0.0101, 0.0119, 0.014, 0.0165, 0.0195, 0.023, 0.0271,
        0.0319, 0.0376, 0.0443, 0.0523, 0.0616, 0.0727, 0.0857,
        0.101, 0.119, 0.14, 0.166, 0.195, 0.23, 0.271, 0.32, 0.377,
        0.445, 0.524, 0.618, 0.729, 0.859, 1.01, 1.19, 1.41, 1.66,
        1.96, 2.31, 2.72, 3.21, 3.78, 4.46, 5.25, 6.2, 7.3, 8.61,
        10.2, 12.0, 14.1, 16.6, 19.6, 23.1, 27.3, 32.1, 37.9, 44.7,
        52.7, 62.1, 73.2, 86.3, 102.0, 120.0
    };

    std::vector<Quantity> sigma_over_delta(10);
    std::vector<Quantity> cs2(10);

    // Set up perturbation computation
    Constants constants;
    constants.m_nu = 0.05;

    std::string class_dir =
        "/home/pettertaule/repos/class_public/output/massive_nu_0.05eV/noFA_lmax_17/";

    int outer_sub_regions = 1000;
    int inner_sub_regions = 20000;
    gsl_integration_workspace* outer_workspace =
        gsl_integration_workspace_alloc(outer_sub_regions);
    gsl_integration_workspace* inner_workspace =
        gsl_integration_workspace_alloc(inner_sub_regions);

    Interpolations interpols(class_dir + "background.dat", constants,
            outer_workspace, outer_sub_regions, 1e-8, 0);

    double k = k_vals[k_index];

    double cutoff = 40;
    double inner_eps_rel = 1e-4;
    double inner_eps_abs = 0;
    double outer_eps_rel = 1e-4;
    double outer_eps_abs = 0;

    for (int i = 0; i < 10; ++i) {
        double tau = gsl_spline_eval(interpols.tau_of_redshift_spline,
                z_vals_subset[i], interpols.tau_of_redshift_acc);

        Perturbations perturbs(tau, k * constants.h, constants, interpols,
                cutoff, outer_workspace, inner_workspace, outer_sub_regions,
                inner_sub_regions, outer_eps_rel, outer_eps_abs, inner_eps_rel,
                inner_eps_abs);
        perturbs.compute();

        sigma_over_delta[i] = (perturbs.sigma / perturbs.delta);
        cs2[i] = perturbs.cs2;
    }

    std::string output_info =
        "_k_" + std::string(argv[1]) + "_z_" + std::string(argv[2]);
    write_results("output/sigma_over_delta" + output_info + ".dat", z_vals_subset,
            sigma_over_delta, k);
    write_results("output/cs2" + output_info + ".dat", z_vals_subset, cs2, k);

    gsl_integration_workspace_free(outer_workspace);
    gsl_integration_workspace_free(inner_workspace);

    return 0;
}
