/*
   main.cpp

   Created by Petter Taule on 29.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <getopt.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

#include "io.hpp"
#include "perturbations.hpp"

#ifndef DEBUG
#define DEBUG 0
#endif

void grid(std::vector<double>& k_grid, std::vector<double>& z_grid,
        size_t& k_index_max, size_t& z_index_max);

int main(int argc, char* argv[]) {

    bool debug = false;
    double cutoff = 40;
    size_t k_index = 0;
    size_t z_index_a = 0;
    size_t z_index_b = 0;
    double z_lambda = 5;
    double m_nu = 0.07;
    std::string m_nu_string = "0.07";
    double rel_tol = 1e-4;

    int c = 0;
    while ((c = getopt(argc, argv, "dk:l:m:r:")) != -1) {
        switch (c) {
            case 'd':
                // Debug mode on
                debug = true;
                break;
            case 'c':
                cutoff = atof(optarg);
                break;
            case 'k':
                k_index = atoi(optarg);
                break;
            case 'l':
                z_lambda = atof(optarg);
                break;
            case 'm':
                m_nu_string = optarg;
                m_nu = atof(optarg);
                break;
            case 'r':
                rel_tol = atof(optarg);
                break;
            case '?':
                if (optopt == 'k' || optopt == 'l' ||
                        optopt == 'm' || optopt == 'r' )
                {
                    std::cerr << "Option " << optopt << " requires a keyword as \
                        argument." << std::endl;
                    return EXIT_FAILURE;
                }
                else {
                    std::cerr << "Unknown option." << std::endl;
                    return EXIT_FAILURE;
                }
            default:
                return EXIT_FAILURE;
        }
    }

    int argOffset = optind;
    if (argc == argOffset + 1) {
        // One z index given: compute perturbations for one redshift grid point
        z_index_a = atoi(argv[argOffset]);
        z_index_b = z_index_a;
    }
    else if (argc == argOffset + 2) {
        // Two z indices given: compute perturbations for range in redshift grid
        z_index_a = atoi(argv[argOffset]);
        z_index_b = atoi(argv[argOffset + 1]);
    }
    else {
        std::cerr << "The number of options and arguments given is not correct."
            << std::endl;
        return EXIT_FAILURE;
    }

    // If DEBUG defined through compilation options or run-time options, do not
    // turn off GSL error handling
#if DEBUG == 0
    if (debug == false) {
        gsl_set_error_handler_off();
    }
#endif

    std::vector<double> k_grid;
    std::vector<double> z_grid;
    size_t k_index_max;
    size_t z_index_max;
    grid(k_grid, z_grid, k_index_max, z_index_max);

    // Check that indices are within range (size_t variables is alwas positive)
    if (k_index > k_index_max) {
        std::cout << "k_index out of range! Exiting." << std::endl;
        return EXIT_FAILURE;
    }
    if (z_index_a > z_index_max) {
        std::cout << "z_index_a out of range! Exiting." << std::endl;
        return EXIT_FAILURE;
    }
    else if (z_index_b > z_index_max) {
        std::cout << "z_index_b out of range! Exiting." << std::endl;
        return EXIT_FAILURE;
    }
    else if (z_index_b < z_index_a) {
        std::cout << "z_index_b is less than z_index_a! Exiting." << std::endl;
        return EXIT_FAILURE;
    }

    double k = k_grid[k_index];

    // Take subset of z_grid depending on arguments
    if (z_index_a == z_index_b) {
        double z = z_grid[z_index_a];
        z_grid.clear();
        z_grid.push_back(z);
    }
    else {
        std::vector<double>::const_iterator first = z_grid.begin() + z_index_a;
        std::vector<double>::const_iterator last = z_grid.begin() + z_index_b;
        z_grid = std::vector<double>(first, last);
    }

    std::vector<Quantity> delta(z_grid.size());
    std::vector<Quantity> sigma(z_grid.size());
    std::vector<Quantity> cs2(z_grid.size());

    // Set up perturbation computation
    Constants constants;
    constants.m_nu = m_nu;

    std::string CLASS_path =
        "/home/pettertaule/repos/class_public/output/massive_nu_" + m_nu_string
        + "eV/noFA_lmax_17/";

    int outer_sub_regions = 10000;
    int inner_sub_regions = 20000;
    gsl_integration_workspace* outer_workspace =
        gsl_integration_workspace_alloc(outer_sub_regions);
    gsl_integration_workspace* inner_workspace =
        gsl_integration_workspace_alloc(inner_sub_regions);

    Interpolations interpols(CLASS_path, constants, outer_workspace,
            outer_sub_regions, 1e-8, 0);

    for (size_t i = 0; i < z_grid.size(); ++i) {
        double tau = gsl_spline_eval(interpols.tau_of_redshift_spline,
                z_grid[i], interpols.tau_of_redshift_acc);

        Perturbations perturbs(tau, k * constants.h, constants, interpols,
                z_lambda, cutoff, outer_workspace, inner_workspace,
                outer_sub_regions, inner_sub_regions, rel_tol, 0.0, rel_tol,
                0.0);
        perturbs.compute();

        delta[i] = perturbs.delta;
        sigma[i] = perturbs.sigma;
        cs2[i] = perturbs.cs2;
    }

    std::string output_info =
        "_k_" + std::to_string(k_index) + "_z_" + std::to_string(z_index_a);
    write_results("output_m_nu_" + m_nu_string + "eV/delta" +
            output_info + ".dat", z_grid, delta, k);
    write_results("output_m_nu_" + m_nu_string + "eV/sigma" +
            output_info + ".dat", z_grid, sigma, k);
    write_results("output_m_nu_" + m_nu_string + "eV/cs2" + output_info +
            ".dat", z_grid, cs2, k);

    gsl_integration_workspace_free(outer_workspace);
    gsl_integration_workspace_free(inner_workspace);

    return 0;
}



void grid(std::vector<double>& k_grid, std::vector<double>& z_grid,
        size_t& k_index_max, size_t& z_index_max) {
    z_grid = { 0.0, 0.54, 1.08, 1.62, 2.16, 2.7, 3.24,
        3.78, 4.32, 4.86, 5.4, 5.94, 6.48, 7.02, 7.56, 8.1, 8.64, 9.18, 9.72,
        10.3, 10.8, 11.3, 11.9, 12.4, 13.0, 13.5, 14.0, 14.6, 15.1, 15.7, 16.2,
        16.7, 17.3, 17.8, 18.4, 18.9, 19.4, 20.0, 20.5, 21.1, 21.6, 22.1, 22.7,
        23.2, 23.8, 24.3, 24.8, 25.4, 25.9, 26.5
    };

    k_grid = {
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

    k_index_max = k_grid.size() - 1;
    z_index_max = z_grid.size() - 1;
}
