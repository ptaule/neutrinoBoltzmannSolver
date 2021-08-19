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

#include "include/io.hpp"
#include "include/measurement.hpp"
#include "include/integrator.hpp"
#include "include/interpolation.hpp"
#include "include/perturbations.hpp"

#ifndef DEBUG
#define DEBUG 0
#endif

using size_t = std::size_t;

template <class T>
using Vec1D = std::vector<T>;
template <class T>
using Vec2D = std::vector<std::vector<T>>;


int main(int argc, char* argv[]) {
    bool debug              = false;
    bool do_interpolate_psi = false;
    double cutoff           = 40;
    double abs_tol          = 0;
    double rel_tol          = 1e-4;
    size_t k_index          = 0;
    size_t z_index_a        = 0;
    size_t z_index_b        = 0;
    double z_lambda         = 5;
    double m_nu             = 0.05;
    std::string k_grid_file = "input/k_grid.dat";
    std::string z_grid_file = "input/z_grid.dat";
    std::string output_path = "output/test/";
    std::string CLASS_path  =
        "/space/ge52sir/class_public/output/m_nu_0.05eV/noFA_lmax_17/";


    int c = 0;
    while ((c = getopt(argc, argv, "a:c:C:dik:K:l:m:o:r:Z:")) != -1) {
        switch (c) {
            case 'a':
                abs_tol = atof(optarg);
                break;
            case 'c':
                cutoff = atof(optarg);
                break;
            case 'C':
                CLASS_path = optarg;
                break;
            case 'd':
                // Debug mode on
                debug = true;
                break;
            case 'i':
                // Interpolate psi2 above k threshold
                do_interpolate_psi = true;
                break;
            case 'k':
                k_index = static_cast<size_t>(atoi(optarg));
                break;
            case 'K':
                k_grid_file = optarg;
                break;
            case 'l':
                z_lambda = atof(optarg);
                break;
            case 'm':
                m_nu = atof(optarg);
                break;
            case 'o':
                output_path = optarg;
                break;
            case 'r':
                rel_tol = atof(optarg);
                break;
            case 'Z':
                z_grid_file = optarg;
                break;
            case '?':
                if (
                        optopt == 'a' ||
                        optopt == 'c' ||
                        optopt == 'C' ||
                        optopt == 'k' ||
                        optopt == 'K' ||
                        optopt == 'l' ||
                        optopt == 'm' ||
                        optopt == 'o' ||
                        optopt == 'r' ||
                        optopt == 'Z')
                {
                    std::cerr << "Option " << static_cast<char>(optopt)
                              << " requires a keyword as argument." << std::endl;
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
        /* One z index given: compute perturbations for one redshift grid point */
        z_index_a = static_cast<size_t>(atoi(argv[argOffset]));
        z_index_b = z_index_a;
    }
    else if (argc == argOffset + 2) {
        // Two z indices given: compute perturbations for range in redshift grid
        z_index_a = static_cast<size_t>(atoi(argv[argOffset]));
        z_index_b = static_cast<size_t>(atoi(argv[argOffset + 1]));
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

    Vec1D<double> k_grid, z_grid;
    Vec2D<double> data;

    read_columns_from_file(k_grid_file, 1, data);
    k_grid = data[0];
    data.clear();
    read_columns_from_file(z_grid_file, 1, data);
    z_grid = data[0];

    // Check that indices are within range (size_t variables are alwas positive)
    if (k_index >= k_grid.size()) {
        std::cout << "k_index out of range! Exiting." << std::endl;
        return EXIT_FAILURE;
    }
    if (z_index_a > z_grid.size()) {
        std::cout << "z_index_a out of range! Exiting." << std::endl;
        return EXIT_FAILURE;
    }
    else if (z_index_b >= z_grid.size()) {
        std::cout << "z_index_b out of range! Exiting." << std::endl;
        return EXIT_FAILURE;
    }
    else if (z_index_b < z_index_a) {
        std::cout << "z_index_b is less than z_index_a! Exiting." << std::endl;
        return EXIT_FAILURE;
    }

    /* k in h/Mpc */
    double k = k_grid[k_index];

    /* Take subset of z_grid depending on arguments */
    if (z_index_a == z_index_b) {
        double z = z_grid[z_index_a];
        z_grid.clear();
        z_grid.push_back(z);
    }
    else {
        Vec1D<double>::const_iterator first =
            z_grid.begin() + static_cast<long int>(z_index_a);
        Vec1D<double>::const_iterator last =
            z_grid.begin() + static_cast<long int>(z_index_b) + 1;
        z_grid = Vec1D<double>(first, last);
    }

    std::vector<Measurement> delta(z_grid.size());
    std::vector<Measurement> sigma(z_grid.size());
    std::vector<Measurement> cs2(z_grid.size());

    /* Set up perturbation computation */
    double hubble = 0.67556;
    double T_ncdm0 = 0.00016773497604334795;
    double tau_ini = 1;

    Integrator outer_integrator(10000, abs_tol, rel_tol);
    Integrator inner_integrator(20000, abs_tol, rel_tol);

    try {
        Background bg(hubble, m_nu, T_ncdm0, tau_ini, CLASS_path +
                "/background.dat");
        Interpolation2D metric_psi;
        interpolate_metric_psi(
            hubble, bg.conf_time_of_redshift, CLASS_path + "/k_grid.dat",
            CLASS_path + "/reverse_z_grid.dat", CLASS_path + "/reverse_psi.dat", metric_psi);

        for (size_t i = 0; i < z_grid.size(); ++i) {
            std::cout << "Iteration " << i << ": z = " << z_grid[i] << ".\t";
            double tau = bg.conf_time_of_redshift(z_grid[i]);

            Perturbations perturbs(tau, k* hubble, z_lambda, cutoff, bg,
                    metric_psi, inner_integrator, outer_integrator,
                    do_interpolate_psi);

            perturbs.compute();

            delta[i] = perturbs.delta;
            sigma[i] = perturbs.sigma;
            cs2[i] = perturbs.cs2;
        }
    }
    catch (const IntegrationException& ex) {
        std::cerr << ex.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return EXIT_FAILURE;
    }

    std::string output_info =
        "_k_" + std::to_string(k_index) + "_z_" + std::to_string(z_index_a);
    write_results(output_path + "/delta" + output_info + ".dat", z_grid, delta, k);
    write_results(output_path + "/sigma" + output_info + ".dat", z_grid, sigma, k);
    write_results(output_path + "/cs2" + output_info + ".dat", z_grid, cs2, k);

    return 0;
}
