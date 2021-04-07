/*
   perturbations.hpp

   Created by Petter Taule on 30.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef PERTURBATIONS_HPP
#define PERTURBATIONS_HPP

#include <cmath>
#include <string>
#include <stdexcept>

#include <gsl/gsl_integration.h>

#include "measurement.hpp"
#include "interpolation.hpp"


class IntegrationException : std::exception {
    private:
        bool critical_;
        std::string ex;
    public:
        IntegrationException(int gsl_status, const std::string& integration_name);
        bool critical() const {return critical_;}
        const char* what() const noexcept {return ex.c_str();}
};


class Background {
    private:
        void compute_y_function();
    public:
        double hubble;
        double m_nu;
        double T_ncdm0;
        double m_nu_over_T_ncdm0_square;
        double tau_ini;

        Interpolation1D scale_factor_of_conf_time;
        Interpolation1D conf_time_of_redshift;
        Interpolation2D y_function;

        Background(
                double h,
                double m_nu,
                double T_ncdm0,
                double tau_ini,
                const std::string& background_file
                );
};


void interpolate_metric_psi(
        double hubble,
        const Interpolation1D& conf_time_of_redshift,
        const std::string& k_grid_file,
        const std::string& z_grid_file,
        const std::string& metric_psi_file,
        Interpolation2D& metric_psi
        );



struct Psi_parameters {
        double tau = 0;
        double k = 0;
        double q = 0;
        double tau_lambda = 0;
        double metric_psi_at_k_and_tau_lambda = 0;
        const Background& bg;
        const Interpolation2D& metric_psi;
        gsl_integration_workspace* workspace = nullptr;
        std::size_t sub_regions;
        double rel_tol;
        double abs_tol;
};


void psi_0(const Psi_parameters& psi_params, Measurement& result);
void psi_1(const Psi_parameters& psi_params, Measurement& result);
void psi_2(const Psi_parameters& psi_params, Measurement& result);


class Perturbations {
    private:
        double tau = 0;
        double k = 0;
        const Background& bg;
        const Interpolation2D& metric_psi;
        double tau_lambda = 0; /* When to turn off EdS approx. */
        double metric_psi_at_k_and_tau_lambda = 0;
        double cutoff = 40;    /* q-integral cutoff */

        /* Interpolate psi_2 if k is above threshold? */
        bool do_interpolate_psi = false;

        gsl_integration_workspace* outer_workspace;
        gsl_integration_workspace* inner_workspace;
        size_t outer_sub_regions;
        size_t inner_sub_regions;
        double outer_rel_tol;
        double outer_abs_tol;
        double inner_rel_tol;
        double inner_abs_tol;

        Measurement integrate_fluid_background(double (*integrand)(double, void*));
        Measurement integrate_fluid_perturbation(
                double (*integrand)(double, void*),
                double q_min = 0,
                const Interpolation1D* psi = nullptr
                );

        void interpolate_psi(
                void (*psi)(
                    const Psi_parameters& psi_params,
                    Measurement& result
                    ),
                double q_min,
                Interpolation1D& psi_spline
                );
    public:
        Measurement rho;
        Measurement pressure;
        Measurement delta;
        Measurement theta;
        Measurement sigma;
        Measurement cs2;

        Perturbations(
                double tau,
                double k,
                const Background& bg,
                const Interpolation2D& metric_psi,
                double z_lambda,
                double cutoff,
                bool do_interpolate_psi,
                gsl_integration_workspace* outer_workspace,
                gsl_integration_workspace* inner_workspace,
                std::size_t outer_sub_regions,
                std::size_t inner_sub_regions,
                double outer_rel_tol,
                double outer_abs_tol,
                double inner_rel_tol,
                double inner_abs_tol
                );

        void compute();
};

#endif /* ifndef PERTURBATIONS_HPP */
