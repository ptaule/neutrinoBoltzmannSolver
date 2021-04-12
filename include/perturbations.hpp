/*
   perturbations.hpp

   Created by Petter Taule on 30.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef PERTURBATIONS_HPP
#define PERTURBATIONS_HPP

#include <cmath>
#include <string>

#include "interpolation.hpp"
#include "measurement.hpp"
#include "integrator.hpp"

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


Measurement compute_psi_l(
        int l,
        double tau,
        double k,
        double q,
        double tau_lambda,
        double metric_psi_at_k_and_tau_lambda,
        const Background& bg,
        const Interpolation2D& metric_psi,
        const Integrator& inner_integrator
        );


class Perturbations {
    private:
        double tau;
        double k;
        double tau_lambda; /* When to turn off EdS approx. */
        double metric_psi_at_k_and_tau_lambda;
        double cutoff;    /* q-integral cutoff */

        /* Interpolate psi_2 if k is above threshold? */
        bool do_interpolate_psi = false;

        const Background& bg;
        const Interpolation2D& metric_psi;

        const Integrator& inner_integrator;
        const Integrator& outer_integrator;

        Measurement integrate_fluid_background(double (*integrand)(double, void*));
        Measurement integrate_fluid_perturbation(
                double (*integrand)(double, void*),
                double q_min = 0,
                const Interpolation1D* psi = nullptr
                );

        void interpolate_psi_l(int l, double q_min, Interpolation1D& psi_spline);
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
                double z_lambda,
                double cutoff,
                const Background& bg,
                const Interpolation2D& metric_psi,
                const Integrator& inner_integrator,
                const Integrator& outer_integrator,
                bool do_interpolate_psi = false
                );

        void compute();
};

#endif /* ifndef PERTURBATIONS_HPP */
