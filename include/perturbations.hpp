/*
   perturbations.hpp

   Created by Petter Taule on 30.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef PERTURBATIONS_HPP
#define PERTURBATIONS_HPP

#include <cmath>

#include <gsl/gsl_integration.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

#include "quantity.hpp"

struct Constants {
    double h = 0.67556;
    double m_nu;
    double T_ncdm0 = 0.00016773497604334795;
    double tau_ini = 1;
};


class Interpolations {
    private:
        void z_function(
                Constants& constants,
                double tau_fin,
                gsl_integration_workspace* workspace,
                int sub_regions,
                double eps_rel,
                double eps_abs
                );
    public:
        gsl_interp_accel* a_of_tau_acc           = nullptr;
        gsl_spline*       a_of_tau_spline        = nullptr;
        gsl_interp_accel* tau_of_redshift_acc    = nullptr;
        gsl_spline*       tau_of_redshift_spline = nullptr;

        gsl_interp_accel* z_spline_x_acc         = nullptr;
        gsl_interp_accel* z_spline_y_acc         = nullptr;
        gsl_spline2d*     z_spline               = nullptr;

        gsl_interp_accel* grav_psi_x_acc         = nullptr;
        gsl_interp_accel* grav_psi_y_acc         = nullptr;
        gsl_spline2d*     grav_psi_spline        = nullptr;

        Interpolations(
                const std::string& CLASS_path,
                Constants& constants,
                gsl_integration_workspace* workspace,
                int sub_regions,
                double eps_rel = 1e-8,
                double eps_abs = 0
                );

        ~Interpolations() {
            gsl_interp_accel_free(a_of_tau_acc);
            gsl_spline_free(a_of_tau_spline);

            gsl_interp_accel_free(tau_of_redshift_acc);
            gsl_spline_free(tau_of_redshift_spline);

            gsl_interp_accel_free(z_spline_x_acc);
            gsl_interp_accel_free(z_spline_y_acc);
            gsl_spline2d_free(z_spline);
        }
};


struct Psi_parameters {
        double tau = 0;
        double k = 0;
        double q = 0;
        double tau_lambda = 0;
        double grav_psi_at_k_and_tau_lambda = 0;
        const Constants& constants;
        const Interpolations& interpols;
        gsl_integration_workspace* workspace = nullptr;
        int sub_regions;
        double eps_rel;
        double eps_abs;
};


void psi_0(const Psi_parameters& psi_params, Quantity& result);
void psi_1(const Psi_parameters& psi_params, Quantity& result);
void psi_2(const Psi_parameters& psi_params, Quantity& result);


class Perturbations {
    private:
        double tau = 0;
        double k = 0;
        const Constants& constants;
        const Interpolations& interpols;
        double tau_lambda = 0; // When to turn off EdS approx.
        double grav_psi_at_k_and_tau_lambda = 0;
        double cutoff = 40;    // q-integral cutoff

        // Interpolate psi_2 if k is above threshold?
        bool do_interpolate_psi = false;

        gsl_integration_workspace* outer_workspace;
        gsl_integration_workspace* inner_workspace;
        int outer_sub_regions;
        int inner_sub_regions;
        double outer_eps_rel;
        double outer_eps_abs;
        double inner_eps_rel;
        double inner_eps_abs;

        Quantity integrate_background(double (*integrand)(double, void*));
        Quantity integrate_perturbations(
                double (*integrand)(double, void*),
                double q_min = 0,
                gsl_interp_accel* psi_acc = nullptr,
                gsl_spline* psi_spline = nullptr
                );

        void interpolate_psi(
                void (*psi)(
                    const Psi_parameters& psi_params,
                    Quantity& result
                    ),
                double q_min,
                gsl_interp_accel** psi_acc,
                gsl_spline** psi_spline
                );
    public:
        Quantity rho;
        Quantity pressure;
        Quantity delta;
        Quantity theta;
        Quantity sigma;
        Quantity cs2;

        Perturbations(
                double tau,
                double k,
                Constants& constants,
                Interpolations& interpols,
                double z_lambda,
                double cutoff,
                bool do_interpolate_psi,
                gsl_integration_workspace* outer_workspace,
                gsl_integration_workspace* inner_workspace,
                int outer_sub_regions,
                int inner_sub_regions,
                double outer_eps_rel,
                double outer_eps_abs,
                double inner_eps_rel,
                double inner_eps_abs
                );

        void compute();
};

#endif /* ifndef PERTURBATIONS_HPP */
