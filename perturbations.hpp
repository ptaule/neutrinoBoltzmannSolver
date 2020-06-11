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

        Interpolations(
                const std::string& filename,
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


class Perturbations {
    private:
        double tau = 0;
        double k = 0;
        const Constants& constants;
        const Interpolations& interpols;
        double cutoff = 40;

        gsl_integration_workspace* outer_workspace;
        gsl_integration_workspace* inner_workspace;
        int outer_sub_regions;
        int inner_sub_regions;
        double outer_eps_rel;
        double outer_eps_abs;
        double inner_eps_rel;
        double inner_eps_abs;

        Quantity integrate_perturbations(double (*integrand)(double, void*));
        Quantity integrate_background(double (*integrand)(double, void*));
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
                double cutoff,
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


void psi_0(
        double tau,
        double k,
        double q,
        const Constants& constants,
        const Interpolations& interpols,
        Quantity& result,
        gsl_integration_workspace* workspace,
        int sub_regions,
        double eps_rel,
        double eps_abs
        );

void psi_1(
        double tau,
        double k,
        double q,
        const Constants& constants,
        const Interpolations& interpols,
        Quantity& result,
        gsl_integration_workspace* workspace,
        int sub_regions,
        double eps_rel,
        double eps_abs
        );

void psi_2(
        double tau,
        double k,
        double q,
        const Constants& constants,
        const Interpolations& interpols,
        Quantity& result,
        gsl_integration_workspace* workspace,
        int sub_regions,
        double eps_rel,
        double eps_abs
        );


#endif /* ifndef PERTURBATIONS_HPP */
