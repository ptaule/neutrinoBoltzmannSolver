/*
   perturbations.cpp

   Created by Petter Taule on 30.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <cmath>
#include <vector>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include "perturbations.hpp"
#include "quantity.hpp"
#include "io.hpp"


void integration_status(
        int status,
        const std::string& integration_info,
        Quantity& result
        )
{
    switch (status) {
        case GSL_EMAXITER:
            std::cerr << "Error during integration of " << integration_info
                << ": Maxiumum number of subdevisions reached.  Aborting."
                << std::endl;
            abort();
        case GSL_EROUND:
            std::cerr << "Warning during integration of " << integration_info
                << ": Rounding error. Result: " << result
                << std::endl;
            return;
        case GSL_ESING:
            std::cerr << "Error during integration of " << integration_info
                << ": Singularity or bad behaviour in integral.  Aborting."
                << std::endl;
            abort();
        case GSL_EDIVERGE:
            std::cerr << "Warning during integration of " << integration_info
                << ": Integrand is divergent or too slowly convergent.  Aborting."
                << std::endl;
            abort();
        case GSL_EDOM:
            std::cerr << "Warning during integration of " << integration_info
                << ": Error in the values of input arguments.  Aborting."
                << std::endl;
            abort();
        default:
                std::cerr << "Unknown status: " << status
                    << " returned from integration of " << integration_info
                    << std::endl;
            abort();
    }
}

Interpolations::Interpolations(
        const std::string& filename,
        Constants& constants,
        gsl_integration_workspace* workspace,
        int sub_regions,
        double eps_rel,
        double eps_abs
        )
{
    read_file_and_interpolate(filename, 2, 0, &a_of_tau_acc,
            &a_of_tau_spline, NULL, &redshift_to_scale_factor);
    read_file_and_interpolate(filename, 0, 2,
            &tau_of_redshift_acc, &tau_of_redshift_spline, NULL,
            NULL);

    constants.tau_ini = 1.0;
    double tau_today = gsl_spline_eval(tau_of_redshift_spline, 0,
            tau_of_redshift_acc);

    z_function(constants, tau_today, workspace, sub_regions, eps_rel, eps_abs);
}

struct z_integration_parameters {
    double q = 0;
    Constants& constants;
    Interpolations& interpols;
};



inline double q_over_eps(double tau, void* parameters) {
    z_integration_parameters* params = (z_integration_parameters*)(parameters);
    double m_nu = params->constants.m_nu;
    double T_ncdm0 = params->constants.T_ncdm0;
    double q = params->q;
    double a = gsl_spline_eval(params->interpols.a_of_tau_spline, tau,
            params->interpols.a_of_tau_acc);

    return pow(1 + a * a * m_nu * m_nu / pow(q * T_ncdm0, 2.0), -0.5);
}

void Interpolations::z_function(
        Constants& constants,
        double tau_fin,
        gsl_integration_workspace* workspace,
        int sub_regions,
        double eps_rel,
        double eps_abs
        )
{
    z_integration_parameters params = {0, constants, *this};

    int n_points = 100;

    std::vector<double> q_vals(n_points, 0);
    std::vector<double> tau_vals(n_points, 0);
    std::vector<double> z_vals(n_points*n_points, 0);

    double q_min = 0.001;
    double q_max = 100;
    // First element of q_vals should remain 0
    for (int i = 1; i < n_points; ++i) {
        q_vals[i] = q_min * pow(q_max/q_min, static_cast<double>(i-1)/(n_points - 2));
    }

    double tau_ini = constants.tau_ini;
    for (int i = 0; i < n_points; ++i) {
        tau_vals[i] = tau_ini * pow(tau_fin/tau_ini,
                static_cast<double>(i)/(n_points - 1));
    }

    gsl_function F;
    F.function = &q_over_eps;
    F.params = &params;

    Quantity result;

    for (int j = 0; j < n_points; ++j) {
        params.q = q_vals[j];
        for (int i = 0; i < n_points; ++i) {
            // Integrate from 1e-6 rather than 0 to avoid extrapolation error
            int status = gsl_integration_qag(&F, 1e-6, tau_vals[i], eps_abs,
                    eps_rel, sub_regions, GSL_INTEG_GAUSS61, workspace,
                    &result.value, &result.error);

            if (status != 0) {
                integration_status(status, "z_function integration", result);
            }
            z_vals[j * n_points + i] = result.value;
        }
    }

    z_spline_x_acc = gsl_interp_accel_alloc();
    z_spline_y_acc = gsl_interp_accel_alloc();
    z_spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, n_points, n_points);
    gsl_spline2d_init(z_spline, tau_vals.data(), q_vals.data(),
            z_vals.data(), n_points, n_points);
}



double eps_over_q(
        double tau,
        double q,
        const Constants& constants,
        const Interpolations& interpols
        )
{
    double m_nu = constants.m_nu;
    double T_ncdm0 = constants.T_ncdm0;
    double a = gsl_spline_eval(interpols.a_of_tau_spline, tau,
            interpols.a_of_tau_acc);

    return pow(1 + a * a * m_nu * m_nu / pow(q * T_ncdm0, 2.0), 0.5);
}



inline double f0(double q) {
    return pow(1 + exp(q), -1);
}



inline double df0_dlnq(double q) {
    return -q * exp(q) * pow(1 + exp(q), -2);
}



struct psi_integration_parameters {
    double q = 0;
    double k = 0;
    double z_func_at_tau = 0;
    const Constants& constants;
    const Interpolations& interpols;
};



double psi_0_integrand(double tau, void* parameters) {
    psi_integration_parameters* params = (psi_integration_parameters*)(parameters);
    double q = params->q;
    double k = params->k;

    double z_func_at_tau = params->z_func_at_tau;
    double z_func_at_tau_prime = gsl_spline2d_eval(params->interpols.z_spline,
            tau, q, params->interpols.z_spline_x_acc,
            params->interpols.z_spline_y_acc);

    return eps_over_q(tau, q, params->constants, params->interpols)
        * gsl_sf_bessel_j1(k * (z_func_at_tau - z_func_at_tau_prime) );
}



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
        )
{
    double z_func_at_tau = gsl_spline2d_eval(interpols.z_spline, tau, q,
            interpols.z_spline_x_acc, interpols.z_spline_y_acc);

    psi_integration_parameters params = {q, k, z_func_at_tau, constants, interpols};

    gsl_function F;
    F.function = &psi_0_integrand;
    F.params = &params;

    int status = gsl_integration_qag(&F, constants.tau_ini, tau, eps_abs,
            eps_rel, sub_regions, GSL_INTEG_GAUSS61, workspace, &result.value,
            &result.error);
    if (status != 0) {
        integration_status(status, "psi_0 integration", result);
    }
}



double psi_1_integrand(double tau, void* parameters) {
    psi_integration_parameters* params = (psi_integration_parameters*)(parameters);
    double q = params->q;
    double k = params->k;

    double z_func_at_tau = params->z_func_at_tau;
    double z_func_at_tau_prime = gsl_spline2d_eval(params->interpols.z_spline,
            tau, q, params->interpols.z_spline_x_acc,
            params->interpols.z_spline_y_acc);

    return eps_over_q(tau, q, params->constants, params->interpols)
        * (
          1.0/3.0 * gsl_sf_bessel_jl(0, k * (z_func_at_tau - z_func_at_tau_prime) )
        - 2.0/3.0 * gsl_sf_bessel_jl(2, k * (z_func_at_tau - z_func_at_tau_prime) )
        );
}



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
        )
{
    double z_func_at_tau = gsl_spline2d_eval(interpols.z_spline, tau, q,
            interpols.z_spline_x_acc, interpols.z_spline_y_acc);

    psi_integration_parameters params = {q, k, z_func_at_tau, constants, interpols};

    gsl_function F;
    F.function = &psi_1_integrand;
    F.params = &params;

    int status = gsl_integration_qag(&F, constants.tau_ini, tau, eps_abs,
            eps_rel, sub_regions, GSL_INTEG_GAUSS61, workspace, &result.value,
            &result.error);
    if (status != 0) {
        integration_status(status, "psi_1 integration", result);
    }

    result *= -1;
}



double psi_2_integrand(double tau, void* parameters) {
    psi_integration_parameters* params = (psi_integration_parameters*)(parameters);
    double q = params->q;
    double k = params->k;

    double z_func_at_tau = params->z_func_at_tau;
    double z_func_at_tau_prime = gsl_spline2d_eval(params->interpols.z_spline,
            tau, q, params->interpols.z_spline_x_acc,
            params->interpols.z_spline_y_acc);

    return eps_over_q(tau, q, params->constants, params->interpols)
        * (
          2.0/5.0 * gsl_sf_bessel_jl(1, k * (z_func_at_tau - z_func_at_tau_prime) )
        - 3.0/5.0 * gsl_sf_bessel_jl(3, k * (z_func_at_tau - z_func_at_tau_prime) )
        );
}



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
        )
{
    double z_func_at_tau = gsl_spline2d_eval(interpols.z_spline, tau, q,
            interpols.z_spline_x_acc, interpols.z_spline_y_acc);

    psi_integration_parameters params = {q, k, z_func_at_tau, constants, interpols};

    gsl_function F;
    F.function = &psi_2_integrand;
    F.params = &params;

    int status = gsl_integration_qag(&F, constants.tau_ini, tau, eps_abs,
            eps_rel, sub_regions, GSL_INTEG_GAUSS61, workspace, &result.value,
            &result.error);

    if (status != 0) {
        integration_status(status, "psi_2 integration", result);
    }

    result *= -1;
}



struct background_integration_parameters {
    double tau;
    const Constants& constants;
    const Interpolations& interpols;
};



inline double rho_integrand(double q, void* parameters) {
    background_integration_parameters* params =
        (background_integration_parameters*)(parameters);
    return pow(q,3) * eps_over_q(params->tau, q, params->constants,
            params->interpols)
        * f0(q);
}



inline double pressure_integrand(double q, void* parameters) {
    background_integration_parameters* params =
        (background_integration_parameters*)(parameters);
    return pow(q,3) * 1.0/eps_over_q(params->tau, q, params->constants,
            params->interpols)
        * f0(q);
}



Quantity Perturbations::integrate_background(double (*integrand)(double, void*)) {
    Quantity result;

    background_integration_parameters params = {tau, constants, interpols};

    gsl_function F;
    F.function = integrand;
    F.params = &params;

    int status = gsl_integration_qag(&F, 0, cutoff, outer_eps_abs,
            outer_eps_rel, outer_sub_regions, GSL_INTEG_GAUSS61,
            outer_workspace, &result.value, &result.error);

    if (status != 0) {
        integration_status(status, "background integration", result);
    }

    return result;
}



struct perturbation_integrand_parameters {
    double tau;
    double k;
    double q_min;
    const Constants& constants;
    const Interpolations& interpols;
    gsl_integration_workspace* inner_workspace;
    int inner_sub_regions;
    double inner_eps_rel;
    double inner_eps_abs;
    gsl_interp_accel* psi_acc = nullptr;
    gsl_spline* psi_spline = nullptr;
};



inline double delta_rho_integrand(double q, void* parameters) {
    perturbation_integrand_parameters* params =
        (perturbation_integrand_parameters*)(parameters);

    Quantity psi_0_result;
    psi_0(params->tau, params->k, q, params->constants, params->interpols,
            psi_0_result, params->inner_workspace, params->inner_sub_regions,
            params->inner_eps_rel, params->inner_eps_abs);

    return pow(q,3) * eps_over_q(params->tau, q, params->constants,
            params->interpols) * df0_dlnq(q) * psi_0_result.value;
}



inline double delta_P_integrand(double q, void* parameters) {
    perturbation_integrand_parameters* params =
        (perturbation_integrand_parameters*)(parameters);

    Quantity psi_0_result;
    psi_0(params->tau, params->k, q, params->constants, params->interpols,
            psi_0_result, params->inner_workspace, params->inner_sub_regions,
            params->inner_eps_rel, params->inner_eps_abs);

    return pow(q,3) * 1.0/eps_over_q(params->tau, q, params->constants,
            params->interpols) * df0_dlnq(q) * psi_0_result.value;
}



inline double theta_integrand(double q, void* parameters) {
    perturbation_integrand_parameters* params =
        (perturbation_integrand_parameters*)(parameters);

    Quantity psi_1_result;
    psi_1(params->tau, params->k, q, params->constants, params->interpols,
            psi_1_result, params->inner_workspace, params->inner_sub_regions,
            params->inner_eps_rel, params->inner_eps_abs);

    return pow(q,3) * df0_dlnq(q) * psi_1_result.value;
}



inline double sigma_integrand(double q, void* parameters) {
    perturbation_integrand_parameters* params =
        (perturbation_integrand_parameters*)(parameters);

    Quantity psi_2_result;
    psi_2(params->tau, params->k, q, params->constants, params->interpols,
            psi_2_result, params->inner_workspace, params->inner_sub_regions,
            params->inner_eps_rel, params->inner_eps_abs);

    return pow(q,3) * 1.0/eps_over_q(params->tau, q, params->constants,
            params->interpols) * df0_dlnq(q) * psi_2_result.value;
}



inline double sigma_integrand_psi_2_interpolated(double q, void* parameters) {
    perturbation_integrand_parameters* params =
        (perturbation_integrand_parameters*)(parameters);

    return pow(q,3) * 1.0/eps_over_q(params->tau, q, params->constants,
            params->interpols) * df0_dlnq(q) *
        gsl_spline_eval(params->psi_spline, q, params->psi_acc);
}



Quantity Perturbations::integrate_perturbations(
        double (*integrand)(double, void*),
        double q_min,
        gsl_interp_accel* psi_acc,
        gsl_spline* psi_spline
        )
{
    Quantity result;

    perturbation_integrand_parameters params = {tau, k, q_min, constants, interpols,
        inner_workspace, inner_sub_regions, inner_eps_rel, inner_eps_abs,
        psi_acc, psi_spline};

    gsl_function F;
    F.function = integrand;
    F.params = &params;

    int status = gsl_integration_qag(&F, 0, cutoff, outer_eps_abs,
            outer_eps_rel, outer_sub_regions, GSL_INTEG_GAUSS61,
            outer_workspace, &result.value, &result.error);

    if (status != 0) {
        integration_status(status, "perturbations integration", result);
    }

    return result;
}



void Perturbations::interpolate_psi(
        void (*psi)(
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
            ),
        double q_min,
        gsl_interp_accel** psi_acc,
        gsl_spline** psi_spline
        )
{
    int n_points = 1000;

    std::vector<double> q_vals(n_points, 0);
    std::vector<double> psi_vals(n_points, 0);

    double q_max = cutoff;
    // First element of q_vals should remain 0
    for (int i = 0; i < n_points; ++i) {
        q_vals[i] = q_min * pow(q_max/q_min, static_cast<double>(i)/(n_points - 1));
    }

    Quantity result;

    for (int i = 0; i < n_points; ++i) {
        psi(tau, k, q_vals[i], constants, interpols, result, inner_workspace,
                inner_sub_regions, inner_eps_rel, inner_eps_abs);

        psi_vals[i] = result.value;
    }

    *psi_acc = gsl_interp_accel_alloc();
    *psi_spline = gsl_spline_alloc(gsl_interp_cspline, n_points);
    gsl_spline_init(*psi_spline, q_vals.data(), psi_vals.data(), n_points);
}



Perturbations::Perturbations(
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
        ) :
    tau(tau), k(k), constants(constants), interpols(interpols), cutoff(cutoff),
    outer_workspace(outer_workspace), inner_workspace(inner_workspace),
    outer_sub_regions(outer_sub_regions), inner_sub_regions(inner_sub_regions),
    outer_eps_rel(outer_eps_rel), outer_eps_abs(outer_eps_abs),
    inner_eps_rel(inner_eps_rel), inner_eps_abs(inner_eps_abs)
{}


void Perturbations::compute() {
    double a = gsl_spline_eval(interpols.a_of_tau_spline, tau,
            interpols.a_of_tau_acc);

    double factor = 4 * M_PI * pow(constants.T_ncdm0/a, 4);

    rho      = integrate_background(rho_integrand);
    pressure = integrate_background(pressure_integrand);

    Quantity delta_rho = integrate_perturbations(delta_rho_integrand);
    Quantity delta_P   = integrate_perturbations(delta_P_integrand);
    /* theta              = integrate_perturbations(theta_integrand); */

    // If k is larger than threshold, interpolate psi_2(q) before integrating
    if (k > 5) {
        double q_min = 1e-5;
        gsl_interp_accel* psi_2_acc = nullptr;
        gsl_spline* psi_2_spline = nullptr;
        interpolate_psi(psi_2, q_min, &psi_2_acc, &psi_2_spline);

        sigma = integrate_perturbations(sigma_integrand_psi_2_interpolated,
                q_min, psi_2_acc, psi_2_spline);
    }
    else {
        sigma = integrate_perturbations(sigma_integrand);
    }

    delta = delta_rho / rho;
    cs2 = delta_P / delta_rho;

    rho *= factor;
    pressure *= factor/3.0;

    theta *= k * factor;
    sigma *= 2.0/3.0 * factor;

    theta /= (rho + pressure);
    sigma /= (rho + pressure);
}
