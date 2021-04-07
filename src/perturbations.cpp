/*
   perturbations.cpp

   Created by Petter Taule on 30.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include "../include/perturbations.hpp"
#include "../include/measurement.hpp"
#include "../include/interpolation.hpp"
#include "../include/io.hpp"

#define SQUARE(x) x*x
#define CUBE(x) x*x

using size_t = std::size_t;

template <class T>
using Vec1D = std::vector<T>;
template <class T>
using Vec2D = std::vector<std::vector<T>>;

IntegrationException::IntegrationException(
                int gsl_status,
                const std::string& integration_name
        )
{
    switch (gsl_status) {
        case GSL_EMAXITER:
            ex = "Maxiumum number of subdevisions reached while integrating " +
                integration_name;
            critical_ = true;
            break;
        case GSL_EROUND:
            ex = "Warning: rounding error during integration of " +
                integration_name;
            critical_ = false;
            break;
        case GSL_ESING:
            ex = "Singularity or bad behaviour in integration of " +
                integration_name;
            critical_ = true;
            break;
        case GSL_EDIVERGE:
            ex = "Divergent or too slowly converging integral: " +
                integration_name;
            critical_ = true;
            break;
        case GSL_EDOM:
            ex = "Error in the values of input arguments in integration of " +
                integration_name;
            critical_ = true;
            break;
        default:
            ex = "Unknown status returned from integration of " + integration_name;
            critical_ = true;
    }
}



struct y_integration_parameters {
    double m_nu_over_q_T_ncdm0_square = 0;
    Interpolation1D& scale_factor_of_conf_time;
};



inline double q_over_eps(double tau, void* parameters) {
    y_integration_parameters& params =
        *static_cast<y_integration_parameters*>(parameters);
    double a = params.scale_factor_of_conf_time(tau);

    return pow(1 + SQUARE(a) * params.m_nu_over_q_T_ncdm0_square, -0.5);
}



void Background::compute_y_function()
{
    size_t sub_regions = 10000;
    double rel_tol  = 1e-8;
    double abs_tol  = 0;

    gsl_integration_workspace* workspace =
        gsl_integration_workspace_alloc(sub_regions);

    y_integration_parameters params = {0, scale_factor_of_conf_time};

    size_t n_points = 1000;
    double q_min = 0.001;
    double q_max = 100;

    Vec1D<double> q_grid(n_points, 0);
    Vec1D<double> tau_grid(n_points, 0);
    Vec1D<double> y(n_points*n_points, 0);

    /* First element of q_vals should remain 0 */
    for (size_t i = 1; i < n_points; ++i) {
        q_grid[i] = q_min * pow(q_max / q_min,
                    static_cast<double>(i - 1) / static_cast<double>(n_points - 2));
    }

    double tau_today = conf_time_of_redshift(0);
    for (size_t i = 0; i < n_points; ++i) {
        tau_grid[i] = tau_ini * pow(tau_today/tau_ini,
                static_cast<double>(i)/static_cast<double>(n_points - 1));
    }

    gsl_function F;
    F.function = &q_over_eps;
    F.params = &params;

    Measurement result;

    for (size_t j = 0; j < n_points; ++j) {
        params.m_nu_over_q_T_ncdm0_square =
            m_nu_over_T_ncdm0_square / SQUARE(q_grid[j]);
        for (size_t i = 0; i < n_points; ++i) {
            /* Integrate from 1e-6 rather than 0 to avoid extrapolation error */
            int status = gsl_integration_qag(&F, 1e-6, tau_grid[i], abs_tol,
                    rel_tol, sub_regions, GSL_INTEG_GAUSS61, workspace,
                    &result.value(), &result.error());

            if (status != 0) {
                IntegrationException ex(status, "y_function");
                if (ex.critical()) {
                    throw ex;
                }
                else {
                    std::cout << ex.what() << std::endl;
                }
            }
            y[j * n_points + i] = result.value();
        }
    }
    y_function = Interpolation2D(tau_grid, q_grid, y);
}



Background::Background(
                double hubble,
                double m_nu,
                double T_ncdm0,
                double tau_ini,
                const std::string& background_file
                ) :
    hubble(hubble), m_nu(m_nu), T_ncdm0(T_ncdm0), tau_ini(tau_ini)
{
    Vec2D<double> data;
    /* Read background (three first columns: redshift, prop. time, conf. time */
    read_columns_from_file(background_file, 3, data);

    /* Copy redshift values and convert to scale_factor */
    Vec1D<double> scale_factor = data[0];
    for (auto& el : scale_factor) {
        el = 1/(el + 1);
    }

    scale_factor_of_conf_time = Interpolation1D(scale_factor, data[2]);

    /* To interpolate tau(z), reverse vectors, since x (redshift) needs to be
     * increasing for interpolation */
    std::reverse(data[0].begin(), data[0].end());
    std::reverse(data[2].begin(), data[2].end());

    conf_time_of_redshift = Interpolation1D(data[0], data[2]);

    m_nu_over_T_ncdm0_square = SQUARE(m_nu) / SQUARE(T_ncdm0);
    compute_y_function();
}



void interpolate_metric_psi(
        double hubble,
        const Interpolation1D& conf_time_of_redshift,
        const std::string& k_grid_file,
        const std::string& z_grid_file,
        const std::string& metric_psi_file,
        Interpolation2D& metric_psi
        )
{
    Vec1D<double> x,y,z;
    Vec2D<double> data;

    /* Interpolate psi of wavenumber and redshift/conformal time */
    data.clear();
    read_columns_from_file(k_grid_file, 1, data);
    x = data[0];
    data.clear();
    read_columns_from_file(z_grid_file, 1, data);
    y = data[0];
    read_data_grid_from_file(metric_psi_file, z, x.size(), y.size());

    /* Convert wavenumber in h/Mpc to 1/Mpc */
    for (auto& el : x) {
        el *= hubble;
    }
    /* Convert redshift (y-grid) to conformal time */
    for (auto& el : y) {
        el = conf_time_of_redshift(el);
    }

    metric_psi = Interpolation2D(x, y, z);
}



double eps_over_q(double tau, double q, const Background& background)
{
    return pow(1 + SQUARE(background.scale_factor_of_conf_time(tau)) *
                       background.m_nu_over_T_ncdm0_square / SQUARE(q),
               0.5);
}



double f0(double q)
{
    return pow(1 + exp(q), -1);
}



double df0_dlnq(double q)
{
    return -q * exp(q) * pow(1 + exp(q), -2);
}



struct psi_integration_parameters {
    double q = 0;
    double k = 0;
    double y_func_at_tau = 0;
    double tau_lambda;
    double metric_psi_at_k_and_tau_lambda;
    const Background& bg;
    const Interpolation2D& metric_psi;
};



double psi_0_integrand(double tau, void* parameters) {
    psi_integration_parameters& params =
        *static_cast<psi_integration_parameters*>(parameters);

    double y_func_at_tau = params.y_func_at_tau;
    double y_func_at_tau_prime = params.bg.y_function(tau, params.q);
    double eps_over_q_val = eps_over_q(tau, params.q, params.bg);

    double result = 0;

    /* Before tau_lambda, assume psi = const */
    if (tau < params.tau_lambda) {
        result = eps_over_q_val
            * gsl_sf_bessel_j1(params.k * (y_func_at_tau - y_func_at_tau_prime));
    }
    /* After, use interpolated psi */
    else {
        result =
            params.metric_psi(params.k, tau) /
            params.metric_psi_at_k_and_tau_lambda *
            (eps_over_q_val + 1 / eps_over_q_val) *
            gsl_sf_bessel_j1(params.k * (y_func_at_tau - y_func_at_tau_prime));
    }

    return result;
}



void psi_0(const Psi_parameters& psi_params, Measurement& result)
{
    double tau                            = psi_params.tau;
    double k                              = psi_params.k;
    double q                              = psi_params.q;
    double tau_lambda                     = psi_params.tau_lambda;
    double metric_psi_at_k_and_tau_lambda = psi_params.metric_psi_at_k_and_tau_lambda;
    const Background& bg                  = psi_params.bg;
    const Interpolation2D& metric_psi     = psi_params.metric_psi;
    gsl_integration_workspace* workspace  = psi_params.workspace;
    size_t sub_regions                    = psi_params.sub_regions;
    double rel_tol                        = psi_params.rel_tol;
    double abs_tol                        = psi_params.abs_tol;

    double y_func_at_tau        = bg.y_function(tau, q);
    double y_func_at_tau_lambda = bg.y_function(tau_lambda, q);

    psi_integration_parameters params = {q, k, y_func_at_tau, tau_lambda,
        metric_psi_at_k_and_tau_lambda, bg, metric_psi};

    gsl_function F;
    F.function = &psi_0_integrand;
    F.params = &params;

    int status = gsl_integration_qag(&F, bg.tau_ini, tau, abs_tol,
            rel_tol, sub_regions, GSL_INTEG_GAUSS61, workspace, &result.value(),
            &result.error());

    if (status != 0) {
        IntegrationException ex(status, "psi_0");
        if (ex.critical()) {
            throw ex;
        }
        else {
            std::cout << ex.what() << std::endl;
        }
    }

    result *= k * metric_psi_at_k_and_tau_lambda;
    result += metric_psi_at_k_and_tau_lambda *
        gsl_sf_bessel_jl(0, k * (y_func_at_tau - y_func_at_tau_lambda))
        - metric_psi(k, tau);
}



double psi_1_integrand(double tau, void* parameters) {
    psi_integration_parameters& params =
        *static_cast<psi_integration_parameters*>(parameters);

    double y_func_at_tau = params.y_func_at_tau;
    double y_func_at_tau_prime = params.bg.y_function(tau, params.q);
    double eps_over_q_val = eps_over_q(tau, params.q, params.bg);

    double result = 0;

    // Before tau_lambda, assume psi = const
    if (tau < params.tau_lambda) {
        result = eps_over_q_val
            * (
              1.0/3.0 * gsl_sf_bessel_jl(0, params.k * (y_func_at_tau - y_func_at_tau_prime))
            - 2.0/3.0 * gsl_sf_bessel_jl(2, params.k * (y_func_at_tau - y_func_at_tau_prime))
            );
    }
    // After, use interpolated psi
    else {
        result = params.metric_psi(params.k, tau)
            / params.metric_psi_at_k_and_tau_lambda
            * (eps_over_q_val + 1/eps_over_q_val)
            * (
              1.0/3.0 * gsl_sf_bessel_jl(0, params.k * (y_func_at_tau - y_func_at_tau_prime))
            - 2.0/3.0 * gsl_sf_bessel_jl(2, params.k * (y_func_at_tau - y_func_at_tau_prime))
        );
    }

    return result;
}



void psi_1(
        const Psi_parameters& psi_params,
        Measurement& result
        )
{
    double tau                            = psi_params.tau;
    double k                              = psi_params.k;
    double q                              = psi_params.q;
    double tau_lambda                     = psi_params.tau_lambda;
    double metric_psi_at_k_and_tau_lambda = psi_params.metric_psi_at_k_and_tau_lambda;
    const Background& bg                  = psi_params.bg;
    const Interpolation2D& metric_psi     = psi_params.metric_psi;
    gsl_integration_workspace* workspace  = psi_params.workspace;
    size_t sub_regions                    = psi_params.sub_regions;
    double rel_tol                        = psi_params.rel_tol;
    double abs_tol                        = psi_params.abs_tol;

    double y_func_at_tau = bg.y_function(tau, q);
    double y_func_at_tau_lambda = bg.y_function(tau_lambda, q);

    psi_integration_parameters params = {q, k, y_func_at_tau, tau_lambda,
        metric_psi_at_k_and_tau_lambda, bg, metric_psi};

    gsl_function F;
    F.function = &psi_1_integrand;
    F.params = &params;

    int status = gsl_integration_qag(&F, bg.tau_ini, tau, abs_tol,
            rel_tol, sub_regions, GSL_INTEG_GAUSS61, workspace, &result.value(),
            &result.error());

    if (status != 0) {
        IntegrationException ex(status, "psi_1");
        if (ex.critical()) {
            throw ex;
        }
        else {
            std::cout << ex.what() << std::endl;
        }
    }

    result *= - k * metric_psi_at_k_and_tau_lambda;
    result += metric_psi(k, tau_lambda)
        * gsl_sf_bessel_jl(1, k * (y_func_at_tau - y_func_at_tau_lambda));
}



double psi_2_integrand(double tau, void* parameters) {
    psi_integration_parameters& params =
        *static_cast<psi_integration_parameters*>(parameters);

    double y_func_at_tau = params.y_func_at_tau;
    double y_func_at_tau_prime = params.bg.y_function(tau, params.q);
    double eps_over_q_val = eps_over_q(tau, params.q, params.bg);

    double result = 0;

    // Before tau_lambda, assume psi = const
    if (tau < params.tau_lambda) {
        result = eps_over_q_val
            * (
                    2.0/5.0 * gsl_sf_bessel_jl(1, params.k * (y_func_at_tau - y_func_at_tau_prime))
                  - 3.0/5.0 * gsl_sf_bessel_jl(3, params.k * (y_func_at_tau - y_func_at_tau_prime))
              );
    }
    // After, use interpolated psi
    else {
        result = params.metric_psi(params.k, tau)
            / params.metric_psi_at_k_and_tau_lambda
            * (eps_over_q_val + 1/eps_over_q_val)
            * (
                    2.0/5.0 * gsl_sf_bessel_jl(1, params.k * (y_func_at_tau - y_func_at_tau_prime))
                  - 3.0/5.0 * gsl_sf_bessel_jl(3, params.k * (y_func_at_tau - y_func_at_tau_prime))
              );
    }

    return result;
}



void psi_2(
        const Psi_parameters& psi_params,
        Measurement& result
        )
{
    double tau                            = psi_params.tau;
    double k                              = psi_params.k;
    double q                              = psi_params.q;
    double tau_lambda                     = psi_params.tau_lambda;
    double metric_psi_at_k_and_tau_lambda = psi_params.metric_psi_at_k_and_tau_lambda;
    const Background& bg                  = psi_params.bg;
    const Interpolation2D& metric_psi     = psi_params.metric_psi;
    gsl_integration_workspace* workspace  = psi_params.workspace;
    size_t sub_regions                    = psi_params.sub_regions;
    double rel_tol                        = psi_params.rel_tol;
    double abs_tol                        = psi_params.abs_tol;

    double y_func_at_tau = bg.y_function(tau, q);
    double y_func_at_tau_lambda = bg.y_function(tau_lambda, q);

    psi_integration_parameters params = {q, k, y_func_at_tau, tau_lambda,
        metric_psi_at_k_and_tau_lambda, bg, metric_psi};

    gsl_function F;
    F.function = &psi_2_integrand;
    F.params = &params;

    int status = gsl_integration_qag(&F, bg.tau_ini, tau, abs_tol,
            rel_tol, sub_regions, GSL_INTEG_GAUSS61, workspace, &result.value(),
            &result.error());

    if (status != 0) {
        IntegrationException ex(status, "psi_2");
        if (ex.critical()) {
            throw ex;
        }
        else {
            std::cout << ex.what() << std::endl;
        }
    }

    result *= - k * metric_psi_at_k_and_tau_lambda;
    result += metric_psi(k, tau_lambda)
        * gsl_sf_bessel_jl(2, k * (y_func_at_tau - y_func_at_tau_lambda));
}



struct fluid_background_integration_parameters {
    double tau;
    const Background& bg;
};



double rho_integrand(double q, void* parameters) {
    fluid_background_integration_parameters& params =
        *static_cast<fluid_background_integration_parameters*>(parameters);
    return CUBE(q) * eps_over_q(params.tau, q, params.bg) * f0(q);
}



double pressure_integrand(double q, void* parameters) {
    fluid_background_integration_parameters& params =
        *static_cast<fluid_background_integration_parameters*>(parameters);
    return CUBE(q) * 1.0/eps_over_q(params.tau, q, params.bg) * f0(q);
}



Measurement Perturbations::integrate_fluid_background(
        double (*integrand)(double, void*)
        )
{
    Measurement result;

    fluid_background_integration_parameters params = {tau, bg};

    gsl_function F;
    F.function = integrand;
    F.params = &params;

    int status = gsl_integration_qag(&F, 0, cutoff, outer_abs_tol,
            outer_rel_tol, outer_sub_regions, GSL_INTEG_GAUSS61,
            outer_workspace, &result.value(), &result.error());

    if (status != 0) {
        IntegrationException ex(status, "fluid_background");
        if (ex.critical()) {
            throw ex;
        }
        else {
            std::cout << ex.what() << std::endl;
        }
    }

    return result;
}



struct fluid_perturbation_integrand_parameters {
    double tau;
    double k;
    double q_min;
    double tau_lambda;
    double metric_psi_at_k_and_tau_lambda;
    const Background& bg;
    const Interpolation2D& metric_psi;
    gsl_integration_workspace* inner_workspace;
    size_t inner_sub_regions;
    double inner_rel_tol;
    double inner_abs_tol;
    const Interpolation1D* psi = nullptr;
};



double delta_rho_integrand(double q, void* parameters) {
    fluid_perturbation_integrand_parameters& params =
        *static_cast<fluid_perturbation_integrand_parameters*>(parameters);

    Measurement psi_0_result;
    const Psi_parameters psi_params = {params.tau, params.k, q,
        params.tau_lambda, params.metric_psi_at_k_and_tau_lambda,
        params.bg, params.metric_psi, params.inner_workspace,
        params.inner_sub_regions, params.inner_rel_tol, params.inner_abs_tol
    };
    psi_0(psi_params, psi_0_result);

    return CUBE(q) * eps_over_q(params.tau, q, params.bg) * df0_dlnq(q) *
        psi_0_result.value();
}



double delta_P_integrand(double q, void* parameters) {
    fluid_perturbation_integrand_parameters& params =
        *static_cast<fluid_perturbation_integrand_parameters*>(parameters);

    Measurement psi_0_result;
    const Psi_parameters psi_params = {params.tau, params.k, q,
        params.tau_lambda, params.metric_psi_at_k_and_tau_lambda,
        params.bg, params.metric_psi, params.inner_workspace,
        params.inner_sub_regions, params.inner_rel_tol, params.inner_abs_tol
    };
    psi_0(psi_params, psi_0_result);

    return CUBE(q) * 1.0/eps_over_q(params.tau, q, params.bg) * df0_dlnq(q) *
        psi_0_result.value();
}



double theta_integrand(double q, void* parameters) {
    fluid_perturbation_integrand_parameters& params =
        *static_cast<fluid_perturbation_integrand_parameters*>(parameters);

    Measurement psi_1_result;
    const Psi_parameters psi_params = {params.tau, params.k, q,
        params.tau_lambda, params.metric_psi_at_k_and_tau_lambda,
        params.bg, params.metric_psi, params.inner_workspace,
        params.inner_sub_regions, params.inner_rel_tol, params.inner_abs_tol
    };
    psi_1(psi_params, psi_1_result);

    return CUBE(q) * df0_dlnq(q) * psi_1_result.value();
}



double sigma_integrand(double q, void* parameters) {
    fluid_perturbation_integrand_parameters& params =
        *static_cast<fluid_perturbation_integrand_parameters*>(parameters);

    Measurement psi_2_result;
    const Psi_parameters psi_params = {params.tau, params.k, q,
        params.tau_lambda, params.metric_psi_at_k_and_tau_lambda,
        params.bg, params.metric_psi, params.inner_workspace,
        params.inner_sub_regions, params.inner_rel_tol, params.inner_abs_tol
    };
    psi_2(psi_params, psi_2_result);

    return CUBE(q) * 1.0/eps_over_q(params.tau, q, params.bg) * df0_dlnq(q) *
        psi_2_result.value();
}



double sigma_integrand_psi_2_interpolated(double q, void* parameters) {
    fluid_perturbation_integrand_parameters& params =
        *static_cast<fluid_perturbation_integrand_parameters*>(parameters);

    return pow(q,3) * 1.0/eps_over_q(params.tau, q, params.bg)
            * df0_dlnq(q) * (*params.psi)(q);
}



Measurement Perturbations::integrate_fluid_perturbation(
        double (*integrand)(double, void*),
        double q_min,
        const Interpolation1D* psi
        )
{
    Measurement result;

    fluid_perturbation_integrand_parameters params = {tau, k, q_min, tau_lambda,
        metric_psi_at_k_and_tau_lambda, bg, metric_psi, inner_workspace,
        inner_sub_regions, inner_rel_tol, inner_abs_tol, psi};

    gsl_function F;
    F.function = integrand;
    F.params = &params;

    int status = gsl_integration_qag(&F, 0, cutoff, outer_abs_tol,
            outer_rel_tol, outer_sub_regions, GSL_INTEG_GAUSS61,
            outer_workspace, &result.value(), &result.error());

    if (status != 0) {
        IntegrationException ex(status, "fluid_perturbation");
        if (ex.critical()) {
            throw ex;
        }
        else {
            std::cout << ex.what() << std::endl;
        }
    }

    return result;
}



void Perturbations::interpolate_psi(
        void (*psi)(
            const Psi_parameters& psi_params,
            Measurement& result
            ),
        double q_min,
        Interpolation1D& psi_spline
        )
{
    size_t n_points = 1000;

    Vec1D<double> q_grid(n_points, 0);
    Vec1D<double> psi_vals(n_points, 0);

    double q_max = cutoff;
    /* First element of q_vals should remain 0 */
    for (size_t i = 0; i < n_points; ++i) {
        q_grid[i] = q_min * pow(q_max / q_min, static_cast<double>(i) /
                static_cast<double>(n_points - 1));
    }

    Psi_parameters psi_params = {tau, k, 0.0,
        tau_lambda, metric_psi_at_k_and_tau_lambda,
        bg, metric_psi, inner_workspace,
        inner_sub_regions, inner_rel_tol, inner_abs_tol
    };
    Measurement result;

    for (size_t i = 0; i < n_points; ++i) {
        psi_params.q = q_grid[i];
        psi(psi_params, result);

        psi_vals[i] = result.value();
    }

    psi_spline = Interpolation1D(q_grid, psi_vals);
}



Perturbations::Perturbations(
        double tau,
        double k,
        const Background& bg,
        const Interpolation2D& metric_psi,
        double z_lambda,
        double cutoff,
        bool do_interpolate_psi,
        gsl_integration_workspace* outer_workspace,
        gsl_integration_workspace* inner_workspace,
        size_t outer_sub_regions,
        size_t inner_sub_regions,
        double outer_rel_tol,
        double outer_abs_tol,
        double inner_rel_tol,
        double inner_abs_tol
        ) :
    tau(tau), k(k), bg(bg), metric_psi(metric_psi), cutoff(cutoff),
    do_interpolate_psi(do_interpolate_psi), outer_workspace(outer_workspace),
    inner_workspace(inner_workspace), outer_sub_regions(outer_sub_regions),
    inner_sub_regions(inner_sub_regions), outer_rel_tol(outer_rel_tol),
    outer_abs_tol(outer_abs_tol), inner_rel_tol(inner_rel_tol),
    inner_abs_tol(inner_abs_tol)
{
    tau_lambda = bg.conf_time_of_redshift(z_lambda);
    metric_psi_at_k_and_tau_lambda = metric_psi(k, tau_lambda);
}



void Perturbations::compute() {
    double a = bg.scale_factor_of_conf_time(tau);
    double factor = 4 * M_PI * pow(bg.T_ncdm0/a, 4);

    try {
        rho      = integrate_fluid_background(rho_integrand);
        pressure = integrate_fluid_background(pressure_integrand);

        std::cout << "Computing delta_rho...";
        std::cout.flush();
        Measurement delta_rho = integrate_fluid_perturbation(delta_rho_integrand);
        std::cout << "done.\nComputing delta_P...";
        std::cout.flush();
        Measurement delta_P   = integrate_fluid_perturbation(delta_P_integrand);
        std::cout << "done.\n";
        std::cout.flush();
        theta              = integrate_fluid_perturbation(theta_integrand);

        std::cout << "Computing sigma...";
        std::cout.flush();
        /* If k is larger than threshold, interpolate psi_2(q) before integrating */
        if (do_interpolate_psi && k > 5) {
            double q_min = 1e-5;
            Interpolation1D psi_2_spline;
            interpolate_psi(psi_2, q_min, psi_2_spline);

            sigma = integrate_fluid_perturbation(sigma_integrand_psi_2_interpolated,
                    q_min, &psi_2_spline);
        }
        else {
            sigma = integrate_fluid_perturbation(sigma_integrand);
        }
        std::cout << "done" << std::endl;

        delta = delta_rho / rho;
        cs2 = delta_P / delta_rho;
        cs2 *= 1.0/3.0;

        rho *= factor;
        pressure *= factor/3.0;

        theta *= k * factor;
        sigma *= 2.0/3.0 * factor;

        theta /= (rho + pressure);
        sigma /= (rho + pressure);
    }
    catch (const IntegrationException& ex) {
        throw ex;
    }
    catch (const std::exception& ex) {
        throw ex;
    }
}
