/*
   perturbations.cpp

   Created by Petter Taule on 30.05.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include <gsl/gsl_sf_bessel.h>

#include "../include/perturbations.hpp"
#include "../include/measurement.hpp"
#include "../include/interpolation.hpp"
#include "../include/integrator.hpp"
#include "../include/io.hpp"

using size_t = std::size_t;

template <class T>
using Vec1D = std::vector<T>;
template <class T>
using Vec2D = std::vector<std::vector<T>>;


struct y_integration_parameters {
    double m_nu_over_q_T_ncdm0_square = 0;
    Interpolation1D& scale_factor_of_conf_time;
};



double q_over_eps(double tau, void* parameters) {
    y_integration_parameters& p =
        *static_cast<y_integration_parameters*>(parameters);
    double a = p.scale_factor_of_conf_time(tau);

    return 1/std::sqrt(1 + a * a * p.m_nu_over_q_T_ncdm0_square);
}



void Background::compute_y_function()
{
    Integrator integrator(10000);

    y_integration_parameters params = {0, scale_factor_of_conf_time};

    size_t n_points = 1000;
    double q_min = 0.001;
    double q_max = 100;

    Vec1D<double> q_grid(n_points, 0);
    Vec1D<double> tau_grid(n_points, 0);
    Vec1D<double> y(n_points*n_points, 0);

    /* First element of q_vals should remain 0 */
    for (size_t i = 1; i < n_points; ++i) {
        q_grid[i] = q_min * std::pow(q_max / q_min,
                    static_cast<double>(i - 1) / static_cast<double>(n_points - 2));
    }

    double tau_today = conf_time_of_redshift(0);
    for (size_t i = 0; i < n_points; ++i) {
        tau_grid[i] = tau_ini * std::pow(tau_today/tau_ini,
                static_cast<double>(i)/static_cast<double>(n_points - 1));
    }

    Measurement result;

    /* for (size_t j = 0; j < n_points; ++j) { */
    for (size_t j = 1; j < n_points; ++j) {
        params.m_nu_over_q_T_ncdm0_square =
            m_nu_over_T_ncdm0_square / (q_grid[j] * q_grid[j]);
        for (size_t i = 0; i < n_points; ++i) {
            /* Integrate from 1e-6 rather than 0 to avoid extrapolation error */
            result = integrator.integrate(q_over_eps, &params, 1e-6, tau_grid[i]);
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

    scale_factor_of_conf_time = Interpolation1D(data[2], scale_factor);

    /* To interpolate tau(z), reverse vectors, since x (redshift) needs to be
     * increasing for interpolation */
    std::reverse(data[0].begin(), data[0].end());
    std::reverse(data[2].begin(), data[2].end());

    conf_time_of_redshift = Interpolation1D(data[0], data[2]);

    m_nu_over_T_ncdm0_square = m_nu * m_nu / (T_ncdm0 * T_ncdm0);
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
    double a = background.scale_factor_of_conf_time(tau);
    return std::sqrt(1 + a * a * background.m_nu_over_T_ncdm0_square / (q * q));
}



double f0(double q)
{
    return 1/(1 + exp(q));
}



double df0_dlnq(double q)
{
    return -q * std::exp(q) * 1/std::pow(1 + std::exp(q), 2);
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
        )
{
    double y_func_at_tau        = bg.y_function(tau, q);
    double y_func_at_tau_lambda = bg.y_function(tau_lambda, q);

    psi_integration_parameters params = {q, k, y_func_at_tau, tau_lambda,
        metric_psi_at_k_and_tau_lambda, bg, metric_psi};

    Measurement result(k * metric_psi_at_k_and_tau_lambda, 0);

    switch (l) {
        case 0:
            result *= inner_integrator.integrate(psi_0_integrand, &params,
                    bg.tau_ini, tau);
            result -= metric_psi(k, tau);
            break;
        case 1:
            result *= - inner_integrator.integrate(psi_1_integrand, &params,
                    bg.tau_ini, tau);
            break;
        case 2:
            result *= - inner_integrator.integrate(psi_2_integrand, &params,
                    bg.tau_ini, tau);
            break;
        default:
            throw std::invalid_argument(
                "Invalid argument l > 2 given to compute_psi().");
    }

    /* Add boundary term from integration split when z < z_lambda */
    if (tau > tau_lambda) {
        result += metric_psi_at_k_and_tau_lambda *
            gsl_sf_bessel_jl(l, k * (y_func_at_tau - y_func_at_tau_lambda));
    }
    return result;
}



struct fluid_background_integration_parameters {
    double tau;
    const Background& bg;
};



double rho_integrand(double q, void* parameters) {
    fluid_background_integration_parameters& p =
        *static_cast<fluid_background_integration_parameters*>(parameters);
    return pow(q,3) * eps_over_q(p.tau, q, p.bg) * f0(q);
}



double pressure_integrand(double q, void* parameters) {
    fluid_background_integration_parameters& p =
        *static_cast<fluid_background_integration_parameters*>(parameters);
    return pow(q,3) * 1.0/eps_over_q(p.tau, q, p.bg) * f0(q);
}



Measurement Perturbations::integrate_fluid_background(
        double (*integrand)(double, void*)
        )
{
    fluid_background_integration_parameters p = {tau, bg};
    Measurement result = outer_integrator.integrate(integrand, &p, 0, cutoff);
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
    const Integrator& inner_integrator;
    const Interpolation1D* psi = nullptr;
};



double delta_rho_integrand(double q, void* parameters) {
    fluid_perturbation_integrand_parameters& p =
        *static_cast<fluid_perturbation_integrand_parameters*>(parameters);

    Measurement psi_0 = compute_psi_l(1, p.tau, p.k, q, p.tau_lambda,
                                      p.metric_psi_at_k_and_tau_lambda, p.bg,
                                      p.metric_psi, p.inner_integrator);

    return pow(q,3) * eps_over_q(p.tau, q, p.bg) * df0_dlnq(q) * psi_0.value();
}



double delta_P_integrand(double q, void* parameters) {
    fluid_perturbation_integrand_parameters& p =
        *static_cast<fluid_perturbation_integrand_parameters*>(parameters);

    Measurement psi_0 = compute_psi_l(0, p.tau, p.k, q, p.tau_lambda,
                                      p.metric_psi_at_k_and_tau_lambda, p.bg,
                                      p.metric_psi, p.inner_integrator);

    return pow(q,3) * 1.0/eps_over_q(p.tau, q, p.bg) * df0_dlnq(q) * psi_0.value();
}



double theta_integrand(double q, void* parameters) {
    fluid_perturbation_integrand_parameters& p =
        *static_cast<fluid_perturbation_integrand_parameters*>(parameters);

    Measurement psi_1 = compute_psi_l(1, p.tau, p.k, q, p.tau_lambda,
                                      p.metric_psi_at_k_and_tau_lambda, p.bg,
                                      p.metric_psi, p.inner_integrator);

    return pow(q,3) * df0_dlnq(q) * psi_1.value();
}



double sigma_integrand(double q, void* parameters) {
    fluid_perturbation_integrand_parameters& p =
        *static_cast<fluid_perturbation_integrand_parameters*>(parameters);

    Measurement psi_2 = compute_psi_l(2, p.tau, p.k, q, p.tau_lambda,
                                      p.metric_psi_at_k_and_tau_lambda, p.bg,
                                      p.metric_psi, p.inner_integrator);

    return pow(q,3) * 1.0/eps_over_q(p.tau, q, p.bg) * df0_dlnq(q) * psi_2.value();
}



double sigma_integrand_psi_2_interpolated(double q, void* parameters) {
    fluid_perturbation_integrand_parameters& p =
        *static_cast<fluid_perturbation_integrand_parameters*>(parameters);
    return pow(q,3) * 1.0/eps_over_q(p.tau, q, p.bg) * df0_dlnq(q) * (*p.psi)(q);
}



Measurement Perturbations::integrate_fluid_perturbation(
        double (*integrand)(double, void*),
        double q_min,
        const Interpolation1D* psi
        )
{
    fluid_perturbation_integrand_parameters p = {tau, k, q_min, tau_lambda,
        metric_psi_at_k_and_tau_lambda, bg, metric_psi, inner_integrator, psi};
    Measurement result = outer_integrator.integrate(integrand, &p, 0, cutoff);
    return result;
}



void Perturbations::interpolate_psi_l(
        int l,
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
        q_grid[i] = q_min * std::pow(q_max / q_min, static_cast<double>(i) /
                static_cast<double>(n_points - 1));
    }

    Measurement result;

    for (size_t i = 0; i < n_points; ++i) {
        result = compute_psi_l(l, tau, k, q_grid[i], tau_lambda,
                metric_psi_at_k_and_tau_lambda, bg,
                metric_psi, inner_integrator);
        psi_vals[i] = result.value();
    }

    psi_spline = Interpolation1D(q_grid, psi_vals);
}



Perturbations::Perturbations(
        double tau,
        double k,
        double z_lambda,
        double cutoff,
        const Background& bg,
        const Interpolation2D& metric_psi,
        const Integrator& inner_integrator,
        const Integrator& outer_integrator,
        bool do_interpolate_psi
        ) :
    tau(tau), k(k), cutoff(cutoff), do_interpolate_psi(do_interpolate_psi),
    bg(bg), metric_psi(metric_psi), inner_integrator(inner_integrator),
    outer_integrator(outer_integrator)
{
    tau_lambda = bg.conf_time_of_redshift(z_lambda);
    metric_psi_at_k_and_tau_lambda = metric_psi(k, tau_lambda);
}



void Perturbations::compute() {
    double a = bg.scale_factor_of_conf_time(tau);
    double factor = 4 * M_PI * std::pow(bg.T_ncdm0/a, 4);

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
            interpolate_psi_l(2, q_min, psi_2_spline);

            sigma = integrate_fluid_perturbation(
                sigma_integrand_psi_2_interpolated, q_min, &psi_2_spline);
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
