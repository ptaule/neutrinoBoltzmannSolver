/*
   integrator.hpp

   Created by Petter Taule on 09.04.2021
   Copyright (c) 2021 Petter Taule. All rights reserved.
*/

#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <stdexcept>

#include <gsl/gsl_integration.h>

class Measurement;

/* Wrapper around gsl_integration quadrature, adaptive general integrand */
class Integrator {
    private:
        gsl_integration_workspace* workspace;
        std::size_t sub_regions;

        double abs_tol;
        double rel_tol;
    public:
        Integrator(
                std::size_t sub_regions = 10000,
                double abs_tol = 0,
                double rel_tol = 1e-8
                ) : sub_regions(sub_regions), abs_tol(abs_tol), rel_tol(rel_tol)
        {
            workspace = gsl_integration_workspace_alloc(sub_regions);
        }
        Integrator(const Integrator&) = delete;
        Integrator& operator=(const Integrator&) = delete;

        Measurement integrate (
                double (*f)(double, void*),
                void* params,
                double x_min,
                double x_max,
                double abs_tol,
                double rel_tol,
                int key = GSL_INTEG_GAUSS61
                ) const;
        Measurement integrate (
                double (*f)(double, void*),
                void* params,
                double x_min,
                double x_max
                ) const;

        ~Integrator() {
            gsl_integration_workspace_free(workspace);
        }
};


class IntegrationException : std::exception {
    private:
        bool critical_;
        std::string ex;
    public:
        IntegrationException(int gsl_status);
        bool critical() const {return critical_;}
        const char* what() const noexcept {return ex.c_str();}
};


#endif /* ifndef INTEGRATOR_HPP */
