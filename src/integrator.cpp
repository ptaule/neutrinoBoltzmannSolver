/*
   integrator.cpp

   Created by Petter Taule on 09.04.2021
   Copyright (c) 2021 Petter Taule. All rights reserved.
*/

#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "../include/integrator.hpp"
#include "../include/measurement.hpp"

struct y_integration_parameters;

Measurement Integrator::integrate (
                double (*f)(double, void*),
                void* params,
                double x_min,
                double x_max,
                double atol,
                double rtol,
                int key
                ) const
{
    Measurement result;

    gsl_function F;
    F.function = f;
    F.params = params;

    int status = gsl_integration_qag(&F, x_min, x_max, atol, rtol, sub_regions,
            key, workspace, &result.value(), &result.error());

    if (status != 0) {
        IntegrationException ex(status);
        if (ex.critical()) {
            throw ex;
        }
        else {
            std::cout << ex.what() << "\t" << "Result: " << result <<
                std::endl;
        }
    }
    return result;
}



Measurement Integrator::integrate (
        double (*f)(double, void*),
        void* params,
        double x_min,
        double x_max
        ) const
{
    return integrate(f, params, x_min, x_max, abs_tol, rel_tol);
}



IntegrationException::IntegrationException(int gsl_status)
{
    switch (gsl_status) {
        case GSL_EMAXITER:
            ex = "Maxiumum number of subdevisions reached.";
            critical_ = true;
            break;
        case GSL_EROUND:
            ex = "Warning: rounding error during integration.";
            critical_ = false;
            break;
        case GSL_ESING:
            ex = "Singularity or bad behaviour in integration.";
            critical_ = true;
            break;
        case GSL_EDIVERGE:
            ex = "Divergent or too slowly converging integral.";
            critical_ = true;
            break;
        case GSL_EDOM:
            ex = "Error in the values of input arguments to the integrand.";
            critical_ = true;
            break;
        default:
            ex = "Unknown status returned from integration.";
    }
}
