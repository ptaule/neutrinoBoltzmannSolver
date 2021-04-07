/*
   measurement.hpp

   Created by Petter Taule on 10.06.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef MEASUREMENT_HPP
#define MEASUREMENT_HPP

#include <iosfwd>

class Measurement {
        double value_ = 0;
        double error_ = 0;
    public:
        Measurement() : value_(0), error_(0) {}
        Measurement(double value, double error) :
            value_(value), error_(error) {}

        double value() const {return value_;};
        double& value() {return value_;};

        double error() const {return error_;};
        double& error() {return error_;};

        Measurement& operator+=(const Measurement& b);
        Measurement operator+(const Measurement& b);

        Measurement& operator-=(const Measurement& b);
        Measurement operator-(const Measurement& b);

        Measurement& operator*=(const Measurement& b);
        Measurement operator*(const Measurement& b);

        Measurement& operator/=(const Measurement& b);
        Measurement operator/(const Measurement& b);

        Measurement& operator+=(double scalar);
        Measurement operator+(double scalar);

        Measurement& operator-=(double scalar);
        Measurement operator-(double scalar);

        Measurement& operator*=(double scalar);
        Measurement operator*(double scalar);
};

std::ostream& operator<<(std::ostream& out, const Measurement& m);


#endif /* ifndef MEASUREMENT_HPP */
