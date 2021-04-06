/*
   measurement.cpp

   Created by Petter Taule on 10.06.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <cmath>
#include <ostream>
#include <iostream>

#include "../include/measurement.hpp"


Measurement& Measurement::operator+=(const Measurement& b) {
    value_ += b.value_;
    /* Gaussian error propagation */
    error_ = sqrt(pow(error_, 2) + pow(b.error_, 2));
    return *this;
}



Measurement Measurement::operator+(const Measurement& b) {
    Measurement result = *this;
    return result += b;
}



Measurement& Measurement::operator-=(const Measurement& b) {
    value_ -= b.value_;
    error_ = sqrt(pow(error_, 2) + pow(b.error_, 2));
    return *this;
}



Measurement Measurement::operator-(const Measurement& b) {
    Measurement result = *this;
    return result -= b;
}



Measurement& Measurement::operator*=(const Measurement& b) {
    /* Gaussian error propagation. Note that value is set after error, so that
     * the appropriate (old) value is used in computation of error */
    error_ = sqrt(
            pow(b.value_ * error_ , 2) +
            pow(value_ * b.error_ , 2)
            );
    value_ *= b.value_;
    return *this;
}



Measurement& Measurement::operator/=(const Measurement& b) {
    /* Gaussian error propagation. Note that value is set after error, so that
     * the appropriate (old) value is used in computation of error */
    error_ = sqrt(
            pow(error_ / b.value_, 2) +
            pow(value_ * b.error_, 2) / pow(b.value_, 4)
            );
    value_ /= b.value_;
    return *this;
}



Measurement Measurement::operator/(const Measurement& b) {
    Measurement result = *this;
    return result /= b;
}



Measurement& Measurement::operator+=(double scalar) {
    value_ += scalar;
    return *this;
}



Measurement Measurement::operator+(double scalar) {
    Measurement result = *this;
    return result += scalar;
}



Measurement& Measurement::operator-=(double scalar) {
    value_ -= scalar;
    return *this;
}



Measurement Measurement::operator-(double scalar) {
    Measurement result = *this;
    return result -= scalar;
}



Measurement& Measurement::operator*=(double scalar) {
    value_ *= scalar;
    error_ *= scalar;
    return *this;
}



Measurement Measurement::operator*(double scalar) {
    Measurement result = *this;
    return result *= scalar;
}



std::ostream & operator<<( std::ostream & out, const Measurement& m) {
    out << m.value() << "\t+- " << m.error();
    return out;
}
