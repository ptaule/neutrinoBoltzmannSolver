/*
   quantity.cpp

   Created by Petter Taule on 10.06.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <cmath>
#include <ostream>
#include <iostream>

#include "quantity.hpp"

Quantity::Quantity() {
    value = 0;
    error = 0;
}



Quantity::Quantity(double value, double error) {
    this->value = value;
    this->error = error;
}



Quantity::Quantity(const Quantity& other) {
    this->value = other.value;
    this->error = other.error;
}



Quantity& Quantity::operator=(const Quantity& other) {
    this->value = other.value;
    this->error = other.error;
    return *this;
}



Quantity& Quantity::operator+=(const Quantity& b) {
    this->value += b.value;
    // Gaussian error propagation
    this->error = sqrt(pow(this->error, 2) + pow(b.error, 2));
    return *this;
}



Quantity Quantity::operator+(const Quantity& b) {
    Quantity result = *this;
    return result += b;
}



Quantity& Quantity::operator-=(const Quantity& b) {
    this->value -= b.value;
    this->error = sqrt(pow(this->error, 2) + pow(b.error, 2));
    return *this;
}



Quantity Quantity::operator-(const Quantity& b) {
    Quantity result = *this;
    return result -= b;
}



Quantity& Quantity::operator/=(const Quantity& b) {
    // Gaussian error propagation
    // Note that value is set after error, so that the appropriate (old) value
    // is used in computation of error
    this->error = sqrt(
            pow(this->error / b.value, 2) +
            pow(this->value * b.error, 2) / pow(b.value, 4)
            );
    this->value /= b.value;
    return *this;
}



Quantity Quantity::operator/(const Quantity& b) {
    Quantity result = *this;
    return result /= b;
}



Quantity& Quantity::operator+=(double scalar) {
    this->value += scalar;
    return *this;
}



Quantity Quantity::operator+(double scalar) {
    Quantity result = *this;
    return result += scalar;
}



Quantity& Quantity::operator-=(double scalar) {
    this->value -= scalar;
    return *this;
}



Quantity Quantity::operator-(double scalar) {
    Quantity result = *this;
    return result -= scalar;
}



Quantity& Quantity::operator*=(double scalar) {
    this->value *= scalar;
    this->error *= scalar;
    return *this;
}



Quantity Quantity::operator*(double scalar) {
    Quantity result = *this;
    return result *= scalar;
}



std::ostream & operator<<( std::ostream & out, const Quantity & quantity ) {
    out << quantity.value << "\t+- " << quantity.error;
    return out;
}
