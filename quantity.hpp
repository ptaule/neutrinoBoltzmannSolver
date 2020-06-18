/*
   quantity.hpp

   Created by Petter Taule on 10.06.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef QUANTITY_HPP
#define QUANTITY_HPP

#include <iosfwd>

class Quantity {
    public:
        double value = 0;
        double error = 0;

        Quantity();
        Quantity(double value, double error);
        Quantity(const Quantity& other);

        Quantity& operator=(const Quantity& other);

        Quantity& operator+=(const Quantity& b);
        Quantity operator+(const Quantity& b);

        Quantity& operator-=(const Quantity& b);
        Quantity operator-(const Quantity& b);

        Quantity& operator/=(const Quantity& b);
        Quantity operator/(const Quantity& b);

        Quantity& operator+=(double scalar);
        Quantity operator+(double scalar);

        Quantity& operator-=(double scalar);
        Quantity operator-(double scalar);

        Quantity& operator*=(double scalar);
        Quantity operator*(double scalar);
};

std::ostream & operator<<( std::ostream & out, const Quantity & quantity );


#endif /* ifndef QUANTITY_HPP */
