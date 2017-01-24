/* 
 MIT License
 
 Copyright (c) 2016 - Jacky Han
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

#ifndef PROC_H
#define PROC_H

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

/*
** NOTE: Not every function in this header file is used. Some of them
** were used in older versions of the program and so they were kept
** in case any comparisons needed to be made or switches were done.
** @Author: Jacky Han, Samuel Foley
** @Date: July 15 2016
*/

const int COL_W = 13; //output column width

//Resting Energy Expenditure (Basal Metabolic Rate) (Watts)
//arguments: mass in kg, height in cm, age in years, and sex (male==true, female==false)
double REE(double m, double h, double a, bool male)
{
    if(male)
    {
        return 0.0484*(10.0*m + 6.2*h - 5.0*a + 5);
    }
    else
    {
        return 0.0484*(10.0*m + 6.2*h - 5.0*a - 161);
    }
}

//convective heat transfer (W)
//arguments: BSA (m^2), skin temp (celcius), ambient temp (C), wind speed (m/s)
double convH(double A, double Ts, double Ta, double v)
{
    //8.0 is the approximate convective coefficient (wheeler)
    return A * 8.0 * (Ts - Ta) * std::sqrt(v);
}

/*
//DuBois Body Surface Area (m^2)
//arguments: mass in kg, height in cm
double BSA(double m, double h)
{
    return 0.007184 * std::pow(m, 0.425) * std::pow(h, 0.725);
}
*/

//Mosteller formula, as recommended by Verbraecken, et al
//arguments: mass in kg, height in cm
double BSA(double m, double h)
{
    return std::sqrt(m*h/3600);
}

//heat lost through evaporation of sweat (Watts)
//arguments: BSA (m^2), wind speed (m/s), ambient vapor pressure,
//saturation vapor pressure at skin temp, ambient air pressure
double Esw(double A, double v, double VPa, double SVPts, double P)
{
    return A * 2416 * 0.00277 * v * 0.662 * (VPa - SVPts)/P;
}

//Black-Body radiation formula
//arguments: BSA, temperature (C), skin temp (C)
double bbRad(double A, double T, double Ts)
{
    return 0.0000000567 * (std::pow(T + 273.15, 4) - std::pow(Ts + 273.15, 4) );
}

//solar radiation heat
//arguments: BSA (m^2), skin reflectivity
double solarRad(double A, double refl)
{
    //the 0.25 term comes from cross-sectional area
    return 0.25 * A * 0 * (1.0 - refl);
}

//sherwood & huber sensible heat
double sherS(double v, double Ts, double T)
{
    return 12.5*v*(T - Ts);
}

//saturation vapor pressure at T (deg C)
double satVP(double T)
{
    return 6.108 * (std::exp((17.27 * T) / (237.3 + T)));
}

//vapor pressure (partial pressure of h2o) at given conditions
double vapP(double T, double Tw)
{
    double Es, Ew;

    Es = satVP(T);
    Ew = satVP(Tw);

    return Ew - (0.00066 * (1 + 0.00115 * Tw) * (T - Tw) * 101.3);
}

//kerslake evaporation heat loss
double evap(double v, double T, double Ts, double Tw)
{
    return 12.4 * std::sqrt(v) * (vapP(T, Tw) - satVP(Ts) );
}

//kerslake convection
double conv(double v, double T, double Ts)
{
    return 8.3 * std::sqrt(v) * (T - Ts);
}

//sherwood & huber latent heat
double sherL(double v, double T, double Tw)
{
    return 12.5*v*(Tw - T);
}

//blood flow in liters/hm^2
double vb(double Tc, double Ts)
{
    // /*Hope 1993*/ return (6.3 + 75*(Tc - 36.6))/(1 + 0.5*(34 - Ts));
    return 0.7*((2.07*Tc)-75.44)*(  (100.0/3.14159)*std::atan(0.75*(Ts-34.7)) + 53);
}

//ratio of shell volume to body volume
double alpha(double vb)
{
    return 0.044 + (0.35/(vb - 0.1386));
}

//core surface area
double coreA(double alpha, double BSA)
{
    return std::pow(1 - alpha, 2.0/3.0) * BSA;
}

//energy exchange between core and shell
//positive is core -> shell
double fc(double Tc, double Ts, double BSA)
{
    double flow = vb(Tc, Ts);
    return coreA(alpha(flow), BSA) * flow * (10e-7) * 1060 * 3860*(Tc - Ts) + bbRad(coreA(alpha(flow), BSA), Tc, Ts);
}

//sweat rate (hoppe 1993)
double SW(double A, double Ts, double Tc)
{
    return A * 8.47 * (10e-5) * ((0.1*Ts + 0.9*Tc) - 36.6);
}

//convection between core and shell
double coreConv(double Tc, double Ts, double BSA)
{
    return coreA(alpha(1), BSA)*10*(Tc - Ts);
}

//output numbers, used in the main loop to output the numerical values
void outputN(std::ostream &o, double _1, double _2, double _3, double _4, double _5, double _6,
                double _7, double _8 = 0, double _9 = 0, double _10 = 0, double _11 = 0)
{
    o << std::left
        << std::setw(COL_W) << _1 << std::setw(COL_W) << _2 << std::setw(COL_W) << _3
        << std::setw(COL_W) << _4 << std::setw(COL_W) << _5 << std::setw(COL_W) << _6
        << std::setw(COL_W) << _7 << std::setw(COL_W) << _8 << std::setw(COL_W) << _9
        << std::setw(COL_W) << _10 << std::setw(COL_W) << _11 << std::endl;
}

//output strings, used in the main loop to output string values
void outputS(std::ostream &o, std::string _1, std::string _2, std::string _3, std::string _4, std::string _5, std::string _6,
                std::string _7, std::string _8 = " ", std::string _9 = " ", std::string _10 = " ", std::string _11 = " ")
{
    o << std::left
        << std::setw(COL_W) << _1 << std::setw(COL_W) << _2 << std::setw(COL_W) << _3
        << std::setw(COL_W) << _4 << std::setw(COL_W) << _5 << std::setw(COL_W) << _6
        << std::setw(COL_W) << _7 << std::setw(COL_W) << _8 << std::setw(COL_W) << _9
        << std::setw(COL_W) << _10 << std::setw(COL_W) << _11 << std::endl;
}

#endif
