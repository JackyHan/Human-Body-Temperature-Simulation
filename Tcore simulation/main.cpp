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

    #include <iostream>
    #include <iomanip>
    #include <fstream>
    #include <string>
    #include <cmath>

    #include "processes.h"

    /*
    ** This program was designed to run in a linux terminal window that
    ** is 143 characters wide, and as such might not display nicely on
    ** other widths unless properly modified. However, the output file
    ** should have nice formating regardless of the running environment.
    **
    ** Most of the computation and ugly details are in the header file
    ** called "processes.h", included above.
    **
    ** The paramters for each run must be specified in a text file called
    ** "config.txt". The order in which the paramters must be listed can
    ** be found below in the code (the several "infile >>" lines).
    ** 
    ** @Author: Jacky Han, Samuel Foley
    ** @Date: July 15 2016
    */

    int main(int argc, char** argv)
    {
        //The program takes as a command line argument the name of a file
        //to which it will save the run's output as text.
        //If no argument is given, the default of "run.txt" is used.
        //The following code implements this
        std::string filename;
        if(argc == 1)
        {
            filename = "run.txt";
        }
        else
        {
            filename = argv[1];
        }
        
        //Create output file stream and open file for output
        std::ofstream outfile;
        outfile.open(filename.c_str() );
        //data output file
        std::ofstream outfile2;
        outfile2.open("data.txt");
        outfile2 << std::left<< std::setw(COL_W) << "Tweb" << std::setw(COL_W) << "Tcore" << std::endl;

        //Create input file stream and open the configuration file
        std::ifstream infile;
        infile.open("config.txt");

        //The following lines set the default values for various parameters
        //before the config values are read in
        const bool MALE   = true;
        const bool FEMALE = false; //for use in REE function argument

        double m = 80;        //mass in kg
        double h = 185;       //height in cm
        double a = 25;        //age in years
        double refl = 0.50;   //skin reflectivity

        double T  = 30;       //dry temp (C)
        double TwL = 22;       //wet bulb temp (C) - low
        double TwH = 35;       //wet bulb temp (C) - high
        double v  = 5.0;       //wind speed (m/s)

        char gend;            //gender ('m' or 'f')

        double met = 0;       //metabolic heat (0 means REE) (in watts)
        
        double L;            //evaporation

        //The next several lines read in the desired parameter values from
        //the file "config.txt". Refer to comments above for their meaning.
        infile >> m;
        infile >> h;
        infile >> a;
        infile >> refl;
        infile >> T;
        infile >> TwL;
        infile >> TwH;
        infile >> v;
        infile >> met;
        infile >> gend;

        infile.close(); //close the config file after reading it

        double A = BSA(m, h); //body surface area (Mosteller)

        if(met == 0)
        {
            if(gend == 'm')
            {
                met = REE(m, h, a, MALE);
            }
            else
            {
                met = REE(m, h, a, FEMALE);
            }
        }
        
        
        //Set up initial values and create vars for later use
        double Tc = 37;       //core temp (initial)
        double Ts = 35;       //skin temp (initial)
        double W = 0;         //necessary water intake (initial)
        double S, M, Rs, Rb, Fs, Fc, coreFlow;
        int i2 = 0;
         
        for (double tempt = TwL; tempt < TwH+0.001; tempt = tempt + 0.02) {
            //Write out the column headers to the output file and console
            outfile << std::endl;
            outfile << "Tweb: " << tempt << std::endl;
            outputS(outfile, "sec", "conv", "evap", "met", "solar", "bbRad", "shell", "core", "Tskin", "Tcore", "water(L)");
            //outputS(std::cout, "sec", "conv", "evap", "met", "solar", "bbRad", "shell", "core", "Tskin", "Tcore", "water(L)");
            
            Fs = 1; //so that the loop initiates
            Fc = 1; //ditto

            //Run the sim for as long as core temp is below 42 or until the core or skin flux is tiny
            for(int i = 0; (Tc < 42)&&((std::abs(Fs) > 1e-4)||(std::abs(Fc) > 1e-4)); i++)
            {
                S  = conv(v, T, Ts);       //sensible heat (convection)
                
                L  = evap(v, T, Ts, tempt);   //Latent heat (evaporation)
                
                M  = met;                  //metabolic heat
                Rs = solarRad(A, refl);    //Solar radiation
                Rb = bbRad(A, T, Ts);      //Blackbody radiation
                coreFlow = fc(Tc, Ts, A);  //Core heat flow

                Fs = 0.1*(S + L + Rs + Rb + coreFlow); //heat flow at skin interface (0.1 corresponding to 1/10 of a sec)
                Fc = 0.1*(M - coreFlow);               //heat flow at core/shell interface

                W  -= L / 2416000;                     //a test value that never got implemented related to water usage
                Ts += Fs / (3874* (alpha(vb(Tc, Ts)) *m));       //skin(shell) temp
                Tc += Fc / (3874* ( (1-alpha(vb(Tc, Ts)) )*m));  //core temp
            
                i2 = i;

            
                //Output a live update every 600 iterations
                if(i % 600 == 0)
                {
                    outputN(outfile, i/10, S, L, M, Rs, Rb, Fs, Fc, Ts, Tc, W);
                }
            
            }
            
            
            //Print out the column headers along the bottom
            outputN(outfile, i2/10, S, L, M, Rs, Rb, Fs, Fc, Ts, Tc, W);
            outfile2 << std::left<< std::setw(COL_W) << tempt << std::setw(COL_W) << Tc << std::endl;
        }

        outfile.close(); //close the output file
        outfile2.close();

        return 0;
    }
