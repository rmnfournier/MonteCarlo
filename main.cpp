#include <iostream>
#include "Interacting_Harmonic_Oscillator.hh"


int main() {
    std::cout << "Hello, World!" << std::endl;
   Interacting_Harmonic_Oscillator sim (64, 0.0025,
                                         10, 10,
                                         0.01, 1,
                                         1000,
                                         100,
                                         "test.csv", 0.25,
                                         20,  1,  225,  1.486*10000,
                                          285, 903, 75,
                                          0.2, 5, 0.0005,
                                          75,  0.001,
                                          64, 1);


   sim.sample();
   return 0;
}