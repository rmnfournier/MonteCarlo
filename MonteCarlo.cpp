/**
 @file  MonteCarlo.cpp
 @brief  Implementation of the methods common to all MonteCarlo simulations
 @author  Romain Fournier
 @date  05.02.2019
**/

#include "MonteCarlo.hh"
#include <random>
#include <ctime>
#include <iostream>
void MonteCarlo::reset (){
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<double> distribution(0.0,initial_variance_);

    for (auto& el : x_){
        el = distribution(generator);
    }
}

bool MonteCarlo::metropolis(){
    bool accepted(false);
    std::random_device rd;
    std::mt19937 generator(rd());

    std::normal_distribution<double> distribution(0.0,step_variance_);
    std::uniform_real_distribution<double> uniform(0,1);
    std::uniform_int_distribution<int> uniint(0,N_-1);

    vector<double>x_new(x_);
    // Propose a new configuration x
    for (unsigned int i(0);i<points_to_update_;i++){
        //select slice to move
        int slice = uniint(generator);
        // generate a move
        double delta_xi = distribution(generator);
        x_new[slice]+=delta_xi;
    }
    // accepts it according to metropolis-hasting rule
    if(p(x_new)/p(x_)>=uniform(generator)){
        x_=std::move(x_new);
        accepted=true;
    }
    return accepted;
}

unsigned int MonteCarlo::warmup(unsigned int nb_steps){
    double ratio(0);
    unsigned int accepted(0);

    do {
        accepted=0;
        for (unsigned int i(0); i < nb_steps; i++) {
            if (metropolis()) accepted++; // Increment accepted if metropolis returns true, i.e. a change has been made
        }
       ratio= accepted / (nb_steps + 0.0);
        //update the step if the ratio does not correspond to the ideal ones.
        if(ratio<0.2) step_variance_*=1.1;
        if(ratio>0.4) step_variance_/=1.15;
        cout<<" New step variance = "<<step_variance_<<"; ratio was : "<<ratio<<endl;
    }while (ratio<0.2 or ratio>0.4);
    return accepted;
}

void MonteCarlo::sample(){
    open_output_file();
    // For each bloc
    for(unsigned int sampling(0);sampling<nb_blocs_;sampling++){
        // Initialize the mean and the variance for the current bloc
        double bloc_mean(0);
        reset();
        warmup(warmup_steps_);

        for(unsigned int bloc(0);bloc<bloc_size_;bloc++) {
            // Perform some updates before sampling again
            warmup(steps_between_samples_);
            // Update the mean value of the term
            bloc_mean+=sampling_term();
        }
        // Devide by the bloc_size to get the mean
        bloc_mean/=(bloc_size_+0.0);
        // Save the results for analysis
        write_bloc_result(bloc_mean);
    }
    close_output_file();
}

void MonteCarlo::write_bloc_result(double bloc_mean){
    outputfile_<<bloc_mean<<",";
    outputfile_.flush();
}