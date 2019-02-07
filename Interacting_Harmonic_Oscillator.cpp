//
// Created by romain on 05/02/19.
//

#include "Interacting_Harmonic_Oscillator.hh"
#include <cmath>
#include <vector>
#include <iostream>
#include <limits>
using namespace std;

Interacting_Harmonic_Oscillator::Interacting_Harmonic_Oscillator(unsigned int N_, double initial_variance_,
                                                                 unsigned int bloc_size_, unsigned int nb_blocs_,
                                                                 double step_variance_, unsigned int points_to_update_,
                                                                 unsigned int warmup_steps_,
                                                                 unsigned int steps_between_samples_,
                                                                 const string &outputfile_name_,
                                                                 double beta_,
                                                                 double omega_0_, double m_, double eps_0_, double a1_,
                                                                 double a2_, double alpha_1_, double alpha_2_,
                                                                 double f_, double t_max_, double dt_,
                                                                 double omega_max_, double domega_,
                                                                 unsigned int F_modes_, unsigned int current_tau_,bool verbose_) : MonteCarlo(N_,
                                                                                                     initial_variance_,
                                                                                                     bloc_size_,
                                                                                                     nb_blocs_,
                                                                                                     step_variance_,
                                                                                                     points_to_update_,
                                                                                                     warmup_steps_,
                                                                                                     steps_between_samples_,
                                                                                                     outputfile_name_),
                                                                                          beta_(beta_),
                                                                                          omega_0_(omega_0_), m_(m_),
                                                                                          eps_0_(eps_0_), a1_(a1_),
                                                                                          a2_(a2_), alpha_1_(alpha_1_),
                                                                                          alpha_2_(alpha_2_), f_(f_),
                                                                                          t_max_(t_max_), dt_(dt_),
                                                                                          omega_max_(omega_max_),
                                                                                          domega_(domega_),
                                                                                          F_modes_(F_modes_),current_tau_(current_tau_),Jb_(ceil(omega_max_/domega_),0),Omega_(F_modes_,0),
                                                                                          x_norm_fourier_squared_(F_modes_,0),Gamma_(N_,0),xi_(ceil(t_max_/dt_),0),verbose_(verbose_)

{
    if (verbose_) cout<<" Start Initialization "<<endl;
    init_Omega();

    init_xi();
    if(verbose_) cout<<"xhi initialized successfully"<<endl;
    init_J();
    if(verbose_) cout<<"J initialized successfully"<<endl;

    init_Gamma();
    if(verbose_) cout<<"Initialization complete "<<endl;
}

void Interacting_Harmonic_Oscillator::init_Omega() {
    for (unsigned int i(0);i<F_modes_;i++){
        Omega_[i]=2*M_PI*i/beta_;
    }
}

void Interacting_Harmonic_Oscillator::init_xi() {
    // Initialize time
    std::ofstream outfile;
    if(verbose_)outfile.open("xhi.csv");
    double t(0);
    for(auto& el : xi_){
        // Compute the term for the current time
        el= eps_0_*(exp (-alpha_1_*pow(f_*t,2))*(1.+a1_*pow(f_*t,4)) +a2_*pow(f_*t,4)*exp (-alpha_2_*pow(f_*t,2))  );
        if(verbose_) {outfile<<el <<",";}
        // update time
        t+=dt_;
    }
    if(verbose_)outfile.close();
}

void Interacting_Harmonic_Oscillator::init_J() {
    std::ifstream infile("Jb.csv");
    // Check if we can load J instead of computing it
    if(infile.is_open()){
        if (verbose_) cout<<"Reading J"<<endl;
        //preprocess
        string line;
        // Non Robust and bad way to read, but it's only an example so who cares :P
        for(auto& el : Jb_){
            getline(infile,line,',');
            el=stod(line);
        }
    }
    else {
        std::ofstream outfile;
        if(verbose_) outfile.open("Jb.csv");
        // Initialize frequency
        double w(0);
        for (auto &el : Jb_) {
            // Perform the cosine transform
            double t(0);
            for (auto it = xi_.begin(); it != xi_.end(); it++) {
                // Count boundary terms once and Inner terms twice
                double factor(0);
                (it == xi_.begin() or it == xi_.end()) ? factor = dt_ / M_PI : factor = 2 * dt_ / M_PI;
                el += *it * cos(w * t) * factor;
                t += dt_;
            }
            if(verbose_) outfile<<el<<',';
            //update the frequency
            w += domega_;
        }
        if(verbose_) outfile.close();
    }
    infile.close();

}

void Interacting_Harmonic_Oscillator::init_Gamma() {
    std::ifstream infile("Gamma.csv");
    // Check if we can load J instead of computing it
    if(infile.is_open()){
        if (verbose_) cout<<"Reading Gamma"<<endl;
        //preprocess
        string line;
        // Non Robust and bad way to read, but it's only an example so who cares :P
        for(auto& el : Gamma_){
            getline(infile,line,',');
            el=stod(line);
        }
    }
    else {
        std::ofstream outfile;
        if (verbose_) outfile.open("Gamma.csv");
        // We start by computing gamma in imaginary time, and then we save its coefficients
        double tau(0), dtau(beta_ /(N_+0.0));

        for (auto &el : Gamma_) {
            double w(domega_ / 100.); // Initialize w close to 0 to avoid dividing by 0
            for (auto it = Jb_.begin(); it != Jb_.end(); it++) {
                double factor(0);
                (it == Jb_.begin() or it == Jb_.end()) ? factor = domega_ * beta_ : factor = 2 * domega_ * beta_;
                el += *it * w * cosh(beta_ * w * (0.5 - tau / beta_)) / sinh(0.5 * beta_ * w) * factor;
                w += domega_;
            }
            if(verbose_)outfile<<el<<",";
            tau += dtau;
        }


        if(verbose_) outfile.close();
    }
}

void Interacting_Harmonic_Oscillator::compute_x_norm_fourier(const vector<double> & x) {
    double beta_inv(1/beta_);
    double dtau(beta_/N_);
    double x_norm_fourier_real(0);
    double x_norm_fourier_im(0);

    for(size_t i(0);i<x_norm_fourier_squared_.size();i++){
        double tau(0);
        double factor(2*beta_inv*dtau);
        x_norm_fourier_im=0;
        x_norm_fourier_real=0;
        for(size_t j(0);j<x_.size();j++){
            x_norm_fourier_real+=factor*x[j]*cos(Omega_[i]*tau);
            x_norm_fourier_im+=factor*x[j]*sin(Omega_[i]*tau);
            tau+=dtau;
        }
        x_norm_fourier_squared_[i]=(pow(x_norm_fourier_im,2)+pow(x_norm_fourier_real,2));
    }
}

double Interacting_Harmonic_Oscillator::action(const vector<double> &x){
    double E_cin(0),E_pot(0), influence(0); // Initialize the three contributions to the action
    double dtau(beta_/(x.size()+0.0));
    for (size_t i(0);i<x.size();i++){
        // Compute Kinetic part
        E_cin+=pow(x[(i+1)%x.size()]-x[i],2);
        // Potential part
        E_pot+=pow(x[i],2);
    }
    E_cin*=m_*0.5/dtau; //only dtau, since I won't multiply the whole action by dtau at the end
    E_pot *= omega_0_*omega_0_*0.5*m_*dtau;
    // influence part
    for(size_t i(0);i<x.size();i++){
        for(size_t j(0);j<x.size();j++){
            influence+=x[i]*x[j]*Gamma_[abs(i-j)];
        }
    }
    influence*=-dtau*dtau*0.5*beta_;

    return E_cin+E_pot+influence;
}

double Interacting_Harmonic_Oscillator::sampling_term(){
    return x_[0]*x_[current_tau_];
}

double Interacting_Harmonic_Oscillator::p(const vector<double> & x){
    return exp(-action(x));
}