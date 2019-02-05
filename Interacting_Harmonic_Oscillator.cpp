//
// Created by romain on 05/02/19.
//

#include "Interacting_Harmonic_Oscillator.hh"
#include <cmath>
#include <vector>
#include <iostream>
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
                                                                 unsigned int F_modes_, unsigned int current_tau_) : MonteCarlo(N_,
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
                                                                                          x_norm_fourier_squared_(F_modes_,0),Gamma_(F_modes_,0),xi_(ceil(t_max_/dt_),0)

{
    cout<<"1"<<endl;
    init_Omega();
    cout<<"2"<<endl;

    init_xi();
    cout<<"3"<<endl;

    init_J();
    cout<<"4"<<endl;

    init_Gamma();
    cout<<"5"<<endl;

}

void Interacting_Harmonic_Oscillator::init_Omega() {
    for (unsigned int i(0);i<F_modes_;i++){
        Omega_[i]=2*M_PI*i/beta_;
    }
}

void Interacting_Harmonic_Oscillator::init_xi() {
    // Initialize time
    double t(0);
    for(auto& el : xi_){
        // Compute the term for the current time
        el= eps_0_*(exp (-alpha_1_*pow(f_*t,2))*(1.+a1_*pow(f_*t,4)) +a2_*pow(f_*t,4)*exp (-alpha_2_*pow(f_*t,2))  );
        // update time
        t+=dt_;
    }
}

void Interacting_Harmonic_Oscillator::init_J() {
    // Initialize frequency
    double w(0);
    cout<<Jb_.size()<<endl<<xi_.size()<<endl;
    for(auto& el : Jb_){
        // Perform the cosine transform
        double t(0);
        for(auto it=xi_.begin();it!=xi_.end();it++){
            // Count boundary terms once and Inner terms twice
            double factor(0);
            (it == xi_.begin() or it ==xi_.end()) ? factor=dt_/M_PI:factor=2*dt_/M_PI;
            el+=*it *cos(w*t)*factor;
            t+=dt_;
        }
        //update the frequency
        w+=domega_;
    }
}

void Interacting_Harmonic_Oscillator::init_Gamma() {
    // We start by computing gamma in imaginary time, and then we save its coefficients
    vector<double> gamma_im_time(ceil(t_max_/dt_),0);
    double tau(0),dtau(0.5*beta_/ceil(t_max_/dt_));

    for(auto& el : gamma_im_time){
        double w(domega_/100.); // Initialize w close to 0 to avoid dividing by 0
        for(auto it=Jb_.begin();it!=Jb_.end();it++){
            double factor(0);
            (it == Jb_.begin() or it ==Jb_.end()) ? factor=domega_*beta_:factor=2*domega_*beta_;
            el+=*it *w*cosh(beta_*w*(0.5-tau/beta_))/sinh(0.5*beta_*w)*factor;
            w+=domega_;
        }
        tau+=dtau;
    }

    // We can now compute the Fourier coefficients
    for(size_t i(0);i<F_modes_;i++){
        tau=0;
        for(auto it=gamma_im_time.begin();it!=gamma_im_time.end();it++){
            double factor(0);
            (it == Jb_.begin() or it ==Jb_.end()) ? factor=dtau/beta_*2:factor=2*dtau/beta_*2;
            Gamma_[i]+=factor* *it *cos(tau*Omega_[i]);
            tau+=dtau;
        }
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

double Interacting_Harmonic_Oscillator::influence_functional(const vector<double> & x){
    compute_x_norm_fourier(x);
    double S(0);

    //First term
    S+=0.5*(m_*omega_0_*omega_0_-Gamma_[0])*x_norm_fourier_squared_[0];
    // other terms
    for (size_t i(1);i<x_norm_fourier_squared_.size();i++){
        S+=(m_*Omega_[i]*Omega_[i]+m_*omega_0_*omega_0_-Gamma_[i])*x_norm_fourier_squared_[i];
    }
    return -beta_*S;
}

double Interacting_Harmonic_Oscillator::sampling_term(){
    return x_[0]*x_[current_tau_];
}

double Interacting_Harmonic_Oscillator::p(const vector<double> & x){
    return exp(-influence_functional(x));
}