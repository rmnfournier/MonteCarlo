/**
 @file  MonteCarlo.hh
 @brief  Abstract class helping to implement Monte Carlo Methods
 @author  Romain Fournier
 @date  05.02.2019
**/

#ifndef MONTECARLO_INTERACTING_HARMONIC_OSCILLATOR_HH
#define MONTECARLO_INTERACTING_HARMONIC_OSCILLATOR_HH

#include "MonteCarlo.hh"
class Interacting_Harmonic_Oscillator : public MonteCarlo{
    public:
        Interacting_Harmonic_Oscillator(unsigned int N_, double initial_variance_, unsigned int bloc_size_,
                                    unsigned int nb_blocs_, double step_variance_, unsigned int points_to_update_,
                                    unsigned int warmup_steps_, unsigned int steps_between_samples_,
                                    const string &outputfile_name_,  double beta_,
                                    double omega_0_, double m_, double eps_0_, double a1_, double a2_, double alpha_1_,
                                    double alpha_2_, double f_, double t_max_, double dt_, double omega_max_,
                                    double domega_, unsigned int F_modes_, unsigned int current_tau_);


        double p(const vector<double> & x) override;
        double sampling_term() override;
    private:
        /**
         * \brief Initialize the friction kernel
         */
            void init_xi();
        /**
        * \brief Initialize the thermal bath
        */
        void init_J();
        /**
         * \brief Initialize the Fourier Spacing
         */
        void init_Omega();
        /**
         * \brief Initialize the modes of the interaction
         */
        void init_Gamma();

        /**
         * \brief Compute the norm of the fourier coefficients of the trajectory x
         * \param x trajectroy to consider
         */
        void compute_x_norm_fourier(const vector<double> & x);
        /**
         * \brief compute the action including the influence functional
         * @param x trajectory to consider
         * @return S
         */
        double influence_functional(const vector<double> & x);


    // Term to consider
        unsigned int current_tau_; /** Term to average */
    // Harmonic Oscillator
        double beta_; /** Inverse temperature */
        double omega_0_; /** Frequency of the harmonic oscillator */
        double m_; /** mass of the harmonic oscillator */
    // Friction Kernel
        double eps_0_;
        double a1_,a2_,alpha_1_,alpha_2_,f_; /** Parameter of the friction kernel */
        double t_max_,dt_; /** Threshold and time step for the integration */
        vector<double> xi_; /** Discretized friction kernel */
    // Thermal Bath
        vector<double> Jb_; /** Discretized thermal bath */
        double omega_max_, domega_;/** Threshold and frequenvy step for the integration */
    // Analytic solution
        unsigned int F_modes_; /** Number of Fourier modes to consider */
        vector<double> Omega_,Gamma_; /** Part of the analytic solution */
        vector<double> x_norm_fourier_squared_; /** Squared Norm of the Fourier coefficients of the trajectory */
};


#endif //MONTECARLO_INTERACTING_HARMONIC_OSCILLATOR_HH
