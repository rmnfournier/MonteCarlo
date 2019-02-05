/**
 @file  MonteCarlo.hh
 @brief  Abstract class helping to implement Monte Carlo Methods
 @author  Romain Fournier
 @date  05.02.2019
**/

#ifndef MONTECARLO_MONTECARLO_HH
#define MONTECARLO_MONTECARLO_HH

#include <vector>
#include <string>
#include <fstream>

using namespace std;

class MonteCarlo {
public:
    MonteCarlo(unsigned int N_, double initial_variance_, unsigned int bloc_size_, unsigned int nb_blocs_,
               double step_variance_, unsigned int points_to_update_, unsigned int warmup_steps_,
               unsigned int steps_between_samples_, const string &outputfile_name_) : N_(
            N_), initial_variance_(initial_variance_), bloc_size_(bloc_size_), nb_blocs_(nb_blocs_), step_variance_(
            step_variance_), points_to_update_(points_to_update_), warmup_steps_(warmup_steps_), steps_between_samples_(
            steps_between_samples_), outputfile_name_(outputfile_name_),x_(N_,0) {}

    /**
     * @brief generate a random configuration x according to the initial_variance_ attribut
     */
    void reset();
    /**
     * @brief Compute the probability of a configuration
     * @param x configuration to test
     * @return double in (0,1) corresponding to the probability of the current configuration
     */
    virtual double p(const vector<double> & x) = 0;

    /**
     * @brief try to make a random modification to the attribut x, and accept it according to the Metropolis-Hasting algorithm
     * @param -
     * @return true if there was a modification, false otherwise
     */
    bool metropolis();
    /**
     * @brief Perform nb_steps updates according to the metropolis-hasting algorithm
     * @param nb_steps number of changes that will be proposed
     * @return Number of accepted changes
     */
    unsigned int warmup(unsigned int nb_steps);

    /**
     * @brief Define the term that has to be sampled
     * @param -
     * @return Value of the term in the current configuration
     */
    virtual double sampling_term()=0;

    /**
     * @brief Perform the MonteCarlo Simulation, save the results in outputfile_
     */
    void sample();

    /**
     * @brief Open the file that saves the results
     */
    void open_output_file() {outputfile_.open(outputfile_name_);}
    /**
    * @brief Open the file that saves the results
    */
    void close_output_file(){outputfile_.close();}

    /**
     * @brief Add a line to the output file
     */
    void write_bloc_result(double);

protected:
    vector<double> x_; /** actual configuration **/
    unsigned int N_; /** Size of the vector x_ **/
    double initial_variance_; /** Initial variance of the distribution x **/

    unsigned int bloc_size_; /** size of a bloc Monte Carlo Bloc **/
    unsigned int nb_blocs_; /** Number of blocs to sample **/

    double step_variance_; /** Variance of the proposed steps **/
    unsigned int points_to_update_ ; /** Number of points to update **/
    unsigned int warmup_steps_; /** Number of updates in the warmup phase **/
    unsigned int steps_between_samples_; /** Number of updates between two consecutive samplings **/

    string outputfile_name_; /** Name of the outputfile */
    ofstream outputfile_; /** Outputfile */
};


#endif //MONTECARLO_MONTECARLO_HH
