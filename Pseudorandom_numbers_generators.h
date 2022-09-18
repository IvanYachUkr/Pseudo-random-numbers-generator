#include <functional>
#include <vector>
#ifndef LAB1_PSEUDORANDOM_NUMBERS_GENERATORS_H
#define LAB1_PSEUDORANDOM_NUMBERS_GENERATORS_H

//CALL FUNCTIONS THAT GENERATE TABLE

void call_table(const std::function<double()>& rand_num_generator, double interval_lower_bound, double interval_upper_bound);
void call_table(const std::function< std::vector<double> () >& multiple_arg_rand_num_gen, double interval_lower_bound, double interval_upper_bound);


//CLASS FOR SAVING CONSTANTS
class Rand_Gen{
private:

    double interval_start{};
    double interval_end{};


public:

    Rand_Gen(double interval_lower_bound, double interval_upper_bound);

    double Get_interval_start() const;
    double Get_interval_end() const;

    void Set_Parameters(double interval_start, double interval_end);

    double start_of_interval;
    double end_of_interval;
    double len_of_interval;
    double end_of_sub_interval;

    static const int num_of_random_numbers = 10000;
    static const int num_of_intervals = 10;

    int distribution_by_intervals [ num_of_intervals ] {};

};

//METHODS FOR GENERATING DISTRIBUTION TABLE

void generate_distribution_info(const std::function<double()>& rand_num_generator, Rand_Gen Rand);

void generate_distribution_info(const std::function< std::vector<double> () >& multiple_arg_rand_num_gen, Rand_Gen Rand);

int find_index(double rand_value, Rand_Gen Rand, double k);

void generate_distribution_table(Rand_Gen Rand);



//METHODS OF GENERATING PSEUDORANDOM NUMBERS

double linear_congruence_method();

double quadratic_congruence_method();

double fibonacci_numbers();

int modInverse(int a, int modulo);
int gcdExtended(int a, int modulo, int* x, int* y);
double inversive_congruential_generator();

double union_method();

double three_sigma_rule();

std::vector<double> polar_coordinates_method();

double ratio_method();

double log_method();

double arends_method();

#endif //LAB1_PSEUDORANDOM_NUMBERS_GENERATORS_H
