//
// Created by Dir on 11.09.2022.
//
#include <functional>
#include <iostream>
#include <cmath>
#include "Pseudorandom_numbers_generators.h"

//Rand_Gen Class DEFINITIONS

Rand_Gen::Rand_Gen(double interval_lower_bound, double interval_upper_bound) {
    Set_Parameters(interval_lower_bound, interval_upper_bound);
    start_of_interval = Rand_Gen::Get_interval_start();
    end_of_interval = Rand_Gen::Get_interval_end();
    len_of_interval = (end_of_interval - start_of_interval) / Rand_Gen::num_of_intervals;
    end_of_sub_interval = start_of_interval + len_of_interval;

}

void Rand_Gen::Set_Parameters(double interval_lower_bound, double interval_upper_bound) {
    interval_start = interval_lower_bound;
    interval_end = interval_upper_bound;
}

double Rand_Gen::Get_interval_start() const {
    return interval_start;
}

double Rand_Gen::Get_interval_end() const {
    return interval_end;
}

//CALL FUNCTIONS THAT GENERATE TABLE


void call_table(const std::function<double()>& rand_num_generator, double interval_lower_bound, double interval_upper_bound){
    Rand_Gen Rand(interval_lower_bound, interval_upper_bound);
    generate_distribution_info(rand_num_generator, Rand);
}
void call_table(const std::function< std::vector<double> () >& multiple_arg_rand_num_gen, double interval_lower_bound, double interval_upper_bound){
    Rand_Gen Rand(interval_lower_bound, interval_upper_bound);
    generate_distribution_info(multiple_arg_rand_num_gen, Rand);
}


//FUNCTIONS THAT GENERATE TABLE

void generate_distribution_info(const std::function<double()>& rand_num_generator, Rand_Gen Rand){

    //ratio of the number of selected unit intervals to the number of unit intervals in the selected interval
    double k = (double)Rand_Gen::num_of_intervals / std::abs(Rand.end_of_interval - Rand.start_of_interval);

    for (int i = 0; i < Rand_Gen::num_of_random_numbers; i++) {

        double rand_num = rand_num_generator();
        if (rand_num < Rand.start_of_interval || rand_num > Rand.end_of_interval){
            continue;
        }
        int index = find_index(rand_num, Rand, k); // upper_bound_of_interval

        Rand.distribution_by_intervals[index]++;
    }
    generate_distribution_table(Rand);

}

void generate_distribution_info(const std::function< std::vector<double> () >& multiple_arg_rand_num_gen, Rand_Gen Rand){

    //ratio of the number of selected unit intervals to the number of unit intervals in the selected interval
    double k = (double)Rand_Gen::num_of_intervals / std::abs(Rand.end_of_interval - Rand.start_of_interval);

    for (int i = 0; i < Rand_Gen::num_of_random_numbers; i++) {
        std::vector<double> random_numbers = multiple_arg_rand_num_gen();
        for (double rand_num: random_numbers) {
            if (rand_num < Rand.start_of_interval || rand_num > Rand.end_of_interval){
                continue;
            }
            int index = find_index(rand_num, Rand, k);

            Rand.distribution_by_intervals[index]++;

        }
    }
    generate_distribution_table(Rand);
}

int find_index(double rand_value, Rand_Gen Rand, double k){
    int index;
    if ( (Rand.start_of_interval>=0 && Rand.end_of_interval>0) or (Rand.start_of_interval<0 and Rand.end_of_interval<=0) ){

        index = (int)(rand_value* k)- (int)Rand.start_of_interval;

    }
    else{
        index = (int) std::abs((rand_value-Rand.start_of_interval)*k);
    }
    return index;
}



void generate_distribution_table(Rand_Gen Rand){

//    int *distribution_by_intervals;
//    distribution_by_intervals = generate_distribution_info(linear_congruence_method);

    std::cout << "Interval" << "      " << "Frequency " << std::endl;

    for (int el_in_interval : Rand.distribution_by_intervals) {

        std::cout << "[" << Rand.start_of_interval << ";" << Rand.end_of_sub_interval << "]";
        std::cout << "  " << (double)el_in_interval / Rand_Gen::num_of_random_numbers << std::endl;
        Rand.start_of_interval += Rand.len_of_interval;
        Rand.end_of_sub_interval += Rand.len_of_interval;
    }
}



//PSEUDORANDOM NUMBERS GENERATORS

double linear_congruence_method(){

    static int x = 9;

    int a = 5, b = 7, m = 9001;
    x = (a * x + b) % m;

    double rand_num = (double)x / m;

    return rand_num;

}

double quadratic_congruence_method(){
    static int x = 17;

    int c = 13, a= 5, d = 8, m = 9001;

    x = (d*x*x + a*x + c)%m;

    double rand_num = (double)x / m;

    return rand_num;
}

double fibonacci_numbers(){
    static int x_1 = 1;
    static int x_2 = 1;
    int m = 9001;

    int temp_x = x_2;

    x_2 = (x_1 + x_2)%m;
    x_1 = temp_x;

    double rand_num = (double)x_2 / m;

    return rand_num;
}


int gcdExtended(int a, int modulo, int* x, int* y)
{

    // Base Case
    if (a == 0)
    {
        *x = 0, *y = 1;
        return modulo;
    }

    // To store results of recursive call
//    And x1 and y1 are results for inputs b%a and a
//    (b%a).x1 + a.y1 = gcd(b%a, a)

    int x1, y1;
    int gcd = gcdExtended(modulo % a, a, &x1, &y1);

    // Update x and y using results of recursive
    // call
    *x = y1 - (modulo / a) * x1;
    *y = x1;

    return gcd;
}

int modInverse(int a, int modulo)
{
    //ax + by = gcd(a, b)
    //mod_inverse_ax +  m* y = gcd(a,m) = 1 modulo m => ax = gcd(a, m) modulo mod_inverse_a
    int x, y;
    int gcd = gcdExtended(a, modulo, &x, &y);
    if (gcd != 1) {
        std::cout << "inv d n ex";
    }

    return (x % modulo + modulo) % modulo;
}

// Function for extended Euclidean Algorithm

double inversive_congruential_generator(){
    static int x = 1;
    int a = 5, c = 6;
    int m = 65536; //2^16

    x = (a * modInverse(x, m) + c)%m;

    double rand_num = (double)x / m;

    return rand_num;
}

double union_method(){
    static int x = 9;

    int a_1 = 19, b_1 = 23, m_1 = 9001;
    x = (a_1 * x + b_1) % m_1;

    static int y = 4;

    int a_2 = 5, b_2 = 7, m_2 = 8821;
    y = (a_2 * y + b_2) % m_2;

    int m = 9473;
    int z = (x - y)%m;

    double rand_num = (double)z / m;

    return rand_num;

}

double three_sigma_rule(){

    double sum = 0;
    for (int i = 0; i < 12; ++i) {
        sum += linear_congruence_method();

    }
    double res = sum - 6;
    return res;
}

std::vector<double> polar_coordinates_method(){

    std::vector<double> result {};
    double u_1 = linear_congruence_method();
    double u_2 = linear_congruence_method();
    double v_1 = 2*u_1 - 1;
    double v_2 = 2*u_2 - 1;
    double s = v_1 * v_1 + v_2 * v_2;
    if ( s >= 1 ){
        result = polar_coordinates_method();

    }else{
        double multiplier = sqrt((-2*log(s))/s );
        double x_1 = v_1 * multiplier;
        double x_2 = v_2 * multiplier;
        result = {x_1, x_2};
    }
    return result;


}

double ratio_method(){
    double u = linear_congruence_method();
    double v = linear_congruence_method();
    double x = sqrt((8.0/exp(1)))*(v-0.5)/u;
    if (x*x > -4*log(u)){
       x = ratio_method();
    }

    return x;
}

double log_method(){
    double u = linear_congruence_method();
    double v = 50.0;
    double x = -v* log(u);
    return x;
}

double arends_method(){

    double u = quadratic_congruence_method();
    double y = std::tan(M_PI*u);
    double a =60.0;
    double x = sqrt(2*a-1)*y+a-1;
    if (x<=0){
        x = arends_method();
    }
    double v = linear_congruence_method();
    if ( v > (1 + y*y)* exp( (a-1)*log(x/(a+1))-sqrt(2*a-1)*y ) ){
        x = arends_method();
    }
return x;

}

