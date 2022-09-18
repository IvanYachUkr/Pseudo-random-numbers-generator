#include<iostream>
#include <set>
#include "Pseudorandom_numbers_generators.h"



int main()


{

    int command;
    std::cout << "Type a number of pseudorandom numbers generator: ";
    std::cin >> command;
    std::cout << std::endl;
    std::cout << "Your number is: " << command << std::endl;

    switch (command) {
        case 1:
            call_table(&linear_congruence_method, 0.0,1.0);
            break;
        case 2:
            call_table(&quadratic_congruence_method, 0.0,1.0);
            break;
        case 3:
            call_table(&fibonacci_numbers, 0.0,1.0);
            break;
        case 4:
            call_table(&inversive_congruential_generator, 0.0,1.0);
            break;
        case 5:
            call_table(&union_method, 0.0,1.0);
            break;
        case 6:
            call_table(&three_sigma_rule, -3.0,3.0);
            break;
        case 7:
            call_table(&polar_coordinates_method, -3.0,3.0);
            break;
        case 8:
            call_table(&ratio_method, -3.0,3.0);
            break;
        case 9:
            call_table(&log_method, 0.0,100.0);
            break;
        case 10:
            call_table(&arends_method, 0.0,100.0);
            break;
        default:
            std::cout << "WRONG INPUT";
    }


    return 0;

}