/* ======================================================================
                            INEM Staff Scheduling
   ====================================================================== */

/*
*	Code for the algorithms belonging to the paper
*	    Vermuyten, H., Namorado Rosa, J., Marques, I., Beliën, J. Barbosa-Póvoa, A. (2018).
*       Integrated staff scheduling at a medical emergency service: An optimisation approach,
*       Expert Systems With Applications, 112, 62-76.
*       https://doi.org/10.1016/j.eswa.2018.06.017
*
*
*	Code author: Hendrik Vermuyten
*/

#include "Problem_Data.h"
#include "Diving_Heuristic.h"
#include "Standard_IP_Algorithm.h"
#include "VNDS.h"
#include "Instance_Generator.h"
#include <iostream>


int main()
{
	// User choice: instance generator vs algorithms
	int option;
	std::cout << "\nTo run the instance generator, enter 1, \n to run one of the algorithms, enter 2: ";
	std::cin >> option;
		
	if (option == 1) // generate dataset
	{
		std::string file_name;
		std::cout << "\nGive the file name: ";
		std::cin >> file_name;
		Instance_Generator::generate_dataset(file_name);
		std::cout << "The dataset was successfully generated!";
	}

	else // run one of the algorithms
	{

		// Input problem data
		{
			std::string input_file_name;
			while (true) {
				std::cout << "\nGive the file name: ";
				std::cin >> input_file_name;
				if (input(input_file_name))
					break;
			}
		}


		// User choice algorithm parameters
		{
			int choice;
			std::cout << "\nDo you want to change the objective function weights? 1 = yes, 0 = no: ";
			std::cin >> choice;
			if (choice == 1)
			{
				std::cout << "Give the objective weight for understaffing in CODU: ";
				std::cin >> obj_weight_Y_REmin_CODU;
				std::cout << "Give the objective weight for understaffing in ambulances: ";
				std::cin >> obj_weight_Y_REmin_amb;
				std::cout << "Give the objective weight for overstaffing: ";
				std::cin >> obj_weight_Y_REplus;
				std::cout << "Give the objective weight for excess hours worked: ";
				std::cin >> obj_weight_Y_Hplus;
				std::cout << "Give the objective weight for shortage hours worked: ";
				std::cin >> obj_weight_Y_Hmin;
				std::cout << "Give the objective weight for full weekends off: ";
				std::cin >> obj_weight_Y_W;
				std::cout << "Give the objective weight for assigning tasks within the proper group for CODU: ";
				std::cin >> obj_weight_Y_G_CODU;
				std::cout << "Give the objective weight for assigning tasks within the proper group for ambulances: ";
				std::cin >> obj_weight_Y_G_ambulances;
			}
		}


		// Algorithm choice
		{
			int choice;
			std::cout << "\nChoose the algorithm. \n1 = Standard LP, \n2 = Standard IP, \n3 = Diving Heuristic, \n4 = VNDS Heuristic: ";
			std::cin >> choice;

			if (choice == 1)
			{
				IP::algorithm_choice = IP::Algorithm_Type::LP;
				IP::Standard_IP_Algorithm algorithm;
				algorithm.run_algorithm();
			}
			else if (choice == 2)
			{
				std::cout << "\nGive the allowed computation time in seconds: ";
				std::cin >> IP::allowed_computation_time;
				IP::algorithm_choice = IP::Algorithm_Type::IP;
				IP::Standard_IP_Algorithm algorithm;
				algorithm.run_algorithm();
			}
			else if (choice == 3)
			{
				std::cout << "\nChoose the column generation method. \n0 = one column per person, \n1 = one column person p and reoptimize: ";
				std::cin >> Diving_Heuristic::column_generation_method;

				std::cout << "\nChoose the branching method. \n0 = largest fractional variable, \n1 = variables with a value above a certain treshold: ";
				std::cin >> Diving_Heuristic::branching_method_diving;
				if (Diving_Heuristic::branching_method_diving == Diving_Heuristic::Branching_Method_Diving::value_above_threshold) {
					std::cout << "\nThreshold: ";
					std::cin >> Diving_Heuristic::branching_threshold_diving;
				}

				std::cout << "\nGive the total allowed computation time (seconds): ";
				std::cin >> Diving_Heuristic::allowed_computation_time;
				std::cout << "\nAllowed computation time root node (seconds): ";
				std::cin >> Diving_Heuristic::allowed_computation_time_root;

				Diving_Heuristic::Diving_Column_Generation algorithm;
				algorithm.run_algorithm();
			}
			else if (choice == 4)
			{
				std::cout << "\nGive the allowed computation time in seconds: ";
				std::cin >> VNDS::allowed_computation_time;
				VNDS::heuristic();
			}
		}
	}

	

	std::cout << "\n\n\nPress enter to exit the program ... ";
	getchar();
	getchar();

	return 0;
}


