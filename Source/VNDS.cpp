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

#include "VNDS.h"
#include "Problem_Data.h"
#include "Solution.h"
#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <memory>


namespace VNDS
{
	// USER CHOICE
	double allowed_computation_time{ 3600 };		// in seconds


	// ALGORITHM SETTINGS
	constexpr bool write_time_evolution_to_file = false;
	constexpr bool write_to_file{ false };			// write IP models to an lp-file


	// STATISTICS
	std::ofstream output_file;
	double elapsed_computation_time{ 0 };			// in seconds
	int nb_shakes{ 0 };								// number of times a shake is performed on the solution
	clock_t start_time;

	struct Performance_Info
	{
		double elapsed_time;
		double objective;
		std::string neighbourhood;
		int nb_shakes;
	};
	std::vector<Performance_Info> time_evolution_best_sol;
	std::vector<Performance_Info> time_evolution_current_sol;


	// RANDOM NUMBER GENERATOR
	std::random_device randdev;
	std::seed_seq seed{ randdev(), randdev(), randdev(), randdev(), randdev(), randdev(), randdev(), randdev() };
	std::mt19937 generator(seed);


	// SOLUTION
	std::unique_ptr<Solution> best_solution;
	std::unique_ptr<Solution> current_solution;





	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// IMPLEMENTATION INITIAL SOLUTION
	
	void Initial_Solution::run_algorithm()
	{
		std::cout << "\n\n\nStarting constructive algorithm ... ";

		// set remaining demands equal to initial demands
		remaining_demands = g_shift_demands;

		// initialize cplex
		initialize_cplex();

		for (int p = 0; p < nb_people; ++p)
		{
			// build and solve subproblem
			build_problem(p);
			solve_problem(p);

			// save partial schedule and update remaining demands
			for (int t = 0; t < nb_tasks; ++t) {
				for (int d = 0; d < nb_days; ++d) {
					for (int s = 0; s < nb_shifts; ++s) {
						if (solution[t*nb_days*nb_shifts + d * nb_shifts + s] > 0.99) {
							current_solution->set_at(p, t, d, s, true);
							--remaining_demands[t*nb_days*nb_shifts + d * nb_shifts + s];
						}
					}
				}
			}
		}

		// free cplex memory
		clear_cplex();
	}

	void Initial_Solution::initialize_cplex()
	{
		// Open the cplex environment
		env = CPXopenCPLEX(&status);

		// Set the output to screen on/off
		status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	}

	void Initial_Solution::build_problem(int person)
	{
		// Local variables
		double obj[1];					// Objective function
		double lb[1];					// Lower bound variables
		double ub[1];					// Upper bound variables
		double rhs[1];					// Right-hand side constraints
		char *colname[1];				// Variable names
		char *rowname[1];				// Constraint names
		char sense[1];					// Sign of constraint
		int matbeg[1];					// Begin position of the constraint
		int matind[20000];				// Position of each element in constraint matrix
		double matval[20000];			// Value of each element in constraint matrix
		char type[1];					// Type of variable (integer, binary, fractional)
		int f{ 0 };						// To calculate number of nonzero coefficients in each constraint

		// Delete the old problem
		status = CPXfreeprob(env, &problem);

		// Create the problem
		{
			std::string name = "INEM_CODU_subproblem_person_" + std::to_string(person + 1);
			problem = CPXcreateprob(env, &status, name.c_str());
		}

		// Problem is minimization
		CPXchgobjsen(env, problem, CPX_MIN);



		// VARIABLES
		// Add the a_tds variables
		for (int t = 0; t < nb_tasks; ++t)
		{
			for (int d = 0; d < nb_days; ++d)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					std::string name = "a_" + std::to_string(t + 1) + "_" + std::to_string(d + 1) + "_" + std::to_string(s + 1);
					colname[0] = const_cast<char*>(name.c_str());

					if (remaining_demands[t*nb_days*nb_shifts + d * nb_shifts + s] >= 1) {
						if (t < nb_tasks_CODU)
							obj[0] = -obj_weight_Y_REmin_CODU;		// objective function coefficient
						else
							obj[0] = -obj_weight_Y_REmin_amb;
					}
					else
						obj[0] = 0;

					lb[0] = 0;		// binary variable
					ub[0] = 1;
					type[0] = 'B';

					status = CPXnewcols(env, problem, 1, obj, lb, ub, type, colname); // Generate columns (the variables) and subsequently add rows (constraints)
				}
			}
		}


		// Add the Y_Wplus_w variables
		for (int w = 0; w < nb_weekends; ++w)
		{
			std::string name = "Y_Wplus_" + std::to_string(w + 1);
			colname[0] = const_cast<char*>(name.c_str());

			obj[0] = obj_weight_Y_W;	// objective function coefficient

			lb[0] = 0;		// integer variable
			type[0] = 'I';

			status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
		}


		// Add the Y_Wmin_w variables
		for (int w = 0; w < nb_weekends; ++w)
		{
			std::string name = "Y_Wmin_" + std::to_string(w + 1);
			colname[0] = const_cast<char*>(name.c_str());

			obj[0] = obj_weight_Y_W;	// objective function coefficient

			lb[0] = 0;		// integer variable
			type[0] = 'I';

			status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
		}


		// Add the Y_Hplus variable
		{
			std::string name = "Y_Hplus";
			colname[0] = const_cast<char*>(name.c_str());
		}

		obj[0] = obj_weight_Y_Hplus;	// objective function coefficient

		lb[0] = 0;		// integer variable
		type[0] = 'I';

		status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);


		// Add the Y_Hmin variable
		{
			std::string name = "Y_Hmin";
			colname[0] = const_cast<char*>(name.c_str());
		}

		obj[0] = obj_weight_Y_Hmin;	// objective function coefficient

		lb[0] = 0;		// integer variable
		type[0] = 'I';

		status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);



		// Add the Y_G_g variables
		for (int g = 0; g < nb_groups; ++g)
		{
			std::string name = "Y_G_" + std::to_string(g + 1);
			colname[0] = const_cast<char*>(name.c_str());

			if (g < nb_groups_CODU)
				obj[0] = obj_weight_Y_G_CODU;	// objective function coefficient
			else
				obj[0] = obj_weight_Y_G_ambulances;

			lb[0] = 0;		// integer variable
			type[0] = 'I';

			status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
		}



		// CONSTRAINTS
		// Constraint set 2: only one shift of three consecutive shifts (11h rest)
		// (2.1) If night (0h-8h), then not morning or afternoon same day
		for (int d = 0; d < nb_days; ++d)
		{
			sense[0] = 'L';
			rhs[0] = 1;

			std::string name = "Min_11_hours_between_night_shift_and_next_shift_on_day_" + std::to_string(d + 1);
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				// night shift day d
				matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::night;
				matval[f] = 1;
				++f;

				// morning shift day d
				matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::morning;
				matval[f] = 1;
				++f;

				// afternoon shift day d
				matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::afternoon;
				matval[f] = 1;
				++f;
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// (2.2) If morning (8h-16h), then not afternoon same day or night next day
		for (int d = 0; d < nb_days - 1; ++d)
		{
			sense[0] = 'L';
			rhs[0] = 1;

			std::string name = "Min_11_hours_between_morning_shift_and_next_shift_on_day_" + std::to_string(d + 1);
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				// morning shift day d
				matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::morning;
				matval[f] = 1;
				++f;

				// afternoon shift day d
				matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::afternoon;
				matval[f] = 1;
				++f;

				// night shift day d+1
				matind[f] = t * nb_days*nb_shifts + (d + 1)*nb_shifts + Shift::night;
				matval[f] = 1;
				++f;
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// (2.3) If afternoon (16h-24h), then not night or morning next day
		for (int d = 0; d < nb_days - 1; ++d)
		{
			sense[0] = 'L';
			rhs[0] = 1;

			std::string name = "Min_11_hours_between_afternoon_shift_and_next_shift_on_day_" + std::to_string(d + 1);
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				// afternoon shift day d
				matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::afternoon;
				matval[f] = 1;
				++f;

				// night shift day d+1
				matind[f] = t * nb_days*nb_shifts + (d + 1)*nb_shifts + Shift::night;
				matval[f] = 1;
				++f;

				// morning shift day d+1
				matind[f] = t * nb_days*nb_shifts + (d + 1)*nb_shifts + Shift::morning;
				matval[f] = 1;
				++f;
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// Constraint set 3: Forbidden to assign a person p to a task t that the person cannot perform
		for (int t = 0; t < nb_tasks; ++t)
		{
			if (!person_task(person, t))
			{
				for (int d = 0; d < nb_days; ++d)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						sense[0] = 'E';
						rhs[0] = 0;

						std::string name = "person_" + std::to_string(person + 1) + "_cannot_do_task_" + std::to_string(t + 1) + "_day_" + std::to_string(d + 1) + "_shift_" + std::to_string(s + 1);
						rowname[0] = const_cast<char*>(name.c_str());

						matbeg[0] = 0;
						f = 0;

						matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
						matval[f] = 1;
						++f;

						status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
					}
				}
			}
		}


		// Constraint set 4: Maximum of 6 consecutive working days 
		for (int r = 0; r < nb_days - 6; ++r)
		{
			sense[0] = 'L';
			rhs[0] = 6;

			std::string name = "maximum_6_consecutive_working_days_from_day_" + std::to_string(r + 1);
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				for (int d = r; d <= r + 6; ++d)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
						matval[f] = 1;
						++f;
					}
				}
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// Constraint set 5: Maximum of 5 consecutive days-off 
		for (int r = 0; r < nb_days - 5; ++r)
		{
			sense[0] = 'G';
			rhs[0] = 1;

			std::string name = "maximum_5_consecutive_days_off_from_day_" + std::to_string(r + 1);
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				for (int d = r; d <= r + 5; ++d)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
						matval[f] = 1;
						++f;
					}
				}
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}

		// Constraint set 6: At least 1 Sunday off in every 4 Sundays 
		sense[0] = 'L';
		rhs[0] = 3 * nb_weekends / 4;
		{
			std::string name = "at_least_one_Sunday_off";
			rowname[0] = const_cast<char*>(name.c_str());
		}
		matbeg[0] = 0;
		f = 0;

		for (int t = 0; t < nb_tasks; ++t)
		{
			// depending on start day of planning horizon
			for (int d = 6 - start_day; d < nb_days; d = d + 7)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
					matval[f] = 1;
					++f;
				}
			}
		}

		status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);


		// Constraint set 7: Ideally Saturday and Sunday off 
		// depending on start day of planning horizon
		int w = 0;
		for (int d = 6 - start_day; d < nb_days; d = d + 7)
		{
			if (d > 0) // if first day is Sunday, this weekend is not included
			{
				sense[0] = 'E';
				rhs[0] = 0;

				std::string name = "whole_weekend_week" + std::to_string(d / 7 + 1);
				rowname[0] = const_cast<char*>(name.c_str());

				matbeg[0] = 0;
				f = 0;

				for (int t = 0; t < nb_tasks; ++t)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						// Sunday
						matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
						matval[f] = 1;
						++f;

						// Saturday
						matind[f] = t * nb_days*nb_shifts + (d - 1)*nb_shifts + s;
						matval[f] = -1;
						++f;
					}
				}

				// Y_Wplus_w
				matind[f] = nb_tasks * nb_days*nb_shifts + w;
				matval[f] = -1;
				++f;

				// Y_Wmin_w
				matind[f] = nb_tasks * nb_days*nb_shifts + nb_weekends + w;
				matval[f] = 1;
				++f;

				status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);

				++w;
			}
		}


		// Constraint set 8: 35 weekly working hours (140 per 28 days)
		sense[0] = 'E';
		rhs[0] = (int)(((double)140 / 28 * nb_days) + 0.5) - 7 * nb_holidays; // depending on the length of the planning horizon and the number of holidays
		{
			std::string name = "140_working_hours";
			rowname[0] = const_cast<char*>(name.c_str());
		}
		matbeg[0] = 0;
		f = 0;

		for (int t = 0; t < nb_tasks; ++t)
		{
			for (int d = 0; d < nb_days; ++d)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
					matval[f] = g_task_durations[t];
					++f;
				}
			}
		}

		// Y_Hplus
		matind[f] = nb_tasks * nb_days*nb_shifts + 2 * nb_weekends;
		matval[f] = -1;
		++f;

		// Y_Hmin
		matind[f] = nb_tasks * nb_days*nb_shifts + 2 * nb_weekends + 1;
		matval[f] = 1;
		++f;

		status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);


		// Constraint set 9: Assign a task t from a group g to a person p from that group 
		for (int g = 0; g < nb_groups; ++g)
		{
			sense[0] = 'E';
			rhs[0] = 0;

			std::string name = "group_" + std::to_string(g + 1) + "_assign_tasks_to_members_of_this_group";
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			if (person_group(person, g))
			{
				for (int t = 0; t < nb_tasks; ++t)
				{
					if (!group_task(g, t))
					{
						for (int d = 0; d < nb_days; ++d)
						{
							for (int s = 0; s < nb_shifts; ++s)
							{
								matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
								matval[f] = 1;
								++f;
							}
						}
					}
				}
			}

			// Y_G
			matind[f] = nb_tasks * nb_days*nb_shifts + 2 * nb_weekends + 2 + g;
			matval[f] = -1;
			++f;

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// (10.1) night shifts
		sense[0] = 'G';
		rhs[0] = 2;
		{
			std::string name = "Person_" + std::to_string(person + 1) + "_at_least_two_night_shifts";
			rowname[0] = const_cast<char*>(name.c_str());
		}
		matbeg[0] = 0;
		f = 0;

		for (int t = 0; t < nb_tasks; ++t)
		{
			if (person_task(person, t))
			{
				for (int d = 0; d < nb_days; ++d)
				{
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::night;
					matval[f] = 1;
					++f;
				}
			}
		}

		status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);


		// (10.2) morning shifts
		sense[0] = 'G';
		rhs[0] = 2;
		{
			std::string name = "Person_" + std::to_string(person + 1) + "_at_least_two_morning_shifts";
			rowname[0] = const_cast<char*>(name.c_str());
		}
		matbeg[0] = 0;
		f = 0;

		for (int t = 0; t < nb_tasks; ++t)
		{
			if (person_task(person, t))
			{
				for (int d = 0; d < nb_days; ++d)
				{
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::morning;
					matval[f] = 1;
					++f;
				}
			}
		}

		status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);


		// (10.3) afternoon shifts
		sense[0] = 'G';
		rhs[0] = 2;
		{
			std::string name = "Person_" + std::to_string(person + 1) + "_at_least_two_afternoon_shifts";
			rowname[0] = const_cast<char*>(name.c_str());
		}
		matbeg[0] = 0;
		f = 0;

		for (int t = 0; t < nb_tasks; ++t)
		{
			if (person_task(person, t))
			{
				for (int d = 0; d < nb_days; ++d)
				{
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::afternoon;
					matval[f] = 1;
					++f;
				}
			}
		}

		status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);


		// Write to the problem to a file
		if (write_to_file && person == 0) {
			status = CPXwriteprob(env, problem, "INEM_subproblem.lp", NULL);
		}
	}

	void Initial_Solution::solve_problem(int person)
	{
		// Allocate memory for solution
		int numcols = CPXgetnumcols(env, problem);	// how many columns (= variables) are there in the problem
		solution = new double[numcols];				// allocate memory/create array for solution


													// Try to optimize and check if the model works
		try
		{
			// Optimize the problem
			std::cout << "\n\nCPLEX is solving the subproblem for person " << person + 1 << " ... ";
			status = CPXmipopt(env, problem);

			// Get the solution
			status = CPXsolution(env, problem, &solstat, &objective, solution, NULL, NULL, NULL);

			if (status != 0)
				throw status;		// if something went wrong, throw exception

			if (!(solstat == 101 || solstat == 102))
				throw solstat;		// if optimal solution not found, throw exception

									// Print solution (only executed when no exception thrown)
			std::cout << "\nObjective value: " << objective;
		} // end try

		catch (int exception)
		{
			// To screen
			std::cout << "\n\n\nCplex didn't find the optimal solution.\n";

			// Give the error reason
			if (exception == 1217) std::cout << "\nError 1217: No solution exists. \nThe requested command cannot be executed because no solution exists for the problem. \nOptimize the problem first.\n\n\n";
			else if (exception == 118) std::cout << "\nProblem is unbounded.\n\n\n";
			else if (exception == 103) std::cout << "\nProblem is infeasible.\n\n\n";
			else if (exception == 119) std::cout << "\nProblem is unbounded or infeasible.\n\n\n";
			else if (exception == 115) std::cout << "\nProblem optimal with unscaled infeasibilities.\n\n\n";
			else if (exception == 107) std::cout << "\nTime limit exceeded, integer solution exists.\n\n\n";
			else if (exception == 108) std::cout << "\nTime limit exceeded, no integer solution.\n\n\n";
			else if (exception == 111) std::cout << "\nTreememory limit, integer solution exists.\n\n\n";
			else if (exception == 112) std::cout << "\nTreememory limit, no integer solution exists.\n\n\n";
			else std::cout << "\nOther reason for termination.\n\n\n";
		} // end catch
	}

	void Initial_Solution::clear_cplex()
	{
		// Free solution 
		delete[] solution;

		// Free the problem
		status = CPXfreeprob(env, &problem);

		// Close the cplex environment
		status = CPXcloseCPLEX(&env);
	}





	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// IMPLEMENTATION COLUMN PERSON

	void Column_Person::find_column(int person)
	{
		std::cout << "\n\n\nStarting shaking phase for person " << person + 1 << " ... ";

		// initialize cplex
		initialize_cplex();

		// build and solve subproblem
		build_problem(person);
		solve_problem(person);

		// save partial schedule 
		current_solution->reset_person(person);
		for (int t = 0; t < nb_tasks; ++t) {
			for (int d = 0; d < nb_days; ++d) {
				for (int s = 0; s < nb_shifts; ++s) {
					if (solution[t*nb_days*nb_shifts + d * nb_shifts + s] > 0.99) {
						current_solution->set_at(person, t, d, s, true);
					}
				}
			}
		}

		// free cplex memory
		clear_cplex();
	}

	void Column_Person::initialize_cplex()
	{
		// Open the cplex environment
		env = CPXopenCPLEX(&status);

		// Set the output to screen on/off
		status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	}

	void Column_Person::build_problem(int person)
	{
		// Local variables
		double obj[1];					// Objective function
		double lb[1];					// Lower bound variables
		double ub[1];					// Upper bound variables
		double rhs[1];					// Right-hand side constraints
		char *colname[1];				// Variable names
		char *rowname[1];				// Constraint names
		char sense[1];					// Sign of constraint
		int matbeg[1];					// Begin position of the constraint
		int matind[20000];				// Position of each element in constraint matrix
		double matval[20000];			// Value of each element in constraint matrix
		char type[1];					// Type of variable (integer, binary, fractional)
		int f{ 0 };						// To calculate number of nonzero coefficients in each constraint

		// Delete the old problem
		status = CPXfreeprob(env, &problem);

		// Create the problem
		{
			std::string name = "INEM_CODU_subproblem_person_" + std::to_string(person + 1);
			problem = CPXcreateprob(env, &status, name.c_str());
		}

		// Problem is minimization
		CPXchgobjsen(env, problem, CPX_MIN);



		// VARIABLES
		// Add the a_tds variables
		std::uniform_int_distribution<int> distribution(1, 10);
		for (int t = 0; t < nb_tasks; ++t)
		{
			for (int d = 0; d < nb_days; ++d)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					std::string name = "a_" + std::to_string(t + 1) + "_" + std::to_string(d + 1) + "_" + std::to_string(s + 1);
					colname[0] = const_cast<char*>(name.c_str());

					obj[0] = distribution(generator) * (-100); 	// objective function coefficient

					lb[0] = 0;		// binary variable
					ub[0] = 1;
					type[0] = 'B';

					status = CPXnewcols(env, problem, 1, obj, lb, ub, type, colname); // Generate columns (the variables) and subsequently add rows (constraints)
				}
			}
		}


		// Add the Y_Wplus_w variables
		for (int w = 0; w < nb_weekends; ++w)
		{
			std::string name = "Y_Wplus_" + std::to_string(w + 1);
			colname[0] = const_cast<char*>(name.c_str());

			obj[0] = obj_weight_Y_W;	// objective function coefficient

			lb[0] = 0;		// integer variable
			type[0] = 'I';

			status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
		}


		// Add the Y_Wmin_w variables
		for (int w = 0; w < nb_weekends; ++w)
		{
			std::string name = "Y_Wmin_" + std::to_string(w + 1);
			colname[0] = const_cast<char*>(name.c_str());

			obj[0] = obj_weight_Y_W;	// objective function coefficient

			lb[0] = 0;		// integer variable
			type[0] = 'I';

			status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
		}


		// Add the Y_Hplus variable
		{
			std::string name = "Y_Hplus";
			colname[0] = const_cast<char*>(name.c_str());
		}

		obj[0] = obj_weight_Y_Hplus;	// objective function coefficient

		lb[0] = 0;		// integer variable
		type[0] = 'I';

		status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);


		// Add the Y_Hmin variable
		{
			std::string name = "Y_Hmin";
			colname[0] = const_cast<char*>(name.c_str());
		}

		obj[0] = obj_weight_Y_Hmin;	// objective function coefficient

		lb[0] = 0;		// integer variable
		type[0] = 'I';

		status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);



		// Add the Y_G_g variables
		for (int g = 0; g < nb_groups; ++g)
		{
			std::string name = "Y_G_" + std::to_string(g + 1);
			colname[0] = const_cast<char*>(name.c_str());

			if (g < nb_groups_CODU)
				obj[0] = obj_weight_Y_G_CODU;	// objective function coefficient
			else
				obj[0] = obj_weight_Y_G_ambulances;

			lb[0] = 0;		// integer variable
			type[0] = 'I';

			status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
		}



		// CONSTRAINTS
		// Constraint set 2: only one shift of three consecutive shifts (11h rest)
		// (2.1) If night (0h-8h), then not morning or afternoon same day
		for (int d = 0; d < nb_days; ++d)
		{
			sense[0] = 'L';
			rhs[0] = 1;

			std::string name = "Min_11_hours_between_night_shift_and_next_shift_on_day_" + std::to_string(d + 1);
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				// night shift day d
				matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::night;
				matval[f] = 1;
				++f;

				// morning shift day d
				matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::morning;
				matval[f] = 1;
				++f;

				// afternoon shift day d
				matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::afternoon;
				matval[f] = 1;
				++f;
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// (2.2) If morning (8h-16h), then not afternoon same day or night next day
		for (int d = 0; d < nb_days - 1; ++d)
		{
			sense[0] = 'L';
			rhs[0] = 1;

			std::string name = "Min_11_hours_between_morning_shift_and_next_shift_on_day_" + std::to_string(d + 1);
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				// morning shift day d
				matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::morning;
				matval[f] = 1;
				++f;

				// afternoon shift day d
				matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::afternoon;
				matval[f] = 1;
				++f;

				// night shift day d+1
				matind[f] = t * nb_days*nb_shifts + (d + 1)*nb_shifts + Shift::night;
				matval[f] = 1;
				++f;
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// (2.3) If afternoon (16h-24h), then not night or morning next day
		for (int d = 0; d < nb_days - 1; ++d)
		{
			sense[0] = 'L';
			rhs[0] = 1;

			std::string name = "Min_11_hours_between_afternoon_shift_and_next_shift_on_day_" + std::to_string(d + 1);
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				// afternoon shift day d
				matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::afternoon;
				matval[f] = 1;
				++f;

				// night shift day d+1
				matind[f] = t * nb_days*nb_shifts + (d + 1)*nb_shifts + Shift::night;
				matval[f] = 1;
				++f;

				// morning shift day d+1
				matind[f] = t * nb_days*nb_shifts + (d + 1)*nb_shifts + Shift::morning;
				matval[f] = 1;
				++f;
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// Constraint set 3: Forbidden to assign a person p to a task t that the person cannot perform
		for (int t = 0; t < nb_tasks; ++t)
		{
			if (!person_task(person, t))
			{
				for (int d = 0; d < nb_days; ++d)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						sense[0] = 'E';
						rhs[0] = 0;

						std::string name = "person_" + std::to_string(person + 1) + "_cannot_do_task_" + std::to_string(t + 1) + "_day_" + std::to_string(d + 1) + "_shift_" + std::to_string(s + 1);
						rowname[0] = const_cast<char*>(name.c_str());

						matbeg[0] = 0;
						f = 0;

						matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
						matval[f] = 1;
						++f;

						status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
					}
				}
			}
		}


		// Constraint set 4: Maximum of 6 consecutive working days 
		for (int r = 0; r < nb_days - 6; ++r)
		{
			sense[0] = 'L';
			rhs[0] = 6;

			std::string name = "maximum_6_consecutive_working_days_from_day_" + std::to_string(r + 1);
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				for (int d = r; d <= r + 6; ++d)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
						matval[f] = 1;
						++f;
					}
				}
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// Constraint set 5: Maximum of 5 consecutive days-off 
		for (int r = 0; r < nb_days - 5; ++r)
		{
			sense[0] = 'G';
			rhs[0] = 1;

			std::string name = "maximum_5_consecutive_days_off_from_day_" + std::to_string(r + 1);
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				for (int d = r; d <= r + 5; ++d)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
						matval[f] = 1;
						++f;
					}
				}
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}

		// Constraint set 6: At least 1 Sunday off in every 4 Sundays 
		sense[0] = 'L';
		rhs[0] = 3 * nb_weekends / 4;
		{
			std::string name = "at_least_one_Sunday_off";
			rowname[0] = const_cast<char*>(name.c_str());
		}
		matbeg[0] = 0;
		f = 0;

		for (int t = 0; t < nb_tasks; ++t)
		{
			// depending on start day of planning horizon
			for (int d = 6 - start_day; d < nb_days; d = d + 7)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
					matval[f] = 1;
					++f;
				}
			}
		}

		status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);


		// Constraint set 7: Ideally Saturday and Sunday off 
		// depending on start day of planning horizon
		int w = 0;
		for (int d = 6 - start_day; d < nb_days; d = d + 7)
		{
			if (d > 0) // if first day is Sunday, this weekend is not included
			{
				sense[0] = 'E';
				rhs[0] = 0;

				std::string name = "whole_weekend_week" + std::to_string(d / 7 + 1);
				rowname[0] = const_cast<char*>(name.c_str());

				matbeg[0] = 0;
				f = 0;

				for (int t = 0; t < nb_tasks; ++t)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						// Sunday
						matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
						matval[f] = 1;
						++f;

						// Saturday
						matind[f] = t * nb_days*nb_shifts + (d - 1)*nb_shifts + s;
						matval[f] = -1;
						++f;
					}
				}

				// Y_Wplus_w
				matind[f] = nb_tasks * nb_days*nb_shifts + w;
				matval[f] = -1;
				++f;

				// Y_Wmin_w
				matind[f] = nb_tasks * nb_days*nb_shifts + nb_weekends + w;
				matval[f] = 1;
				++f;

				status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);

				++w;
			}
		}


		// Constraint set 8: 35 weekly working hours (140 per 28 days)
		sense[0] = 'E';
		rhs[0] = (int)(((double)140 / 28 * nb_days) + 0.5) - 7 * nb_holidays; // depending on the length of the planning horizon and the number of holidays
		{
			std::string name = "140_working_hours";
			rowname[0] = const_cast<char*>(name.c_str());
		}
		matbeg[0] = 0;
		f = 0;

		for (int t = 0; t < nb_tasks; ++t)
		{
			for (int d = 0; d < nb_days; ++d)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
					matval[f] = g_task_durations[t];
					++f;
				}
			}
		}

		// Y_Hplus
		matind[f] = nb_tasks * nb_days*nb_shifts + 2 * nb_weekends;
		matval[f] = -1;
		++f;

		// Y_Hmin
		matind[f] = nb_tasks * nb_days*nb_shifts + 2 * nb_weekends + 1;
		matval[f] = 1;
		++f;

		status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);


		// Constraint set 9: Assign a task t from a group g to a person p from that group 
		for (int g = 0; g < nb_groups; ++g)
		{
			sense[0] = 'E';
			rhs[0] = 0;

			std::string name = "group_" + std::to_string(g + 1) + "_assign_tasks_to_members_of_this_group";
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			if (person_group(person, g))
			{
				for (int t = 0; t < nb_tasks; ++t)
				{
					if (!group_task(g, t))
					{
						for (int d = 0; d < nb_days; ++d)
						{
							for (int s = 0; s < nb_shifts; ++s)
							{
								matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
								matval[f] = 1;
								++f;
							}
						}
					}
				}
			}

			// Y_G
			matind[f] = nb_tasks * nb_days*nb_shifts + 2 * nb_weekends + 2 + g;
			matval[f] = -1;
			++f;

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// (10.1) night shifts
		sense[0] = 'G';
		rhs[0] = 2;
		{
			std::string name = "Person_" + std::to_string(person + 1) + "_at_least_two_night_shifts";
			rowname[0] = const_cast<char*>(name.c_str());
		}
		matbeg[0] = 0;
		f = 0;

		for (int t = 0; t < nb_tasks; ++t)
		{
			if (person_task(person, t))
			{
				for (int d = 0; d < nb_days; ++d)
				{
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::night;
					matval[f] = 1;
					++f;
				}
			}
		}

		status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);


		// (10.2) morning shifts
		sense[0] = 'G';
		rhs[0] = 2;
		{
			std::string name = "Person_" + std::to_string(person + 1) + "_at_least_two_morning_shifts";
			rowname[0] = const_cast<char*>(name.c_str());
		}
		matbeg[0] = 0;
		f = 0;

		for (int t = 0; t < nb_tasks; ++t)
		{
			if (person_task(person, t))
			{
				for (int d = 0; d < nb_days; ++d)
				{
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::morning;
					matval[f] = 1;
					++f;
				}
			}
		}

		status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);


		// (10.3) afternoon shifts
		sense[0] = 'G';
		rhs[0] = 2;
		{
			std::string name = "Person_" + std::to_string(person + 1) + "_at_least_two_afternoon_shifts";
			rowname[0] = const_cast<char*>(name.c_str());
		}
		matbeg[0] = 0;
		f = 0;

		for (int t = 0; t < nb_tasks; ++t)
		{
			if (person_task(person, t))
			{
				for (int d = 0; d < nb_days; ++d)
				{
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + Shift::afternoon;
					matval[f] = 1;
					++f;
				}
			}
		}

		status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);

	}

	void Column_Person::solve_problem(int person)
	{
		// Allocate memory for solution
		int numcols = CPXgetnumcols(env, problem);	// how many columns (= variables) are there in the problem
		solution = new double[numcols];				// allocate memory/create array for solution

													// Try to optimize and check if the model works
		try
		{
			// Optimize the problem
			std::cout << "\n\nCPLEX is solving the subproblem for person " << person + 1 << " ... ";
			status = CPXmipopt(env, problem);

			// Get the solution
			status = CPXsolution(env, problem, &solstat, &objective, solution, NULL, NULL, NULL);

			if (status != 0)
				throw status;		// if something went wrong, throw exception

			if (!(solstat == 101 || solstat == 102))
				throw solstat;		// if optimal solution not found, throw exception

									// Print solution (only executed when no exception thrown)
			std::cout << "\nObjective value: " << objective;
		} // end try

		catch (int exception)
		{
			// To screen
			std::cout << "\n\n\nCplex didn't find the optimal solution.\n";

			// Give the error reason
			if (exception == 1217) std::cout << "\nError 1217: No solution exists. \nThe requested command cannot be executed because no solution exists for the problem. \nOptimize the problem first.\n\n\n";
			else if (exception == 118) std::cout << "\nProblem is unbounded.\n\n\n";
			else if (exception == 103) std::cout << "\nProblem is infeasible.\n\n\n";
			else if (exception == 119) std::cout << "\nProblem is unbounded or infeasible.\n\n\n";
			else if (exception == 115) std::cout << "\nProblem optimal with unscaled infeasibilities.\n\n\n";
			else if (exception == 107) std::cout << "\nTime limit exceeded, integer solution exists.\n\n\n";
			else if (exception == 108) std::cout << "\nTime limit exceeded, no integer solution.\n\n\n";
			else if (exception == 111) std::cout << "\nTreememory limit, integer solution exists.\n\n\n";
			else if (exception == 112) std::cout << "\nTreememory limit, no integer solution exists.\n\n\n";
			else std::cout << "\nOther reason for termination.\n\n\n";
		} // end catch
	}

	void Column_Person::clear_cplex()
	{
		// Free solution 
		delete[] solution;

		// Free the problem
		status = CPXfreeprob(env, &problem);

		// Close the cplex environment
		status = CPXcloseCPLEX(&env);
	}





	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// IMPLEMENTATION MIP_HEURISTIC

	void MIP_Heuristic::run_algorithm()
	{
		// initialize cplex
		initialize_cplex();

		// build problem
		build_problem();

		// set initial solution from constructive method
		set_initial_solution();

		// set time limit subproblems
		status = CPXsetdblparam(env, CPXPARAM_TimeLimit, time_limit_subproblem);

		// improvement heuristic
		std::cout << "\n\n\n\nStarting improvement heuristic ... ";
		probabilistic_VNDS();

		// clear memory cplex
		clear_cplex();
	}

	void MIP_Heuristic::initialize_cplex()
	{
		// Open the cplex environment
		env = CPXopenCPLEX(&status);

		// Set the output to screen on/off
		status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	}

	void MIP_Heuristic::build_problem()
	{
		// Local variables
		double obj[1];					// Objective function
		double lb[1];					// Lower bound variables
		double ub[1];					// Upper bound variables
		double rhs[1];					// Right-hand side constraints
		char *colname[1];				// Variable names
		char *rowname[1];				// Constraint names
		char sense[1];					// Sign of constraint
		int matbeg[1];					// Begin position of the constraint
		int* matind;					// Position of each element in constraint matrix
		double* matval;					// Value of each element in constraint matrix
		char type[1];					// Type of variable (integer, binary, fractional)
		int f{ 0 };						// To calculate number of nonzero coefficients in each constraint

		// Allocate memory
		matind = new int[nb_people*nb_tasks*nb_days*nb_shifts];
		matval = new double[nb_people*nb_tasks*nb_days*nb_shifts];

		// Create the problem
		problem = CPXcreateprob(env, &status, "INEM_CODU_scheduling_problem");

		// Problem is minimization
		CPXchgobjsen(env, problem, CPX_MIN);


		// VARIABLES
		// Add the X_ptds variables
		for (int p = 0; p < nb_people; ++p)
		{
			for (int t = 0; t < nb_tasks; ++t)
			{
				for (int d = 0; d < nb_days; ++d)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						std::string name = "x_" + std::to_string(p + 1) + "_" + std::to_string(t + 1) + "_" + std::to_string(d + 1) + "_" + std::to_string(s + 1);
						colname[0] = const_cast<char*>(name.c_str());

						obj[0] = 0;	// objective function coefficient

						lb[0] = 0;		// binary variable
						ub[0] = 1;
						type[0] = 'B';

						status = CPXnewcols(env, problem, 1, obj, lb, ub, type, colname); // Generate columns (the variables) and subsequently add rows (constraints)
					}
				}
			}
		}

		// Add the Y_REplus variables
		for (int t = 0; t < nb_tasks; ++t)
		{
			for (int d = 0; d < nb_days; ++d)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					std::string name = "Y_REplus_" + std::to_string(t + 1) + "_" + std::to_string(d + 1) + "_" + std::to_string(s + 1);
					colname[0] = const_cast<char*>(name.c_str());

					if (t < nb_tasks_CODU)
						obj[0] = obj_weight_Y_REplus;	// objective function coefficient
					else
						obj[0] = 10000;

					lb[0] = 0;		// integer variable
					type[0] = 'I';

					status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
				}
			}
		}

		// Add the Y_REmin variables
		for (int t = 0; t < nb_tasks; ++t)
		{
			for (int d = 0; d < nb_days; ++d)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					std::string name = "Y_REmin_" + std::to_string(t + 1) + "_" + std::to_string(d + 1) + "_" + std::to_string(s + 1);
					colname[0] = const_cast<char*>(name.c_str());

					if (t < nb_tasks_CODU)
						obj[0] = obj_weight_Y_REmin_CODU;	// objective function coefficient
					else
						obj[0] = obj_weight_Y_REmin_amb;

					lb[0] = 0;		// integer variable
					type[0] = 'I';

					status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
				}
			}
		}

		// Add the Y_Wplus_pw variables
		for (int p = 0; p < nb_people; ++p)
		{
			for (int w = 0; w < nb_weekends; ++w)
			{
				std::string name = "Y_Wplus_" + std::to_string(p + 1) + "_" + std::to_string(w + 1);
				colname[0] = const_cast<char*>(name.c_str());

				obj[0] = obj_weight_Y_W;	// objective function coefficient

				lb[0] = 0;		// integer variable
				type[0] = 'I';

				status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
			}
		}

		// Add the Y_Wmin_pw variables
		for (int p = 0; p < nb_people; ++p)
		{
			for (int w = 0; w < nb_weekends; ++w)
			{
				std::string name = "Y_Wmin_" + std::to_string(p + 1) + "_" + std::to_string(w + 1);
				colname[0] = const_cast<char*>(name.c_str());

				obj[0] = obj_weight_Y_W;	// objective function coefficient

				lb[0] = 0;		// integer variable
				type[0] = 'I';

				status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
			}
		}

		// Add the Y_Hplus_p variables
		for (int p = 0; p < nb_people; ++p)
		{
			std::string name = "Y_Hplus_" + std::to_string(p + 1);
			colname[0] = const_cast<char*>(name.c_str());

			obj[0] = obj_weight_Y_Hplus;	// objective function coefficient

			lb[0] = 0;		// integer variable
			type[0] = 'I';

			status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
		}

		// Add the Y_Hmin_p variables
		for (int p = 0; p < nb_people; ++p)
		{
			std::string name = "Y_Hmin_" + std::to_string(p + 1);
			colname[0] = const_cast<char*>(name.c_str());

			obj[0] = obj_weight_Y_Hmin;	// objective function coefficient

			lb[0] = 0;		// integer variable
			type[0] = 'I';

			status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
		}

		// Add the Y_G_g variables
		for (int g = 0; g < nb_groups; ++g)
		{
			std::string name = "Y_G_" + std::to_string(g + 1);
			colname[0] = const_cast<char*>(name.c_str());

			if (g < nb_groups_CODU)
				obj[0] = obj_weight_Y_G_CODU;	// objective function coefficient
			else
				obj[0] = obj_weight_Y_G_ambulances;

			lb[0] = 0;		// integer variable
			type[0] = 'I';

			status = CPXnewcols(env, problem, 1, obj, lb, NULL, type, colname);
		}


		// CONSTRAINTS
		// Constraint set 1: coverage requirement
		// (1.1) CODU
		for (int t = 0; t < nb_tasks_CODU; ++t)
		{
			for (int d = 0; d < nb_days; ++d)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					sense[0] = 'E';
					rhs[0] = shift_demands(t, d, s);

					std::string name = "Coverage_constraint_task_" + std::to_string(t + 1) + "_day_" + std::to_string(d + 1) + "_shift_" + std::to_string(s + 1);
					rowname[0] = const_cast<char*>(name.c_str());

					matbeg[0] = 0;
					f = 0;

					for (int p = 0; p < nb_people; ++p)
					{
						if (person_task(p, t))
						{
							matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
							matval[f] = 1;
							++f;
						}
					}

					// REplus
					matind[f] = nb_people * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
					matval[f] = -1;
					++f;

					// REmin
					matind[f] = nb_people * nb_tasks*nb_days*nb_shifts + nb_tasks * nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
					matval[f] = 1;
					++f;

					status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
				}
			}
		}

		// (1.2) Ambulances
		for (int t = nb_tasks_CODU; t < nb_tasks; ++t)
		{
			for (int d = 0; d < nb_days; ++d)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					sense[0] = 'E';
					rhs[0] = shift_demands(t, d, s);

					std::string name = "Coverage_constraint_task_" + std::to_string(t + 1) + "_day_" + std::to_string(d + 1) + "_shift_" + std::to_string(s + 1);
					rowname[0] = const_cast<char*>(name.c_str());

					matbeg[0] = 0;
					f = 0;

					for (int p = 0; p < nb_people; ++p)
					{
						if (person_task(p, t))
						{
							matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
							matval[f] = 1;
							++f;
						}
					}

					// REplus
					matind[f] = nb_people * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
					matval[f] = -1;
					++f;

					// REmin
					matind[f] = nb_people * nb_tasks*nb_days*nb_shifts + nb_tasks * nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
					matval[f] = 1;
					++f;

					status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
				}
			}
		}

		// Constraint set 2: only one shift of three consecutive shifts (11h rest)
		// (2.1) If night (0h-8h), then not morning or afternoon same day
		for (int p = 0; p < nb_people; ++p)
		{
			for (int d = 0; d < nb_days; ++d)
			{
				sense[0] = 'L';
				rhs[0] = 1;

				std::string name = "Min_11_hours_between_night_shift_and_next_shift_for_person_" + std::to_string(p + 1) + "_on_day_" + std::to_string(d + 1);
				rowname[0] = const_cast<char*>(name.c_str());

				matbeg[0] = 0;
				f = 0;

				for (int t = 0; t < nb_tasks; ++t)
				{
					// night shift day d
					matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + Shift::night;
					matval[f] = 1;
					++f;

					// morning shift day d
					matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + Shift::morning;
					matval[f] = 1;
					++f;

					// afternoon shift day d
					matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + Shift::afternoon;
					matval[f] = 1;
					++f;
				}

				status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
			}
		}

		// (2.2) If morning (8h-16h), then not afternoon same day or night next day
		for (int p = 0; p < nb_people; ++p)
		{
			for (int d = 0; d < nb_days - 1; ++d)
			{
				sense[0] = 'L';
				rhs[0] = 1;

				std::string name = "Min_11_hours_between_morning_shift_and_next_shift_for_person_" + std::to_string(p + 1) + "_on_day_" + std::to_string(d + 1);
				rowname[0] = const_cast<char*>(name.c_str());

				matbeg[0] = 0;
				f = 0;

				for (int t = 0; t < nb_tasks; ++t)
				{
					// morning shift day d
					matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + Shift::morning;
					matval[f] = 1;
					++f;

					// afternoon shift day d
					matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + Shift::afternoon;
					matval[f] = 1;
					++f;

					// night shift day d+1
					matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + (d + 1)*nb_shifts + Shift::night;
					matval[f] = 1;
					++f;
				}

				status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
			}
		}

		// (2.3) If afternoon (16h-24h), then not night or morning next day
		for (int p = 0; p < nb_people; ++p)
		{
			for (int d = 0; d < nb_days - 1; ++d)
			{
				sense[0] = 'L';
				rhs[0] = 1;

				std::string name = "Min_11_hours_between_afternoon_shift_and_next_shift_for_person_" + std::to_string(p + 1) + "_on_day_" + std::to_string(d + 1);
				rowname[0] = const_cast<char*>(name.c_str());

				matbeg[0] = 0;
				f = 0;

				for (int t = 0; t < nb_tasks; ++t)
				{
					// afternoon shift day d
					matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + Shift::afternoon;
					matval[f] = 1;
					++f;

					// night shift day d+1
					matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + (d + 1)*nb_shifts + Shift::night;
					matval[f] = 1;
					++f;

					// morning shift day d+1
					matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + (d + 1)*nb_shifts + Shift::morning;
					matval[f] = 1;
					++f;
				}

				status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
			}
		}

		// Constraint set 3: Forbidden to assign a person p to a task t that the person cannot perform
		for (int p = 0; p < nb_people; ++p)
		{
			for (int t = 0; t < nb_tasks; ++t)
			{
				if (!person_task(p, t))
				{
					for (int d = 0; d < nb_days; ++d)
					{
						for (int s = 0; s < nb_shifts; ++s)
						{
							sense[0] = 'E';
							rhs[0] = 0;

							std::string name = "person_" + std::to_string(p + 1) + "_cannot_do_task_" + std::to_string(t + 1) + "_day_" + std::to_string(d + 1) + "_shift_" + std::to_string(s + 1);
							rowname[0] = const_cast<char*>(name.c_str());

							matbeg[0] = 0;
							f = 0;

							matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
							matval[f] = 1;
							++f;

							status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
						}
					}
				}
			}
		}

		// Constraint set 4: Maximum of 6 consecutive working days 
		for (int p = 0; p < nb_people; ++p)
		{
			for (int r = 0; r < nb_days - 6; ++r)
			{
				sense[0] = 'L';
				rhs[0] = 6;

				std::string name = "person_" + std::to_string(p + 1) + "_maximum_6_consecutive_working_days_from_day_" + std::to_string(r + 1);
				rowname[0] = const_cast<char*>(name.c_str());

				matbeg[0] = 0;
				f = 0;

				for (int t = 0; t < nb_tasks; ++t)
				{
					for (int d = r; d <= r + 6; ++d)
					{
						for (int s = 0; s < nb_shifts; ++s)
						{
							matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
							matval[f] = 1;
							++f;
						}
					}
				}

				status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
			}
		}

		// Constraint set 5: Maximum of 5 consecutive days-off 
		for (int p = 0; p < nb_people; ++p)
		{
			for (int r = 0; r < nb_days - 5; ++r)
			{
				sense[0] = 'G';
				rhs[0] = 1;

				std::string name = "person_" + std::to_string(p + 1) + "_maximum_5_consecutive_days_off_from_day_" + std::to_string(r + 1);
				rowname[0] = const_cast<char*>(name.c_str());

				matbeg[0] = 0;
				f = 0;

				for (int t = 0; t < nb_tasks; ++t)
				{
					for (int d = r; d <= r + 5; ++d)
					{
						for (int s = 0; s < nb_shifts; ++s)
						{
							matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
							matval[f] = 1;
							++f;
						}
					}
				}

				status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
			}
		}

		// Constraint set 6: At least 1 Sunday off in every 4 Sundays 
		for (int p = 0; p < nb_people; ++p)
		{
			sense[0] = 'L';
			rhs[0] = 3 * nb_weekends / 4;

			std::string name = "person_" + std::to_string(p + 1) + "_at_least_one_Sunday_off";
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				// depending on start day of planning horizon
				for (int d = 6 - start_day; d < nb_days; d = d + 7)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
						matval[f] = 1;
						++f;
					}
				}
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}

		// Constraint set 7: Ideally Saturday and Sunday off 
		for (int p = 0; p < nb_people; ++p)
		{
			int w = 0;
			// depending on start day of planning horizon
			for (int d = 6 - start_day; d < nb_days; d = d + 7)
			{
				if (d > 0)
				{
					sense[0] = 'E';
					rhs[0] = 0;

					std::string name = "person_" + std::to_string(p + 1) + "_whole_weekend_week_" + std::to_string(d % 7 + 1);
					rowname[0] = const_cast<char*>(name.c_str());

					matbeg[0] = 0;
					f = 0;

					for (int t = 0; t < nb_tasks; ++t)
					{
						for (int s = 0; s < nb_shifts; ++s)
						{
							// Sunday
							matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
							matval[f] = 1;
							++f;

							// Saturday
							matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + (d - 1)*nb_shifts + s;
							matval[f] = -1;
							++f;
						}
					}

					// Y_Wplus
					matind[f] = nb_people * nb_tasks*nb_days*nb_shifts + 2 * nb_tasks*nb_days*nb_shifts + p * nb_weekends + w;
					matval[f] = -1;
					++f;

					// Y_Wmin
					matind[f] = nb_people * nb_tasks*nb_days*nb_shifts + 2 * nb_tasks*nb_days*nb_shifts + nb_people * nb_weekends + p * nb_weekends + w;
					matval[f] = 1;
					++f;

					status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);

					++w;
				}
			}
		}

		// Constraint set 8: 35 weekly working hours (140 per 28 days)
		for (int p = 0; p < nb_people; ++p)
		{
			sense[0] = 'E';
			rhs[0] = (int)(((double)140 / 28 * nb_days) + 0.5) - 7 * nb_holidays; // depending on the length of the planning horizon and the number of holidays

			std::string name = "person_" + std::to_string(p + 1) + "_140_working_hours";
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				for (int d = 0; d < nb_days; ++d)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
						matval[f] = g_task_durations[t];
						++f;
					}
				}
			}

			// Y_Hplus
			matind[f] = nb_people * nb_tasks*nb_days*nb_shifts + 2 * nb_tasks*nb_days*nb_shifts + 2 * nb_people*nb_weekends + p;
			matval[f] = -1;
			++f;

			// Y_Hmin
			matind[f] = nb_people * nb_tasks*nb_days*nb_shifts + 2 * nb_tasks*nb_days*nb_shifts + 2 * nb_people*nb_weekends + nb_people + p;
			matval[f] = 1;
			++f;

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}

		// Constraint set 9: Assign a task t from a group g to a person p from that group 
		for (int g = 0; g < nb_groups; ++g)
		{
			sense[0] = 'E';
			rhs[0] = 0;

			std::string name = "group_" + std::to_string(g + 1) + "_assign_tasks_to_members_of_this_group";
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int p = 0; p < nb_people; ++p)
			{
				if (person_group(p, g))
				{
					for (int t = 0; t < nb_tasks; ++t)
					{
						if (!group_task(g, t))
						{
							for (int d = 0; d < nb_days; ++d)
							{
								for (int s = 0; s < nb_shifts; ++s)
								{
									matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
									matval[f] = 1;
									++f;
								}
							}
						}
					}
				}
			}

			// Y_G
			matind[f] = nb_people * nb_tasks*nb_days*nb_shifts + 2 * nb_tasks*nb_days*nb_shifts + 2 * nb_people*nb_weekends + 2 * nb_people + g;
			matval[f] = -1;
			++f;

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}

		// Constraint set 10: Every person at least two shifts of every type (morning, afternoon, and night)
		// (10.1) night shifts
		for (int p = 0; p < nb_people; ++p)
		{
			sense[0] = 'G';
			rhs[0] = 2;

			std::string name = "person_" + std::to_string(p + 1) + "_at_least_two_night_shifts";
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				if (person_task(p, t))
				{
					for (int d = 0; d < nb_days; ++d)
					{
						matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + Shift::night;
						matval[f] = 1;
						++f;
					}
				}
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}

		// (10.2) morning shifts
		for (int p = 0; p < nb_people; ++p)
		{
			sense[0] = 'G';
			rhs[0] = 2;

			std::string name = "person_" + std::to_string(p + 1) + "_at_least_two_morning_shifts";
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				if (person_task(p, t))
				{
					for (int d = 0; d < nb_days; ++d)
					{
						matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + Shift::morning;
						matval[f] = 1;
						++f;
					}
				}
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}

		// (10.3) afternoon shifts
		for (int p = 0; p < nb_people; ++p)
		{
			sense[0] = 'G';
			rhs[0] = 2;

			std::string name = "person_" + std::to_string(p + 1) + "_at_least_two_afternoon_shifts";
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			for (int t = 0; t < nb_tasks; ++t)
			{
				if (person_task(p, t))
				{
					for (int d = 0; d < nb_days; ++d)
					{
						matind[f] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + Shift::afternoon;
						matval[f] = 1;
						++f;
					}
				}
			}

			status = CPXaddrows(env, problem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}



		// Write to the problem to a file
		if (write_to_file)
			status = CPXwriteprob(env, problem, "INEM_CODU_scheduling_problem.lp", NULL);


		// Free memory
		delete[] matind;
		delete[] matval;


		int numcols = CPXgetnumcols(env, problem);
		solution = new double[numcols];

		nb_main_constraints = CPXgetnumrows(env, problem);
	}

	void MIP_Heuristic::solve_problem()
	{
		// Try to optimize and check if the model works
		try
		{
			// Optimize the problem
			status = CPXmipopt(env, problem);

			// Get the solution
			status = CPXsolution(env, problem, &solstat, &candidate_objective, solution, NULL, NULL, NULL);

			if (status != 0)
				throw status;		// if something went wrong, throw exception

			if (!(solstat == 101 || solstat == 102))
				throw solstat;		// if optimal solution not found, throw exception
		} // end try

		catch (int exception)
		{
			// To screen
			std::cout << "\n\n\nCplex didn't find the optimal solution.\n";

			// Give the error reason
			if (exception == 1217) std::cout << "\nError 1217: No solution exists. \nThe requested command cannot be executed because no solution exists for the problem. \nOptimize the problem first.\n\n\n";
			else if (exception == 118) std::cout << "\nProblem is unbounded.\n\n\n";
			else if (exception == 103) std::cout << "\nProblem is infeasible.\n\n\n";
			else if (exception == 119) std::cout << "\nProblem is unbounded or infeasible.\n\n\n";
			else if (exception == 115) std::cout << "\nProblem optimal with unscaled infeasibilities.\n\n\n";
			else if (exception == 107) std::cout << "\nTime limit exceeded, integer solution exists.\n\n\n";
			else if (exception == 108) std::cout << "\nTime limit exceeded, no integer solution.\n\n\n";
			else if (exception == 111) std::cout << "\nTreememory limit, integer solution exists.\n\n\n";
			else if (exception == 112) std::cout << "\nTreememory limit, no integer solution exists.\n\n\n";
			else std::cout << "\nOther reason for termination.\n\n\n";
		} // end catch
	}

	void MIP_Heuristic::clear_cplex()
	{
		// Free solution 
		delete[] solution;

		// Free the problem
		status = CPXfreeprob(env, &problem);

		// Close the cplex environment
		status = CPXcloseCPLEX(&env);
	}


	void MIP_Heuristic::set_initial_solution()
	{
		// set solution
		for (int p = 0; p < nb_people; ++p)
		{
			for (int t = 0; t < nb_tasks; ++t)
			{
				for (int d = 0; d < nb_days; ++d)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						if (current_solution->at(p, t, d, s))
						{
							double rhs[1];					// Right-hand side constraints
							char sense[1];					// Sign of constraint
							int matbeg[1];					// Begin position of the constraint
							int matind[1];					// Position of each element in constraint matrix
							double matval[1];				// Value of each element in constraint matrix

							rhs[0] = 1; // X_ptds = 1
							sense[0] = 'E';
							matbeg[0] = 0;
							matind[0] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
							matval[0] = 1;

							status = CPXaddrows(env, problem, 0, 1, 1, rhs, sense, matbeg, matind, matval, NULL, NULL);
						}
					}
				}
			}
		}

		// solve problem and save objective values
		solve_problem();
		current_solution->objective = candidate_objective;
		best_solution->objective = candidate_objective;
		std::cout << "\n\n\nInitial solution: objective value = " << best_solution->objective;
	}


	void MIP_Heuristic::probabilistic_VNDS()
	{
		// PARAMETERS
		int nb_people_shake = 1;
		constexpr int max_iter_without_improvement = 10;
		constexpr double probabilities[] = { 0.40, 0.30, 0.20, 0.08, 0.02 };		// D-2, D-3, D-4, S-28, T-10


		std::uniform_int_distribution<int> dist_people(0, nb_people - 1);
		std::uniform_real_distribution<double> dist_neighbourhood(0, 1);

		int previous_neighbourhood = -1;	// 0 = D-2, 1 = D-3, 2 = D-4, 3 = S, 4 = T
		int previous_start_day = -1;
		int previous_end_day = -1;
		int iterations_without_improvement = 0;

		// VNDS
		// 1. Local search
		// 2. Move or not 
		// 3. Shake
		Solution base_solution;
		base_solution.save(*current_solution);

		while (true)
		{
			// A. Randomly select neighbourhoods until iterations_without_improvement exceeds threshold, then go to B.
			iterations_without_improvement = 0;
			while (iterations_without_improvement < max_iter_without_improvement)
			{
				// A.1. choose neighbourhood
				int neighbourhood = -1;
				do {
					double prob = dist_neighbourhood(generator);
					double cumul_prob = 0;
					while (cumul_prob < prob) {
						++neighbourhood;
						cumul_prob += probabilities[neighbourhood];
					}
					// if previously neighbourhood S or T, choose other neighbourhood
				} while (neighbourhood >= 3 && neighbourhood == previous_neighbourhood);
				previous_neighbourhood = neighbourhood;


				// A.2. solve subproblem
				if (neighbourhood == 0) // D-2
				{
					std::uniform_int_distribution<int> dist_start_day(0, nb_days - 2);

					int start_day, end_day;
					do {
						start_day = dist_start_day(generator);
						end_day = start_day + 1;
						// new interval of days should have different start and end day than previous interval
					} while (start_day == previous_start_day || end_day == previous_end_day);
					previous_start_day = start_day; previous_end_day = end_day;

					if (neighbourhood_days(start_day, end_day))
						iterations_without_improvement = 0;
					else
						++iterations_without_improvement;
				}
				else if (neighbourhood == 1) // D-3
				{
					std::uniform_int_distribution<int> dist_start_day(0, nb_days - 3);

					int start_day, end_day;
					do {
						start_day = dist_start_day(generator);
						end_day = start_day + 2;
						// new interval of days should have different start and end day than previous interval
					} while (start_day == previous_start_day || end_day == previous_end_day);
					previous_start_day = start_day; previous_end_day = end_day;

					if (neighbourhood_days(start_day, end_day))
						iterations_without_improvement = 0;
					else
						++iterations_without_improvement;
				}
				else if (neighbourhood == 2) // D-4
				{
					std::uniform_int_distribution<int> dist_start_day(0, nb_days - 4);

					int start_day, end_day;
					do {
						start_day = dist_start_day(generator);
						end_day = start_day + 3;
						// new interval of days should have different start and end day than previous interval
					} while (start_day == previous_start_day || end_day == previous_end_day);
					previous_start_day = start_day; previous_end_day = end_day;

					if (neighbourhood_days(start_day, end_day))
						iterations_without_improvement = 0;
					else
						++iterations_without_improvement;
				}
				else if (neighbourhood == 3) // S-28
				{
					// for every individual shift, solve a seperate subproblem
					bool solution_improved = false;
					for (int s = 0; s < nb_shifts; ++s) {
						if (neighbourhood_shifts(s)) {
							solution_improved = true;
						}
						// TIME
						elapsed_computation_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
						if (elapsed_computation_time >= allowed_computation_time)
							return;
					}

					// check if solution has been improved
					if (solution_improved)
						iterations_without_improvement = 0;
					else
						++iterations_without_improvement;
				}
				else  // T-10
				{
					// choose tasks at random
					std::bernoulli_distribution dist_choose_task(10.0 / nb_tasks);
					std::vector<int> tasks; tasks.reserve(10);
					for (int t = 0; t < nb_tasks; ++t) {
						if (dist_choose_task(generator))
							tasks.push_back(t);
					}

					// for every individual task, solve a seperate subproblem
					bool solution_improved = false;
					for (int i = 0; i < tasks.size(); ++i) {
						if (neighbourhood_tasks(tasks[i])) {
							solution_improved = true;
						}
						// TIME
						elapsed_computation_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
						if (elapsed_computation_time >= allowed_computation_time)
							return;
					}

					// check if solution has been improved
					if (solution_improved)
						iterations_without_improvement = 0;
					else
						++iterations_without_improvement;
				}

				// TIME
				elapsed_computation_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
				if (elapsed_computation_time >= allowed_computation_time)
					return;
			}



			// B. Move or not
			if (current_solution->objective < base_solution.objective) {
				base_solution.save(*current_solution);
				nb_people_shake = 1;
			}
			else
				++nb_people_shake;


			// C. Shake
			// start search from base solution
			current_solution->save(base_solution);

			std::vector<int> people_to_shake;
			int person = dist_people(generator);
			people_to_shake.push_back(person);
			for (int pp = 1; pp < nb_people_shake; ++pp) {
				do {
					person = dist_people(generator);
				} while (person == people_to_shake.back());
				people_to_shake.push_back(person);
			}
			shaking(people_to_shake);

			// TIME
			elapsed_computation_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
			if (elapsed_computation_time >= allowed_computation_time)
				return;
		}
	}


	bool MIP_Heuristic::neighbourhood_days(int first_day, int last_day)
	{
		// delete old constraints
		int nb_rows_to_delete = CPXgetnumrows(env, problem) - nb_main_constraints;
		if (nb_rows_to_delete > 0)
			status = CPXdelrows(env, problem, nb_main_constraints, CPXgetnumrows(env, problem) - 1);

		// add new constraints
		for (int p = 0; p < nb_people; ++p)
		{
			for (auto&& ass : current_solution->get_tasks_person(p))
			{
				if (ass.day < first_day || ass.day > last_day)
				{
					double rhs[1];					// Right-hand side constraints
					char sense[1];					// Sign of constraint
					int matbeg[1];					// Begin position of the constraint
					int matind[1];					// Position of each element in constraint matrix
					double matval[1];				// Value of each element in constraint matrix

					rhs[0] = 1;
					sense[0] = 'E';
					matbeg[0] = 0;
					matind[0] = p * nb_tasks*nb_days*nb_shifts + ass.task*nb_days*nb_shifts + ass.day*nb_shifts + ass.shift;
					matval[0] = 1;

					status = CPXaddrows(env, problem, 0, 1, 1, rhs, sense, matbeg, matind, matval, NULL, NULL);
				}
			}
		}


		// solve the problem
		elapsed_computation_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
		if (allowed_computation_time - elapsed_computation_time < time_limit_subproblem)
			status = CPXsetdblparam(env, CPXPARAM_TimeLimit, allowed_computation_time - elapsed_computation_time);
		std::cout << "\n\n\nCPLEX is solving the subproblem from day " << first_day + 1 << " till day " << last_day + 1 << " ... ";
		solve_problem();

		// update the solution
		// if better than current solution, save it; else do nothing
		if (candidate_objective < current_solution->objective)
		{
			// candidate -> current
			current_solution->objective = candidate_objective;
			current_solution->reset();
			for (int p = 0; p < nb_people; ++p) {
				for (int t = 0; t < nb_tasks; ++t) {
					for (int d = 0; d < nb_days; ++d) {
						for (int s = 0; s < nb_shifts; ++s) {
							if (solution[p*nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s] > 0.99) {
								current_solution->set_at(p, t, d, s, true);
							}
						}
					}
				}
			}
			// save algorithm progress
			time_evolution_current_sol.push_back(Performance_Info());
			time_evolution_current_sol.back().elapsed_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
			time_evolution_current_sol.back().neighbourhood.assign("D-");
			time_evolution_current_sol.back().neighbourhood.append(std::to_string(first_day));
			time_evolution_current_sol.back().neighbourhood.push_back('-');
			time_evolution_current_sol.back().neighbourhood.append(std::to_string(last_day));
			time_evolution_current_sol.back().objective = current_solution->objective;
			time_evolution_current_sol.back().nb_shakes = nb_shakes;

			// best solution
			if (current_solution->objective < best_solution->objective)
			{
				best_solution->save(*current_solution);

				// save algorithm progress
				time_evolution_best_sol.push_back(Performance_Info());
				time_evolution_best_sol.back().elapsed_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
				time_evolution_best_sol.back().neighbourhood.assign("D-");
				time_evolution_best_sol.back().neighbourhood.append(std::to_string(first_day));
				time_evolution_best_sol.back().neighbourhood.push_back('-');
				time_evolution_best_sol.back().neighbourhood.append(std::to_string(last_day));
				time_evolution_best_sol.back().objective = best_solution->objective;
				time_evolution_best_sol.back().nb_shakes = nb_shakes;
			}



			std::cout << "\nCandidate objective value = " << candidate_objective;
			std::cout << "\nCurrent objectiva value = " << current_solution->objective;
			std::cout << "\nBest objectiva value = " << best_solution->objective;

			return true;
		}




		std::cout << "\nCandidate objective value = " << candidate_objective;
		std::cout << "\nCurrent objectiva value = " << current_solution->objective;
		std::cout << "\nBest objectiva value = " << best_solution->objective;

		return false;
	}

	bool MIP_Heuristic::neighbourhood_tasks(int task)
	{
		// delete old constraints
		int nb_rows_to_delete = CPXgetnumrows(env, problem) - nb_main_constraints;
		if (nb_rows_to_delete > 0)
			status = CPXdelrows(env, problem, nb_main_constraints, CPXgetnumrows(env, problem) - 1);

		// add new constraints
		for (int p = 0; p < nb_people; ++p)
		{
			for (auto&& ass : current_solution->get_tasks_person(p))
			{
				if (ass.task != task) // not for the unfixed task
				{
					double rhs[1];					// Right-hand side constraints
					char sense[1];					// Sign of constraint
					int matbeg[1];					// Begin position of the constraint
					int matind[1];					// Position of each element in constraint matrix
					double matval[1];				// Value of each element in constraint matrix

					rhs[0] = 1;
					sense[0] = 'E';
					matbeg[0] = 0;
					matind[0] = p * nb_tasks*nb_days*nb_shifts + ass.task*nb_days*nb_shifts + ass.day*nb_shifts + ass.shift;
					matval[0] = 1;

					status = CPXaddrows(env, problem, 0, 1, 1, rhs, sense, matbeg, matind, matval, NULL, NULL);
				}
			}
		}

		// solve the problem
		elapsed_computation_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
		if (allowed_computation_time - elapsed_computation_time < time_limit_subproblem)
			status = CPXsetdblparam(env, CPXPARAM_TimeLimit, allowed_computation_time - elapsed_computation_time);
		std::cout << "\n\n\nCPLEX is solving the subproblem for task " << task << " ... ";
		solve_problem();

		// update the solution
		// if better than current solution, save it; else do nothing
		if (candidate_objective < current_solution->objective)
		{
			// candidate -> current
			current_solution->objective = candidate_objective;
			current_solution->reset();
			for (int p = 0; p < nb_people; ++p) {
				for (int t = 0; t < nb_tasks; ++t) {
					for (int d = 0; d < nb_days; ++d) {
						for (int s = 0; s < nb_shifts; ++s) {
							if (solution[p*nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s] > 0.99) {
								current_solution->set_at(p, t, d, s, true);
							}
						}
					}
				}
			}
			// save algorithm progress
			time_evolution_current_sol.push_back(Performance_Info());
			time_evolution_current_sol.back().elapsed_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
			time_evolution_current_sol.back().neighbourhood.assign("T-individual-");
			time_evolution_current_sol.back().neighbourhood.append(std::to_string(task));
			time_evolution_current_sol.back().objective = current_solution->objective;
			time_evolution_current_sol.back().nb_shakes = nb_shakes;

			// best solution
			if (current_solution->objective < best_solution->objective)
			{
				best_solution->save(*current_solution);

				// save algorithm progress
				time_evolution_best_sol.push_back(Performance_Info());
				time_evolution_best_sol.back().elapsed_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
				time_evolution_best_sol.back().neighbourhood.assign("T-individual-");
				time_evolution_best_sol.back().neighbourhood.append(std::to_string(task));
				time_evolution_best_sol.back().objective = best_solution->objective;
				time_evolution_best_sol.back().nb_shakes = nb_shakes;
			}


			std::cout << "\nCandidate objective value = " << candidate_objective;
			std::cout << "\nCurrent objectiva value = " << current_solution->objective;
			std::cout << "\nBest objectiva value = " << best_solution->objective;

			return true;
		}



		std::cout << "\nCandidate objective value = " << candidate_objective;
		std::cout << "\nCurrent objectiva value = " << current_solution->objective;
		std::cout << "\nBest objectiva value = " << best_solution->objective;

		return false;
	}

	bool MIP_Heuristic::neighbourhood_shifts(int shift)
	{
		// delete old constraints
		int nb_rows_to_delete = CPXgetnumrows(env, problem) - nb_main_constraints;
		if (nb_rows_to_delete > 0)
			status = CPXdelrows(env, problem, nb_main_constraints, CPXgetnumrows(env, problem) - 1);

		// add new constraints
		for (int p = 0; p < nb_people; ++p)
		{
			for (auto&& ass : current_solution->get_tasks_person(p))
			{
				if (ass.shift != shift) // not for the unfixed task
				{
					double rhs[1];					// Right-hand side constraints
					char sense[1];					// Sign of constraint
					int matbeg[1];					// Begin position of the constraint
					int matind[1];					// Position of each element in constraint matrix
					double matval[1];				// Value of each element in constraint matrix

					rhs[0] = 1;
					sense[0] = 'E';
					matbeg[0] = 0;
					matind[0] = p * nb_tasks*nb_days*nb_shifts + ass.task*nb_days*nb_shifts + ass.day*nb_shifts + ass.shift;
					matval[0] = 1;

					status = CPXaddrows(env, problem, 0, 1, 1, rhs, sense, matbeg, matind, matval, NULL, NULL);

				}
			}
		}

		// solve the problem
		elapsed_computation_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
		if (allowed_computation_time - elapsed_computation_time < time_limit_subproblem)
			status = CPXsetdblparam(env, CPXPARAM_TimeLimit, allowed_computation_time - elapsed_computation_time);
		std::cout << "\n\n\nCPLEX is solving the subproblem for shift " << shift << " ... ";
		solve_problem();

		// update the solution
		// if better than current solution, save it; else do nothing
		if (candidate_objective < current_solution->objective)
		{
			// candidate -> current
			current_solution->objective = candidate_objective;
			current_solution->reset();
			for (int p = 0; p < nb_people; ++p) {
				for (int t = 0; t < nb_tasks; ++t) {
					for (int d = 0; d < nb_days; ++d) {
						for (int s = 0; s < nb_shifts; ++s) {
							if (solution[p*nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s] > 0.99) {
								current_solution->set_at(p, t, d, s, true);
							}
						}
					}
				}
			}
			// save algorithm progress
			time_evolution_current_sol.push_back(Performance_Info());
			time_evolution_current_sol.back().elapsed_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
			time_evolution_current_sol.back().neighbourhood.assign("S-");
			time_evolution_current_sol.back().neighbourhood.append(std::to_string(shift));
			time_evolution_current_sol.back().objective = current_solution->objective;
			time_evolution_current_sol.back().nb_shakes = nb_shakes;

			// best solution
			if (current_solution->objective < best_solution->objective)
			{
				best_solution->save(*current_solution);

				// save algorithm progress
				time_evolution_best_sol.push_back(Performance_Info());
				time_evolution_best_sol.back().elapsed_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
				time_evolution_best_sol.back().neighbourhood.assign("S-");
				time_evolution_best_sol.back().neighbourhood.append(std::to_string(shift));
				time_evolution_best_sol.back().objective = best_solution->objective;
				time_evolution_best_sol.back().nb_shakes = nb_shakes;
			}


			std::cout << "\nCandidate objective value = " << candidate_objective;
			std::cout << "\nCurrent objectiva value = " << current_solution->objective;
			std::cout << "\nBest objectiva value = " << best_solution->objective;

			return true;
		}



		std::cout << "\nCandidate objective value = " << candidate_objective;
		std::cout << "\nCurrent objectiva value = " << current_solution->objective;
		std::cout << "\nBest objectiva value = " << best_solution->objective;

		return false;
	}

	void MIP_Heuristic::shaking(const std::vector<int>& people_to_shake)
	{
		++nb_shakes;

		// Generate new column(s)
		Column_Person shake_method;
		for (auto&& p : people_to_shake)
			shake_method.find_column(p);



		// Calculate new objective value
		// delete old constraints
		int nb_rows_to_delete = CPXgetnumrows(env, problem) - nb_main_constraints;
		if (nb_rows_to_delete > 0)
			status = CPXdelrows(env, problem, nb_main_constraints, CPXgetnumrows(env, problem) - 1);

		// set solution
		for (int p = 0; p < nb_people; ++p)
		{
			for (int t = 0; t < nb_tasks; ++t)
			{
				for (int d = 0; d < nb_days; ++d)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						if (current_solution->at(p, t, d, s))
						{
							double rhs[1];					// Right-hand side constraints
							char sense[1];					// Sign of constraint
							int matbeg[1];					// Begin position of the constraint
							int matind[1];					// Position of each element in constraint matrix
							double matval[1];				// Value of each element in constraint matrix

							rhs[0] = 1;
							sense[0] = 'E';
							matbeg[0] = 0;
							matind[0] = p * nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
							matval[0] = 1;

							status = CPXaddrows(env, problem, 0, 1, 1, rhs, sense, matbeg, matind, matval, NULL, NULL);
						}
					}
				}
			}
		}

		// solve problem and save objective values
		elapsed_computation_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
		if (allowed_computation_time - elapsed_computation_time < time_limit_subproblem)
			status = CPXsetdblparam(env, CPXPARAM_TimeLimit, allowed_computation_time - elapsed_computation_time);
		solve_problem();
		current_solution->objective = candidate_objective;
		std::cout << "\nCurrent objective value after shake = " << current_solution->objective;
	}





	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// MAIN HEURISTIC 

	void heuristic()
	{
		// time
		start_time = clock();

		// open output file
		output_file.open("algorithm_output.txt");

		// constructive method
		current_solution = std::make_unique<Solution>();
		Initial_Solution initial_algorithm;
		initial_algorithm.run_algorithm();

		// improvement heuristic
		best_solution = std::make_unique<Solution>();
		best_solution->save(*current_solution);
		MIP_Heuristic algorithm;
		algorithm.run_algorithm();

		// Print output
		// Calculate unmet demands
		int unmet_demand = 0;
		int excess_supply = 0;
		for (int d = 0; d < nb_days; ++d)
		{
			for (int s = 0; s < nb_shifts; ++s)
			{
				for (int t = 0; t < nb_tasks; ++t)
				{
					int supply = 0;
					for (int p = 0; p < nb_people; ++p) {
						if (best_solution->at(p, t, d, s))
							++supply;
					}

					if (supply < shift_demands(t, d, s))
						unmet_demand += (shift_demands(t, d, s) - supply);
					else
						excess_supply += (supply - shift_demands(t, d, s));
				}
			}
		}


		// To screen
		std::cout << "\n\n\n\nAlgorithm has terminated.";
		std::cout << "\n\nTime evolution best solution: ";
		std::cout << "\nTime \t\tBest obj. \tNeighbourhood \tNb shakes";
		for (auto&& e : time_evolution_best_sol) {
			std::cout << "\n" << e.elapsed_time;
			std::cout << "\t\t" << e.objective;
			std::cout << "\t\t" << e.neighbourhood;
			std::cout << "\t\t" << e.nb_shakes;
		}
		std::cout << "\n\nTime evolution current solution: ";
		std::cout << "\nTime \t\tCurrent obj. \tNeighbourhood \tNb shakes";
		for (auto&& e : time_evolution_current_sol) {
			std::cout << "\n" << e.elapsed_time;
			std::cout << "\t\t" << e.objective;
			std::cout << "\t\t" << e.neighbourhood;
			std::cout << "\t\t" << e.nb_shakes;
		}
		std::cout << "\n\n\nBest solution value: " << best_solution->objective;
		std::cout << "\nUnmet demand = " << unmet_demand;
		std::cout << "\nExcess supply = " << excess_supply;

		// To file
		output_file << "\nAlgorithm choice: Probabilistic VNDS";
		output_file << "\nAllowed computation time (s): " << allowed_computation_time;
		output_file << "\n\nBest found solution: " << best_solution->objective;
		output_file << "\nUnmet demand = " << unmet_demand;
		output_file << "\nExcess supply = " << excess_supply;
		if (write_time_evolution_to_file)
		{
			output_file << "\n\nTime \t\tBest obj. \tNeighbourhood \tNb shakes";
			for (auto&& e : time_evolution_best_sol) {
				output_file << "\n" << e.elapsed_time;
				output_file << "\t" << e.objective;
				output_file << "\t" << e.neighbourhood;
				output_file << "\t" << e.nb_shakes;
			}
			output_file << "\n\nTime \t\tCurrent obj. \tNeighbourhood \tNb shakes";
			for (auto&& e : time_evolution_current_sol) {
				output_file << "\n" << e.elapsed_time;
				output_file << "\t" << e.objective;
				output_file << "\t" << e.neighbourhood;
				output_file << "\t" << e.nb_shakes;
			}
		}

		output_file.flush();


		// free solutions
		current_solution.reset();
		best_solution.reset();

		// reset statistics
		elapsed_computation_time = 0;
		nb_shakes = 0;
		time_evolution_best_sol.clear();
		time_evolution_current_sol.clear();

		// close output file
		output_file.close();
	}

}