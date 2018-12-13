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

#include "Standard_IP_Algorithm.h"
#include "Problem_Data.h"
#include <iostream>
#include <ctime>


namespace IP
{
	// USER CHOICE
	int algorithm_choice{ Algorithm_Type::LP };
	double allowed_computation_time{ 3600 };		// in seconds


	// ALGORITHM STATISTICS
	double elapsed_computation_time{ 0 };	// in seconds
	constexpr bool write_to_file = true;




	// CPLEX IP IMPLEMENTATION
	void Standard_IP_Algorithm::run_algorithm()
	{
		initialize_cplex();
		build_problem();

		// set time limit
		status = CPXsetdblparam(env, CPXPARAM_TimeLimit, allowed_computation_time);

		clock_t start_time = clock();
		if (algorithm_choice == Algorithm_Type::IP)
			solve_problem_as_IP();
		else if (algorithm_choice == Algorithm_Type::LP)
			solve_problem_as_LP();
		elapsed_computation_time = (double)(clock() - start_time) / CLK_TCK;
		std::cout << "\n\nElapsed time: " << elapsed_computation_time;
		std::cout << "\nObjective: " << objective;

		double unmet_demand_CODU = 0, unmet_demand_amb = 0;
		for (int i = nb_people * nb_tasks*nb_days*nb_shifts + nb_tasks * nb_days*nb_shifts; i < nb_people*nb_tasks*nb_days*nb_shifts + 2 * nb_tasks*nb_days*nb_shifts; ++i) {
			if (solution[i] > 0.001) {
				if (i < nb_people*nb_tasks*nb_days*nb_shifts + nb_tasks * nb_days*nb_shifts + nb_tasks_CODU * nb_days*nb_shifts) {
					unmet_demand_CODU += solution[i];
				}
				else {
					unmet_demand_amb += solution[i];
				}
			}
		}
		std::cout << "\n\nUnmet demand CODU = " << unmet_demand_CODU;
		std::cout << "\nUnmet demand EV   = " << unmet_demand_amb;

		clear_cplex();
	}

	void Standard_IP_Algorithm::initialize_cplex()
	{
		// Open the cplex environment
		env = CPXopenCPLEX(&status);

		// Set the output to screen on/off
		status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	}

	void Standard_IP_Algorithm::build_problem()
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
	}

	void Standard_IP_Algorithm::solve_problem_as_IP()
	{
		// Allocate memory for solution
		numcols = CPXgetnumcols(env, problem);		// how many columns (= variables) are there in the problem
		solution = new double[numcols];				// allocate memory/create array for solution

		// Try to optimize and check if the model works
		try
		{
			// Optimize the problem
			std::cout << "\n\nCPLEX is solving the IP model ... \n\n";
			status = CPXmipopt(env, problem);

			// Get the solution
			status = CPXsolution(env, problem, &solstat, &objective, solution, NULL, NULL, NULL);

			if (status != 0)
				throw status;		// if something went wrong, throw exception

			if (!(solstat == 101 || solstat == 102))
				throw solstat;		// if optimal solution not found, throw exception

			// Print solution (only executed when no exception thrown)
			std::cout << "\nCLPEX has found the optimal solution!";
			std::cout << "\nObjective value : " << objective;
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
			else if (exception == 107) {
				std::cout << "\nTime limit exceeded, integer solution exists.\n";
				std::cout << "\nObjective value : " << objective << "\n\n\n";
			}
			else if (exception == 108) std::cout << "\nTime limit exceeded, no integer solution.\n\n\n";
			else if (exception == 111) std::cout << "\nTreememory limit, integer solution exists.\n\n\n";
			else if (exception == 112) std::cout << "\nTreememory limit, no integer solution exists.\n\n\n";
			else std::cout << "\nOther reason for termination.\n\n\n";
		} // end catch
	}

	void Standard_IP_Algorithm::solve_problem_as_LP()
	{
		// Allocate memory for solution
		numcols = CPXgetnumcols(env, problem);		// how many columns (= variables) are there in the problem
		solution = new double[numcols];				// allocate memory/create array for solution

		// Make LP
		status = CPXchgprobtype(env, problem, CPXPROB_LP);

		// Try to optimize and check if the model works
		// Optimize the problem
		std::cout << "\n\nCPLEX is solving the LP model ... \n\n";
		status = CPXlpopt(env, problem);

		// Get the solution
		status = CPXsolution(env, problem, &solstat, &objective, solution, NULL, NULL, NULL);

		std::cout << "\n\nOptimal solution value = " << objective;
	}

	void Standard_IP_Algorithm::clear_cplex()
	{
		// Free solution 
		delete[] solution;

		// Free the problem
		status = CPXfreeprob(env, &problem);

		// Close the cplex environment
		status = CPXcloseCPLEX(&env);
	}

}