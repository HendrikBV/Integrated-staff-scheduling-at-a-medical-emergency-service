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

#include "Diving_Heuristic.h"
#include "Problem_Data.h"
#include <iostream>
#include <fstream>


namespace Diving_Heuristic
{
	// USER CHOICE
	int column_generation_method{ Column_Generation_Method::one_column_person_p_and_reoptimize };
	int branching_method_diving{ Branching_Method_Diving::value_above_threshold };
	double branching_threshold_diving{ 0.6 };

	double allowed_computation_time{ 3600 };			// in seconds
	double allowed_computation_time_root{ 1800 };		// in seconds
	double allowed_computation_time_node{ 10 };			// in seconds


	// FIXED ALGORITHM SETTINGS
	constexpr double reduced_cost_tolerance_LP{ 0.001 };		// the tolerance for deciding whether a reduced cost is smaller than 0; also used by CPLEX when solving an LP to determine whether optimum has been found (CPLEX default value = 0.000001 (1e-6))
	constexpr double fractionality_tolerance_master{ 0.001 };	// the tolerance for deciding whether the master problem is integer
	
	constexpr bool write_to_file{ false };						// write models to an lp-file
	constexpr bool write_to_file_each_iteration{ false };		// write newly added columns to file during each iteration
	constexpr bool write_to_file_each_branch{ false };			// write new branch to file
	constexpr bool print_each_CG_iter{ true };					// print the values of each iteration during the column generation to screen
	constexpr bool print_diving_info{ true };

	
	// OUTPUT / STATISTICS
	double computation_time{ 0 };	// computation time used by the algorithm in seconds
	std::ofstream output_file;





	// COLUMN GENERATION
	void Diving_Column_Generation::run_algorithm()
	{
		// Time
		start_time_total_diving = clock();

		// Open output file
		output_file.open("solution.txt");

		// Initialize cplex
		initialize_cplex();

		// Build the problems (master is initialized with |P| supercolumns)
		//std::cout << "\n\nCPLEX is building the masterproblem ... ";
		build_master_problem();
		for (int p = 0; p < nb_people; ++p) {
			//std::cout << "\nCPLEX is building subproblem " << p+1 << " ... ";
			build_subproblem_problem(p);
		}

		// Create arrays for CPLEX
		create_cplex_arrays();

		// Column generation + diving
		diving_heuristic();

		// Time
		elapsed_time_total = (double)(clock() - start_time_total_diving) / CLOCKS_PER_SEC;

		// Print solution
		print_solution();

		// Clear arrays for CPLEX and CPLEX models
		clear_cplex();

		// Close output file
		output_file.close();
	}



	double Diving_Column_Generation::column_generation()
	{
		bool LP_optimum_found{ false };			// Do we need a new iteration for the column generation?
		double reduced_cost;					// Reduced cost of the current subproblem
		double obj_master_LP = 1e+100;			// Current value of the restricted master LP (optimum will be lower when all optimal columns present)
		clock_t start_time = clock();


		// Stop after fixed computation time
		int iteration = 0;

		// DIFFERENT CG METHODS
		// CG METHOD 1
		if (column_generation_method == Column_Generation_Method::one_column_per_person)
		{
			while (!LP_optimum_found)
			{
				++iterations_CG;
				++iteration;

				LP_optimum_found = true;

				// solve master problem
				{
					clock_t st_ = clock();
					obj_master_LP = solve_masterproblem();
					elapsed_time_CG_masterproblem += (double)(clock() - st_) / CLK_TCK;
				}


				if (print_each_CG_iter) {
					std::cout << "\n\nIteration: " << iterations_CG;
					std::cout << "\nObjective master: " << obj_master_LP;
				}


				// allowed computation time
				double elapsedtime = (double)(clock() - start_time) / CLK_TCK;
				if (first_call_CG) {
					if (elapsedtime > allowed_computation_time_root)
						break;
				}
				else {
					if (elapsedtime > allowed_computation_time_node)
						break;
				}



				for (int p = 0; p < nb_people; ++p)
				{
					if (solve_subproblem_p[p])
					{
						solve_subproblem_p[p] = false;

						change_coefficients_subproblem(p);

						// solve subproblem
						clock_t st_ = clock();
						reduced_cost = solve_subproblem(p);
						elapsed_time_CG_subproblem += (double)(clock() - st_) / CLK_TCK;

						if (print_each_CG_iter) {
							std::cout << "\nReduced cost problem " << p + 1 << ": " << reduced_cost;
						}

						if (reduced_cost < -reduced_cost_tolerance_LP) // rounding errors
						{
							++nb_columns_added_total;
							add_column_to_master(p);
							LP_optimum_found = false;
							solve_subproblem_p[p] = true;
						}
					}
				}
			}
		}


		// CG METHOD 2
		else if (column_generation_method == Column_Generation_Method::one_column_person_p_and_reoptimize)
		{
			int person = 0;
			int last_person_improved = 0;

			while (!LP_optimum_found)
			{
				++iterations_CG;
				++iteration;

				// allowed computation time
				double elapsedtime = (double)(clock() - start_time) / CLK_TCK;
				if (first_call_CG) {
					if (elapsedtime > allowed_computation_time_root)
						break;
				}
				else {
					if (elapsedtime > allowed_computation_time_node)
						break;
				}

				// solve master problem
				{
					clock_t st_ = clock();
					obj_master_LP = solve_masterproblem();
					elapsed_time_CG_masterproblem += (double)(clock() - st_) / CLK_TCK;
				}

				if (print_each_CG_iter) {
					std::cout << "\n\nIteration: " << iterations_CG;
					std::cout << "\nObjective master: " << obj_master_LP;
				}

				reduced_cost = 999999; // to prevent solving master if last RC < 0, but then solve_subproblem_p = false for next person 

				do {
					if (solve_subproblem_p[person]) {
						solve_subproblem_p[person] = false;
						change_coefficients_subproblem(person);

						// solve subproblem
						clock_t st_ = clock();
						reduced_cost = solve_subproblem(person);
						elapsed_time_CG_subproblem += (double)(clock() - st_) / CLK_TCK;

						if (print_each_CG_iter) {
							std::cout << "\nReduced cost person " << person + 1 << " = " << reduced_cost;
						}


						if (reduced_cost < -reduced_cost_tolerance_LP) {
							++nb_columns_added_total;
							add_column_to_master(person);
							last_person_improved = person;
							solve_subproblem_p[person] = true;
						}
					}

					++person;
					if (person >= nb_people)
						person = 0;
					if (person == last_person_improved) {
						LP_optimum_found = true;
						break;
					}
				} while (reduced_cost > -reduced_cost_tolerance_LP);
			}
		}






		// statistics
		if (first_call_CG) {
			iterations_CG_root_node = iterations_CG;
			nb_columns_added_root_node = nb_columns_added_total;
			elapsed_time_CG_root_node = (double)(clock() - start_time) / CLK_TCK;
			lowerbound_root_node = obj_master_LP;
			first_call_CG = false;
		}

		elapsed_time_CG_total += (double)(clock() - start_time) / CLK_TCK;


		return obj_master_LP;
	}



	void Diving_Column_Generation::initialize_cplex()
	{
		// Open the cplex environment
		env = CPXopenCPLEX(&status);

		// Set the output to screen on/off
		status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);

		// Initialize the vector for subproblem pointers
		subproblems.reserve(nb_people);
		for (int p = 0; p < nb_people; ++p)
			subproblems.push_back(nullptr);
	}



	void Diving_Column_Generation::build_master_problem()
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
		int *matind;					// Position of each element in constraint matrix
		double *matval;					// Value of each element in constraint matrix
		int f{ 0 };						// To calculate number of nonzero coefficients in each constraint

										// create arrays
		matind = new int[100000];
		matval = new double[100000];


		// Create the problem
		masterproblem = CPXcreateprob(env, &status, "INEM_CODU_master_problem");

		// Problem is minimization
		CPXchgobjsen(env, masterproblem, CPX_MIN);



		// VARIABLES
		// Add the Y_REplus variables
		for (int t = 0; t < nb_tasks; ++t)
		{
			for (int d = 0; d < nb_days; ++d)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					std::string name = "Y_REplus_" + std::to_string(t + 1) + "_" + std::to_string(d + 1) + "_" + std::to_string(s + 1);
					colname[0] = const_cast<char*>(name.c_str());


					obj[0] = obj_weight_Y_REplus;	// objective function coefficient

					lb[0] = 0;

					status = CPXnewcols(env, masterproblem, 1, obj, lb, NULL, NULL, colname);
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
									//type[0] = 'I';

					status = CPXnewcols(env, masterproblem, 1, obj, lb, NULL, NULL, colname);
				}
			}
		}

		// Add the Z_pk variables (start with one supercolumn for every person)
		for (int p = 0; p < nb_people; ++p)
		{
			std::string name = "z_" + std::to_string(p + 1) + "_" + std::to_string(0);
			colname[0] = const_cast<char*>(name.c_str());

			obj[0] = obj_value_super_column;	// objective function coefficient

			lb[0] = 0;		// binary variable
			ub[0] = 1;
			//type[0] = 'B';

			status = CPXnewcols(env, masterproblem, 1, obj, lb, ub, NULL, colname); // Generate columns (the variables) and subsequently add rows (constraints)
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

					std::string name = "Coverage_constraint_CODU_task_" + std::to_string(t + 1) + "_day_" + std::to_string(d + 1) + "_shift_" + std::to_string(s + 1);
					rowname[0] = const_cast<char*>(name.c_str());

					matbeg[0] = 0;
					f = 0;

					for (int p = 0; p < nb_people; ++p)
					{
						matind[f] = 2 * nb_tasks*nb_days*nb_shifts + p;
						//matval[f] = 1; // supercolumn: every person can do every task
						matval[f] = 0; // supercolumn are empty
						++f;
					}

					// REplus
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
					matval[f] = -1;
					++f;

					// REmin
					matind[f] = nb_tasks * nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
					matval[f] = 1;
					++f;

					status = CPXaddrows(env, masterproblem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
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

					std::string name = "Coverage_constraint_ambulances_task_" + std::to_string(t + 1) + "_day_" + std::to_string(d + 1) + "_shift_" + std::to_string(s + 1);
					rowname[0] = const_cast<char*>(name.c_str());

					matbeg[0] = 0;
					f = 0;

					for (int p = 0; p < nb_people; ++p)
					{
						matind[f] = 2 * nb_tasks*nb_days*nb_shifts + p;
						//matval[f] = 1; // supercolumn: every person can do every tasks
						matval[f] = 0; // supercolumn are empty
						++f;
					}

					// REplus
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
					matval[f] = -1;
					++f;

					// REmin
					matind[f] = nb_tasks * nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s;
					matval[f] = 1;
					++f;

					status = CPXaddrows(env, masterproblem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
				}
			}
		}

		// Constraint set 2: one column per person
		for (int p = 0; p < nb_people; ++p)
		{
			sense[0] = 'E';
			rhs[0] = 1;

			std::string name = "One_column_for_person_" + std::to_string(p + 1);
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			matind[f] = 2 * nb_tasks*nb_days*nb_shifts + p;
			matval[f] = 1;
			++f;

			status = CPXaddrows(env, masterproblem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}



		// Write to the problem to a file
		if (write_to_file)
			status = CPXwriteprob(env, masterproblem, "INEM_master_problem.lp", NULL);


		// delete memory
		delete[] matind;
		delete[] matval;
	}



	void Diving_Column_Generation::build_subproblem_problem(int person)
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


										// Create the problem
		{
			std::string name = "INEM_CODUE_subproblem_person_" + std::to_string(person + 1);
			subproblems[person] = CPXcreateprob(env, &status, name.c_str());
		}

		// Problem is minimization
		CPXchgobjsen(env, subproblems[person], CPX_MIN);



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

					obj[0] = 0;	// objective function coefficient, will be changed later at each iteration

					lb[0] = 0;		// binary variable
					ub[0] = 1;
					type[0] = 'B';

					status = CPXnewcols(env, subproblems[person], 1, obj, lb, ub, type, colname); // Generate columns (the variables) and subsequently add rows (constraints)
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

			status = CPXnewcols(env, subproblems[person], 1, obj, lb, NULL, type, colname);
		}


		// Add the Y_Wmin_w variables
		for (int w = 0; w < nb_weekends; ++w)
		{
			std::string name = "Y_Wmin_" + std::to_string(w + 1);
			colname[0] = const_cast<char*>(name.c_str());

			obj[0] = obj_weight_Y_W;	// objective function coefficient

			lb[0] = 0;		// integer variable
			type[0] = 'I';

			status = CPXnewcols(env, subproblems[person], 1, obj, lb, NULL, type, colname);
		}


		// Add the Y_Hplus variable
		{
			std::string name = "Y_Hplus";
			colname[0] = const_cast<char*>(name.c_str());
		}

		obj[0] = obj_weight_Y_Hplus;	// objective function coefficient

		lb[0] = 0;		// integer variable
		type[0] = 'I';

		status = CPXnewcols(env, subproblems[person], 1, obj, lb, NULL, type, colname);


		// Add the Y_Hmin variable
		{
			std::string name = "Y_Hmin";
			colname[0] = const_cast<char*>(name.c_str());
		}

		obj[0] = obj_weight_Y_Hmin;	// objective function coefficient

		lb[0] = 0;		// integer variable
		type[0] = 'I';

		status = CPXnewcols(env, subproblems[person], 1, obj, lb, NULL, type, colname);



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

			status = CPXnewcols(env, subproblems[person], 1, obj, lb, NULL, type, colname);
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

			status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
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

			status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
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

			status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
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

						std::string name = "Person_" + std::to_string(person + 1) + "_cannot_do_task" + std::to_string(t + 1) + "_day_" + std::to_string(d + 1) + "_shift_" + std::to_string(s + 1);
						rowname[0] = const_cast<char*>(name.c_str());

						matbeg[0] = 0;
						f = 0;

						matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s;
						matval[f] = 1;
						++f;

						status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
					}
				}
			}
		}


		// Constraint set 4: Maximum of 6 consecutive working days 
		for (int r = 0; r < nb_days - 6; ++r)
		{
			sense[0] = 'L';
			rhs[0] = 6;

			std::string name = "Maximum_6_consecutive_working_days_from_day_" + std::to_string(person + 1);
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

			status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// Constraint set 5: Maximum of 5 consecutive days-off 
		for (int r = 0; r < nb_days - 5; ++r)
		{
			sense[0] = 'G';
			rhs[0] = 1;

			std::string name = "maximum_5_consecutive_days_off_from_day_" + std::to_string(person + 1);
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

			status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// Constraint set 7: At least 1 Sunday off in every 4 Sundays 
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

		status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);


		// Constraint set 8: Ideally Saturday and Sunday off 
		// depending on start day of planning horizon
		int w = 0;
		for (int d = 6 - start_day; d < nb_days; d = d + 7)
		{
			if (d > 0)
			{
				sense[0] = 'E';
				rhs[0] = 0;

				std::string name = "whole_weekend_week_" + std::to_string(d / 7 + 1);
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

				status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);

				++w;
			}
		}


		// Constraint set 9: 35 weekly working hours (140 per 28 days)
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

		status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);



		// Constraint set 12: Assign a task t from a group g to a person p from that group 
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

			status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}

		// (13.1) night shifts
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

		status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);


		// (13.2) morning shifts
		sense[0] = 'G';
		rhs[0] = 2;

		std::string name = "Person_" + std::to_string(person + 1) + "_at_least_two_morning_shifts";
		rowname[0] = const_cast<char*>(name.c_str());

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

		status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);


		// (13.3) afternoon shifts
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

		status = CPXaddrows(env, subproblems[person], 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);





		// Write to the problem to a file
		if (write_to_file && person == 0) {
			status = CPXwriteprob(env, subproblems[person], "INEM_subproblem.lp", NULL);
		}
	}



	double Diving_Column_Generation::solve_masterproblem()
	{
		// Create memory for solution
		numrows_master = CPXgetnumrows(env, masterproblem);
		delete[] dual_var_master;
		dual_var_master = new double[numrows_master];

		numcols_master = CPXgetnumcols(env, masterproblem);
		delete[] solution_master;
		solution_master = new double[numcols_master];


		// Optimize the problem
		status = CPXlpopt(env, masterproblem);

		// Get the solution
		status = CPXsolution(env, masterproblem, &solstat, &objective_master, solution_master, dual_var_master, NULL, NULL);

		// Check the solution status
		check_solution_status(false, status, solstat);

		// Return the solution value
		return objective_master;
	}



	double Diving_Column_Generation::solve_subproblem(int person)
	{
		// Upperbound is equal to mu_p, i.e. if no solution can be found with obj < mu_p => obj-mu_p > 0, i.e. no column feasible column exists with RC < 0
		status = CPXsetdblparam(env, CPX_PARAM_CUTLO, mu_p);

		// Optimize the problem
		status = CPXmipopt(env, subproblems[person]);

		// Get the solution
		status = CPXsolution(env, subproblems[person], &solstat, &objective_subproblem, solution_subproblem, NULL, NULL, NULL);

		// Check the solution status
		if (!check_solution_status(true, status, solstat))
			objective_subproblem = 1000000000000; // no column with reduced cost found

												  // Calculate the cost of the current column
		cost_column_k_person_p = 0;
		for (int w = 0; w < 2 * nb_weekends; ++w)  // Y_wplus and Y_wmin
			cost_column_k_person_p += obj_weight_Y_W * solution_subproblem[nb_tasks*nb_days*nb_shifts + w];
		cost_column_k_person_p += obj_weight_Y_Hplus * solution_subproblem[nb_tasks*nb_days*nb_shifts + 2 * nb_weekends]; // Y_Hplus
		cost_column_k_person_p += obj_weight_Y_Hmin * solution_subproblem[nb_tasks*nb_days*nb_shifts + 2 * nb_weekends + 1]; // Y_Hmin
		for (int g = 0; g < nb_groups_CODU; ++g) // Y_g_CODU
			cost_column_k_person_p += obj_weight_Y_G_CODU * solution_subproblem[nb_tasks*nb_days*nb_shifts + 2 * nb_weekends + 2 + g];
		for (int g = nb_groups_CODU; g < nb_groups; ++g) // Y_g_ambulances
			cost_column_k_person_p += obj_weight_Y_G_ambulances * solution_subproblem[nb_tasks*nb_days*nb_shifts + 2 * nb_weekends + 2 + g];

		// Return the reduced cost of the current column
		return objective_subproblem - mu_p;
	}



	void Diving_Column_Generation::change_coefficients_subproblem(int person)
	{
		int cnt = nb_tasks * nb_days*nb_shifts;
		mu_p = dual_var_master[cnt + person];
		for (int i = 0; i < cnt; ++i)
			values_coef[i] = -dual_var_master[i];

		status = CPXchgobj(env, subproblems[person], cnt, indices_sub, values_coef);

		// Write to the problem to a file
		if (write_to_file_each_iteration && person == 0)
			status = CPXwriteprob(env, subproblems[person], "INEM_subproblem.lp", NULL);
	}



	void Diving_Column_Generation::add_column_to_master(int person)
	{
		// Local variables
		double obj[1];					// Objective function
		double lb[1];					// Lower bound variables
		double ub[1];					// Upper bound variables
		char *colname[1];				// Variable names
		int matbeg[1];					// Begin position of the constraint
		int matind[20000];				// Position of each element in constraint matrix
		double matval[20000];			// Value of each element in constraint matrix
		int f{ 0 };						// To calculate number of nonzero coefficients in each constraint

		++nb_calls_CG;

		// Variable z_pk
		std::string name = "z_" + std::to_string(person + 1) + "_" + std::to_string(nb_calls_CG);
		colname[0] = const_cast<char*>(name.c_str());

		obj[0] = cost_column_k_person_p;	// objective function coefficient

		lb[0] = 0;		// binary variable
		ub[0] = 1;
		//type[0] = 'B';


		// Create new column object to references columns in the future
		Column new_column;
		new_column.person = person;
		new_column.name.assign(name);
		new_column.cost = cost_column_k_person_p;
		new_column.a_tds.reserve(nb_tasks*nb_days*nb_shifts);


		// Coefficients in the constraints
		// Constraint set 1: coverage requirement
		for (int t = 0; t < nb_tasks; ++t)
		{
			for (int d = 0; d < nb_days; ++d)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					matind[f] = t * nb_days*nb_shifts + d * nb_shifts + s; // index of row
					matval[f] = solution_subproblem[t*nb_days*nb_shifts + d * nb_shifts + s]; // a_tds
					++f;

					// store a_pktds in the columns vector
					new_column.a_tds.push_back((int)(solution_subproblem[t*nb_days*nb_shifts + d * nb_shifts + s] + 0.1));
				}
			}
		}

		// Constraint set 2: one column per person
		matind[f] = nb_tasks * nb_days*nb_shifts + person;
		matval[f] = 1;
		++f;

		matbeg[0] = 0;
		status = CPXaddcols(env, masterproblem, 1, f, obj, matbeg, matind, matval, lb, ub, colname);


		// Add the column object to the vector of columns
		current_columns.push_back(std::move(new_column));


		// Write to the problem to a file
		if (write_to_file_each_iteration)
			status = CPXwriteprob(env, masterproblem, "INEM_master_problem.lp", NULL);
	}



	void Diving_Column_Generation::save_solution()
	{
		// solution value
		objective_best_solution = (int)(objective_master + 0.5);

		// the solution itself
		{
			std::vector<int> temp_sol;
			temp_sol.reserve(nb_people*nb_tasks*nb_days*nb_shifts);
			for (int i = 0; i < nb_people*nb_tasks*nb_days*nb_shifts; ++i)
				temp_sol.push_back(0);

			for (auto&& col : fixed_columns) {
				for (int tds = 0; tds < nb_tasks*nb_days*nb_shifts; ++tds) {
					if (col.a_tds[tds] == 1) {
						temp_sol[col.person*nb_tasks*nb_days*nb_shifts + tds] = 1;
					}
				}
			}



			solution.clear();
			solution.reserve(nb_people*nb_days*nb_shifts);
			for (int i = 0; i < nb_people*nb_days*nb_shifts; ++i)
				solution.push_back(-1);

			for (int p = 0; p < nb_people; ++p) {
				for (int d = 0; d < nb_days; ++d) {
					for (int s = 0; s < nb_shifts; ++s) {
						for (int t = 0; t < nb_tasks; ++t) {
							if (temp_sol[p*nb_tasks*nb_days*nb_shifts + t * nb_days*nb_shifts + d * nb_shifts + s] == 1) {
								solution[p*nb_days*nb_shifts + d * nb_shifts + s] = t;
							}
						}
					}
				}
			}
		}
	}



	void Diving_Column_Generation::print_solution()
	{
		// TO SCREEN
		std::cout << "\nAllowed computation time column generation root node: " << allowed_computation_time_root;
		std::cout << "\nAllowed computation time column generation overall diving: " << 3600 - allowed_computation_time_root;

		/*if (column_generation_method == Column_Generation_Method::one_column_per_person)
		std::cout << "\n\n\n\nColumn generation method: one column per person";
		else if (column_generation_method == Column_Generation_Method::one_column_person_p_and_reoptimize)
		std::cout << "\n\n\n\nColumn generation method: one column for person p and reoptimize";

		std::cout << "\nBranching method: variables above threshold. Beta = " << branching_threshold_diving;

		std::cout << "\nAllowed computation time column generation root node: " << allowed_computation_time_root;
		std::cout << "\nAllowed computation time column generation during each dive: " << allowed_computation_time_node;*/


		//std::cout << "\n\n\n\nAlgorithm has found a solution!";
		std::cout << "\n\nObjective value:                                   " << objective_best_solution;
		std::cout << "\nLowerbound root node:                              " << lowerbound_root_node;
		std::cout << "\nMaximal optimality gap:                            " << (double)(objective_best_solution - lowerbound_root_node) / (double)objective_best_solution;
		std::cout << "\nTotal elapsed time (seconds):                      " << elapsed_time_total;
		std::cout << "\nTotal elapsed time column generation:              " << elapsed_time_CG_total;
		std::cout << "\nElapsed time column generation root node:          " << elapsed_time_CG_root_node;
		std::cout << "\nElapsed time column generation master problem:     " << elapsed_time_CG_masterproblem;
		std::cout << "\nElapsed time column generation subproblem:         " << elapsed_time_CG_subproblem;
		std::cout << "\nNumber of total iterations column generation:      " << iterations_CG;
		std::cout << "\nNumber of iterations column generation root node:  " << iterations_CG_root_node;
		std::cout << "\nNumber of columns added total:                     " << nb_columns_added_total;
		std::cout << "\nNumber of columns added root node:                 " << nb_columns_added_root_node;


		double unmet_demand = 0, excess_supply = 0;
		for (int i = nb_tasks * nb_days*nb_shifts; i < 2 * nb_tasks*nb_days*nb_shifts; ++i) {
			if (solution_master[i] > fractionality_tolerance_master) {
				unmet_demand += solution_master[i];
			}
		}
		for (int i = 0; i < nb_tasks*nb_days*nb_shifts; ++i) {
			if (solution_master[i] > fractionality_tolerance_master) {
				excess_supply += solution_master[i];
			}
		}
		std::cout << "\n\nUnmet demand = " << unmet_demand;
		std::cout << "\nExcess supply  = " << excess_supply;




		// TO FILE
		// first time: info and headers of table
		static bool first_call_print = true;
		if (first_call_print)
		{
			first_call_print = false;

			if (column_generation_method == Column_Generation_Method::one_column_per_person)
				output_file << "Column generation method: one column per person";
			else if (column_generation_method == Column_Generation_Method::one_column_person_p_and_reoptimize)
				output_file << "Column generation method: one column for person p and reoptimize";

			if (branching_method_diving == Branching_Method_Diving::largest_fractional_variable)
				output_file << "\nBranching method: largest fractional variable";
			else if (branching_method_diving == Branching_Method_Diving::value_above_threshold) {
				output_file << "\nBranching method: largest fractional variable";
				output_file << "\nBeta = " << branching_threshold_diving;
			}

			output_file << "\n\n\n";
			output_file << "Instance \tObj value \tComp time (s) \tTime master (s) \tTime subproblem (s) "
				"\tCols added total \tUnmet demand \tExcess supply";
		}

		// each dataset
		output_file << "\n";
		output_file << objective_best_solution << "\t" << elapsed_time_total << "\t"
			<< elapsed_time_CG_masterproblem << "\t" << elapsed_time_CG_subproblem << "\t" << nb_columns_added_total
			<< "\t" << unmet_demand << "\t" << excess_supply;

		output_file.flush();
	}



	bool Diving_Column_Generation::check_solution_status(bool integer, int status_, int solstat_)
	{
		// general failures 
		if (status_ != 0) {
			algorithm_failed = true;
			algorithm_solution_status = "\nError 1217: No solution exists. \nThe requested command cannot be executed because no solution exists for the problem. \nOptimize the problem first.";
			return false;
		}


		// LP optimizations
		if (!integer)
		{
			if (solstat_ == CPX_STAT_OPTIMAL) {
				algorithm_failed = false;
				algorithm_solution_status = "CPLEX has found the optimal solution.";
			}
			else if (solstat_ == CPX_STAT_UNBOUNDED) {
				algorithm_failed = true;
				algorithm_solution_status = "\nProblem is unbounded.";
			}
			else if (solstat_ == CPX_STAT_INFEASIBLE) {
				algorithm_failed = true;
				algorithm_solution_status = "\nProblem is infeasible.";
			}
			else if (solstat_ == CPX_STAT_INForUNBD) {
				algorithm_failed = true;
				algorithm_solution_status = "\nProblem is unbounded or infeasible.";
			}
		}


		// MIP optimizations
		else
		{
			// algorithm succesful
			if (solstat_ == CPXMIP_OPTIMAL || solstat_ == CPXMIP_OPTIMAL_TOL) {
				algorithm_failed = false;
				algorithm_solution_status = "CPLEX has found the optimal solution.";
			}
			/*else if (solstat_ == CPXMIP_ABORT_FEAS && m_solution_exists && elapsed_computation_time > AS::allowed_computation_time) {
			algorithm_failed = false;
			algorithm_solution_status = QStringLiteral("Search has been terminated due to the user-specified time limit; but a solution exists.");
			}
			else if (solstat_ == CPXMIP_ABORT_FEAS && m_solution_exists) {
			algorithm_failed = false;
			algorithm_solution_status = QStringLiteral("CPLEX has found a solution within the user-specified optimality gap.");
			}*/

			// algorithm failure
			else if (solstat_ == CPXMIP_UNBOUNDED) {
				algorithm_failed = true;
				algorithm_solution_status = "\nProblem is unbounded.";
			}
			else if (solstat_ == CPXMIP_INFEASIBLE) {
				algorithm_failed = true;
				algorithm_solution_status = "\nProblem is infeasible.";
			}
			else if (solstat_ == CPXMIP_INForUNBD) {
				algorithm_failed = true;
				algorithm_solution_status = "\nProblem is unbounded or infeasible.";
			}
			else if (solstat_ == CPXMIP_OPTIMAL_INFEAS) {
				algorithm_failed = true;
				algorithm_solution_status = "\nProblem optimal with unscaled infeasibilities.";
			}
			else if (solstat_ == CPXMIP_TIME_LIM_FEAS) {
				algorithm_failed = true;
				algorithm_solution_status = "\nTime limit exceeded, integer solution exists.";
			}
			else if (solstat_ == CPXMIP_TIME_LIM_INFEAS) {
				algorithm_failed = true;
				algorithm_solution_status = "\nTime limit exceeded, no integer solution.";
			}
			else if (solstat_ == CPXMIP_MEM_LIM_FEAS) {
				algorithm_failed = true;
				algorithm_solution_status = "\nTreememory limit, integer solution exists.";
			}
			else if (solstat_ == CPXMIP_MEM_LIM_INFEAS) {
				algorithm_failed = true;
				algorithm_solution_status = "\nTreememory limit, no integer solution exists.";
			}
			else {
				algorithm_failed = true;
				algorithm_solution_status = "\nOther reason for failure.";
			}
		}

		if (algorithm_failed)
			return false;
		else
			return true;
	}



	void Diving_Column_Generation::create_cplex_arrays()
	{
		// Master problem (dual_var_master & solution_master) are created
		// in 'solve_masterproblem()' because they change in dimension each iteration

		// Subproblems (solution_subproblem)
		numcols_subproblem = CPXgetnumcols(env, subproblems[0]);
		solution_subproblem = new double[numcols_subproblem];

		// In 'change_coefficients_subproblem()' (indices_sub & values_coef)
		indices_sub = new int[nb_tasks*nb_days*nb_shifts];
		for (int i = 0; i < nb_tasks*nb_days*nb_shifts; ++i)
			indices_sub[i] = i;

		values_coef = new double[nb_tasks*nb_days*nb_shifts];
	}



	void Diving_Column_Generation::clear_cplex()
	{
		// Free the CPLEX arrays
		delete[] solution_master;
		delete[] dual_var_master;
		delete[] values_coef;
		delete[] indices_sub;
		delete[] solution_subproblem;

		// Free the problems
		status = CPXfreeprob(env, &masterproblem);
		for (int p = 0; p < nb_people; ++p)
			status = CPXfreeprob(env, &subproblems[p]);

		// Close the cplex environment
		status = CPXcloseCPLEX(&env);
	}





	// DIVING
	void Diving_Column_Generation::diving_heuristic()
	{
		int level_tree = 0;
		double lowerbound = 0;

		// initialize
		solve_subproblem_p.reserve(nb_people);
		for (int i = 0; i < nb_people; ++i)
			solve_subproblem_p.push_back(true);


		// dive until integer solution
		while (true)
		{
			// (1) do column generation
			lowerbound = column_generation();

			// (2) check if solution is fractional, if integer save solution
			if (!is_solution_fractional()) {
				fix_residual_schedule();
				save_solution();
				break;
			}

			// (3) else dive
			++level_tree;
			if (print_diving_info)
				std::cout << "\n\n\nLevel tree = " << level_tree;

			// (4) find variable(s) (= person/people) to fix
			find_branching_variables(level_tree);

			// (5) set these variables to one in cplex master
			add_branching_restrictions(level_tree);

			// (6) remove columns z_pk' from current_columns and put them in removed_columns
			remove_columns_that_violate_branching_restrictions(level_tree);

			// (7) set solve_subproblem_p to correct value for new column generation
			for (int p = 0; p < nb_people; ++p) {
				solve_subproblem_p[p] = true;
				for (auto&& pp : fixed_people) {
					if (p == pp) {
						solve_subproblem_p[p] = false;
						break;
					}
				}
			}

			// (8) check remaining time per node for diving
			elapsed_time_total = (double)(clock() - start_time_total_diving) / CLOCKS_PER_SEC;
			allowed_computation_time_node = (allowed_computation_time - elapsed_time_total)
				/ (nb_people - fixed_people.size());
			std::cout << "\nElapsed time for this dataset: " << elapsed_time_total;
		};
	}



	void Diving_Column_Generation::fix_residual_schedule()
	{
		char *cur_colname[1];
		char cur_colnamestore[200];
		int surplus;
		std::string column_name;
		int column_index;

		fixed_columns.clear();

		for (int i = 2 * nb_tasks*nb_days*nb_shifts + nb_people; i < CPXgetnumcols(env, masterproblem); ++i)
		{
			if (solution_master[i] >= 1 - fractionality_tolerance_master)
			{
				status = CPXgetcolname(env, masterproblem, cur_colname, cur_colnamestore, 200, &surplus, i, i);

				column_name.assign(cur_colnamestore);

				for (int j = 0; j < current_columns.size(); ++j) {
					if (current_columns[j].name == column_name) {
						column_index = j;
						break;
					}
				}

				fixed_columns.push_back(current_columns[column_index]);
			}
		}
	}



	bool Diving_Column_Generation::is_solution_fractional()
	{
		for (int i = 0; i < CPXgetnumcols(env, masterproblem); ++i) {
			if (solution_master[i] > fractionality_tolerance_master && solution_master[i] < 1 - fractionality_tolerance_master) {
				return true;
			}
		}

		return false;
	}



	void Diving_Column_Generation::find_branching_variables(int level)
	{
		char *cur_colname[1];
		char cur_colnamestore[200];
		int surplus;
		std::string column_name;
		int cplex_index, column_index;
		int branching_person;
		double z_val_current_highest_frac_col = 0;


		// reset
		fixed_indices_cplex.clear();
		nb_people_fixed_current_iteration = 0;


		// find branching variable
		// method 1: largest fractional variable
		if (branching_method_diving == Branching_Method_Diving::largest_fractional_variable)
		{
			for (int i = 2 * nb_tasks*nb_days*nb_shifts + nb_people; i < CPXgetnumcols(env, masterproblem); ++i)
			{
				if (solution_master[i] < 1 - fractionality_tolerance_master && solution_master[i] > z_val_current_highest_frac_col)
				{
					status = CPXgetcolname(env, masterproblem, cur_colname, cur_colnamestore, 200, &surplus, i, i);

					column_name.assign(cur_colnamestore);
					cplex_index = i;

					for (int j = 0; j < current_columns.size(); ++j) {
						if (current_columns[j].name == column_name) {
							branching_person = current_columns[j].person;
							column_index = j;
							break;
						}
					}

					z_val_current_highest_frac_col = solution_master[i];
				}
			}

			fixed_indices_cplex.push_back(cplex_index);
			fixed_people.push_back(branching_person);
			++nb_people_fixed_current_iteration;
			fixed_columns.push_back(current_columns[column_index]);
		}

		// method 2: value above a certain threshold
		else /*if(branching_method_diving == Branching_Method_Diving::value_above_threshold)*/
		{
			int variables_fixed = 0;

			for (int i = 2 * nb_tasks*nb_days*nb_shifts + nb_people; i < CPXgetnumcols(env, masterproblem); ++i)
			{
				if (solution_master[i] < 1 - fractionality_tolerance_master && solution_master[i] > branching_threshold_diving)
				{
					status = CPXgetcolname(env, masterproblem, cur_colname, cur_colnamestore, 200, &surplus, i, i);

					column_name.assign(cur_colnamestore);
					cplex_index = i;

					for (int j = 0; j < current_columns.size(); ++j) {
						if (current_columns[j].name == column_name) {
							column_index = j;
							branching_person = current_columns[j].person;
							break;
						}
					}

					fixed_indices_cplex.push_back(cplex_index);
					fixed_people.push_back(branching_person);
					++nb_people_fixed_current_iteration;
					fixed_columns.push_back(current_columns[column_index]);

					++variables_fixed;
				}
			}

			// if no branching variables have been found, branch on value with highest fractional value
			if (variables_fixed == 0)
			{
				for (int i = 2 * nb_tasks*nb_days*nb_shifts + nb_people; i < CPXgetnumcols(env, masterproblem); ++i)
				{
					if (solution_master[i] < 1 - fractionality_tolerance_master && solution_master[i] > z_val_current_highest_frac_col)
					{
						status = CPXgetcolname(env, masterproblem, cur_colname, cur_colnamestore, 200, &surplus, i, i);

						column_name.assign(cur_colnamestore);
						cplex_index = i;

						for (int j = 0; j < current_columns.size(); ++j) {
							if (current_columns[j].name == column_name) {
								branching_person = current_columns[j].person;
								column_index = j;
								break;
							}
						}

						z_val_current_highest_frac_col = solution_master[i];
					}
				}

				fixed_indices_cplex.push_back(cplex_index);
				fixed_people.push_back(branching_person);
				++nb_people_fixed_current_iteration;
				fixed_columns.push_back(current_columns[column_index]);
			}
		}




		if (print_diving_info) {
			std::cout << "\nColumns fixed for person: ";
			for (auto&& i : fixed_people)
				std::cout << i << ", ";
		}
	}



	void Diving_Column_Generation::add_branching_restrictions(int level)
	{
		// Local variables
		double rhs[1];					// Right-hand side constraints
		char *rowname[1];				// Constraint names
		char jchar[6], name[200];		// Name storage
		char sense[1];					// Sign of constraint
		int matbeg[1];					// Begin position of the constraint
		int matind[1];					// Position of each element in constraint matrix
		double matval[1];				// Value of each element in constraint matrix
		int f{ 0 };						// To calculate number of nonzero coefficients in each constraint

		for (auto i : fixed_indices_cplex)
		{
			++nb_branching_dec;

			rhs[0] = 1; // one
			sense[0] = 'E';

			std::string name = "Branching_Restriction_" + std::to_string(nb_branching_dec);
			rowname[0] = const_cast<char*>(name.c_str());

			matbeg[0] = 0;
			f = 0;

			matind[f] = i;
			matval[f] = 1;
			++f;

			status = CPXaddrows(env, masterproblem, 0, 1, f, rhs, sense, matbeg, matind, matval, NULL, rowname);
		}


		// Write to the problem to a file
		if (write_to_file_each_branch)
			status = CPXwriteprob(env, masterproblem, "INEM_masterproblem.lp", NULL);
	}



	void Diving_Column_Generation::remove_columns_that_violate_branching_restrictions(int level)
	{
		//std::vector<int> indices_to_remove;

		for (int i = 0; i < current_columns.size(); ++i)
		{
			for (int j = 0; j < nb_people_fixed_current_iteration; ++j)
				//for (auto&& j: fixed_columns)
			{
				if (current_columns[i].person == fixed_columns[fixed_columns.size() - 1 - j].person && current_columns[i].name != fixed_columns[fixed_columns.size() - 1 - j].name)
				{
					// put in removed vector
					//removed_columns.push_back(current_columns[i]);
					//removed_columns.back().level_removed = level;

					// delete from cplex
					int col_index = -1;
					status = CPXgetcolindex(env, masterproblem, current_columns[i].name.data(), &col_index);

					if (col_index != -1) {
						status = CPXdelcols(env, masterproblem, col_index, col_index);
					}

					// which columns will need to be removed
					//indices_to_remove.push_back(i);
				}
			}
		}

		//for (int i = indices_to_remove.size() - 1; i >= 0; --i)
		//	current_columns.erase(current_columns.begin() + indices_to_remove[i]);


		// Write to the problem to a file
		if (write_to_file_each_iteration)
			status = CPXwriteprob(env, masterproblem, "INEM_master_problem.lp", NULL);
	}


}