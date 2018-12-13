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

#ifndef VNDS_H
#define VNDS_H


#include <vector>
#include "ilcplex\cplex.h"


namespace VNDS
{
	// USER CHOICE
	extern double allowed_computation_time;		// in seconds

	// MAIN FUNCTION
	extern void heuristic();



	class Initial_Solution
	{
		// CPLEX VARIABLES
		CPXENVptr env{ nullptr };			// Cplex environment
		CPXLPptr problem{ nullptr };		// Cplex problem pointer
		int status{ 0 };					// Cplex status parameter for output
		int solstat{ 0 };					// Solution status problem
		double objective;					// Objective value of the problem's solution
		double* solution{ nullptr };		// Solution of the problem

		std::vector<int> remaining_demands;

	public:
		void run_algorithm();

	private:
		void initialize_cplex();
		void build_problem(int person);
		void solve_problem(int person);
		void clear_cplex();
	};



	class Column_Person
	{
		// CPLEX VARIABLES
		CPXENVptr env{ nullptr };			// Cplex environment
		CPXLPptr problem{ nullptr };		// Cplex problem pointer
		int status{ 0 };					// Cplex status parameter for output
		int solstat{ 0 };					// Solution status problem
		double objective;					// Objective value of the problem's solution
		double* solution{ nullptr };		// Solution of the problem


	public:
		void find_column(int person);

	private:
		void initialize_cplex();
		void build_problem(int person);
		void solve_problem(int person);
		void clear_cplex();
	};



	class MIP_Heuristic
	{
		// CPLEX VARIABLES
		CPXENVptr env{ nullptr };				// Cplex environment
		CPXLPptr problem{ nullptr };			// Cplex problem pointer
		int status{ 0 };						// Cplex status parameter for output
		int solstat{ 0 };						// Solution status problem
		double candidate_objective{ 1e20 };		// Objective value of the candidate solution 
		double* solution{ nullptr };			// Solution of the problem
		int nb_main_constraints;

		// ALGORITHM PARAMETERS
		int time_limit_subproblem = 120;		// in seconds



	public:
		void run_algorithm();

	private:
		void initialize_cplex();
		void build_problem();
		void solve_problem();
		void clear_cplex();

		void set_initial_solution();

		void probabilistic_VNDS();

		bool neighbourhood_days(int first_day, int last_day);
		bool neighbourhood_tasks(int task);
		bool neighbourhood_shifts(int shift);
		void shaking(const std::vector<int>& people_to_shake);
	};


}


#endif // !VNDS_H

