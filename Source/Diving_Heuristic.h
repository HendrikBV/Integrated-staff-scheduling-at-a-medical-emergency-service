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

#ifndef DIVING_HEURISTIC_H
#define DIVING_HEURISTIC_H


#include <string>
#include <vector>
#include <list>
#include <ctime>
#include "ilcplex\cplex.h"


namespace Diving_Heuristic
{
	// DEFINITIONS
	namespace Column_Generation_Method {
		enum {
			one_column_per_person,					// solve the master and the subproblems for every person, then add all columns with negative reduced costs
			one_column_person_p_and_reoptimize,		// solve the master and the subproblem for person p, then add the column if RC < 0, reoptimize the master and go to person p+1
		};
	}

	namespace Branching_Method_Diving {
		enum {
			largest_fractional_variable,
			value_above_threshold
		};
	}


	// USER CHOICE
	extern int column_generation_method;
	extern int branching_method_diving;
	extern double branching_threshold_diving;

	extern double allowed_computation_time;			// in seconds
	extern double allowed_computation_time_root;	// in seconds
	extern double allowed_computation_time_node;	// in seconds


	// CLASS DEFINITION ALGORITHM
	class Diving_Column_Generation
	{
		struct Column {
			int person{ -1 };			// for which person is this column z_pk
			std::string name;			// name of the column in CPLEX
			double cost;				// c_pk of the column
			std::vector<int> a_tds;		// the (t,d,s) coefficients of the column
			Column(const Column& other) { person = other.person; name = other.name; cost = other.cost; a_tds = other.a_tds; }
			Column() {}
		};


		// CPLEX VARIABLES
		CPXENVptr env{ nullptr };						// Cplex environment
		CPXLPptr masterproblem{ nullptr };				// Cplex masterproblem pointer
		std::vector<CPXLPptr> subproblems;				// Cplex pointers for the subproblems
		int status{ 0 };								// Cplex status parameter for output
		int solstat{ 0 };								// Solution status problem
		int numcols_master{ 0 };						// Number of columns (variables) in the masterproblem
		int numrows_master{ 0 };						// Number of rows (constraints) in the masterproblem
		int numcols_subproblem{ 0 };					// Number of columns in the subproblem
		double objective_master;						// Objective value of the masterproblem's solution
		double objective_subproblem;					// Objective value of the subproblem's solution
		double* solution_master{ nullptr };				// Solution master problem
		double* dual_var_master{ nullptr };				// Dual variables master problem's solution
		double* values_coef{ nullptr };					// Coefficients in the subproblem (-1 x dual var)
		int* indices_sub{ nullptr };					// Indices of the coefficients in the subproblem that will be changed
		double cost_column_k_person_p;					// C_pk for the newly found column for person p
		double mu_p;									// Dual variable for constraint set (2) for person p
		double* solution_subproblem{ nullptr };			// Solution subproblem
		const int obj_value_super_column{ 10000 };		// Objective function values for the super columns (startvariables)

		// COLUMNS
		std::vector<bool> solve_subproblem_p;			// does the pricing problem need to be solved for person p in the next iteration?
		std::vector<Column> current_columns;			// The columns currently in the master
		std::vector<Column> removed_columns;			// The columns removed from the master due to branching 
		std::vector<Column> fixed_columns;
		std::list<int> fixed_indices_cplex;
		std::list<int> fixed_people;
		int nb_people_fixed_current_iteration;

		// MISCELLANEOUS
		clock_t start_time_total_diving;
		bool first_call_CG{ true };
		int nb_calls_CG{ 0 };
		int nb_branching_dec{ 0 };

		// ALGORITHM STATUS AND SOLUTION
		bool algorithm_failed;
		std::string algorithm_solution_status;
		std::vector<int> solution;	// (p,d,s) = t 

		// STATISTICS
		int iterations_CG{ 0 };
		int iterations_CG_root_node{ 0 };
		int nb_columns_added_total{ 0 };
		int nb_columns_added_root_node{ 0 };
		double elapsed_time_CG_root_node{ 0 };
		double elapsed_time_CG_total{ 0 };
		double elapsed_time_CG_masterproblem{ 0 };
		double elapsed_time_CG_subproblem{ 0 };
		double elapsed_time_total{ 0 };
		double lowerbound_root_node;
		double objective_best_solution{ 1e20 };



	public:
		void run_algorithm();

	private:
		double column_generation();
		void build_master_problem();
		void build_subproblem_problem(int person);
		double solve_masterproblem();
		double solve_subproblem(int person);
		void change_coefficients_subproblem(int person);
		void add_column_to_master(int person);
		bool check_solution_status(bool integer, int status_, int solstat_);

		void initialize_cplex();
		void create_cplex_arrays();
		void clear_cplex();

		void fix_residual_schedule();
		void save_solution();
		void print_solution();

		void diving_heuristic();
		bool is_solution_fractional();
		void find_branching_variables(int level);
		void add_branching_restrictions(int level);
		void remove_columns_that_violate_branching_restrictions(int level);
	};
}


#endif // !DIVING_HEURISTIC_H

