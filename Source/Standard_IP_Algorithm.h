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

#ifndef STANDARD_IP_ALGORITHM_H
#define STANDARD_IP_ALGORITHM_H


#include "ilcplex\cplex.h"


namespace IP
{
	// DEFINITIONS
	namespace Algorithm_Type {
		enum {
			LP,
			IP
		};
	}

	
	// USER CHOICE
	extern int algorithm_choice;
	extern double allowed_computation_time;	// in seconds


	// CLASS DEFINITION ALGORITHM
	class Standard_IP_Algorithm
	{
		// CPLEX VARIABLES
		CPXENVptr env{ nullptr };			// Cplex environment
		CPXLPptr problem{ nullptr };		// Cplex problem pointer
		int status{ 0 };					// Cplex status parameter for output
		int solstat{ 0 };					// Solution status problem
		int numcols;						// Number of columns (variables) in the problem
		double objective;					// Objective value of the problem's solution
		double* solution{ nullptr };		// Solution of the problem


	public:
		void run_algorithm();

	private:
		void initialize_cplex();
		void build_problem();
		void solve_problem_as_IP();
		void solve_problem_as_LP();
		void clear_cplex();
	};

}


#endif // !STANDARD_IP_ALGORITHM_H

