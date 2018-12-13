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

#ifndef SOLUTION_H
#define SOLUTION_H


#include <vector>


namespace VNDS
{
	struct Task_Assignment
	{
		int task;
		int day;
		int shift;
	};

	class Solution
	{
		std::vector<std::vector<Task_Assignment>> people_task_assignments;

	public:
		Solution();
		bool at(int person, int task, int day, int shift);
		void set_at(int person, int task, int day, int shift, bool value);
		void reset();
		void reset_person(int person);
		void save(const Solution& othersol);
		const std::vector<Task_Assignment>& get_tasks_person(int person);

		double objective{ 1e20 };
	};

}


#endif // !SOLUTION_H

