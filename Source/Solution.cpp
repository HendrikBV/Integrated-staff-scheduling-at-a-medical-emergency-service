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

#include "Solution.h"
#include "Problem_Data.h"

namespace VNDS
{
	Solution::Solution()
	{
		people_task_assignments.reserve(nb_people);
		for (int p = 0; p < nb_people; ++p) {
			std::vector<Task_Assignment> vec;
			vec.reserve(nb_days*0.7);
			people_task_assignments.push_back(vec);
		}
	}

	bool Solution::at(int person, int task, int day, int shift)
	{
		for (auto&& ass : people_task_assignments[person]) {
			if (ass.task == task && ass.day == day && ass.shift == shift) {
				return true;
			}
		}

		return false;
	}

	void Solution::set_at(int person, int task, int day, int shift, bool value)
	{
		if (value)
		{
			people_task_assignments[person].push_back(Task_Assignment());
			people_task_assignments[person].back().task = task;
			people_task_assignments[person].back().day = day;
			people_task_assignments[person].back().shift = shift;
		}
	}

	void Solution::reset()
	{
		for (int p = 0; p < nb_people; ++p)
			people_task_assignments[p].clear();
	}

	void Solution::reset_person(int person)
	{
		people_task_assignments[person].clear();
	}

	void Solution::save(const Solution& othersol)
	{
		this->people_task_assignments = othersol.people_task_assignments;
		this->objective = othersol.objective;
	}

	const std::vector<Task_Assignment>& Solution::get_tasks_person(int person)
	{
		return people_task_assignments[person];
	}
}