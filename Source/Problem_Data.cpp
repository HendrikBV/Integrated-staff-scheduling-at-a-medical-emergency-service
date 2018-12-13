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
#include <fstream>
#include <iostream>


// PROBLEM DATA
const int nb_shifts{ 3 };
int nb_people;
int nb_tasks;
int nb_tasks_CODU;
int nb_days;
int start_day;	// the first day of the planning horizon (0 == Monday, 1 == Tuesday, ..., 6 == Sunday)
int nb_groups;
int nb_groups_CODU;
int nb_weekends;
int nb_holidays;

std::vector<bool> g_people_group;			// [p][g] == true if person p belongs to group g
std::vector<bool> g_people_task;			// [p][t] == true if person p can to task t
std::vector<bool> g_group_task;				// [g][t] == true if task t belongs to group g
std::vector<int> g_shift_demands;			// [t][d][s] == number of people required for task t on shift s of day d
std::vector<int> g_task_durations;			// [t] == duration (in hours) of task t (8 or 12 hours)


// OBJECTIVE FUNCTION WEIGHTS
int obj_weight_Y_REplus{ 10 };			// objective function weight for overstaffing
int obj_weight_Y_REmin_CODU{ 100 };		// objective function weight for understaffing in CODU
int obj_weight_Y_REmin_amb{ 1000 };		// objective function weight for understaffing in ambulances
int obj_weight_Y_Hplus{ 1 };			// objective function weight for overtime of a person
int obj_weight_Y_Hmin{ 1 };				// objective function weight for undertime of a person
int obj_weight_Y_W{ 10 };				// objective function weight for not full weekends off
int obj_weight_Y_G_CODU{ 10 };			// objective function weight for assigning tasks to other groups for CODU
int obj_weight_Y_G_ambulances{ 20 };	// objective function weight for assigning tasks to other groups for ambulances


// AUXILIARY FUNCTIONS
bool person_group(int person, int group) { return g_people_group[person*nb_groups + group]; }
bool person_task(int person, int task) { return g_people_task[person*nb_tasks + task]; }
bool group_task(int group, int task) { return g_group_task[group*nb_tasks + task]; }
int shift_demands(int task, int day, int shift) { return g_shift_demands[task*nb_days*nb_shifts + day * nb_shifts + shift]; }


// INPUT
bool input(const std::string& file_name)
{
	std::ifstream myfile;
	myfile.open(file_name);
	if (!myfile.is_open()) {
		std::cout << "\nError: couldn't open input file.\n\n";
		return false;
	}
	else
	{
		int value;

		// basic data
		myfile >> nb_people;
		myfile >> nb_groups;
		myfile >> nb_groups_CODU;
		myfile >> nb_tasks;
		myfile >> nb_tasks_CODU;
		myfile >> nb_days;
		myfile >> nb_holidays;
		myfile >> start_day;



		// people_group
		g_people_group.reserve(nb_people*nb_groups);
		for (int p = 0; p < nb_people; ++p) {
			for (int g = 0; g < nb_groups; ++g) {
				myfile >> value;
				g_people_group.push_back(value);
			}
		}

		// people_task
		g_people_task.reserve(nb_people*nb_tasks);
		for (int p = 0; p < nb_people; ++p) {
			for (int t = 0; t < nb_tasks; ++t) {
				myfile >> value;
				g_people_task.push_back(value);
			}
		}

		// group_task
		g_group_task.reserve(nb_groups*nb_tasks);
		for (int g = 0; g < nb_groups; ++g) {
			for (int t = 0; t < nb_tasks; ++t) {
				myfile >> value;
				g_group_task.push_back(value);
			}
		}

		// requirements
		g_shift_demands.reserve(nb_tasks*nb_days*nb_shifts);
		for (int t = 0; t < nb_tasks; ++t) {
			for (int d = 0; d < nb_days; ++d) {
				for (int s = 0; s < nb_shifts; ++s) {
					myfile >> value;
					g_shift_demands.push_back(value);
				}
			}
		}

		// task durations
		g_task_durations.reserve(nb_tasks);
		for (int t = 0; t < nb_tasks; ++t) {
			myfile >> value;
			g_task_durations.push_back(value);
		}


		// CALCULATIONS
		// Calculate nb_weekends
		nb_weekends = 0;
		for (int d = 6 - start_day; d < nb_days; d = d + 7)
			++nb_weekends;

		return true;
	}
}

void clear_data()
{
	// reset instance data
	nb_people = 0;
	nb_tasks = 0;
	nb_tasks_CODU = 0;
	nb_days = 0;
	start_day = 0;
	nb_groups = 0;
	nb_groups_CODU = 0;
	nb_weekends = 0;

	g_group_task.clear();
	g_people_group.clear();
	g_people_task.clear();
	g_shift_demands.clear();
	g_task_durations.clear();
}


