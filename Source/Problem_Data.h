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

#ifndef PROBLEM_DATA_H
#define PROBLEM_DATA_H


#include <string>
#include <vector>


// DATA DEFINITIONS
namespace Day {
	enum {
		Monday,
		Tuesday,
		Wednesday,
		Thursday,
		Friday,
		Saturday,
		Sunday
	};
}

namespace Shift {
	enum {
		night,
		morning,
		afternoon
	};
}


// PROBLEM DATA
extern const int nb_shifts;
extern int nb_people;
extern int nb_tasks;
extern int nb_tasks_CODU;
extern int nb_days;
extern int start_day;	// the first day of the planning horizon (0 == Monday, 1 == Tuesday, ..., 6 == Sunday)
extern int nb_groups;
extern int nb_groups_CODU;
extern int nb_weekends;
extern int nb_holidays;

extern std::vector<bool> g_people_group;		// [p][g] == true if person p belongs to group g
extern std::vector<bool> g_people_task;			// [p][t] == true if person p can to task t
extern std::vector<bool> g_group_task;			// [g][t] == true if task t belongs to group g
extern std::vector<int> g_shift_demands;		// [t][d][s] == number of people required for task t on shift s of day d
extern std::vector<int> g_task_durations;		// [t] == duration (in hours) of task t (8 or 12 hours)


// OBJECTIVE FUNCTION WEIGHTS
extern int obj_weight_Y_REplus;			// objective function weight for overstaffing
extern int obj_weight_Y_REmin_CODU;		// objective function weight for understaffing in CODU
extern int obj_weight_Y_REmin_amb;		// objective function weight for understaffing in ambulances
extern int obj_weight_Y_Hplus;			// objective function weight for overtime of a person
extern int obj_weight_Y_Hmin;			// objective function weight for undertime of a person
extern int obj_weight_Y_W;				// objective function weight for not full weekends off
extern int obj_weight_Y_G_CODU;			// objective function weight for assigning tasks to other groups for CODU
extern int obj_weight_Y_G_ambulances;	// objective function weight for assigning tasks to other groups for ambulances


// AUXILIARY FUNCTIONS
extern bool person_group(int person, int group);
extern bool person_task(int person, int task);
extern bool group_task(int group, int task);
extern int shift_demands(int task, int day, int shift);


// INPUT
extern bool input(const std::string& file_name);
extern void clear_data();


#endif // !PROBLEM_DATA_H

