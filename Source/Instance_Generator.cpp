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

#include "Instance_Generator.h"
#include <iostream>
#include <fstream>
#include <random>

/*
*   This code generates a random data set for the INEM scheduling problem.
*	The number of people is randomly chosen from a uniform distribution in(min_nb_people, max_nb_people).
*	The number of groups is randomly chosen from a uniform distribution in(min_nb_groups, max_nb_groups).
*	The number of tasks is randomly chosen from a uniform distribution in(min_nb_tasks, max_nb_tasks).
*	The number of days is 28, with a probability of 'prob_28_days' and 56 days with a probability of 1 - prob_28_days.
*
*	People and tasks are randomly assigned to groups.
*
*   The skills coverage is randomly chosen from a uniform distribution in(min_skills_coverage, max_skills_coverage).
*	Then each person is assigned the skills for each respective task with the chosen skills coverage probability.
*
*   Task-day-shift requirements are calculated based on the number of people and their skills
*   as well as the demand coverage.This demand coverage is the ratio of demand to the supply of people and their skills.
*   It is randomly chosen from a uniform distribution in (min_demand_coverage, max_demand_coverage).
*/

namespace Instance_Generator
{
	// DEFINITIONS
	constexpr int min_nb_people = 250;
	constexpr int max_nb_people = 350;
	constexpr int min_nb_groups = 5;
	constexpr int max_nb_groups = 15;
	constexpr int min_nb_tasks = 40;
	constexpr int max_nb_tasks = 80;
	constexpr double prob_28_days = 0.8;
	constexpr int min_nb_holidays = 0;
	constexpr int max_nb_holidays = 2;
	constexpr int min_start_day = 0;
	constexpr int max_start_day = 6;
	constexpr double min_skills_coverage = 0.2;
	constexpr double max_skills_coverage = 1.0;
	constexpr double min_demand_coverage = 0.6;
	constexpr double max_demand_coverage = 1.2;
	constexpr int max_nb_variables = 3000000;		// x_ptds variables, i.e. nb_people*nb_tasks*nb_days*nb_shifts


	// RANDOM NUMBER GENERATOR
	std::random_device randdev;
	std::seed_seq seed{ randdev(), randdev(), randdev(), randdev(), randdev(), randdev(), randdev(), randdev() };
	std::mt19937 generator(seed);


	// DATA
	constexpr int nb_shifts = 3;
	int nb_people;
	int nb_groups;
	int nb_groups_CODU;
	int nb_tasks;
	int nb_tasks_CODU;
	int nb_days;
	int nb_holidays;
	int start_day;	// the first day of the planning horizon (0 == Monday, 1 == Tuesday, ..., 6 == Sunday)

	std::vector<bool> people_group;			// [p][g] == true if person p belongs to group g
	std::vector<bool> people_task;			// [p][t] == true if person p can to task t
	std::vector<bool> group_task;			// [g][t] == true if task t belongs to group g
	std::vector<int> shift_demands;			// [t][d][s] == number of people required for task t on shift s of day d
	std::vector<int> task_durations;		// [t] == duration (in hours) of task t (8 or 12 hours)



	void generate_basic_data()
	{
		std::uniform_int_distribution<int> dist_people(min_nb_people, max_nb_people);
		std::uniform_int_distribution<int> dist_groups(min_nb_groups, max_nb_groups);
		std::uniform_int_distribution<int> dist_tasks(min_nb_tasks, max_nb_tasks);
		std::bernoulli_distribution dist_days(prob_28_days);
		std::uniform_int_distribution<int> dist_holidays(min_nb_holidays, max_nb_holidays);
		std::uniform_int_distribution<int> dist_startday(min_start_day, max_start_day);

		// first generate number of days
		if (dist_days(generator))
			nb_days = 28;
		else
			nb_days = 56;

		do {
			nb_people = dist_people(generator);

			nb_groups = dist_groups(generator);
			nb_groups_CODU = nb_groups / 4;

			nb_tasks = dist_tasks(generator);
			nb_tasks_CODU = nb_tasks / 4;
		} while (nb_people*nb_tasks*nb_days*nb_shifts > max_nb_variables);

		nb_holidays = dist_holidays(generator);
		start_day = dist_startday(generator);
	}



	void assign_people_group()
	{
		std::uniform_int_distribution<int> dist(0, nb_groups - 1);

		people_group.reserve(nb_people*nb_groups);
		for (int p = 0; p < nb_people; ++p)
		{
			int group = dist(generator);
			for (int g = 0; g < nb_groups; ++g)
			{
				if (g == group)
					people_group.push_back(true);
				else
					people_group.push_back(false);
			}
		}
	}



	void assign_people_task()
	{
		std::uniform_real_distribution<double> dist_skills_coverage(min_skills_coverage, max_skills_coverage);
		double skills_coverage = dist_skills_coverage(generator);
		std::bernoulli_distribution dist_skills(skills_coverage);

		people_task.reserve(nb_people*nb_tasks);
		for (int p = 0; p < nb_people; ++p)
		{
			for (int t = 0; t < nb_tasks; ++t)
			{
				if (dist_skills(generator))
					people_task.push_back(true);
				else
					people_task.push_back(false);
			}
		}
	}



	void assign_group_task()
	{
		std::uniform_int_distribution<int> dist(0, nb_groups - 1);

		std::vector<bool> task_group;
		task_group.reserve(nb_tasks*nb_groups);
		group_task.reserve(nb_groups*nb_tasks);

		for (int t = 0; t < nb_tasks; ++t)
		{
			int group = dist(generator);
			for (int g = 0; g < nb_groups; ++g)
			{
				// task group
				if (g == group)
					task_group.push_back(true);
				else
					task_group.push_back(false);

				// group task
				group_task.push_back(false);
			}
		}


		for (int g = 0; g < nb_groups; ++g)
		{
			for (int t = 0; t < nb_tasks; ++t)
			{
				group_task[g*nb_tasks + t] = task_group[t*nb_groups + g];
			}
		}

	}



	void assign_shift_demands()
	{
		std::uniform_real_distribution<double> dist_coverage(min_demand_coverage, max_demand_coverage);
		double demand_coverage = dist_coverage(generator);

		// Each person contributes 1/total_tasks_person_p to task t if they have the required skill, 0 otherwise (ex post)
		// Furthermore people can only work 17.5/28 days and 1/3 shifts
		// Finally we adjust for the demand_coverage (utilisation, i.e. demand/supply) we want
		// Then based on the available supply we generate a demand from a Poisson distribution ~ P(supply)

		// 1. count number of tasks for every person p
		std::vector<int> number_tasks_person_p;
		number_tasks_person_p.reserve(nb_people);
		for (int p = 0; p < nb_people; ++p)
		{
			number_tasks_person_p.push_back(0);

			for (int t = 0; t < nb_tasks; ++t)
			{
				if (people_task[p*nb_tasks + t])
					++number_tasks_person_p[p];
			}
		}

		// 2. for every task, count the available supply (contribution of each person)
		//    then adjust for contract hours (17.5/28), one shift per day (1/3), and demand_coverage
		//	  finally generate demands for each day and each shift
		for (int t = 0; t < nb_tasks; ++t)
		{
			double supply_task_t = 0; // per day per shift for task t
			for (int p = 0; p < nb_people; ++p)
			{
				if (people_task[p*nb_tasks + t])
					supply_task_t += 1.0 / number_tasks_person_p[p];
			}
			supply_task_t *= (17.5 / 28.0) * (1.0 / 3.0) * demand_coverage;


			// now we generate a demand from a Poisson distribution ~ P(supply)
			std::poisson_distribution<int> dist_demand(supply_task_t);
			for (int d = 0; d < nb_days; ++d)
			{
				for (int s = 0; s < nb_shifts; ++s)
				{
					int demand = dist_demand(generator);
					shift_demands.push_back(demand);
				}
			}
		}
	}



	void assign_task_durations()
	{
		task_durations.reserve(nb_tasks);
		for (int t = 0; t < nb_tasks; ++t)
			task_durations.push_back(8);
	}



	void clear_data()
	{
		people_group.clear();
		people_task.clear();
		group_task.clear();
		shift_demands.clear();
		task_durations.clear();
	}



	void write_to_file(const std::string& file_name)
	{
		std::ofstream file;
		file.open(file_name);
		if (!file.is_open()) {
			std::cerr << "\nError: couldn't open output file.\n\n";
			return;
		}
		else
		{
			// basic data
			file << nb_people;
			file << "\n" << nb_groups;
			file << "\n" << nb_groups_CODU;
			file << "\n" << nb_tasks;
			file << "\n" << nb_tasks_CODU;
			file << "\n" << nb_days;
			file << "\n" << nb_holidays;
			file << "\n" << start_day;

			// people_group;	
			file << "\n";
			for (int p = 0; p < nb_people; ++p)
			{
				file << "\n";
				for (int g = 0; g < nb_groups; ++g)
				{
					if (g != 0)
						file << "\t";
					file << people_group[p*nb_groups + g];
				}
			}

			// people_task;	
			file << "\n";
			for (int p = 0; p < nb_people; ++p)
			{
				file << "\n";
				for (int t = 0; t < nb_tasks; ++t)
				{
					if (t != 0)
						file << "\t";
					file << people_task[p*nb_tasks + t];
				}
			}

			// group_task;	
			file << "\n";
			for (int g = 0; g < nb_groups; ++g)
			{
				file << "\n";
				for (int t = 0; t < nb_tasks; ++t)
				{
					if (t != 0)
						file << "\t";
					file << group_task[g*nb_tasks + t];
				}
			}

			// shift_demands;	
			// [t][d][s] == number of people required for task t on shift s of day d
			file << "\n";
			for (int t = 0; t < nb_tasks; ++t)
			{
				file << "\n";
				for (int d = 0; d < nb_days; ++d)
				{
					for (int s = 0; s < nb_shifts; ++s)
					{
						if (!(d == 0 && s == 0))
							file << "\t";
						file << shift_demands[t*nb_days*nb_shifts + d * nb_shifts + s];
					}
				}
			}

			// task_durations;		
			file << "\n\n";
			for (int t = 0; t < nb_tasks; ++t)
			{
				if (t != 0)
					file << "\t";
				file << task_durations[t];
			}

			file.close();
		}
	}


	void generate_dataset(const std::string& instance_name)
	{
		generate_basic_data();
		assign_people_group();
		assign_people_task();
		assign_group_task();
		assign_shift_demands();
		assign_task_durations();

		write_to_file(instance_name);

		clear_data();
	}
}