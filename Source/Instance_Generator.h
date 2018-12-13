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

#ifndef INSTANCE_GENERATOR_H
#define INSTANCE_GENERATOR_H


#include <string>


namespace Instance_Generator
{
	extern void generate_dataset(const std::string& instance_name);
}


#endif // !INSTANCE_GENERATOR_H

