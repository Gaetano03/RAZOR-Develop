/* ------------------------------------------------------------------------------
\file low_dimensional_solution.cpp
* \brief Subroutines and functions to generate a low dimensional solution.
*
* Copyright 2016-2020, Aerospace Centre of Excellence University of Strathclyde
*
* RAZOR is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* RAZOR is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with RAZOR. If not, see <http://www.gnu.org/licenses/>.
* ------------------------------------------------------------------------------*/
//
#include "low_dimensional_solution.hpp"
#include "read_data.hpp"
//
// -------------------------------------------------------------------------------------------
int compute_low_dimensional_solution ( const std::string filedata ) {
// -------------------------------------------------------------------------------------------
//
//  Variables declaration
    target_set_data target_data;
    aero_data aerodynamic_data;
//
//  Reading configuration file
    read_soldata ( filedata, target_data, aerodynamic_data);        
    //print_soldata ( lowdim_data, training_data, error_data, 
    //    modal_data, manifold_data ); 
//
    return 0;
}
//
//
// -------------------------------------------------------------------------------------------
smartuq::surrogate::RBF_FUNCTION get_key_rbf ( const std::string &key_string ) {
// -------------------------------------------------------------------------------------------
//
    if ( key_string == "LINEAR" )
        return smartuq::surrogate::LINEAR;
    else if ( key_string == "CUBIC" )
        return smartuq::surrogate::CUBIC;
    else if ( key_string == "GAUSSIAN" )
        return smartuq::surrogate::GAUSSIAN;
    else if ( key_string == "THIN_PLATE" )
        return smartuq::surrogate::THIN_PLATE;
    else if ( key_string == "MULTIQUADRATICS" )
        return smartuq::surrogate::MULTIQUADRATICS;
//
    else {
//
        std::cout << " " << std::endl;
        std::cout << "-> Error: Unknown REDUCTION STRATEGY type" << std::endl;
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//     
    }
//
}