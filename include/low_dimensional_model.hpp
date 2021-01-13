/* ------------------------------------------------------------------------------
\file generate_low_dimensional_model.hpp
* \brief Subroutines and functions to generate a low dimensional model.
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
#ifndef low_dimensional_model_hpp
#define low_dimensional_model_hpp
//
#include "read_data.hpp"
#include "modal_identification.hpp"
#include "acquire_training_set.hpp"
#include "write_data.hpp"
#include <eigen3-hdf5.hpp>
//
//------------------------------------------------------------------------------
//  Data structure
//------------------------------------------------------------------------------
//
enum lowdim_strat { POD, SPOD, DMD, RDMD, MRDMD, ISOMAP };
//
//------------------------------------------------------------------------------            
// Functions
//------------------------------------------------------------------------------            
//
//!< \brief Compare keywords with input string
lowdim_strat read_ldmethod ( const std::string &ldstrat_string );
//
//!< \brief Generate low dimensional model and store it into a database
int generate_low_dimensional_model ( const std::string filedata );
//
//!< \brief Save low dimensional model info into a database 
void save_lowdim_model_info ( ld_model_data &lowdim_data, training_set_data &training_data, 
	ld_error_data &error_data );
//
#endif