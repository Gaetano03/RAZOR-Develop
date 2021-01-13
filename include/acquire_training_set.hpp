/* ------------------------------------------------------------------------------
\file generate_snset.hpp
* \brief Subroutines and functions to acquire a training set of observations.
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
#ifndef acquire_training_set_hpp
#define acquire_training_set_hpp
//
#include "read_data.hpp"
//
//------------------------------------------------------------------------------            
// Functions
//------------------------------------------------------------------------------            
//
//!< \brief Generate snapshot matrix
Eigen::MatrixXd generate_snap_matrix( const int fld_id, const int Ni, std::ifstream 
	filestream[], ld_model_data lowdim_data, training_set_data training_data );
//
#endif