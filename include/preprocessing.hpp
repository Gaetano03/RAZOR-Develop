/* ------------------------------------------------------------------------------
\file preprocessing.hpp
* \brief Pre-processing snapshot and mesh data, saving to HDF5 file
*
* Copyright 2016-2021, Aerospace Centre of Excellence University of Strathclyde
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
#ifndef preprocessing_hpp
#define preprocessing_hpp
//
#include "read_data.hpp"
#include "acquire_training_set.hpp"
#include "write_data.hpp"
//
//------------------------------------------------------------------------------            
// Functions
//------------------------------------------------------------------------------            
//
//!< \brief 
int preprocess_data ( const std::string filedata );
//
//!< \brief Create H5 file and save fields' and trainingset info as attributes
void save_attributes ( flds_data &fields_data, training_set_data &training_data );
//
//!< \brief Save snapshots into database (created by function above)
void save_snap_matrix ( Eigen::MatrixXd &snap_matrix, flds_data &fields_data, int &field_id );
//
//!< \brief Read and save mesh into h5 database
Eigen::VectorXi save_mesh (flds_data &fields_data, training_set_data &training_data );                                      
//
#endif
