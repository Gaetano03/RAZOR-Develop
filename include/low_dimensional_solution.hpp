/* ------------------------------------------------------------------------------
\file low_dimensional_solution.hpp
* \brief Functions to reconstruct a solution at untried/untested conditions.
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
#ifndef low_dimensional_solution_hpp
#define low_dimensional_solution_hpp
//
#include "Surrogates/rbf.h"
//
using namespace smartuq::surrogate;
//
//------------------------------------------------------------------------------            
// Functions
//------------------------------------------------------------------------------            
//
//!< \brief Seelction of RBF method
smartuq::surrogate::RBF_FUNCTION get_key_rbf ( const std::string &key_string ); 
//
//!< \brief Compute low dimensional solution(s)
int compute_low_dimensional_solution ( const std::string filedata );
//
#endif

