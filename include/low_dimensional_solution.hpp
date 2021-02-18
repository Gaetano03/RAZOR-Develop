/* ------------------------------------------------------------------------------
\file low_dimensional_solution.hpp
* \brief Functions to reconstruct a solution at untried/untested conditions.
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
#ifndef low_dimensional_solution_hpp
#define low_dimensional_solution_hpp
//
#include "Surrogates/rbf.h"
#include "read_data.hpp"
//
using namespace smartuq::surrogate;
//

class CModalReconstruction {
//
//changed coeff eigens to Xd from Xcd
protected:
//
    Eigen::MatrixXcd m_Phi_flds; //!< \brief Matrix of flow features computed through a generic linear algorithm for all fields
    Eigen::MatrixXd m_Coefs_flds;
    Eigen::MatrixXcd m_Phi; //!< \brief Matrix of flow features computed through a generic linear algorithm for a specific field
    Eigen::MatrixXd m_Coefs; //!< \brief Training coefficients with each column representing evolution for a specific mode over the parameter space
    MatrixXd m_Rec_field; //!< \brief Matrix containing high-dimensional solution at each target point (per column)
    target_set_data m_target_data;
    training_set_data m_training_data;
    flds_data m_flds_data;
    std::vector<double> m_mean_params; //!< \brief vector of means (needed for feature scaling)
    std::vector<double> m_std_params; //!< \brief vector of standard deviations (needed for feature scaling)
    std::vector<std::vector<double> > m_train_pnts_scaled; //!< \brief scaled list of parameters for interpolation purposes
    std::vector<std::vector<double> > m_targ_pnts_scaled;
//
public:
//
//Constructor
    CModalReconstruction() = default;

    CModalReconstruction(target_set_data &target_data, training_set_data &training_data, flds_data &fields_data) :
        m_target_data(target_data), m_training_data(training_data), m_flds_data(fields_data)
    {
    }

//
    void features_scaling(); //!< \brief scale parameters if they are on different ranges
    virtual void load_lowdim(const std::string filename) {}; //!< \brief load infos about the low dimensional model
    virtual void compute_lowdim_surrogate(const int field_id) {}; //!< \brief Compute a low dimensional surrogate
    virtual void compute_lowdim_sol(const int field_id) {}; //!< \brief Compute a low dimensional solution
    virtual void compute_highdim_sol(const int field_id) {}; //!< \brief Project low dimensional solution back to the high-dim space
    void save_solutions();
//
};
//
//!< \brief Class for online computation of low dimensional solutions through interpolation
class CInterpolation : public CModalReconstruction {
//
protected:
//
    MatrixXd m_coefs_interp; //!< \brief Matrix containing low dimensional solution at each target point (columns represent lowdim sol at each target point)
    std::vector<std::vector<rbf> > surrogates; //!< \brief vector containing surrogates for each mode coefficient
//
public:
//
//Constructor
    CInterpolation( target_set_data &target_data, training_set_data &training_data, flds_data &fields_data ) :
        CModalReconstruction(target_data, training_data, fields_data)
        {
            surrogates = std::vector<std::vector<rbf> >(fields_data.n_flds);
        }
//
    void load_lowdim(const std::string filename);
    void compute_lowdim_surrogate(const int field_id);
    void compute_lowdim_sol(const int field_id);
    void compute_highdim_sol(const int field_id);
//
};
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

