/* ------------------------------------------------------------------------------
\file write_data.hpp
* \brief Functions to write information and data
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
#ifndef write_data_hpp
#define write_data_hpp
//
//#include "low_dimensional_solution.hpp"
#include "read_data.hpp"
//
//------------------------------------------------------------------------------            
// Data structure
//------------------------------------------------------------------------------            
//
const int CGNS_STRING_SIZE = 33;
//
//------------------------------------------------------------------------------            
// Functions                                                                     
//------------------------------------------------------------------------------            
//
//!< \brief Print to screen data needed for low-dimensional model generation
void print_gendata ( ld_model_data &lowdim_data, training_set_data &training_data, ld_error_data 
    &error_data, modal_identification_data &modal_data, manifold_learning_data &manifold_data );        //devnote: again I guess I will need to split print_gendata into two exactly as done with read_gendata
//
//!< \brief Print to screen data needed for low-dimensional solution computation
void print_soldata ( target_set_data &target_data, aero_data &aerodynamic_data ); 
//
/*
void write_modes_sPOD ( const Eigen::MatrixXd &Phi_cut, 
                        const Eigen::MatrixXd &Coords,
                        std::string flag_prob );

void write_modes_DMD ( const Eigen::MatrixXcd &Phi_cut,
                    const Eigen::MatrixXd &Coords, 
                    std::string flag_prob );

void write_alfa_lam_DMD( const Eigen::VectorXcd Alfa,
                        const Eigen::VectorXcd Lambda);

void write_coeffs_sPOD ( const Eigen::MatrixXd &Coeffs,
                        const std::vector<double> &t_vec,
                        const Eigen::VectorXd &lam );

void write_TimeDynamics_DMD ( const Eigen::VectorXcd omega,
                            const Eigen::VectorXcd alfa,
                            const Eigen::VectorXd t);

void write_CoefsDynamics_mrDMD( std::vector<node_mrDMD> &nodes, 
                                const int level, 
                                const int ns, 
                                const int max_levels );

void write_Reconstructed_fields ( const Eigen::MatrixXd Rec,
                                    const Eigen::MatrixXd &Coords,
                                    std::string filename,
                                    std::string flag_prob,
                                    const int nt );

void Write_Restart_Cons_Time ( const Eigen::MatrixXd &Rec,
                                    const Eigen::MatrixXd &Coords,
                                    std::string filename,
                                    const int nt,
                                    const int nC,
                                    const double alpha,
                                    const double beta,
                                    const std::string flag = "NO" );

void write_modes ( const Eigen::MatrixXd &Phi_cut ); //write modes RDMD

void write_coefs ( const Eigen::MatrixXd &Coefs ); //write coefs RDMD

void write_err_j ( const Eigen::MatrixXd data, std::string filename ); //write error/jaccard surface for RBM method

void Write_Plot3d_Modes( Eigen::MatrixXd Phi,           //Modes have to be scalar functions
                        std::string filename, 
                        plot3d_info Info );
*/
#endif
