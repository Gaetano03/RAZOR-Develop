/* ------------------------------------------------------------------------------
\file read_data.hpp
* \brief All the information about the definition of the physical problem.
*        The subroutines and functions are in the <i>read_Inputs.cpp</i> file.
*
* Copyright 2016-2020, Aerospace Centre of Excellence University of Strathclyde
*
* MODES is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* MODES is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with MODES. If not, see <http://www.gnu.org/licenses/>.
* ------------------------------------------------------------------------------*/
//
#ifndef read_data_hpp
#define read_data_hpp
//
#include <iostream>
#include <sstream>
#include <stdio.h> 
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
//
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
//
using namespace Eigen;
//
//------------------------------------------------------------------------------
//  Defining class and function for matrix indexing
//  Mimic matlab capability to select specific columns and rows from a matrix
//------------------------------------------------------------------------------
//
template<class ArgType, class RowIndexType, class ColIndexType>
class indexing_functor {
  const ArgType &m_arg;
  const RowIndexType &m_rowIndices;
  const ColIndexType &m_colIndices;
public:
  typedef Eigen::Matrix<typename ArgType::Scalar,
                 RowIndexType::SizeAtCompileTime,
                 ColIndexType::SizeAtCompileTime,
                 ArgType::Flags&RowMajorBit?RowMajor:ColMajor,
                 RowIndexType::MaxSizeAtCompileTime,
                 ColIndexType::MaxSizeAtCompileTime> MatrixType;
  indexing_functor(const ArgType& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
    : m_arg(arg), m_rowIndices(row_indices), m_colIndices(col_indices)
  {}
  const typename ArgType::Scalar& operator() (Eigen::Index row, Eigen::Index col) const {
    return m_arg(m_rowIndices[row], m_colIndices[col]);
  }
};
//
template <class ArgType, class RowIndexType, class ColIndexType>
CwiseNullaryOp<indexing_functor<ArgType,RowIndexType,ColIndexType>, typename indexing_functor<ArgType,RowIndexType,ColIndexType>::MatrixType>
indexing(const MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
{
  typedef indexing_functor<ArgType,RowIndexType,ColIndexType> Func;
  typedef typename Func::MatrixType MatrixType;
  return MatrixType::NullaryExpr(row_indices.size(), col_indices.size(), Func(arg.derived(), row_indices, col_indices));
}
//
//  Offset for reading fortran binary files
const int RECORD_DELIMITER_LENGTH = 4;
//
//------------------------------------------------------------------------------
//  Data structure
//------------------------------------------------------------------------------
//
struct training_set_data {
    int n_snap;                      //!< \brief Number of snapshots

};
//
struct modal_extraction_data {

};
//
struct manifold_learning_data {

};
//
struct adaptive_reduction_data {

};
//
struct target_set_data {
  
};
//
struct prob_settings {
//
//  Definition of general variables
    int n_targ;                      //!< \brief Number of target solution files in reconstruction
    int n_parm;                      //!< \brief Number of parameters
    int n_flds;                      //!< \brief Number of fields to process
    int n_eval_pnts_res;             //!< \brief Number of points where residual error will be evaluated during adaptive solution
    int n_eval_pnts_dir;             //!< \brief Number of points where direct error will be evaluated during adaptive solution    
    int n_extr_meth;                 //!< \brief Number of extraction methods to consider (>1 only if adaptive)
    int n_modes;                     //!< \brief Number of modes used in reconstruction formula        
    int SPOD_fs;                     //!< \brief Filter size for SPOD
    int DMD_rank;                    //!< \brief User-defined rank (can be used for POD/SPOD/DMD) if r=0 SVHT is used    
    int MRDMD_max_cycles;            //!< \brief Max number of levels for the mrDMD
    int MRDMD_max_levels;            //!< \brief Minimum number of samples per time window is 4 and can be obtained with MAX_CYCLES=1
    int RDMD_rank;                   //!< \brief Number of modes to extract for recursive DMD 
//    
    std::vector<int> fields;         //!< \brief Fields columns to process
//
    double en;                       //!< \brief Energetic content desired for reconstruction
    double SPOD_sigma;               //!< \brief Sigma value in case of Gaussian filter
    double Mach;                     //!< \brief Freestream Mach number
    double Reynolds;                 //!< \brief Freestream Reynolds number
    double temp_inf;                 //!< \brief Freestream Temperature
    double pres_inf;                 //!< \brief Freestream Pressure
    double dens_inf;                 //!< \brief Freestream Density
    double alpha;                    //!< \brief Angle of attack of the problem
    double beta;                     //!< \brief Angle of sideslip
    double mu;                       //!< \brief Freestream dynamic viscosity
    double dt_res;                   //!< \brief Delta t used for residual evaluation
//    
    std::vector<std::vector<double>> parm;             //!< \brief Array of problem parameters 
    std::vector<std::vector<double>> snap_pnts;        //!< \brief Vector of training points in parameter space
    std::vector<std::vector<double>> targ_pnts;        //!< \brief Vector of target points in parameter space
    std::vector<std::vector<double>> res_eval_pnts;    //!< \brief Points where residual error will be evaluated in adaptive approach
    std::vector<std::vector<double>> dir_eval_pnts;    //!< \brief Points where direct error will be evaluated in adaptive approach    
//
    std::string prob_type;           //!< \brief Problem to be solved (Identification, reconstruction)
    std::string lowdim_approach;     //!< \brief Flag to set adaptive solution or single approach
    std::string flow_type;           //!< \brief Flow under investigation (Steady, unsteady)
    std::string snap_fmt;            //!< \brief Format of snapshot files    
    std::string targ_fmt;            //!< \brief Format of target files for visualisation
    std::string refsol_fmt;          //!< \brief Format of reference solution files    
    std::string field_type;          //!< \brief Defines if scalar or vector field to be processed
    std::string err_fmla;            //!< \brief Flag to set type of error formula, i.e. direct or residual
    std::string lowdim_model_name;   //!< \brief Flag write database basis extraction    
    std::string use_lowdim_model;    //!< \brief Flag to use RB database in reconstruction
    std::string SPOD_ft;             //!< \brief SPOD filter type
    std::string SPOD_flag_bc;        //!< \brief Type of boundary condition for correlation matrix, flag
    std::string ref_field;           //!< \brief Reference field to be subtracted
    std::string DMD_coef_flag;       //!< \brief Method for coefficients calculation
    std::string recon_fmla;          //!< \brief Reconstruction formula
    std::string coeff_fmla;          //!< \brief Formula to compute modal coefficients
    std::string modes_sel_meth;      //!< \brief Method to select number of terms in reconstruction
    std::string mesh_file;           //!< \brief Name of file containing the mesh for visualisation
    std::string mesh_fmt;            //!< \brief Format of the mesh
//
    std::vector<std::string> modes_extraction;  //!< \brief Method to extract basis from snapshots
    std::vector<std::string> parm_name;    //!< \brief Name of parameters
    std::vector<std::string> flds_name;    //!< \brief Name of fields to process
    std::vector<std::string> snap_list;    //!< \brief List of snapshot files
    std::vector<std::string> refsol_list;  //!< \brief List of reference solution files for direct error evaluation
    std::vector<std::string> target_list;  //!< \brief List of target solution files for reconstruction
//
    bool generate_lowdim_model;    //!< \brief Enable basis generation w or w/o saving a database
    bool compute_lowdim_solution;  //!< \brief Enable flow reconstruction w or w/o using a database
//
//------------------------------------------------------------------------------    
//  TO BE REMOVED FROM HERE AND PUT SOMEWHERE ELSE (MF) kept only temporarily
//------------------------------------------------------------------------------    
    int ndim;                       
    int Ds;
    int nstart;
    std::string in_file;
    std::vector<int> cols_coords;    //!< \brief Columns with coordinates    
    double tol;                 //!< \brief Error tolerance for adaptive reconstruction
//
};
//
struct plot3d_info {
//
    int nblocks;            //!< \brief Number of blocks
    std::vector<int> ni;    //!< \brief Number of elements along X for each block
    std::vector<int> nj;    //!< \brief Number of elements along Y for each block
    std::vector<int> nk;    //!< \brief Number of elements along Z for each block
//    
    float Ma, alpha, Re, time;  //!< \brief Mach, angle of attack, Reynolds, Time non-dimensional
};
//
//------------------------------------------------------------------------------
// List of keywords in config file
//------------------------------------------------------------------------------
//
enum keywords { PROBLEM,               // General problem keywords
                LOW_DIMENSIONAL_APPROACH,
                FLOW_TYPE,       
                NAME_PARM,
                SNAPSHOT_LIST,        // Training set data
                SNAPSHOT_FMT,
                FLAG_FIELDS,
                FIELDS_NAME,                
                MODES_EXTRACTION,      // Modal extraction keywords 
                LOW_DIMENSIONAL_MODEL,
                FLAG_REF_FIELD,
                SPOD_FILTER_SIZE,
                SPOD_FILTER_TYPE,
                SPOD_SIGMA,
                SPOD_FLAG_BC,
                DMD_RANK,
                DMD_COEFF_FLAG,
                MR_DMD_MAX_LEVELS,
                MR_DMD_MAX_CYCLES,
                RDMD_RANK,
                MANIFOLD_IDENTIFICATION,  // Manifold learning keywords
                ERROR_FMLA,               // Adaptive framework keywords
                RESID_EVAL_POINTS,
                RESID_DT_BDF,
                REF_SOLUTION_LIST,
                REF_SOLUTION_FMT,                                
                TARGET_LIST,           // Reconstruction keywords
                TARGET_FMT,
                USE_LOW_DIMENSIONAL_MODEL,
                FLAG_METHOD_RECONSTRUCTION,
                FLAG_METHOD_COEFFICIENTS,
                MODES_SELECTION,
                ENERGY_TRESHOLD,
                NUMBER_MODES,                
                FREESTREAM_DENSITY,    // Plugin for aerodynamics 
                FREESTREAM_PRESSURE,
                FREESTREAM_TEMPERATURE,                 
                ALPHA, 
                BETA, 
                MACH, 
                REYNOLDS,                 
                VISCOSITY,
                MESH_FILE,
                MESH_FMT              
};
//
//------------------------------------------------------------------------------
// List of snapshot formats
//------------------------------------------------------------------------------
//
enum snapformat { CSV_ASCII,
                  SU2_CSV,
                  SU2_BINARY,
                  PLOT3D_BINARY,
                  HDF5_BINARY,
                  CGNS 
};
//
//------------------------------------------------------------------------------
// List of field types
//------------------------------------------------------------------------------
//
enum fieldtype { SCALAR,
                 VECTOR,
                 CONSERVATIVE,  // rho, rhoU, rhoV, (rhoW), rhoE
                 PRIMITIVE,     // P, T, U, V, (W)
                 VELOCITY,      // U, V, (W)
                 AEROLOADS,     // P, shearX, shearY, (shearZ)
                 QCRITERION
};
//
//------------------------------------------------------------------------------
// List of reference solutions for subtraction
//------------------------------------------------------------------------------
//
enum refsol { MEAN_FIELD,
              INITIAL_FIELD,
              CONSTANT_FIELD,
              READ_FROM_FILE
};
//
//------------------------------------------------------------------------------            
// Functions
//------------------------------------------------------------------------------            
//
//!< \brief Check if a string contains any non-number character
bool is_number( const std::string& s );
//
//!< \brief Compare keywords with input string
keywords read_keyword_type ( const std::string &key_string );
//
//!< \brief Check snapshot format
snapformat read_snapshot_format ( const std::string &snapfmt );
//
//!< \brief Check field types
fieldtype read_field_type ( const std::string &fldtyp );
//
//!< \brief Check reference solution type for subtraction
refsol read_refsol_type ( const std::string &reftyp );
//
//!< \brief Read config file and store info in prob_settings
void read_config ( std::string filename, prob_settings &settings );
//
//!< \brief Get size of the snapshots. If snapshots are CFD solutions this is the number of grid points in the CFD mesh
int snapshot_size ( const std::string snapshot_file, std::string snapshot_format );
//
//!< \brief Get fields on the specified columns 
Eigen::MatrixXd read_field ( int snp, int Ni, const std::vector<int> &fields, std::ifstream filestream[], std::string filefmt );
//Eigen::MatrixXd read_col ( std::string filename, int Nr, std::vector<int> Cols );
//
//!< \brief Read modes for Recursive DMD
Eigen::MatrixXd read_modes ( std::string filename, int Nr, int r_RDMD );
//
//!< \brief Read coeffs for Recursive DMD
Eigen::MatrixXd read_coefs( std::string filename, int Ns, int r_RDMD );
//
//!< \brief Read Errors and Jaccard index for adaptive reconstruction
Eigen::MatrixXd read_err_j ( std::string filename, int Ns );
//
//!< \brief Read .q (plot3d) file for visualization info
plot3d_info read_plot3d_info (std::string filename);
//
//!< \brief Read .q (plot3d) file, fortran binary for flow field information
std::vector<Eigen::VectorXd> read_plot3d ( std::string filename, plot3d_info Info );
//
//!< \brief Read and change su2 file
void modify_su2_cfg ( std::string file_in, std::string file_out, double dt_res );
//
#endif
