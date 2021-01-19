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
#include "H5Cpp.h"
//
using namespace Eigen;
using namespace H5;
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
struct ld_model_data {
    int n_flds;                      //!< \brief Number of fields to process  
    int n_meth;                      //!< \brief Number of reduction methods to consider    
    std::vector<int> fields;         //!< \brief Fields ID to process
    std::string field_type;          //!< \brief Defines if scalar or vector field to be processed    
    std::string field_attr;          //!< \brief Defines if field is a surface or volume field wrt to the mesh   
    std::string lowdim_model_name;   //!< \brief Flag write database basis extraction        
    std::string ref_field;           //!< \brief Reference field to be subtracted
    std::vector<std::string> flds_name;           //!< \brief Name of fields to process    
    std::vector<std::string> reduction_strategy;  //!< \brief Flag to set adaptive solution or single method 
//
    std::string inpstn_lowdim_model_name;   //!<brief Configuration input string for lowdim_model_name keyword
    std::string inpstn_flag_fields;         //!<brief Configuration input string for flag_fields keyword
    std::string inpstn_fields_name;         //!<brief Configuration input string for field_name keyword
    std::string inpstn_fields_attribute;    //!<brief Configuration input string for field_attr keyword
    std::string inpstn_ref_field;           //!<brief Configuration input string for ref_field keyword
    std::string inpstn_reduction_strategy;  //!<brief Configuration input string for reduction strategy keyword
};
//
struct training_set_data {
    int n_snap;                      //!< \brief Number of snapshots
    int n_parm;                      //!< \brief Number of parameters
    int time_ID;                     //!< \brief For unsteady flows indicate what parameter is time  
    std::vector<std::vector<double>> snap_pnts;   //!< \brief Vector of training points in parameter space    
    std::string flow_type;           //!< \brief Flow under investigation (Steady, unsteady)
    std::string snap_fmt;            //!< \brief Format of snapshot files    
    std::string mesh_file;                 //!< \brief Name of file containing the mesh for visualisation
    std::string mesh_fmt;                  //!< \brief Format of the mesh    
    std::vector<std::string> parm_name;    //!< \brief Name of parameters
    std::vector<std::string> snap_list;    //!< \brief List of snapshot files
//
    std::string inpstn_flow_type;       //!<brief Configuration input string for flow_type keyword
    std::string inpstn_name_parm;       //!<brief Configuration input string for name_parm keyword
    std::string inpstn_snapshot_list;   //!<brief Configuration input string for snapshot_list keyword
    std::string inpstn_snapshot_fmt;    //!<brief Configuration input string for snapshot_fmt keyword
    std::string inpstn_mesh_file;       //!<brief Configuration input string for mesh_file keyword
    std::string inpstn_mesh_fmt;        //!<brief Configuration input string for mesh_fmt keyword
};
//
struct ld_error_data {
    int n_err_pnts;                  //!< \brief Number of points where  error will be evaluated during error solution
    double dtBDF;                    //!< \brief Value of dt for residual erro evaluation through BDF formula
    std::vector<std::vector<double>> err_pnts;    //!< \brief Points where error will be evaluated in error approach    
    std::string err_fmla;                   //!< \brief Flag to set type of error formula, i.e. direct or residual
    std::string err_sol_fmt;                //!< \brief Format of reference solution files    
    std::string error_file;                 //!<brief Configuration input string for error_file keyword
    std::vector<std::string> err_sol_list;  //!< \brief List of solution files for error evaluation
//
    std::string inpstn_error_file;      //!<brief Configuration input string for error_file keyword
    std::string inpstn_error_fmla;      //!<brief Configuration input string for error_fmla keyword    
    std::string inpstn_error_points;    //!<brief Configuration input string for error_points keyword
    std::string inpstn_error_sol_fmt;   //!<brief Configuration input string for err_solution_fmt keyword
    std::string inpstn_resid_dt_bdf;    //!<brief Configuration input string for resid_dt_bdf keyword
};
//
struct modal_identification_data {
    int SPOD_fs;                     //!< \brief Filter size for SPOD
    int DMD_rank;                    //!< \brief User-defined rank (can be used for POD/SPOD/DMD) if r=0 SVHT is used    
    int MRDMD_max_cycles;            //!< \brief Max number of levels for the mrDMD
    int MRDMD_max_levels;            //!< \brief Minimum number of samples per time window is 4 and can be obtained with MAX_CYCLES=1
    int RDMD_rank;                   //!< \brief Number of modes to extract for recursive DMD     
    double SPOD_sigma;               //!< \brief Sigma value in case of Gaussian filter
    std::string SPOD_ft;             //!< \brief SPOD filter type
    std::string SPOD_flag_bc;        //!< \brief Type of boundary condition for correlation matrix, flag
    std::string DMD_coef_flag;       //!< \brief Method for coefficients calculation
//
    std::string inpstn_spod_filter_size;   //!<brief Configuration input string for spod_filter_size keyword
    std::string inpstn_spod_filter_type;   //!<brief Configuration input string for spod_filter_type keyword
    std::string inpstn_spod_sigma;         //!<brief Configuration input string for spod_sigma keyword
    std::string inpstn_spod_flag_bc;       //!<brief Configuration input string for spod_flag_bc keyword
    std::string inpstn_dmd_rank;           //!<brief Configuration input string for dmd_rank keyword
    std::string inpstn_dmd_coef_flag;      //!<brief Configuration input string for dmd_coef_flag keyword
    std::string inpstn_mr_dmd_max_levels;  //!<brief Configuration input string for mr_dmd_max_levels keyword 
    std::string inpstn_mr_dmd_max_cycles;  //!<brief Configuration input string for mr_dmd_max_cycles keyword
    std::string inpstn_rdmd_rank;          //!<brief Configuration input string for rdmd_rank keyword
};
//
struct manifold_learning_data {
    std::string manifold_meth;
    std::string inpstn_manifold_meth;
};
//
struct target_set_data {
    int n_targ;      //!< \brief Number of target solution files in reconstruction
    int n_flds;      //!< \brief Number of fields to reconstruct
    int n_parm;      //!< \brief Number of parameters
    int n_modes;     //!< \brief Number of modes used in reconstruction formula        
    double en;       //!< \brief Energetic content desired for reconstruction
    std::vector<std::vector<double>> target_pnts;    //!< \brief Vector of target points in parameter space      
    std::string target_fmt;                  //!< \brief Format of target files for visualisation
    std::string recon_fmla;                  //!< \brief Reconstruction formula
    std::string coeff_fmla;                  //!< \brief Formula to compute modal coefficients
    std::string modes_sel_meth;              //!< \brief Method to select number of terms in reconstruction
    std::vector<std::string> target_list;    //!< \brief List of target solution files for reconstruction
    std::vector<std::string> target_ldmodel; //!< \brief List of low dimensional models to compute solution files for reconstruction
    std::vector<std::string> target_lderror; //!< \brief File containing the error
//
    std::string inpstn_target_list;                 //!<brief Configuration input string for target_list keyword
    std::string inpstn_target_fmt;                  //!<brief Configuration input string for target_fmt keyword
    std::string inpstn_flag_method_coefficients;    //!<brief Configuration input string for flag_method_coefficients keyword
    std::string inpstn_modes_selection;             //!<brief Configuration input string for modes_selection keyword
    std::string inpstn_energy_treshold;             //!<brief Configuration input string for energy_treshold keyword
    std::string inpstn_number_modes;                //!<brief Configuration input string for number_modes keyword
};
//
struct aero_data {
    int Nfom;               //!< \brief Number of figure of merits to evaluate
    double Mach;            //!< \brief Freestream Mach number
    double Reynolds;        //!< \brief Freestream Reynolds number
    double temp_inf;        //!< \brief Freestream Temperature
    double pres_inf;        //!< \brief Freestream Pressure
    double dens_inf;        //!< \brief Freestream Density
    double alpha;           //!< \brief Angle of attack of the problem
    double beta;            //!< \brief Angle of sideslip
    double mu;              //!< \brief Freestream dynamic viscosity
    double refA;            //!< \brief Reference area
    double refL;            //!< \brief Reference lenght
    std::vector<double> cmcoord;         //!< \brief Coordinates for moment
    std::string cmaxis;                  //!< \brief Axis around which the moment is to be computed
    std::string target_profile_file;     //!< \brief Name of file containing the reference profile for inverse design
    std::string target_profile_fmt;      //!< \brief Format of the reference profile for inverse design
    std::vector<std::string> aeroutput;  //!< \brief Type of aerodynamic output needed
//
    std::string inpstn_aeroutput;               //!<brief Configuration input string for the type of aerodynamic output required
    std::string inpstn_freestream_density;      //!<brief Configuration input string for freestream_density keyword
    std::string inpstn_freestream_pressure;     //!<brief Configuration input string for freestream_pressure keyword
    std::string inpstn_freestream_temperature;  //!<brief Configuration input string for freestream_temperature keyword
    std::string inpstn_alpha;       //!<brief Configuration input string for alpha keyword
    std::string inpstn_beta;        //!<brief Configuration input string for beta keyword
    std::string inpstn_mach;        //!<brief Configuration input string for mach keyword
    std::string inpstn_reynolds;    //!<brief Configuration input string for reynolds keyword
    std::string inpstn_viscosity;   //!<brief Configuration input string for viscosity keyword
//
    std::string inpstn_reference_area;         //!<brief Configuration input string for the reference area
    std::string inpstn_reference_lenght;       //!<brief Configuration input string for the reference lenght
    std::string inpstn_cmcoordinates;          //!<brief Configuration input string for the reference coordinates for moment
    std::string inpstn_cmaxis;                 //!<brief Configuration input string for the axis around which compute the moment
    std::string inpstn_target_profile;         //!<brief Configuration input string for the target profile for inverse design
    std::string inpstn_target_profile_fmt;     //!<brief Configuration input string for the format of the target profile
};
//
struct plot3d_info {
    int nblocks;            //!< \brief Number of blocks
    std::vector<int> ni;    //!< \brief Number of elements along X for each block
    std::vector<int> nj;    //!< \brief Number of elements along Y for each block
    std::vector<int> nk;    //!< \brief Number of elements along Z for each block
    float Ma, alpha, Re, time;  //!< \brief Mach, angle of attack, Reynolds, Time non-dimensional
};
//
// List of keywords in config file
enum keywords { LOW_DIMENSIONAL_MODEL_NAME,  // Low dimensional model keywords
                REDUCTION_STRATEGY,   
                FIELDS_ID,
                FIELDS_NAME,
                FIELDS_ATTRIBUTE,
                FLAG_REF_FIELD,
                SNAPSHOT_LIST,  // Training set data
                SNAPSHOT_FMT,
                NAME_PARM,
                FLOW_TYPE,
                MESH_FILE,
                MESH_FMT,                                              
                ERROR_FMLA,  // Error evaluation
                ERROR_EVAL_POINTS,
                ERROR_SOLUTION_LIST,
                ERROR_SOLUTION_FMT,
                RESID_DT_BDF,
                SPOD_FILTER_SIZE,   // Modal identification keywords                                 
                SPOD_FILTER_TYPE,
                SPOD_SIGMA,
                SPOD_FLAG_BC,
                DMD_RANK,
                DMD_COEFF_FLAG,
                MR_DMD_MAX_LEVELS,
                MR_DMD_MAX_CYCLES,
                RDMD_RANK,                
                TARGET_LIST,           // Reconstruction keywords
                TARGET_FMT,                
                MODAL_METHOD_COEFFICIENTS,
                MODAL_MODES_SELECTION,
                MODAL_ENERGY_TRESHOLD,
                MODAL_NUMBER_MODES,                
                AERODYNAMIC_OUTPUT,    // Plugin for aerodynamics 
                FREESTREAM_DENSITY,
                FREESTREAM_PRESSURE,
                FREESTREAM_TEMPERATURE,                 
                ALPHA, 
                BETA, 
                MACH, 
                REYNOLDS,                 
                VISCOSITY,
                REFERENCE_AREA,
                REFERENCE_LENGHT,
                CM_AXIS,
                CM_COORDINATES,
                TARGET_PROFILE,
                TARGET_PROFILE_FMT };
//
// List of snapshot formats
enum snapformat { CSV_ASCII, SU2_CSV, SU2_BINARY, PLOT3D_BINARY, HDF5_BINARY,
                  CGNS };
//
// List of field types
enum fieldtype { SCALAR, VECTOR,
                 CONSERVATIVE,  // rho, rhoU, rhoV, (rhoW), rhoE
                 PRIMITIVE,     // P, T, U, V, (W)
                 VELOCITY,      // U, V, (W)
                 AEROLOADS,     // P, shearX, shearY, (shearZ)
                 QCRITERION };
//
// List of reference solutions for subtraction
enum refsol { NONE, MEAN_FIELD, INITIAL_FIELD, CONSTANT_FIELD, READ_FROM_FILE };
//
// List of file format. Applies only to data for low-dimensional model generation
enum datfmt { RAZR, HDF5 };
//
// List of mesh formats
enum mshfmt { SU2, CGNSM, VTK };
//
// List of mesh formats
enum su2key { NDIME, NELEM, NPOIN, NMARK, MARKER_TAG, MARKER_ELEMS, NOKEYW };
//
//  /*!< \brief Main class for object mesh
class CMesh {
//
protected:
    ld_model_data m_lowdim_data;
//
    int sdim;
    int nelements;
    int nboundaries;
//
    int npoints;
    std::vector<std::vector<double>> coordinates;
    Eigen::MatrixXd Coords;
    Eigen::MatrixXi Connectivity;
    Eigen::VectorXi PointID;
//
//  1D elements    
    int nsegm;
    std::vector<std::vector<int>> segm;
//
//  2D elements
    int ntria; 
    std::vector<std::vector<int>> tria;
    int nquad;
    std::vector<std::vector<int>> quad;
//
//  3D elements
    int ntetra; 
    std::vector<std::vector<int>> tetra;
    int nhexa; 
    std::vector<std::vector<int>> hexa;
    int nprism; 
    std::vector<std::vector<int>> prism;
    int npyra;    
    std::vector<std::vector<int>> pyra;
//
public:
//
//  Constructor
    CMesh( ld_model_data lowdim_data) {
        m_lowdim_data = lowdim_data;
        npoints = -1; nelements = -1; nboundaries = -1;
        nsegm = -1; ntria= -1; nquad = -1; ntetra = -1; nhexa = -1;
        nprism = -1; npyra = -1; };
//
    void read( std::ifstream &fstream, const std::string &format );
    void saveH5_volume();
    void saveH5_surface();

    Eigen::MatrixXd GetCoords();
    Eigen::MatrixXi GetConnectivity();
//
};
//
//------------------------------------------------------------------------------            
// Functions
//------------------------------------------------------------------------------            
//
//!< \brief Check if a string contains any non-number character
bool is_number( const std::string &s );
//
//!< \brief Remove leading spaces in a string
std::string ltrim(const std::string &s);
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
//!< \brief Check mesh format to be read
mshfmt read_mesh_format ( const std::string &format );
//
//!< \brief Check SU@ mesh keyword
su2key su2mesh_keyword( const std::string &key );
//
//!< \brief Read data needed for low-dimensional model generation
void read_gendata ( const std::string filename, const std::string filefmt, ld_model_data 
    &lowdim_data, training_set_data &training_data, ld_error_data &error_data, 
    modal_identification_data &modal_data, manifold_learning_data &manifold_data );
//
//!< \brief Read data needed for low-dimensional solution computation
void read_soldata ( const std::string filename, target_set_data &target_data, 
    aero_data &aerodynamic_data );
//
//!< \brief Read low dimensional model data
void read_lowdim_model_data( ld_model_data &lowdim_data );
//
//!< \brief Read training set data
void read_training_set_data( training_set_data &training_data );
//
//!< \brief Read modal extraction data
void read_modal_identification_data( modal_identification_data &modal_data );
//
//!< \brief Read manifold learning data
void read_manifold_learning_data( manifold_learning_data &manifold_data );
//
//!< \brief Read error reduction
void read_error_data( ld_error_data &error_data );
//
//!< \brief Read target set data
void read_target_set_data( target_set_data &target_data );
//
//!< \brief Read aerodynamics
void read_aerodynamic_data( aero_data &aerodynamic_data );
//
//!< \brief Get size of the snapshots. If snapshots are CFD solutions this is the number of grid points in the CFD mesh
int snapshot_size ( const std::string snapshot_file, std::string snapshot_format );
//
//!< \brief Get fields on the specified columns 
Eigen::MatrixXd read_field ( int snp, int Ni, const std::vector<int> &fields, std::ifstream filestream[], 
  std::string filefmt );
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
