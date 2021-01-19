/* ------------------------------------------------------------------------------
\file read_data.cpp
* \brief Subroutines and functions to read configuration file.
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
#include "read_data.hpp"
//
// -------------------------------------------------------------------------------------------
bool is_number( const std::string &s ) {
// -------------------------------------------------------------------------------------------
    long double ld;
    return( (std::istringstream(s) >> ld >> std::ws).eof() );
}
//
// -------------------------------------------------------------------------------------------
std::string ltrim(const std::string& s) {
// -------------------------------------------------------------------------------------------
//    
    const std::string WHITESPACE = " \n\r\t\f\v";
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
//
// -------------------------------------------------------------------------------------------
inline void wait_on_enter () {
// -------------------------------------------------------------------------------------------    
    std::string dummy;
    std::cout << "Enter to continue..." << std::endl;
    std::getline(std::cin, dummy);
}
//
// -------------------------------------------------------------------------------------------
keywords read_keyword_type( const std::string &key_string ) {
// -------------------------------------------------------------------------------------------    
//
    if ( key_string == "LOW_DIMENSIONAL_MODEL_NAME" )
        return LOW_DIMENSIONAL_MODEL_NAME;
    else if ( key_string == "REDUCTION_STRATEGY" )
        return REDUCTION_STRATEGY;
    else if( key_string == "FIELDS_ID" )
        return FIELDS_ID;
    else if( key_string == "FIELDS_NAME" )
        return FIELDS_NAME;    
    else if( key_string == "FLAG_REF_FIELD" )
        return FLAG_REF_FIELD;    
    else if( key_string == "FIELDS_ATTRIBUTE" )
        return FIELDS_ATTRIBUTE;        
    else if( key_string == "SNAPSHOT_LIST" )
        return SNAPSHOT_LIST;    
    else if( key_string == "SNAPSHOT_FMT" )
        return SNAPSHOT_FMT;
    else if( key_string == "NAME_PARM" )
        return NAME_PARM;
    else if( key_string == "FLOW_TYPE" )
        return FLOW_TYPE;        
    else if( key_string == "MESH_FILE" )
        return MESH_FILE;
    else if( key_string == "MESH_FMT" )
        return MESH_FMT;    
    else if( key_string == "SPOD_FILTER_SIZE" )
        return SPOD_FILTER_SIZE;
    else if( key_string == "SPOD_FILTER_TYPE" )
        return SPOD_FILTER_TYPE;
    else if( key_string == "SPOD_SIGMA" )
        return SPOD_SIGMA;
    else if( key_string == "SPOD_FLAG_BC" )
        return SPOD_FLAG_BC;
    else if( key_string == "DMD_RANK" )
        return DMD_RANK;
    else if( key_string == "DMD_COEFF_FLAG" )
        return DMD_COEFF_FLAG;
    else if( key_string == "MR_DMD_MAX_CYCLES" )
        return MR_DMD_MAX_CYCLES;
    else if( key_string == "MR_DMD_MAX_LEVELS" )
        return MR_DMD_MAX_LEVELS;
    else if( key_string == "RDMD_RANK" )
        return RDMD_RANK;
    else if( key_string == "ERROR_FMLA" )
        return ERROR_FMLA;
    else if( key_string == "ERROR_EVAL_POINTS" )
        return ERROR_EVAL_POINTS;
    else if( key_string == "ERROR_SOLUTION_FMT" )
        return ERROR_SOLUTION_FMT;    
    else if( key_string == "RESID_DT_BDF" )
        return RESID_DT_BDF;
    else if( key_string == "TARGET_LIST" )
        return TARGET_LIST;
    else if( key_string == "TARGET_FMT" )
        return TARGET_FMT;
    else if( key_string == "MODAL_METHOD_COEFFICIENTS" )
        return MODAL_METHOD_COEFFICIENTS;
    else if( key_string == "MODAL_MODES_SELECTION" )
        return MODAL_MODES_SELECTION;
    else if( key_string == "MODAL_ENERGY_TRESHOLD" )
        return MODAL_ENERGY_TRESHOLD;
    else if( key_string == "MODAL_NUMBER_MODES" )
        return MODAL_NUMBER_MODES;
    else if( key_string == "AERODYNAMIC_OUTPUT" )
        return AERODYNAMIC_OUTPUT;    
    else if( key_string == "FREESTREAM_DENSITY" )
        return FREESTREAM_DENSITY;
    else if( key_string == "FREESTREAM_PRESSURE" )
        return FREESTREAM_PRESSURE;
    else if( key_string == "FREESTREAM_TEMPERATURE" )
        return FREESTREAM_TEMPERATURE;
    else if( key_string == "ALPHA" )
        return ALPHA;
    else if( key_string == "BETA" )
        return BETA;
    else if( key_string == "MACH" )
        return MACH;
    else if( key_string == "REYNOLDS" )
        return REYNOLDS;
    else if( key_string == "VISCOSITY" )
        return VISCOSITY;
    else if( key_string == "REFERENCE_AREA" )
        return REFERENCE_AREA;
    else if( key_string == "REFERENCE_LENGHT" )
        return REFERENCE_LENGHT;
    else if( key_string == "CM_AXIS" )
        return CM_AXIS;
    else if( key_string == "CM_COORDINATES" )
        return CM_COORDINATES;
    else if( key_string == "TARGET_PROFILE" )
        return TARGET_PROFILE;
    else if( key_string == "TARGET_PROFILE_FMT" )
        return TARGET_PROFILE_FMT;
//    
    else {
//
        std::cout << " " << std::endl;
        std::cout << "-> Error: Unknown keyword in configuration file" << std::endl;
        std::cout << key_string << std::endl; 
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//     
    }
}
//
// -------------------------------------------------------------------------------------------
snapformat read_snapshot_format ( const std::string &snapfmt ) {
// -------------------------------------------------------------------------------------------
//
    if ( snapfmt == "CSV_ASCII" )
        return CSV_ASCII;
    else if ( snapfmt == "SU2_CSV" )
        return SU2_CSV;
    else if ( snapfmt == "SU2_BINARY" )
        return SU2_BINARY;
    else if ( snapfmt == "PLOT3D_BINARY" )
        return PLOT3D_BINARY;
    else if ( snapfmt == "HDF5_BINARY" )
        return HDF5_BINARY;
    else if ( snapfmt == "CGNS" )
        return CGNS; 
//
    else {
//
        std::cout << " " << std::endl;
        std::cout << "-> Error: Unknown snapshot format" << std::endl;
        std::cout << snapfmt << std::endl; 
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//     
    }
}
//
// -------------------------------------------------------------------------------------------
fieldtype read_field_type ( const std::string &fldtyp ) {
// -------------------------------------------------------------------------------------------
//
    if ( fldtyp == "SCALAR" )
        return SCALAR;
    else if ( fldtyp == "VECTOR" )
        return VECTOR;
    else if ( fldtyp == "CONSERVATIVE" )
        return CONSERVATIVE;
    else if ( fldtyp == "PRIMITIVE" )
        return PRIMITIVE;
    else if ( fldtyp == "AEROLOADS" )
        return AEROLOADS;
    else if ( fldtyp == "VELOCITY" )
        return VELOCITY;
    else if ( fldtyp == "QCRITERION" )
        return QCRITERION; 
//
    else {
//
        std::cout << " " << std::endl;
        std::cout << "-> Error: Unknown field type" << std::endl;
        std::cout << fldtyp << std::endl; 
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//     
    }
}
//
// -------------------------------------------------------------------------------------------
refsol read_refsol_type ( const std::string &reftyp ) {
// -------------------------------------------------------------------------------------------
//
    if ( reftyp == "NONE" )
        return NONE;
    else if ( reftyp == "MEAN_FIELD" )
        return MEAN_FIELD;
    else if ( reftyp == "INITIAL_FIELD" )
        return INITIAL_FIELD;
    else if ( reftyp == "CONSTANT_FIELD" )
        return CONSTANT_FIELD;
    else if ( reftyp == "READ_FROM_FILE" )
        return READ_FROM_FILE;
//
    else {
//
        std::cout << " " << std::endl;
        std::cout << "-> Error: Unknown FLAG_REF_FIELD in configuration file" << std::endl;
        std::cout << reftyp << std::endl; 
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//     
    }
}
//
// -------------------------------------------------------------------------------------------
datfmt read_data_format ( const std::string &datfmt ) {
// -------------------------------------------------------------------------------------------
//
    if ( datfmt == "RAZR" )
        return RAZR;
    else if ( datfmt == "HDF5" )
        return HDF5;
    else {
//
        std::cout << std::endl;
        std::cout << "-> Error: Unknown data file format" << std::endl;        
        std::cout << std::endl;
        exit (EXIT_FAILURE);
    }
}
//
// -------------------------------------------------------------------------------------------
mshfmt read_mesh_format ( const std::string &format ) {
// -------------------------------------------------------------------------------------------
//
    if ( format == "SU2" )
        return SU2;
    else if ( format == "CGNSM" )
        return CGNSM;
    else if ( format == "VTK" )
        return VTK;
    else {
//
        std::cout << std::endl;
        std::cout << "-> Error: Unknown mesh file format" << std::endl;        
        std::cout << std::endl;
        exit (EXIT_FAILURE);
    }
}
//
// -------------------------------------------------------------------------------------------
su2key su2mesh_keyword( const std::string &key ) {
// -------------------------------------------------------------------------------------------
//
    if ( key == "NDIME" )
        return NDIME;
    else if ( key == "NELEM" )
        return NELEM;
    else if ( key == "NPOIN" )
        return NPOIN;
    else if ( key == "NMARK" )
        return NMARK;
    else if ( key == "MARKER_TAG" )
        return MARKER_TAG;
    else if ( key == "MARKER_ELEMS" )
        return MARKER_ELEMS;
    else 
        return NOKEYW;
//
}
//
// ----------------------------------------------------------------------------------------
void read_gendata ( const std::string filename, const std::string filefmt, ld_model_data 
    &lowdim_data, training_set_data &training_data, ld_error_data &error_data, 
    modal_identification_data &modal_data, manifold_learning_data &manifold_data ) {
// ----------------------------------------------------------------------------------------
//
//  Low dimensional data
    lowdim_data.inpstn_reduction_strategy = "POD";
    lowdim_data.inpstn_lowdim_model_name = {};
    lowdim_data.inpstn_flag_fields = {};      
    lowdim_data.inpstn_fields_name = {};      
    lowdim_data.inpstn_fields_attribute = "VOLUME";
    lowdim_data.inpstn_ref_field = "NONE";            
//
//  Training set data
    training_data.inpstn_snapshot_list = {};
    training_data.inpstn_snapshot_fmt = "SU2_CSV";    
    training_data.inpstn_name_parm = {};
    training_data.inpstn_flow_type = "STEADY";  
    training_data.inpstn_mesh_file = "NONE";
    training_data.inpstn_mesh_fmt = {};    
//
//  Error data
    error_data.inpstn_error_file = {};
    error_data.inpstn_error_fmla = "RESIDUAL";
    error_data.inpstn_error_points = {};
    error_data.inpstn_error_sol_fmt = {};
    error_data.inpstn_resid_dt_bdf = "1.e-3";
//
//  Modal extraction data
    modal_data.inpstn_spod_filter_size = "-1";
    modal_data.inpstn_spod_filter_type = "BOX";
    modal_data.inpstn_spod_sigma = "1.0";
    modal_data.inpstn_spod_flag_bc = "ZERO";
    modal_data.inpstn_dmd_rank = "-2";
    modal_data.inpstn_dmd_coef_flag = "OPT";
    modal_data.inpstn_mr_dmd_max_levels = "-1";
    modal_data.inpstn_mr_dmd_max_cycles = "-1";
    modal_data.inpstn_rdmd_rank = "-1";    
//
//  Manifold learning data
//  ...
//
// ----------------------------------------------------------------------------------------------------------------------        
//  Reading each line of the input file into a string that will be parsed later on to avoid issues with dependencies
//  or have to introduce multiple times the same information. In this way the information in the configuration file
//  can be entered in whatsoever order. The file is the ASCII configuration file if generationg low dimensinoal model
//  otherwise, in the case of solution reconstruction the file is the HDF5 containing the low dimensional model
// ----------------------------------------------------------------------------------------------------------------------
//
    switch ( read_data_format (filefmt) ) {
//        
        case RAZR: {
//
            std::ifstream cFile ( filename );
            if ( cFile.is_open() ) {
//
                size_t delimiterPos; 
                std::string name, value;
                std::string line;
//
                while ( getline(cFile, line) ) {
//            
                    line.erase( remove_if(line.begin(), line.end(), isspace), line.end() );  //include ::  in front of isspace if using namespace std
                    if ( line[0] == '#' || line.empty() )
                        continue;
                    delimiterPos = line.find("=");
                    name = line.substr(0, delimiterPos);
                    value = line.substr(delimiterPos + 1);
//
                    switch (read_keyword_type(name)) {                
                        case LOW_DIMENSIONAL_MODEL_NAME: { lowdim_data.inpstn_lowdim_model_name = value; break; }
                        case REDUCTION_STRATEGY: { lowdim_data.inpstn_reduction_strategy = value; break; }
                        case FIELDS_ID: { lowdim_data.inpstn_flag_fields = value; break; }
                        case FIELDS_NAME: { lowdim_data.inpstn_fields_name = value; break; }
                        case FIELDS_ATTRIBUTE: { lowdim_data.inpstn_fields_attribute = value; break; }
                        case FLAG_REF_FIELD: { lowdim_data.inpstn_ref_field = value; break; }
                        case SNAPSHOT_LIST: { training_data.inpstn_snapshot_list = value; break; }
                        case SNAPSHOT_FMT: { training_data.inpstn_snapshot_fmt = value; break; }
                        case NAME_PARM: { training_data.inpstn_name_parm = value; break; }                 
                        case FLOW_TYPE: { training_data.inpstn_flow_type = value; break; }                
                        case MESH_FILE: { training_data.inpstn_mesh_file = value; break; }
                        case MESH_FMT: { training_data.inpstn_mesh_fmt = value; break; }                
                        case SPOD_FILTER_SIZE: { modal_data.inpstn_spod_filter_size = value; break; }
                        case SPOD_FILTER_TYPE: { modal_data.inpstn_spod_filter_type = value; break; }
                        case SPOD_SIGMA: { modal_data.inpstn_spod_sigma = value; break; }
                        case SPOD_FLAG_BC: { modal_data.inpstn_spod_flag_bc = value; break; }
                        case DMD_RANK: { modal_data.inpstn_dmd_rank = value; break; }
                        case DMD_COEFF_FLAG: { modal_data.inpstn_dmd_coef_flag = value; break; }
                        case MR_DMD_MAX_LEVELS: { modal_data.inpstn_mr_dmd_max_levels = value; break; }
                        case MR_DMD_MAX_CYCLES: { modal_data.inpstn_mr_dmd_max_cycles = value; break; }
                        case RDMD_RANK: { modal_data.inpstn_rdmd_rank = value; break; }
                        case ERROR_FMLA: { error_data.inpstn_error_fmla = value; break; }
                        case ERROR_EVAL_POINTS: { error_data.inpstn_error_points = value; break; }
                        case ERROR_SOLUTION_FMT: { error_data.inpstn_error_sol_fmt = value; break; }
                        case RESID_DT_BDF: { error_data.inpstn_resid_dt_bdf = value; break; }                                
                        default: { break; } 
                    }
                }
                cFile.close();                
//
            } else {
                std::cout << std::endl;
                std::cout << "-> Error: Unable to open configuration file." << std::endl;
                std::cout << filename  << std::endl;
                std::cout << std::endl;
                exit (EXIT_FAILURE);
            } break; }
//
        case HDF5: { std::cout << "coming ..." << std::endl; break; }
//
    }
//
//  Parsing each input string and assign values to variables and arrays
    read_lowdim_model_data( lowdim_data );
    read_training_set_data( training_data );
    read_error_data( error_data );
//
    if ( std::any_of( lowdim_data.reduction_strategy.begin(), lowdim_data.reduction_strategy.end(), [](std::string a){return a == "POD";} ) ||
        std::any_of( lowdim_data.reduction_strategy.begin(), lowdim_data.reduction_strategy.end(), [](std::string a){return a == "SPOD";} ) ||
        std::any_of( lowdim_data.reduction_strategy.begin(), lowdim_data.reduction_strategy.end(), [](std::string a){return a == "DMD";} ) || 
        std::any_of( lowdim_data.reduction_strategy.begin(), lowdim_data.reduction_strategy.end(), [](std::string a){return a == "MRDMD";} ) ||                   
        std::any_of( lowdim_data.reduction_strategy.begin(), lowdim_data.reduction_strategy.end(), [](std::string a){return a == "RDMD";} ) ) 
        read_modal_identification_data( modal_data );
//        
    if ( std::any_of( lowdim_data.reduction_strategy.begin(), lowdim_data.reduction_strategy.end(), [](std::string a){return a == "ISOMAP";} ) )
        read_manifold_learning_data( manifold_data );
//
}
//
//
// ----------------------------------------------------------------------------------------
void read_soldata ( const std::string filename, target_set_data &target_data, 
    aero_data &aerodynamic_data ) {
// ----------------------------------------------------------------------------------------
//
//  Target set data  
    target_data.inpstn_target_list = {};
    target_data.inpstn_target_fmt = "SU2_CSV";
    target_data.inpstn_flag_method_coefficients = "RBF_LINEAR";
    target_data.inpstn_modes_selection = "ENERGY";
    target_data.inpstn_energy_treshold = "0.95";
    target_data.inpstn_number_modes = "3";
//
//  Aerodynamic plugin data
    aerodynamic_data.inpstn_aeroutput = {};    
    aerodynamic_data.inpstn_freestream_density = "1.225";
    aerodynamic_data.inpstn_freestream_pressure = "101325";
    aerodynamic_data.inpstn_freestream_temperature = "288.15";                 
    aerodynamic_data.inpstn_alpha = "0.0"; 
    aerodynamic_data.inpstn_beta = "0.0"; 
    aerodynamic_data.inpstn_mach = "0.3"; 
    aerodynamic_data.inpstn_reynolds = "100";                 
    aerodynamic_data.inpstn_viscosity = "1.78e-5";
    aerodynamic_data.inpstn_reference_area = "1";
    aerodynamic_data.inpstn_reference_lenght = "1";
    aerodynamic_data.inpstn_cmaxis = "y-axis";
    aerodynamic_data.inpstn_cmcoordinates = {};
    aerodynamic_data.inpstn_target_profile = {};
    aerodynamic_data.inpstn_target_profile_fmt = {};
//
// ----------------------------------------------------------------------------------------------------------------------        
//  Reading each line of the input file (ASCII) into a string that will be parsed later on to avoid issues 
//  with dependencies or have to introduce multiple times the same information. In this way the information 
//  in the configuration file can be entered in whatsoever order.
// ----------------------------------------------------------------------------------------------------------------------
//
    std::ifstream cFile ( filename );
    if ( cFile.is_open() ) {
//
        size_t delimiterPos; 
        std::string name, value;
        std::string line;
//
        while ( getline(cFile, line) ) {
//            
            line.erase( remove_if(line.begin(), line.end(), isspace), line.end() );  //include ::  in front of isspace if using namespace std
            if ( line[0] == '#' || line.empty() )
                continue;
            delimiterPos = line.find("=");
            name = line.substr(0, delimiterPos);
            value = line.substr(delimiterPos + 1);
//
            switch (read_keyword_type(name)) {        
                case TARGET_LIST: { target_data.inpstn_target_list = value; break; }
                case TARGET_FMT: { target_data.inpstn_target_fmt = value; break; }
                case MODAL_METHOD_COEFFICIENTS: { target_data.inpstn_flag_method_coefficients = value; break; }
                case MODAL_MODES_SELECTION: { target_data.inpstn_modes_selection = value; break; }
                case MODAL_ENERGY_TRESHOLD: { target_data.inpstn_energy_treshold = value; break; }
                case MODAL_NUMBER_MODES: { target_data.inpstn_number_modes = value; break; }
                case AERODYNAMIC_OUTPUT: { aerodynamic_data.inpstn_aeroutput = value; break; }
                case FREESTREAM_DENSITY: { aerodynamic_data.inpstn_freestream_density = value; break; }
                case FREESTREAM_PRESSURE: { aerodynamic_data.inpstn_freestream_pressure = value; break; }
                case FREESTREAM_TEMPERATURE: { aerodynamic_data.inpstn_freestream_temperature = value; break; }
                case ALPHA: { aerodynamic_data.inpstn_alpha = value; break; }
                case BETA: { aerodynamic_data.inpstn_beta = value; break; }
                case MACH: { aerodynamic_data.inpstn_mach = value; break; }
                case REYNOLDS: { aerodynamic_data.inpstn_reynolds = value; break; }     
                case VISCOSITY: { aerodynamic_data.inpstn_viscosity = value; break; }
                case REFERENCE_AREA: { aerodynamic_data.inpstn_reference_area = value; break; }
                case REFERENCE_LENGHT: { aerodynamic_data.inpstn_reference_lenght = value; break; }
                case CM_AXIS: { aerodynamic_data.inpstn_cmaxis = value; break; }
                case CM_COORDINATES: { aerodynamic_data.inpstn_cmcoordinates = value; break; }
                case TARGET_PROFILE: { aerodynamic_data.inpstn_target_profile = value; break; }
                case TARGET_PROFILE_FMT: { aerodynamic_data.inpstn_target_profile_fmt = value; break; }
                default: { break; }
            }
        }
        cFile.close();                
//
    } else {
//
        std::cout << std::endl;
        std::cout << "-> Error: Unable to open configuration file." << std::endl;
        std::cout << filename  << std::endl;
        std::cout << std::endl;
        exit (EXIT_FAILURE);
    }
//
    read_target_set_data( target_data );
    if ( !aerodynamic_data.inpstn_aeroutput.empty() ) {        
    read_aerodynamic_data( aerodynamic_data ); }
}
//
//-------------------------------------------------------------------------------------------------------
void read_lowdim_model_data( ld_model_data &lowdim_data ) {
//-------------------------------------------------------------------------------------------------------
//
    int counter;
    int pcount;
    std::string substr;
    std::vector<double> tmp;
//
//  Low dimensional model data
    lowdim_data.lowdim_model_name = {};
    lowdim_data.reduction_strategy = {};
    lowdim_data.field_type = {};
    lowdim_data.field_attr = {};
    lowdim_data.fields = {};   
    lowdim_data.n_flds = 0;  
    lowdim_data.n_meth = 0;    
    lowdim_data.flds_name = {};
    lowdim_data.ref_field = {};
//
    if ( lowdim_data.inpstn_lowdim_model_name.empty() ) {
//
        std::cout << " " << std::endl;
        std::cout << "-> Error: The keyword LOW_DIMENSIONAL_MODEL_NAME is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
    } else {
        lowdim_data.lowdim_model_name = lowdim_data.inpstn_lowdim_model_name;
    }
//
    if ( lowdim_data.inpstn_reduction_strategy.empty() ) {
//
        std::cout << " " << std::endl;
        std::cout << "-> Error: The keyword REDUCTION_STRATEGY is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
    } else {
//        
        std::stringstream sstream(lowdim_data.inpstn_reduction_strategy);
        while(sstream.good()) {
            getline(sstream, substr, ',');
            lowdim_data.reduction_strategy.push_back(substr);
        }
        lowdim_data.n_meth = lowdim_data.reduction_strategy.size();
//
    }
//
    if ( lowdim_data.inpstn_flag_fields.empty() ) {
//
        std::cout << " " << std::endl;
        std::cout << "-> Error: The keyword FLAG_FIELDS is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
    } else {
//
        std::stringstream sstream(lowdim_data.inpstn_flag_fields);
        counter = 0;                
        while ( sstream.good() ) {
            getline(sstream, substr, ','); 
//
            if ( counter == 0 ) {
                lowdim_data.field_type = substr;
                counter ++;
            } else {
                lowdim_data.fields.push_back(std::stoi(substr));
            }
//
        }
        sort(lowdim_data.fields.begin(), lowdim_data.fields.end());
//        
        if ( lowdim_data.field_type == "VECTOR" || lowdim_data.field_type == "CONSERVATIVE" || 
            lowdim_data.field_type == "PRIMITIVE" || lowdim_data.field_type == "VELOCITY" || 
            lowdim_data.field_type == "AEROLOADS") {
            lowdim_data.n_flds = 1;
        } else {
            lowdim_data.n_flds = lowdim_data.fields.size();
        }
    }
//
    if ( !lowdim_data.inpstn_fields_name.empty() ) {
//
        std::stringstream sstream(lowdim_data.inpstn_fields_name);
        while ( sstream.good() ) {
            getline(sstream, substr, ',');
            lowdim_data.flds_name.push_back(substr);
        }
    }
//
    lowdim_data.field_attr = lowdim_data.inpstn_fields_attribute;
    lowdim_data.ref_field = lowdim_data.inpstn_ref_field;
//
}
//
//-------------------------------------------------------------------------------------------------------
void read_training_set_data( training_set_data &training_data ) {
//-------------------------------------------------------------------------------------------------------
//
    int counter;
    int pcount;
    std::string substr;
    std::vector<double> tmp;
//
//  Training set data
    training_data.snap_list = {};
    training_data.snap_fmt = {};    
    training_data.n_snap = 0;
    training_data.n_parm = 0;
    training_data.snap_pnts = {}; 
    training_data.parm_name = {};    
    training_data.flow_type = {};     
    training_data.time_ID = 0;    
    training_data.mesh_file = {};
    training_data.mesh_fmt = {};
//
    if ( training_data.inpstn_flow_type != "STEADY" ) {
//
        std::stringstream sstream(training_data.inpstn_flow_type);
        counter = 0;
        while ( sstream.good() ) {
            getline(sstream, substr, ',');
            if ( counter == 0 ) {
                training_data.flow_type = substr;
                counter++; 
            } else if ( counter != 0 && training_data.flow_type == "UNSTEADY" ) {
                training_data.time_ID = std::stoi(substr);
            }
        }
//        
    } else {
        training_data.flow_type = "STEADY";
    }
//
    if ( training_data.inpstn_snapshot_list.empty() ) {
            std::cout << " " << std::endl;
            std::cout << "-> Error: The keyword SNAPSHOT_LIST is missing in configuration file." << std::endl;        
            std::cout << " " << std::endl;
            exit (EXIT_FAILURE);
    } else {
//
        std::stringstream sstream(training_data.inpstn_snapshot_list);
        counter = 0;
        while ( sstream.good() ) {
            getline(sstream, substr, ',');
//
            if ( is_number(substr) ) {
                tmp.push_back(std::stod(substr));
                counter++;
            } else {
                training_data.n_parm = counter;
                training_data.snap_pnts.push_back(tmp);            
                training_data.snap_list.push_back(substr);
                tmp.clear();
                counter = 0;
            }
        }
        training_data.n_snap = training_data.snap_list.size();
    }
//
    training_data.snap_fmt = training_data.inpstn_snapshot_fmt;
//
    if ( training_data.inpstn_name_parm.empty() ) {
        for ( int i = 0; i < training_data.n_parm; i++ ){
            substr = "Field ";
            training_data.parm_name.push_back(substr);
        }
//
    } else {
//
        std::stringstream sstream(training_data.inpstn_name_parm);
        while(sstream.good()) {
            getline(sstream, substr, ',');
            training_data.parm_name.push_back(substr);
        }
//
    }
//
    if ( training_data.inpstn_mesh_file != "NONE" ) {
        training_data.mesh_file = training_data.inpstn_mesh_file;
//
        if ( training_data.inpstn_mesh_fmt.empty() ) {
//
            std::cout << " " << std::endl;
            std::cout << "-> Error: The keyword MESH_FMT is missing in configuration file." << std::endl;        
            std::cout << " " << std::endl;
            exit (EXIT_FAILURE);
//
        } else { training_data.mesh_fmt = training_data.inpstn_mesh_fmt; }

    }
//
}
//
//-------------------------------------------------------------------------------------------------------
void read_modal_identification_data( modal_identification_data &modal_data ) {
//-------------------------------------------------------------------------------------------------------
//
    std::string substr;
//
//  Modal identification data
    modal_data.SPOD_fs = 0;                    
    modal_data.DMD_rank = 0; 
    modal_data.MRDMD_max_cycles = 0;
    modal_data.MRDMD_max_levels = 0;
    modal_data.RDMD_rank = 0;
    modal_data.SPOD_sigma = 0;               
    modal_data.SPOD_ft = {};             
    modal_data.SPOD_flag_bc = {};            
    modal_data.DMD_coef_flag = {};
//
    modal_data.SPOD_fs = std::stoi(modal_data.inpstn_spod_filter_size);
    modal_data.SPOD_ft = modal_data.inpstn_spod_filter_type;
    modal_data.SPOD_sigma = std::stod(modal_data.inpstn_spod_sigma);        
    modal_data.SPOD_flag_bc = modal_data.inpstn_spod_flag_bc;
//
    modal_data.DMD_rank = std::stoi(modal_data.inpstn_dmd_rank);
    modal_data.DMD_coef_flag = modal_data.inpstn_dmd_coef_flag;
//    
    modal_data.MRDMD_max_levels = std::stoi(modal_data.inpstn_mr_dmd_max_levels);        
    modal_data.MRDMD_max_cycles = std::stoi(modal_data.inpstn_mr_dmd_max_cycles);
//    
    modal_data.RDMD_rank = std::stoi(modal_data.inpstn_rdmd_rank);
//
}
//
//-------------------------------------------------------------------------------------------------------
void read_manifold_learning_data( manifold_learning_data &manifold_data ) {
//-------------------------------------------------------------------------------------------------------
//
//  Manifold learning data
    manifold_data.manifold_meth = {};    
//
    manifold_data.manifold_meth = manifold_data.inpstn_manifold_meth;
//
}
//
//-------------------------------------------------------------------------------------------------------
void read_error_data( ld_error_data &error_data ) {
//-------------------------------------------------------------------------------------------------------
//
    int counter;
    int pcount;
    std::string substr;
    std::vector<double> tmp;
//
//  Adaptive reduction data
    error_data.n_err_pnts = 0;
    error_data.dtBDF = 0.0;    
    error_data.err_pnts = {};
    error_data.err_fmla = {};
    error_data.err_sol_fmt = {};
    error_data.err_sol_list = {};
//
    error_data.err_fmla = error_data.inpstn_error_fmla;
//
    if ( error_data.inpstn_error_points.empty() ) {
        std::cout << " " << std::endl;
        std::cout << "-> Error: The keyword ERROR_EVAL_POINTS is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
    } else {
//
        std::stringstream sstream(error_data.inpstn_error_points);
        counter = 0;
        while ( sstream.good() ) {
            getline(sstream, substr, ',');
//
            if ( is_number(substr) ) {
                tmp.push_back(std::stod(substr));
                counter++;
            } else {
                error_data.n_err_pnts = counter;
                error_data.err_pnts.push_back(tmp);            
                error_data.err_sol_list.push_back(substr);
                tmp.clear();
                counter = 0;
            }
        }
        error_data.n_err_pnts = error_data.err_sol_list.size();
    }
//
    error_data.err_sol_fmt = error_data.inpstn_error_sol_fmt;
    if ( error_data.err_fmla == "RESIDUAL" ) error_data.dtBDF = std::stod(error_data.inpstn_resid_dt_bdf);
//
}
//
//-------------------------------------------------------------------------------------------------------
void read_target_set_data( target_set_data &target_data ) {
//-------------------------------------------------------------------------------------------------------
//
    int counter;
    int pcount;
    std::string substr;
    std::vector<double> tmp;
//
//  Target set data
    target_data.n_targ = 0;
    target_data.n_modes = 0;
    target_data.n_flds = 0;    
    target_data.n_parm = 0; 
    target_data.en = 0.0;                 
    target_data.target_pnts = {};
    target_data.target_fmt = {};    
    target_data.recon_fmla = {};
    target_data.coeff_fmla = {};
    target_data.modes_sel_meth = {};
    target_data.target_list = {};
    target_data.target_ldmodel = {};
//
    if ( target_data.inpstn_target_list.empty() ) {
        std::cout << " " << std::endl;
        std::cout << "-> Error: The keyword TARGET_LIST is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
    } else {
//
        std::stringstream sstream(target_data.inpstn_target_list);
        counter = 0;
        pcount = 0;
        while ( sstream.good() ) {
            getline(sstream, substr, ',');
//
            if ( is_number(substr) ) {
                tmp.push_back(std::stod(substr));
                counter++;
            } else {
//
                if ( pcount == 0 ) { 
                    target_data.n_parm = counter;
                    target_data.target_pnts.push_back(tmp);            
                    target_data.target_list.push_back(substr);
                    tmp.clear();
                    counter = 0;
                    pcount++;
                } else {
                    target_data.target_ldmodel.push_back(substr);
                    pcount = 0;
                }
//                
            }
        }
        target_data.n_targ = target_data.target_list.size();
    }
//
    { std::stringstream sstream(target_data.inpstn_target_fmt);
    getline(sstream, substr, '#');
    target_data.target_fmt = substr; }
//    
    { std::stringstream sstream(target_data.inpstn_flag_method_coefficients);
    getline(sstream, substr, '#');
    target_data.coeff_fmla = substr; }
//   
    { std::stringstream sstream(target_data.inpstn_modes_selection);
    getline(sstream, substr, '#');
    target_data.modes_sel_meth = substr; }
//
    { std::stringstream sstream(target_data.inpstn_energy_treshold);
    getline(sstream, substr, '#');
    target_data.en = std::stod(substr); }
//   
    { std::stringstream sstream(target_data.inpstn_number_modes);
    getline(sstream, substr, '#');
    target_data.n_modes = std::stoi(substr); }
//
}
//
//-------------------------------------------------------------------------------------------------------
void read_aerodynamic_data( aero_data &aerodynamic_data ) {
//-------------------------------------------------------------------------------------------------------
//
//  Aerodynamic plugin data
    aerodynamic_data.aeroutput = {};
    aerodynamic_data.Mach = 0.0;                  
    aerodynamic_data.Reynolds = 0.0;  
    aerodynamic_data.temp_inf = 0.0;              
    aerodynamic_data.pres_inf = 0.0;
    aerodynamic_data.dens_inf = 0.0;
    aerodynamic_data.alpha = 0.0;                 
    aerodynamic_data.beta = 0.0;                  
    aerodynamic_data.mu = 0.0;
    aerodynamic_data.refA = 1.0;
    aerodynamic_data.refL = 1.0;
    aerodynamic_data.cmaxis = {};
    aerodynamic_data.cmcoord ={};
    aerodynamic_data.target_profile_file = {};
    aerodynamic_data.target_profile_fmt = {};
    aerodynamic_data.Nfom = 0;
//
    std::string substr;    
    std::stringstream sstream(aerodynamic_data.inpstn_aeroutput);
    while(sstream.good()) {
        getline(sstream, substr, ',');
        aerodynamic_data.aeroutput.push_back(substr);
    }
    aerodynamic_data.Nfom = aerodynamic_data.aeroutput.size();
//
    aerodynamic_data.Mach = std::stod(aerodynamic_data.inpstn_mach); 
    aerodynamic_data.Reynolds = std::stod(aerodynamic_data.inpstn_reynolds); 
    aerodynamic_data.temp_inf = std::stod(aerodynamic_data.inpstn_freestream_temperature);
    aerodynamic_data.pres_inf = std::stod(aerodynamic_data.inpstn_freestream_pressure);
    aerodynamic_data.dens_inf = std::stod(aerodynamic_data.inpstn_freestream_density);
    aerodynamic_data.alpha = std::stod(aerodynamic_data.inpstn_alpha);
    aerodynamic_data.beta = std::stod(aerodynamic_data.inpstn_beta);
    aerodynamic_data.mu = std::stod(aerodynamic_data.inpstn_viscosity);
    aerodynamic_data.refA = std::stod(aerodynamic_data.inpstn_reference_area);
//
    if ( std::any_of(aerodynamic_data.aeroutput.begin(), aerodynamic_data.aeroutput.end(), [](std::string a){return a == "CM";}) ) {
//        
        aerodynamic_data.refL = std::stod(aerodynamic_data.inpstn_reference_lenght);
        aerodynamic_data.cmaxis = aerodynamic_data.inpstn_cmaxis;
//        
        std::stringstream sstream(aerodynamic_data.inpstn_cmcoordinates);
        while(sstream.good()) {
            getline(sstream, substr, ',');
            aerodynamic_data.cmcoord.push_back(std::stod(substr));
        }
//
    }
//
    if ( std::any_of(aerodynamic_data.aeroutput.begin(), aerodynamic_data.aeroutput.end(), [](std::string a){return a == "INVDES";}) ) {
//
        if ( aerodynamic_data.inpstn_target_profile.empty() ) {
//    
                std::cout << " " << std::endl;
                std::cout << "-> Error: The keyword TARGET_PROFILE is missing in configuration file." << std::endl;        
                std::cout << " " << std::endl;
                exit (EXIT_FAILURE);
//    
            } else { aerodynamic_data.target_profile_file = aerodynamic_data.inpstn_target_profile; }
//    
            if ( aerodynamic_data.inpstn_target_profile_fmt.empty() ) {
//    
                std::cout << " " << std::endl;
                std::cout << "-> Error: The keyword TARGET_PROFILE_FMT is missing in configuration file." << std::endl;        
                std::cout << " " << std::endl;
                exit (EXIT_FAILURE);
//    
        } else { aerodynamic_data.target_profile_fmt = aerodynamic_data.inpstn_target_profile_fmt; }
//        
    }
//
}
//
// -------------------------------------------------------------------------------------------
int snapshot_size ( const std::string snapshot_file, std::string snapshot_format ) {
// -------------------------------------------------------------------------------------------
//
    std::string line_flow_data;
    std::ifstream flow_data;
//
    flow_data.open(snapshot_file.c_str());
//
    if( !flow_data.is_open() ) {
        std::cout << " " << std::endl;
        std::cout << " ERROR. File: " << snapshot_file << " not found " << std::endl;
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
    }
//
    int snap_size = 0;
//
    switch ( read_snapshot_format(snapshot_format) ) {
//
        case CSV_ASCII: std::cout << "CSV_ASCII not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        case SU2_CSV:
//
            while(getline( flow_data, line_flow_data )) {
                if ( line_flow_data.compare(0,1,"E") == 0 || line_flow_data.compare(0,1,"A") == 0 )
                    break;            
                    snap_size++;
            }; break;
//        
        case SU2_BINARY: std::cout << "-> Error: SU2_BINARY not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        case PLOT3D_BINARY: std::cout << "-> Error: PLOT3D_BINARY not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        case HDF5_BINARY: std::cout << "-> Error: HDF5_BINARY not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        case CGNS: std::cout << "-> Error: CGNS not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        default: std::cout << "-> Error: Unknown snapshot format" << std::endl; exit (EXIT_FAILURE); break;
//
    }
//
    flow_data.close();
    return (snap_size-1);
//
}
//
// --------------------------------------------------------------------------------------------------------------
Eigen::MatrixXd read_field( int snp, int Ni, const std::vector<int> &fields, std::ifstream filestream[], std::string filefmt ) {
// --------------------------------------------------------------------------------------------------------------
//
    std::string line_flow_data;
    std::string substr;
    Eigen::MatrixXd snap_field;
    snap_field = Eigen::MatrixXd::Zero(Ni,fields.size());
//    
    switch ( read_snapshot_format(filefmt) ) {
//
        case CSV_ASCII: { std::cout << "CSV_ASCII snapshot format not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
//        
        case SU2_CSV: {
//
//          Reading header of SU2 file
            getline( filestream[snp], line_flow_data );
            int n = 0;
            while(getline( filestream[snp], line_flow_data )) {
//
                std::stringstream sstream(line_flow_data);
                int counter = 0;            
                for ( int k = 0; k < fields.size(); k++ ) {
//                 
                    while ( sstream.good() ) {
                        getline(sstream, substr, ','); 
                        counter++;
                        if ( counter == fields[k] ) {
                            snap_field(n,k) = std::stod(substr);
                            break;
                        }
                    }
                }
//                         
               n++;
               if ( n == Ni ) {
                   filestream[snp].seekg(0, std::ios::beg);                   
                   break;
               }
            }; break; }
//        
        case SU2_BINARY: std::cout << "-> Error: SU2_BINARY snapshot format not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        case PLOT3D_BINARY: std::cout << "-> Error: PLOT3D_BINARY snapshot format not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        case HDF5_BINARY: std::cout << "-> Error: HDF5_BINARY snapshot format not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        case CGNS: std::cout << "-> Error: CGNS snapshot format not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        default: std::cout << "-> Error: Unknown snapshot format" << std::endl; exit (EXIT_FAILURE); break;
//
    }
    return snap_field;
}
//
// --------------------------------------------------------------------------------------------------------------
Eigen::MatrixXd CMesh::GetCoords(){
// --------------------------------------------------------------------------------------------------------------
    return Coords;
}
//
// --------------------------------------------------------------------------------------------------------------
Eigen::MatrixXi CMesh::GetConnectivity(){
// --------------------------------------------------------------------------------------------------------------
    return Connectivity;
}
//
//// --------------------------------------------------------------------------------------------------------------
//void CMesh::read( std::ifstream &fstream, const std::string &format ) {
//// --------------------------------------------------------------------------------------------------------------
////
//    int ln = 0;
//    std::string line;
//    std::string value;
//    std::string substr;
//    size_t delimiterPos;
////
//    std::cout << std::endl;
//    std::cout << " Reading " << format << " mesh ..." << std::endl;
////
//    switch ( read_mesh_format(format) ) {
////
//        case SU2: {
////
//            while ( getline( fstream, line ) ) {
////
////              Remove leading spaces, ignore hashtag and blank lines
//                line = ltrim(line);
//                if ( line[0] == '#' || line.empty() )
//                    continue;
////
//                delimiterPos = line.find("=");
//                substr = line.substr(0, delimiterPos);
//                value = line.substr(delimiterPos + 1);
////
//                switch ( su2mesh_keyword(substr) ) {
////
//                    case NDIME: { sdim = std::stoi(value); break; }
////
//                    case NELEM: {
////
//                        nelements = std::stoi(value);
//                        std::vector<int> nodes;
////
//                        for ( int i = 0; i < nelements; i++ ) {
////
//                            int p = 0;
//                            getline( fstream, line );
//                            line = ltrim(line);
////
//                            delimiterPos = line.find(" ");
//                            value = line.substr(0, delimiterPos);
//                            std::stringstream sstream(line.substr(delimiterPos + 1));
////
//                            if ( std::stoi(value) == 3 ) {           // Line element
////
//                                while ( sstream.good() ) {
//                                    getline(sstream, substr, ' ');
//                                    nodes.push_back(std::stoi(substr));
//                                    p++;
//                                    if ( p == 2 ) {
//                                        segm.push_back(nodes);
//                                        nodes.clear();
//                                        break;
//                                    }
//                                }
////
//                            } else if ( std::stoi(value) == 5 ) {    // Triangular element
////
//                                while ( sstream.good() ) {
//                                    getline(sstream, substr, ' ');
//                                    nodes.push_back(std::stoi(substr));
//                                    p++;
//                                    if ( p == 3 ) {
//                                        tria.push_back(nodes);
//                                        nodes.clear();
//                                        break;
//                                    }
//                                }
////
//                            } else if ( std::stoi(value) == 9 ) {    // Quadrilateral element
////
//                                while ( sstream.good() ) {
//                                    getline(sstream, substr, ' ');
//                                    nodes.push_back(std::stoi(substr));
//                                    p++;
//                                    if ( p == 4 ) {
//                                        quad.push_back(nodes);
//                                        nodes.clear();
//                                        break;
//                                    }
//                                }
////
//                            } else if ( std::stoi(value) == 10 ) {   // Tetrahedron
////
//                                while ( sstream.good() ) {
//                                    getline(sstream, substr, ' ');
//                                    nodes.push_back(std::stoi(substr));
//                                    p++;
//                                    if ( p == 4 ) {
//                                        tetra.push_back(nodes);
//                                        nodes.clear();
//                                        break;
//                                    }
//                                }
////
//                            } else if ( std::stoi(value) == 12 ) {   // Hexahedron
////
//                                while ( sstream.good() ) {
//                                    getline(sstream, substr, ' ');
//                                    nodes.push_back(std::stoi(substr));
//                                    p++;
//                                    if ( p == 8 ) {
//                                        hexa.push_back(nodes);
//                                        nodes.clear();
//                                        break;
//                                    }
//                                }
////
//                            } else if ( std::stoi(value) == 13 ) {   // Prism
////
//                                while ( sstream.good() ) {
//                                    getline(sstream, substr, ' ');
//                                    nodes.push_back(std::stoi(substr));
//                                    p++;
//                                    if ( p == 6 ) {
//                                        prism.push_back(nodes);
//                                        nodes.clear();
//                                        break;
//                                    }
//                                }
////
//                            } else if ( std::stoi(value) == 14 ) {   // Pyramid
////
//                                while ( sstream.good() ) {
//                                    getline(sstream, substr, ' ');
//                                    nodes.push_back(std::stoi(substr));
//                                    p++;
//                                    if ( p == 5 ) {
//                                        pyra.push_back(nodes);
//                                        nodes.clear();
//                                        break;
//                                    }
//                                }
////
//                            }
////
//                        }
//                        break; }
////
//                    case NPOIN: {
////
//                        npoints = std::stoi(value);
//                        std::vector<double> coords;
//                        int p;
////
//                        for ( int i = 0; i < npoints; i++ ) {
////
//                            getline( fstream, line );
//                            line = ltrim(line);
//                            std::stringstream sstream(line);
////
//                            p = 0;
//                            while ( sstream.good() ) {
//                                getline(sstream, substr, ' ');
//                                coords.push_back(std::stod(substr));
//                                p++;
//                                if ( p == sdim ) {
//                                    coordinates.push_back(coords);
//                                    coords.clear();
//                                    break;
//                                }
//                            }
//                        }
////
//                        break; }
////
//                    case NMARK: { nboundaries = std::stoi(value); break; }
////
//                    case MARKER_TAG: { break; }
//                    case MARKER_ELEMS: { break; }
//                    default: { continue; break; }
//                }
//            }
//            break; }
////
//        case CGNSM: { std::cout << "-> Error: CGNS mesh format not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
//        case VTK:  { std::cout << "-> Error: VTK mesh format not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
////
//    }
//
//std::cout << npoints << " " << nelements << " " << nboundaries << " " << tetra.size() << " " << hexa.size() << " " << pyra.size() << std::endl;
//
//}
////
// --------------------------------------------------------------------------------------------------------------
void CMesh::read( std::ifstream &fstream, const std::string &format ) {
// --------------------------------------------------------------------------------------------------------------
//
    int ln = 0;
    std::string line;
    std::string value;
    std::string substr;
    size_t delimiterPos;
//
    std::cout << std::endl;
    std::cout << " Reading " << format << " mesh ..." << std::endl;
//
    switch ( read_mesh_format(format) ) {
//
        case SU2: {
//
            while ( getline( fstream, line ) ) {
//
//              Remove leading spaces, ignore hashtag and blank lines
                line = ltrim(line);
                if ( line[0] == '#' || line.empty() )
                    continue;
//
                delimiterPos = line.find("=");
                substr = line.substr(0, delimiterPos);
                value = line.substr(delimiterPos + 1);
//
                switch ( su2mesh_keyword(substr) ) {
//
                    case NDIME: { sdim = std::stoi(value); break; }
//
                    case NELEM: {
//
                        nelements = std::stoi(value);
                        if ( sdim == 2 ) Connectivity = (-1)*Eigen::MatrixXi::Ones(nelements,6);
                        else Connectivity = (-1)*Eigen::MatrixXi::Ones(nelements,10);
//
                        for ( int i = 0; i < nelements; i++ ) {
//
                            getline( fstream, line );
                            line = ltrim(line);
                            std::stringstream sstream(line);
//
                            int p = 0;
                            while ( sstream.good() ) {
                                getline(sstream, substr, ' ');
                                Connectivity(i,p) = std::stod(substr);
                                p++;
                            }
                        }
//
                        break; }
//
                    case NPOIN: {
//
                        npoints = std::stoi(value);
                        Coords = Eigen::MatrixXd::Zero(npoints,sdim);
//
                        for ( int i = 0; i < npoints; i++ ) {
//
                            getline( fstream, line );
                            line = ltrim(line);
                            std::stringstream sstream(line);
//
                            for ( int j = 0; j < sdim; j++ ) {
                                getline(sstream, substr, ' ');
                                Coords(i,j) = std::stod(substr);
                            }
                        }
//
                        break; }
//
                    case NMARK: { nboundaries = std::stoi(value); break; }
//
                    case MARKER_TAG: { break; }
                    case MARKER_ELEMS: { break; }
                    default: { continue; break; }
                }
            }
            break; }
//
        case CGNSM: { std::cout << "-> Error: CGNS mesh format not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
        case VTK:  { std::cout << "-> Error: VTK mesh format not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
//
    }
}
//
// --------------------------------------------------------------------------------------------------------------
Eigen::MatrixXd read_modes( std::string filename, int Nr, int r_RDMD ) {
// --------------------------------------------------------------------------------------------------------------
//    
    Eigen::MatrixXd f(Nr, r_RDMD);

    std::ifstream Modes_data;
    Modes_data.open( filename );

        if ( !Modes_data.is_open() )
    {
        std::cout << "File : " << filename << " not found" << std::endl;    
        exit (EXIT_FAILURE);
    }

    std::string line_flow_data ;

    // Read row of headers
    getline( Modes_data, line_flow_data );

    int n_row = 0;

    while ( getline( Modes_data, line_flow_data ) )
    {

        Eigen::RowVectorXd point(r_RDMD);
        std::istringstream iss(line_flow_data);
        std::string token;
        double rubbish;
        int count = 0, c = 0; 

        while( getline( iss, token, ' ') )
        {
            rubbish = std::stod(token);
            if ( count < r_RDMD )
                point(count) = rubbish;
            
            count ++;
        } 

        f.row(n_row) = point; 
        n_row++;

    }

    Modes_data.close();

    return f;

}


Eigen::MatrixXd read_coefs( std::string filename, int Ns, int r_RDMD )
{

    Eigen::MatrixXd f(r_RDMD, Ns);

    std::ifstream Coefs_data;
    Coefs_data.open( filename );

    if ( !Coefs_data.is_open() )
    {
        std::cout << "File : " << filename << " not found" << std::endl;    
        exit (EXIT_FAILURE);
    }

    std::string line_flow_data ;

    // Read row of headers
    getline( Coefs_data, line_flow_data );

    int n_row = 0;

    int count = 0;
    while ( getline( Coefs_data, line_flow_data ) && n_row < r_RDMD )
    {

        Eigen::RowVectorXd point(Ns);
        std::istringstream iss(line_flow_data);
        std::string token;
        double rubbish;
        int count = 0, c = 0; 

        while( getline( iss, token, ' ') && count < Ns )
        {
            rubbish = std::stod(token);
            point(count) = rubbish;
            
            count ++;
        } 

        f.row(n_row) = point; 
        n_row++;

    }

    Coefs_data.close();

    return f;
}



Eigen::MatrixXd read_err_j ( std::string filename, int Ns )
{
        std::ifstream file_data;
        file_data.open( filename );

            if ( !file_data.is_open() )
        {
            std::cout << "File : " << filename << " not found" << std::endl;    
            exit (EXIT_FAILURE);
        }

        std::string line_flow_data ;

        int n_row = 0, count = 0;
        Eigen::MatrixXd Err_map = Eigen::MatrixXd::Zero(Ns, Ns); 

        while ( getline( file_data, line_flow_data ) )
        {

            
            std::istringstream iss(line_flow_data);
            std::string token;
            double err;
            count = 0; 

            while( getline( iss, token, '\t') )
            {
                err = std::stod(token);
                Err_map(n_row, count) = err;

                count ++;
            } 
 
            n_row++;

        }

        file_data.close();

        return Err_map;

}
//
plot3d_info read_plot3d_info (std::string filename)
{

    plot3d_info Info;
    int nC = 5; // Density, Momentum (3 Components), Total Energy

    // open file in binary mode
    std::ifstream input(filename, std::ios::binary);

    if (input.good())
    {
     // read number of points
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        input.read((char*) &Info.nblocks, sizeof(int));
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);

        Info.ni.resize(Info.nblocks);
        Info.nj.resize(Info.nblocks);
        Info.nk.resize(Info.nblocks);
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        
        for ( int i = 0; i < Info.nblocks; i++)
        {
            input.read((char*) &Info.ni[i], sizeof(int));
            input.read((char*) &Info.nj[i], sizeof(int));
            input.read((char*) &Info.nk[i], sizeof(int));
        }

        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        input.read((char*) &Info.Ma, sizeof(float));
        input.read((char*) &Info.alpha, sizeof(float));
        input.read((char*) &Info.Re, sizeof(float));
        input.read((char*) &Info.time, sizeof(float));

        input.close();
        // std::cout << "Successfully read plot3d Info" << std::endl;
    }
    else
    {
     std::cout << "Unable to open file .q" << std::endl;
    }

    return Info;
}


std::vector<Eigen::VectorXd> read_plot3d (std::string filename, plot3d_info Info)
{

    int nC = 5; // Density, Momentum (3 Components), Total Energy
    // Read .q file fortran binary

    float dum_f;
    int dum_i;
    // open file in binary mode
    std::ifstream input(filename, std::ios::binary);

    std::vector<Eigen::VectorXd> flow_data(Info.nblocks);

    if (input.good())
    {
     // read number of points
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        input.read((char*) &dum_i, sizeof(int));
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);

        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        
        for ( int i = 0; i < Info.nblocks; i++)
        {

            input.read((char*) &dum_i, sizeof(int));
            input.read((char*) &dum_i, sizeof(int));
            input.read((char*) &dum_i, sizeof(int));
        
        }

        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);

        for ( int i = 0; i < Info.nblocks; i++ )
        {
            Eigen::VectorXf data_i(nC*Info.ni[i]*Info.nj[i]*Info.nk[i]);
            input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
            input.read((char*) &dum_f, sizeof(float));
            input.read((char*) &dum_f, sizeof(float));
            input.read((char*) &dum_f, sizeof(float));
            input.read((char*) &dum_f, sizeof(float));
            input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);

            input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
            
            flow_data[i] = Eigen::VectorXd::Zero(nC*Info.ni[i]*Info.nj[i]*Info.nk[i]);
            for ( int j = 0; j < nC*Info.ni[i]*Info.nj[i]*Info.nk[i]; j++ )
                input.read((char*) &data_i(j), sizeof(float));

            flow_data[i] = data_i.cast<double>();
            input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);

        }

        input.close();
        // std::cout << "Successfully read .q file" << std::endl;
        return flow_data;

    }
    else
    {
     std::cout << "Unable to open .q file" << std::endl;
     std::cout << "Terminating ... " << std::endl;
     exit (EXIT_FAILURE);
    }


}
//
// -------------------------------------------------------------------------------------------
void modify_su2_cfg ( std::string file_in, std::string file_out, double dt_res ) {
// -------------------------------------------------------------------------------------------
//
    std::ifstream inFile ( file_in );
    std::ofstream outFile ( file_out );
//
    if ( inFile.is_open() && outFile.is_open() ) {
//
        size_t delimiterPos;
        std::string name, value;
        std::string line;
//
        while (getline(inFile, line)) {
//
            line.erase(remove_if(line.begin(), line.end(), isspace),
                       line.end());  //include ::  in front of isspace if using namespace std
            if (line[0] == '#' || line.empty())
                continue;
//
            delimiterPos = line.find("=");
            name = line.substr(0, delimiterPos);
//
            if ( name == "UNST_TIMESTEP" )
                outFile << "UNST_TIMESTEP=" << std::setprecision(16) << dt_res;
            else
                outFile << line;

            outFile << std::endl;
//
        }
//
        outFile.close();
        inFile.close();    
//
    } else {
//        
        std::cout << " " << std::endl;
        std::cout << "ERROR: Unable to open SU2 configuration file." << std::endl;
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
    }
}

