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
bool is_number(const std::string& s) {
// -------------------------------------------------------------------------------------------
    long double ld;
    return((std::istringstream(s) >> ld >> std::ws).eof());
}
//
inline void wait_on_enter()
{
    std::string dummy;
    std::cout << "Enter to continue..." << std::endl;
    std::getline(std::cin, dummy);
}
//
// -------------------------------------------------------------------------------------------
keywords read_keyword_type( const std::string &key_string ) {
// -------------------------------------------------------------------------------------------    
//
    if( key_string == "PROBLEM" )
        return PROBLEM;
    else if( key_string == "LOW_DIMENSIONAL_APPROACH" )
        return LOW_DIMENSIONAL_APPROACH;    
    else if( key_string == "FLOW_TYPE" )
        return FLOW_TYPE;
    else if( key_string == "NAME_PARM" )
        return NAME_PARM;
    else if( key_string == "SNAPSHOT_LIST" )
        return SNAPSHOT_LIST;    
    else if( key_string == "SNAPSHOT_FMT" )
        return SNAPSHOT_FMT;
    else if( key_string == "FLAG_FIELDS" )
        return FLAG_FIELDS;
    else if( key_string == "FIELDS_NAME" )
        return FIELDS_NAME;    
    else if( key_string == "ERROR_FMLA" )
        return ERROR_FMLA;
    else if( key_string == "RESID_EVAL_POINTS" )
        return RESID_EVAL_POINTS;
    else if( key_string == "RESID_DT_BDF" )
        return RESID_DT_BDF;
    else if( key_string == "REF_SOLUTION_LIST" )
        return REF_SOLUTION_LIST;
    else if( key_string == "REF_SOLUTION_FMT" )
        return REF_SOLUTION_FMT;
    else if( key_string == "ERROR_FMLA" )
        return ERROR_FMLA;
    else if( key_string == "MODES_EXTRACTION" )
        return MODES_EXTRACTION;
    else if( key_string == "LOW_DIMENSIONAL_MODEL" )
        return LOW_DIMENSIONAL_MODEL;
    else if( key_string == "FLAG_REF_FIELD" )
        return FLAG_REF_FIELD;
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
    else if( key_string == "TARGET_LIST" )
        return TARGET_LIST;
    else if( key_string == "TARGET_FMT" )
        return TARGET_FMT;
    else if( key_string == "USE_LOW_DIMENSIONAL_MODEL" )
        return USE_LOW_DIMENSIONAL_MODEL;
    else if( key_string == "FLAG_METHOD_RECONSTRUCTION" )
        return FLAG_METHOD_RECONSTRUCTION;
    else if( key_string == "FLAG_METHOD_COEFFICIENTS" )
        return FLAG_METHOD_COEFFICIENTS;
    else if( key_string == "MODES_SELECTION" )
        return MODES_SELECTION;
    else if( key_string == "ENERGY_TRESHOLD" )
        return ENERGY_TRESHOLD;
    else if( key_string == "NUMBER_MODES" )
        return NUMBER_MODES;
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
    else if( key_string == "MESH_FILE" )
        return MESH_FILE;
    else if( key_string == "MESH_FMT" )
        return MESH_FMT;
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
    if ( reftyp == "MEAN_FIELD" )
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
        std::cout << "-> Error: Unknown reference field type" << std::endl;
        std::cout << reftyp << std::endl; 
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//     
    }
}
//
// ----------------------------------------------------------------------------------------
void read_config ( const std::string filename, prob_settings &settings ) {
// ----------------------------------------------------------------------------------------
//
//  Initializing values in settings to NULL or default values    
    std::string inpstn_problem = {};
    std::string inpstn_flow_type = "STEADY";       
    std::string inpstn_name_parm = {};
    std::string inpstn_lowdim_approach = {};
    std::string inpstn_error_fmla = "RESIDUAL";
    std::string inpstn_resid_eval_points = {};
    std::string inpstn_resid_dt_bdf = "1.e-3";
    std::string inpstn_ref_solution_list = {};
    std::string inpstn_ref_solution_fmt = "CSV_ASCII";  
    std::string inpstn_snapshot_list = {};
    std::string inpstn_snapshot_fmt = "CSV_ASCII";
    std::string inpstn_flag_fields = {};                  
    std::string inpstn_fields_name = {};                  
    std::string inpstn_modes_extraction = "POD";
    std::string inpstn_lowdim_model_name = {};
    std::string inpstn_ref_field = "NONE";
    std::string inpstn_spod_filter_size = "-1";
    std::string inpstn_spod_filter_type = "BOX";
    std::string inpstn_spod_sigma = "1.0";
    std::string inpstn_spod_flag_bc = "ZERO";
    std::string inpstn_dmd_rank = "-2";
    std::string inpstn_dmd_coef_flag = "OPT";
    std::string inpstn_mr_dmd_max_levels = "-1";
    std::string inpstn_mr_dmd_max_cycles = "-1";
    std::string inpstn_rdmd_rank = "-1";
    std::string inpstn_target_list = {};
    std::string inpstn_target_fmt = "CSV_ASCII";
    std::string inpstn_use_lowdim_model = {};
    std::string inpstn_flag_method_reconstruction = "LINEAR_COMBINATION";
    std::string inpstn_flag_method_coefficients = "RBF_LINEAR";
    std::string inpstn_modes_selection = "ENERGY";
    std::string inpstn_energy_treshold = "0.95";
    std::string inpstn_number_modes = "1";                
    std::string inpstn_freestream_density = "1.225";
    std::string inpstn_freestream_pressure = "101325";
    std::string inpstn_freestream_temperature = "288.15";                 
    std::string inpstn_alpha = "0.0"; 
    std::string inpstn_beta = "0.0"; 
    std::string inpstn_mach = "0.3"; 
    std::string inpstn_reynolds = "100";                 
    std::string inpstn_viscosity = "1.78e-5";
    std::string inpstn_mesh_file = {};
    std::string inpstn_mesh_fmt = {};
//
    settings.n_snap = 0;
    settings.n_targ = 0;
    settings.n_parm = 0;
    settings.n_flds = 0;    
    settings.SPOD_fs = 0;                    
    settings.DMD_rank = 0; 
    settings.MRDMD_max_cycles = 0;
    settings.MRDMD_max_levels = 0;
    settings.RDMD_rank = 0;    
    settings.n_modes = 0;    
//
    settings.fields = {};
//
    settings.en = 0;                 
    settings.SPOD_sigma = 0;               
    settings.Mach = 0.0;                  
    settings.Reynolds = 0.0;  
    settings.temp_inf = 0.0;              
    settings.pres_inf = 0.0;
    settings.dens_inf = 0.0;
    settings.alpha = 0.0;                 
    settings.beta = 0.0;                  
    settings.mu = 0.0;
//
    settings.dt_res = {};     
    settings.snap_pnts = {}; 
    settings.targ_pnts = {}; 
    settings.res_eval_pnts = {};     
    settings.dir_eval_pnts = {};         
//
    settings.prob_type = {};
    settings.flow_type = {};
    settings.snap_fmt = {};
    settings.targ_fmt = {};
    settings.refsol_fmt = {};
    settings.field_type = {};
    settings.lowdim_approach = {};
    settings.err_fmla = {};
    settings.lowdim_model_name = {};
    settings.use_lowdim_model = {};
    settings.SPOD_ft = {};             
    settings.SPOD_flag_bc = {};
    settings.ref_field = {};          
    settings.DMD_coef_flag = {};
    settings.recon_fmla = {};
    settings.coeff_fmla = {};
    settings.modes_sel_meth = {};
    settings.mesh_file = {};
    settings.mesh_fmt = {};
//
    settings.modes_extraction = {};
    settings.parm_name = {};  
    settings.flds_name = {};  
    settings.snap_list = {};  
    settings.refsol_list = {};
    settings.target_list = {};
//
//  Parsing configuration file
    std::ifstream cFile ( filename );
    if ( cFile.is_open() ) {
//
        size_t delimiterPos; 
        std::string name, value;
        std::string line;
//
//  Reading each line of the config file into a string that will be parsed later on to avoid issues with dependencies
//  or have to introduce multiple times the same information. In this way the information in the configuration file
//  can be entered in whatsoever order
        while( getline(cFile, line) ) {
//            
            line.erase( remove_if(line.begin(), line.end(), isspace), line.end() );  //include ::  in front of isspace if using namespace std
            if ( line[0] == '#' || line.empty() )
                continue;
            delimiterPos = line.find("=");
            name = line.substr(0, delimiterPos);
            value = line.substr(delimiterPos + 1);
//
            switch (read_keyword_type(name)) {
                case PROBLEM: { inpstn_problem = value; break; }
                case LOW_DIMENSIONAL_APPROACH: { inpstn_lowdim_approach = value; break; }
                case FLOW_TYPE: { inpstn_flow_type = value; break; }                
                case NAME_PARM: { inpstn_name_parm = value; break; }                 
                case SNAPSHOT_LIST: { inpstn_snapshot_list = value; break; }
                case SNAPSHOT_FMT: { inpstn_snapshot_fmt = value; break; }
                case FLAG_FIELDS: { inpstn_flag_fields = value; break; }
                case FIELDS_NAME: { inpstn_fields_name = value; break; }
                case ERROR_FMLA: { inpstn_error_fmla = value; break; }
                case RESID_EVAL_POINTS: { inpstn_resid_eval_points = value; break; }
                case RESID_DT_BDF: { inpstn_resid_dt_bdf = value; break; }
                case REF_SOLUTION_LIST: { inpstn_ref_solution_list = value; break; }
                case REF_SOLUTION_FMT: { inpstn_ref_solution_fmt = value; break; }
                case MODES_EXTRACTION: { inpstn_modes_extraction = value; break; }
                case LOW_DIMENSIONAL_MODEL: { inpstn_lowdim_model_name = value; break; }
                case FLAG_REF_FIELD: { inpstn_ref_field = value; break; }
                case SPOD_FILTER_SIZE: { inpstn_spod_filter_size = value; break; }
                case SPOD_FILTER_TYPE: { inpstn_spod_filter_type = value; break; }
                case SPOD_SIGMA: { inpstn_spod_sigma = value; break; }
                case SPOD_FLAG_BC: { inpstn_spod_flag_bc = value; break; }
                case DMD_RANK: { inpstn_dmd_rank = value; break; }
                case DMD_COEFF_FLAG: { inpstn_dmd_coef_flag = value; break; }
                case MR_DMD_MAX_LEVELS: { inpstn_mr_dmd_max_levels = value; break; }
                case MR_DMD_MAX_CYCLES: { inpstn_mr_dmd_max_cycles = value; break; }
                case RDMD_RANK: { inpstn_rdmd_rank = value; break; }
                case TARGET_LIST: { inpstn_target_list = value; break; }
                case TARGET_FMT: { inpstn_target_fmt = value; break; }
                case USE_LOW_DIMENSIONAL_MODEL: { inpstn_use_lowdim_model = value; break; }
                case FLAG_METHOD_RECONSTRUCTION: { inpstn_flag_method_reconstruction = value; break; }
                case FLAG_METHOD_COEFFICIENTS: { inpstn_flag_method_coefficients = value; break; }
                case MODES_SELECTION: { inpstn_modes_selection = value; break; }
                case ENERGY_TRESHOLD: { inpstn_energy_treshold = value; break; }
                case NUMBER_MODES: { inpstn_number_modes = value; break; }
                case FREESTREAM_DENSITY: { inpstn_freestream_density = value; break; }
                case FREESTREAM_PRESSURE: { inpstn_freestream_pressure = value; break; }
                case FREESTREAM_TEMPERATURE: { inpstn_freestream_temperature = value; break; }
                case ALPHA: { inpstn_alpha = value; break; }
                case BETA: { inpstn_beta = value; break; }
                case MACH: { inpstn_mach = value; break; }
                case REYNOLDS: { inpstn_reynolds = value; break; }     
                case VISCOSITY: { inpstn_viscosity = value; break; }
                case MESH_FILE: { inpstn_mesh_file = value; break; }
                case MESH_FMT: { inpstn_mesh_fmt = value; break; }
                default: { break; }
            }
        }
        cFile.close();
//
    } else {
//
        std::cout << " " << std::endl;
        std::cout << "ERROR: Unable to open configuration file." << std::endl;
        std::cout << filename  << std::endl;
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
    }
//
//  Parsing each input string and assign values to variables and arrays
    int counter;
    int pcount;
    std::string substr;
    std::vector<double> tmp;
//
//  General problem definitions ---------------------------------------------------------
    if ( inpstn_problem.empty() ) {
//
        std::cout << " " << std::endl;
        std::cout << "ERROR: The keyword PROBLEM is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
    } else { settings.prob_type = inpstn_problem; }
//
    settings.flow_type = inpstn_flow_type;
//
    if ( inpstn_name_parm.empty() ) {
//
        std::cout << " " << std::endl;
        std::cout << "ERROR: The keyword NAME_PARAM is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
    } else {
//
        std::stringstream sstream(inpstn_name_parm);
        while(sstream.good()) {
            getline(sstream, substr, ',');
            settings.parm_name.push_back(substr);
        }
        settings.n_parm = settings.parm_name.size();
//
    }
//
//  Adaptive framework parameters --------------------------------------------------------- 
    settings.lowdim_approach = inpstn_lowdim_approach;
    if ( settings.lowdim_approach == "ADAPTIVE" ) {
//
        settings.err_fmla = inpstn_error_fmla;
//
        if ( settings.err_fmla == "RESIDUAL" ) {
//            
            std::stringstream sstream(inpstn_resid_eval_points);
            counter = 0;
            pcount = 0;
            while ( sstream.good() ) {
                getline(sstream, substr, ',');
//
                if ( counter == 0 ) {
                    settings.n_eval_pnts_res = std::stoi(substr);
                    counter++; 
                } else {
                    if ( pcount < settings.n_parm-1 ) {
                        tmp.push_back(std::stod(substr));
                        pcount++;
                    } else if ( pcount == settings.n_parm-1 ) {
                        tmp.push_back(std::stod(substr));
                        settings.res_eval_pnts.push_back(tmp);
                        pcount ++;
                        tmp.clear();
                    }
                }
            }
            settings.dt_res = std::stod(inpstn_resid_dt_bdf);
//
        } else if ( settings.err_fmla == "DIRECT" ) {
//
            std::stringstream sstream(inpstn_ref_solution_list);
            counter = 0;
            pcount = 0;
            while ( sstream.good() ) {
                getline(sstream, substr, ',');
//
                if ( counter == 0 ) {
                    settings.n_eval_pnts_dir = std::stoi(substr);
                    counter++; 
                } else {
                    if ( pcount < settings.n_parm-1 ) {
                        tmp.push_back(std::stod(substr));
                        pcount++;
                    } else if ( pcount == settings.n_parm-1 ) {
                        tmp.push_back(std::stod(substr));
                        settings.dir_eval_pnts.push_back(tmp);
                        pcount ++;
                        tmp.clear();
                    } else {
                        settings.refsol_list.push_back(substr);
                        pcount = 0;
                    }
                }
            }
            settings.refsol_fmt = inpstn_ref_solution_fmt;
//            
        }    
//    
    }    
//
//  Basis extraction parameters ---------------------------------------------------------
    if ( inpstn_snapshot_list.empty() ) {
//
        std::cout << " " << std::endl;
        std::cout << "ERROR: The keyword SNAPSHOT_LIST is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
    } else {
//
        std::stringstream sstream(inpstn_snapshot_list);
        counter = 0;
        pcount = 0;
        while ( sstream.good() ) {
            getline(sstream, substr, ',');
//
            if ( counter == 0 ) {
                settings.n_snap = std::stoi(substr);
                counter++; 
            } else {
                if ( pcount < settings.n_parm-1 ) {
                    tmp.push_back(std::stod(substr));
                    pcount++;
                } else if ( pcount == settings.n_parm-1 ) {
                    tmp.push_back(std::stod(substr));
                    settings.parm.push_back(tmp);
                    pcount ++;
                    tmp.clear();
                } else {
                    settings.snap_list.push_back(substr);
                    pcount = 0;
                }
            }
        }
    }
//
    settings.snap_fmt = inpstn_snapshot_fmt;
//
    if ( inpstn_flag_fields.empty() ) {
//
        std::cout << " " << std::endl;
        std::cout << "ERROR: The keyword FLAG_FIELDS is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
    } else {
//
        std::stringstream sstream(inpstn_flag_fields);
        counter = 0;                
        while ( sstream.good() ) {
            getline(sstream, substr, ','); 
//
            if ( counter == 0 ) {
                settings.field_type = substr;
                counter ++;
            } else {
                settings.fields.push_back(std::stoi(substr));
            }
//
        }
        sort(settings.fields.begin(), settings.fields.end());
//        
        if ( settings.field_type == "VECTOR") {
            settings.n_flds = 1;
        } else {
            settings.n_flds = settings.fields.size();
        }
    }
//
    if ( !inpstn_fields_name.empty() ) {
//
        std::stringstream sstream(inpstn_fields_name);
        while ( sstream.good() ) {
            getline(sstream, substr, ',');
            settings.flds_name.push_back(substr);
        }
    }
//
    if ( inpstn_modes_extraction.empty() ) {
//
        std::cout << " " << std::endl;
        std::cout << "ERROR: The keyword MODES_EXTRACTION is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
    } else {
//
        std::stringstream sstream(inpstn_modes_extraction);
        while ( sstream.good() ) {
            getline(sstream, substr, ',');
            settings.modes_extraction.push_back(substr);
        }
        settings.n_extr_meth = settings.modes_extraction.size();
    }
//
    if ( inpstn_lowdim_model_name.empty() ) {
//
        std::cout << " " << std::endl;
        std::cout << "ERROR: The keyword LOW_DIMENSIONAL_MODEL is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
    } else {
        settings.lowdim_model_name = inpstn_lowdim_model_name;
    }
//
    settings.ref_field = inpstn_ref_field;
//
    if ( std::any_of(settings.modes_extraction.begin(), settings.modes_extraction.end(), 
                     [](std::string c){ return ( c == "SPOD" );} ) ) {
        settings.SPOD_fs = std::stoi(inpstn_spod_filter_size);
        settings.SPOD_ft = inpstn_spod_filter_type;
        settings.SPOD_sigma = std::stod(inpstn_spod_sigma);        
        settings.SPOD_flag_bc = inpstn_spod_flag_bc;
    }
//
    if ( std::any_of(settings.modes_extraction.begin(), settings.modes_extraction.end(), 
                     [](std::string c){ return ( c == "DMD" );} ) ) {
        settings.DMD_rank = std::stoi(inpstn_dmd_rank);
        settings.DMD_coef_flag = inpstn_dmd_coef_flag;
        settings.MRDMD_max_levels = std::stoi(inpstn_mr_dmd_max_levels);        
        settings.MRDMD_max_cycles = std::stoi(inpstn_mr_dmd_max_cycles);
    }
//
    if ( std::any_of(settings.modes_extraction.begin(), settings.modes_extraction.end(), 
                     [](std::string c){ return ( c == "RDMD" );} ) ) {
        settings.RDMD_rank = std::stoi(inpstn_rdmd_rank);
    }
//
//  Solution reconstruction parameters ------------------------------------------------------
    if ( inpstn_target_list.empty() && inpstn_problem == "RECONSTRUCTION" ) {
//
        std::cout << " " << std::endl;
        std::cout << "ERROR: The keyword TARGET_LIST is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
        } else {
//
        std::stringstream sstream(inpstn_target_list);
        counter = 0;
        pcount = 0;
        while ( sstream.good() ) {
            getline(sstream, substr, ',');
//
            if ( counter == 0 ) {
                settings.n_targ = std::stoi(substr);
                counter++; 
            } else {
                if ( pcount < settings.n_parm-1 ) {
                    tmp.push_back(std::stod(substr));
                    pcount++;
                } else if ( pcount == settings.n_parm-1 ) {
                    tmp.push_back(std::stod(substr));
                    settings.targ_pnts.push_back(tmp);
                    pcount ++;
                    tmp.clear();
                } else {
                    settings.target_list.push_back(substr);
                    pcount = 0;
                }
            }
        }
    }
//
    settings.targ_fmt = inpstn_target_fmt; 
//
    if ( inpstn_use_lowdim_model.empty() ) {
//
        std::cout << " " << std::endl;
        std::cout << "ERROR: The keyword USE_LOW_DIMENSIONAL_MODEL is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
    } else {
        settings.use_lowdim_model = inpstn_use_lowdim_model;
    }
//
    settings.recon_fmla = inpstn_flag_method_reconstruction;
    settings.coeff_fmla = inpstn_flag_method_coefficients;
    settings.modes_sel_meth = inpstn_modes_selection;
    settings.en = std::stod(inpstn_energy_treshold);
    settings.n_modes = std::stoi(inpstn_number_modes); 
//
    settings.Mach = std::stod(inpstn_mach); 
    settings.Reynolds = std::stod(inpstn_reynolds); 
    settings.temp_inf = std::stod(inpstn_freestream_temperature);
    settings.pres_inf = std::stod(inpstn_freestream_pressure);
    settings.dens_inf = std::stod(inpstn_freestream_density);
    settings.alpha = std::stod(inpstn_alpha);
    settings.beta = std::stod(inpstn_beta);
    settings.mu = std::stod(inpstn_viscosity);
//
    if ( inpstn_mesh_file.empty() && (settings.targ_fmt == "PARAVIEW" || settings.targ_fmt == "TECPLOT") ) {
//
        std::cout << " " << std::endl;
        std::cout << "ERROR: The keyword MESH_FILE is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
    } else { settings.mesh_file = inpstn_mesh_file; }
//
    if ( inpstn_mesh_fmt.empty() && (settings.targ_fmt == "PARAVIEW" || settings.targ_fmt == "TECPLOT") ) {
//
        std::cout << " " << std::endl;
        std::cout << "ERROR: The keyword MESH_FMT is missing in configuration file." << std::endl;        
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//
    } else { settings.mesh_fmt = inpstn_mesh_fmt; }
//
//---------------------------------------------------------------------------------
//  Selection of execution mode:
//  IDENTIFICATION = acquire snapshots and generate basis functions or manifold
//  RECONSTRUCTION = reconstruct a specific set of fields at untried conditions
//                   this step may require an exisitng database of modes/manifold
//                   or may require a generation of basis/manifold on the fly
//---------------------------------------------------------------------------------
//
    if ( settings.prob_type == "GENERATE_LOW_DIMENSIONAL_MODEL" ) {
        settings.generate_lowdim_model = true;
        settings.compute_lowdim_solution = false;
    } else if ( settings.prob_type == "COMPUTE_LOW_DIMENSIONAL_SOLUTION" ) {
        settings.generate_lowdim_model = false;
        settings.compute_lowdim_solution = true;
    }
//
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
        case CSV_ASCII: { std::cout << "CSV_ASCII not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
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
        case SU2_BINARY: std::cout << "-> Error: SU2_BINARY not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        case PLOT3D_BINARY: std::cout << "-> Error: PLOT3D_BINARY not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        case HDF5_BINARY: std::cout << "-> Error: HDF5_BINARY not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        case CGNS: std::cout << "-> Error: CGNS not supported yet" << std::endl; exit (EXIT_FAILURE); break;
        default: std::cout << "-> Error: Unknown snapshot format" << std::endl; exit (EXIT_FAILURE); break;
//
    }
    return snap_field;
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

