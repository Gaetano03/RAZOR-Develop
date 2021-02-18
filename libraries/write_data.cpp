/* ------------------------------------------------------------------------------
\file write_data.cpp
* \brief Print and write defintions.
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
#include "read_data.hpp"
#include "write_data.hpp"
//

/*

// ----------------------------------------------------------------------------------------
void print_gendata ( ld_model_data &lowdim_data, training_set_data &training_data, ld_error_data 
    &error_data, modal_identification_data &modal_data, manifold_learning_data &manifold_data ) {
// ----------------------------------------------------------------------------------------
//
    std::cout << " Generating low dimensional model: " << lowdim_data.lowdim_model_name << std::endl;    
    std::cout << " Reduction strategy: "; 
        for ( int i=0; i<lowdim_data.reduction_strategy.size() ; i++ ) {
            std::cout << lowdim_data.reduction_strategy[i] << " "; }
            std::cout << std::endl;
    std::cout << " Flow: " << training_data.flow_type << std::endl;
    std::cout << " Parameters: ";
        for ( int i=0; i<training_data.parm_name.size() ; i++ ) {
            std::cout << training_data.parm_name[i] << " "; }
            std::cout << std::endl;    
    std::cout << " Snapshots number and format: " << training_data.n_snap << ", " << training_data.snap_fmt << std::endl;
    std::cout << " Error method: " << error_data.err_fmla << std::endl;
//        
    if ( !lowdim_data.flds_name.empty() ) {
        std::cout << " Fields name: ";
        for ( int i=0; i<lowdim_data.flds_name.size() ; i++ ) {
            std::cout << lowdim_data.flds_name[i] << " "; }
            std::cout << std::endl;    
    }
    std::cout << " Field type (IDs): " << lowdim_data.field_type << " ( " ;
        for ( int i=0; i<lowdim_data.fields.size() ; i++ ) {
            std::cout << lowdim_data.fields[i] << " "; }
            std::cout << ")" << std::endl;    
//
    for ( int i = 0; i < lowdim_data.reduction_strategy.size(); i++ ) { 
//
        if ( lowdim_data.reduction_strategy[i] == "SPOD" ) {
//                
            std::cout << " SPOD filter size: " << modal_data.SPOD_fs << std::endl;
            std::cout << " SPOD filter type: " << modal_data.SPOD_ft << std::endl;
            std::cout << " SPOD sigma: " << modal_data.SPOD_sigma << std::endl;        
            std::cout << " SPOD flag BC: " << modal_data.SPOD_flag_bc << std::endl;
//                
        } else if ( lowdim_data.reduction_strategy[i] == "DMD" || 
            lowdim_data.reduction_strategy[i] == "MRDMD") {
//                
            std::cout << " DMD rank: " << modal_data.DMD_rank << std::endl;
            std::cout << " DMD coefficient flag: " << modal_data.DMD_coef_flag << std::endl;
            std::cout << " MRDMD max levels: " << modal_data.MRDMD_max_levels << std::endl;        
            std::cout << " MRDMD max cycles: " << modal_data.MRDMD_max_cycles << std::endl;
//
        } else if ( lowdim_data.reduction_strategy[i] == "RDMD" ) {
            std::cout << " DMD rank: " << modal_data.DMD_rank << std::endl;
            std::cout << " DMD coefficient flag: " << modal_data.DMD_coef_flag << std::endl;
            std::cout << " RDMD rank: " << modal_data.RDMD_rank << std::endl;
        }
    }
//
//  Saving list of snapshots in csv format
    std::string snpsfilesum = "list-snapshots.csv";
    std::ofstream snps_data;
    snps_data.open(snpsfilesum.c_str());
//
    for ( int i = 0; i < training_data.n_parm; i++ ) {
        snps_data << training_data.parm_name[i] << ", ";
    }
    snps_data << " File name" << std::endl;
//
    for ( int i = 0; i < training_data.n_snap; i++ ) {
        for ( int j = 0; j < training_data.n_parm; j++) {
            snps_data << std::setprecision(4) << std::scientific << training_data.snap_pnts[i][j] << ", ";
        }
        snps_data << training_data.snap_list[i] << std::endl;
    }
    snps_data << std::endl;
    snps_data.close();
//
//  Saving list of error locations in csv format
    std::string errfilesum = "list-error-points.csv";
    std::ofstream err_data;
    err_data.open(errfilesum.c_str());
//
    for ( int i = 0; i < training_data.n_parm; i++ ) {
        err_data << training_data.parm_name[i] << ", ";
    }
    err_data << " File name" << std::endl;
//
    for ( int i = 0; i < error_data.n_err_pnts; i++ ) {
        for ( int j = 0; j < training_data.n_parm; j++) {
            err_data << std::setprecision(4) << std::scientific << error_data.err_pnts[i][j] << ", ";
        }
        err_data << error_data.err_sol_list[i] << std::endl;
    }
    err_data << std::endl;
    err_data.close();
}
//

*/

// ----------------------------------------------------------------------------------------
void print_soldata ( target_set_data &target_data, aero_data &aerodynamic_data ) {
// ----------------------------------------------------------------------------------------
//
    std::cout << " Target number and visualization format: " << target_data.n_targ << ", " << target_data.target_fmt << std::endl;      
    std::cout << " Modal coefficients formula: " << target_data.coeff_fmla << std::endl;
//
    if ( target_data.modes_sel_meth == "ENERGY" ) { 
        std::cout << " Modes selection method: " << target_data.modes_sel_meth << ", " << target_data.en << std::endl; 
    } else if ( target_data.modes_sel_meth == "USERDEF" ) { 
        std::cout << " Modes selection method: " << target_data.modes_sel_meth << ", "  << target_data.n_modes << std::endl; 
    }
//
//  Saving list of targets in csv format        
    std::string trgtfilesum = "list-targets.csv";
    std::ofstream trgt_data;
    trgt_data.open(trgtfilesum.c_str());
//
    for ( int i = 0; i < target_data.n_targ; i++ ) {
        for ( int j = 0; j < target_data.n_parm; j++) {
            trgt_data << std::setprecision(4) << std::scientific << target_data.target_pnts[i][j] << ", ";
        }
        trgt_data << target_data.target_list[i] << std::endl;
    }
    trgt_data << std::endl;
    trgt_data.close();
//
}
