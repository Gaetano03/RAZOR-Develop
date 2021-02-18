/* ------------------------------------------------------------------------------
\file low_dimensional_model.cpp
* \brief Main to use a single RBM technique to reduce the complexity of an
*  aerodynamics problem
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
#include "low_dimensional_model.hpp"
//
// -------------------------------------------------------------------------------------------
int generate_low_dimensional_model ( const std::string filedata ) {
// -------------------------------------------------------------------------------------------    
//
//  Variable declaration.
    ld_model_data lowdim_data;
    ld_error_data error_data;
    modal_identification_data modal_data;
    manifold_learning_data manifold_data;
//
//  Extended (with all fields) and regular snapshots matrix
    Eigen::MatrixXd snap_flds;
    Eigen::MatrixXd snap_matrix;  
//
    H5File file;
//
//  Read config file
    read_gendata ( filedata, lowdim_data, error_data, modal_data, manifold_data );
//
    const std::string filename = { lowdim_data.database_name };
//
//  Open and load H5 database.
    file.openFile ( filename, H5F_ACC_RDWR, FileAccPropList::DEFAULT );
    EigenHDF5::load(file, "TrainingSet/Snapshots", snap_flds);
//
//  Save LD model attributes, read number of fields
    Group group_root = file.openGroup("/");
    Group group_ld = file.openGroup("/LDmodels");
//
    create_h5_attr (group_ld, "REDUCTION STRATEGY", lowdim_data.inpstn_reduction_strategy );
    create_h5_attr (group_ld, "NUMBER OF METHODS", std::to_string(lowdim_data.n_meth) );
//
    std::string n_flds_str = read_h5_attr ( group_root, "NUMBER OF FIELDS" );
    lowdim_data.n_flds = std::stoi(n_flds_str);
//
    group_ld.close();
    group_root.close();
//
//  Init snap_matrix
    snap_matrix.setZero(snap_flds.rows()/lowdim_data.n_flds, snap_flds.cols());
//
//  Loop over the different fields to be reduced
    for ( int i = 0; i < lowdim_data.n_flds; i++ ) {
//
//      Select sub-matrix. NOTE: it'd be good to check the values
        snap_matrix = snap_flds.block(i*snap_matrix.rows(), 0, snap_matrix.rows(), snap_matrix.cols());
//
//-----------------------------------------------------------------------------------------------------------------------
//  Perform order reduction by identifying a reduced space. Two high-level approaches are implemented: the first one is
//  MODAL IDENTIFICATION and allows for a variety of methods for basis (aka modes, aka primitives) identification
//  the second one is based on MANIFOLD LEARNING and aims at identifying the manifold onto which the solutions lie
//-----------------------------------------------------------------------------------------------------------------------
//
        for ( int m = 0; m < lowdim_data.n_meth; m++ ) {
//
            CModalIdentification *lowmodel_MI;
            //CManifoldLearning *lowmodel_ML;
//            
            std::string lowdim_strategy;
            switch ( read_ldmethod( lowdim_data.reduction_strategy[m] ) ) {
//    
                case POD: { lowdim_strategy = "MODAL_IDENTIFICATION";
                    lowmodel_MI = new CModalIdentification_POD ; break; }                            
                case SPOD: { lowdim_strategy = "MODAL_IDENTIFICATION";
                    /*CModalIdentification_SPOD *lowmodel;*/ break; }             
                case DMD: { lowdim_strategy = "MODAL_IDENTIFICATION";
                    /*CModalIdentification_DMD *lowmodel;*/ break; } 
                case RDMD: { lowdim_strategy = "MODAL_IDENTIFICATION"; 
                    /*CModalIdentification_rDMD *lowmodel;*/ break; } 
                case MRDMD: { lowdim_strategy = "MODAL_IDENTIFICATION";  
                    /*CModalIdentification_mrDMD *lowmodel;*/ break; } 
                case ISOMAP: { lowdim_strategy = "MANIFOLD_LEARNING";
                    /*CManifoldLearning_ISOMAP *lowmodel;*/ break; } 
            }
//
            if ( lowdim_strategy == "MODAL_IDENTIFICATION" ) {
//
                lowmodel_MI->compute_correlation_matrix( snap_matrix, modal_data );           
                lowmodel_MI->compute_modes( snap_matrix );            
                lowmodel_MI->get_modal_coefficients();
                lowmodel_MI->save_modal_representation( lowdim_data, i);
//
            }
        }
    }
//
    file.close();
    return 0;
}
//
// -------------------------------------------------------------------------------------------
lowdim_strat read_ldmethod ( const std::string &ldstrat_string ) {
// -------------------------------------------------------------------------------------------
//
    if ( ldstrat_string == "POD" )
        return POD;
    else if ( ldstrat_string == "SPOD" )
        return SPOD;
    else if ( ldstrat_string == "DMD" )
        return DMD;
    else if ( ldstrat_string == "RDMD" )
        return RDMD;
    else if ( ldstrat_string == "MRDMD" )
        return MRDMD;
    else if ( ldstrat_string == "ISOMAP" )
        return ISOMAP;
//
    else {
//
        std::cout << " " << std::endl;
        std::cout << "-> Error: Unknown REDUCTION STRATEGY" << std::endl;
        std::cout << "   " << ldstrat_string << std::endl;
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//     
    }
}
//
