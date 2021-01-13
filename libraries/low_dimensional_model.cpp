/* ------------------------------------------------------------------------------
\file generate_low_dimansional_model.cpp
* \brief Main to use a single RBM technique to reduce the complexity of an
*  aerodynamics problem
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
#include "low_dimensional_model.hpp"
//
// -------------------------------------------------------------------------------------------
int generate_low_dimensional_model ( const std::string filedata ) {
// -------------------------------------------------------------------------------------------    
//
//  Variables declaration
    ld_model_data lowdim_data;
    training_set_data training_data;
    ld_error_data error_data;
    modal_identification_data modal_data;
    manifold_learning_data manifold_data;
//
//  Reading configuration file
    std::string filefmt = "RAZR"; 
    read_gendata ( filedata, filefmt, lowdim_data, training_data, error_data, 
        modal_data, manifold_data );        
    print_gendata ( lowdim_data, training_data, error_data, 
        modal_data, manifold_data );        
//
//  Compute size of snapshots
    int snap_size = snapshot_size ( training_data.snap_list[1], training_data.snap_fmt );
//
//  Opening snapshots stream for reading fields
    std::ifstream snap_file_stream [training_data.n_snap];
    for ( int i = 0; i < training_data.n_snap; i++ ) {
//
        snap_file_stream[i].open( training_data.snap_list[i] );
//
        if ( !snap_file_stream[i].is_open() ) {
            std::cout << " " << std::endl;            
            std::cout << "-> Error: Snapshot file: " << training_data.snap_list[i] << " not found" << std::endl;            
            std::cout << std::endl;
            exit (EXIT_FAILURE);
        }
    }
//
//  Save the low dimensional model info into h5 database
    save_lowdim_model_info ( lowdim_data, training_data, error_data );
//
//  Matrix storing the snapshots
    Eigen::MatrixXd snap_matrix;
//
//  Loop over the different fields to be reduced
    for ( int i = 0; i < lowdim_data.n_flds; i++ ) { 
//
//---------------------------------------------------------------------------------------------------------------------
//  Read snapshots and create matrix of snapshots. This operation is done within a loop over the fields to be
//  processed. This is to manage memory requirements if many large snapshot files need to be processed. The exception
//  is in the case of field types such as VECTOR or CONSERVATIVE or PRIMITIVE etc. where the field are processed 
//  all togehter in the attempt to ensure consistency in capturing the time dynamics for unsteady problems. This 
//  is implemented in sucha  way that when a VECTOR or CONSERVATIVE etc. type is specified the variable setting.n_flds
//  is set to 1 field only so that the components are treated in a single snapshot matrix   
//---------------------------------------------------------------------------------------------------------------------
//
        snap_matrix = generate_snap_matrix( i, snap_size, snap_file_stream, lowdim_data, training_data );
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
                //lowmodel_MI->save_modal_representation( lowdim_data.lowdim_model_name );
//
            }
        }
    }
//
//  Closing all snapshots filestreams
    for ( int i = 0; i < training_data.n_snap; i++ )
        snap_file_stream[i].close();
    //
    return 0;
}
//
// ----------------------------------------------------------------------------------------
void save_lowdim_model_info ( ld_model_data &lowdim_data, training_set_data &training_data, 
    ld_error_data &error_data ) {
// ----------------------------------------------------------------------------------------
//
//  Creating the h5 database for the low dimensional model
    const H5std_string FILE_NAME(lowdim_data.lowdim_model_name);
    H5File file(FILE_NAME, H5F_ACC_TRUNC);
//
//  Create the baseline structure of the low dimensional database
    Group group01(file.createGroup("/Common"));
    Group group02(file.createGroup("/LDmodels"));
    Group group03(file.createGroup("/TrainingSet"));
//
//  Checking the existence of a mesh file and possibly read it
    CMesh mesh;
    if ( training_data.mesh_file != "NONE" ) {
//
        std::ifstream mesh_file_stream;
        mesh_file_stream.open( training_data.mesh_file );
//
        if ( !mesh_file_stream.is_open() ) {
                std::cout << " " << std::endl;
                std::cout << "-> Error: Mesh file: " << training_data.mesh_file << " not found" << std::endl;
                std::cout << "   If mesh file is not available set the keyword MESH_FILE to NONE" << std::endl;
                std::cout << std::endl;
                exit (EXIT_FAILURE); }
//
        mesh.read( mesh_file_stream, training_data.mesh_fmt );
        mesh_file_stream.close();
//
    }
//
    const int nstring = 9;
//
//  List of low dimensional model key info
    std::string header01 = "<REDUCTION STRATEGY> "+lowdim_data.inpstn_reduction_strategy;    
    std::string header02 = "<FIELD(S) NAME> "+lowdim_data.inpstn_fields_name;    
    std::string header03 = "<FIELD(S) ID> "+lowdim_data.inpstn_flag_fields;
    std::string header04 = "<REFERENCE FIELD> "+lowdim_data.inpstn_ref_field;
    std::string header05 = "<SNAPSHOTS FORMAT AND NUMBER> "+training_data.inpstn_snapshot_fmt+", "
      +std::to_string(training_data.n_snap);
    std::string header06 = "<PARAMETERS> "+training_data.inpstn_name_parm;
    std::string header07 = "<FLOW TYPE, TIME FIELD ID (for unsteady)> "+training_data.inpstn_flow_type;
    std::string header08 = "<MESH FILE> "+training_data.inpstn_mesh_file;
    std::string header09 = "<ERROR FORMULA> "+error_data.inpstn_error_fmla;
//
    const char *wdata[nstring]= { header01.c_str(), header02.c_str(), header03.c_str(), header04.c_str(), 
        header05.c_str(), header06.c_str(), header07.c_str(), header08.c_str(), header09.c_str() };
//
//  Low dimensional model specifications
    const H5std_string DATASET_COMMON_INFO("/Common/LDinfo");
//            
//  Dimsnsions of the dataset
    hsize_t dimspec[1] = {nstring};
    DataSpace specs_dims(1, dimspec);
    StrType stype(0, H5T_VARIABLE);
//
//  Writing the content on the dataset
    DataSet specs_data = file.createDataSet(DATASET_COMMON_INFO, stype, specs_dims);
    specs_data.write((void*)wdata, stype);
//
    specs_dims.close();
    specs_data.close();

    EigenHDF5::save(file, "TrainingSet/Coords", mesh.GetCoords());
    EigenHDF5::save(file, "TrainingSet/Connectivity", mesh.GetConnectivity());

    file.close();
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
