/* ------------------------------------------------------------------------------
\file preprocessing.cpp
* \brief Preprocessing snapshot and mesh data, saving to HDF5 file
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
#include "preprocessing.hpp"
//
// ------------------------------------------------------------------------------------------- 
int preprocess_data ( const std::string filedata ) {
// ------------------------------------------------------------------------------------------- 
//
//  Variables declaration
//
//  Structures to store input config data
    flds_data fields_data;
    training_set_data training_data;
//
//  Matrix storing the snapshots and vector storing ID of required (sub)set of nodes
    Eigen::MatrixXd snap_matrix_full;
    Eigen::MatrixXd snap_matrix;
    Eigen::VectorXi snap_PointID;
//
//  Reading configuration file
    read_procdata ( filedata, fields_data, training_data );       
    //print_gendata ( lowdim_data, training_data, error_data, modal_data, manifold_data );
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
//  Create H5 file, save data as root/group attributes
    save_attributes ( fields_data , training_data );
//
//  Read and save mesh to H5. Also returns vector with PointID to select subsets within a snapshot.                            
    snap_PointID = save_mesh ( fields_data, training_data );                                                                         
//
//  Generating matrix of snapshots and saving to H5 database while looping over* [long comment at bottom of code]
    for ( int i = 0; i < fields_data.n_flds; i++ ) { 
        snap_matrix_full = generate_snap_matrix ( i, snap_size, snap_file_stream, fields_data, training_data );
        snap_matrix.setZero(snap_PointID.rows(), snap_matrix_full.cols());
//
        for ( int j = 0; j < snap_PointID.rows(); j++ )
            snap_matrix.row(j) = snap_matrix_full.row(snap_PointID(j, 1));
//
        save_snap_matrix ( snap_matrix, fields_data, i);
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
void save_attributes ( flds_data &fields_data, training_set_data &training_data ) {       
// ----------------------------------------------------------------------------------------
//
//  Creating the h5 database for the low dimensional model
    const std::string filename = { fields_data.database_name };
    H5File file(filename, H5F_ACC_TRUNC);
//
//  Create the baseline structure of the low dimensional database
    Group group_root = file.openGroup("/");
    Group group01(file.createGroup("/LDmodels"));
    Group group02(file.createGroup("/TrainingSet"));
//
//  Write attributes to root group
    create_h5_attr ( group_root, "FIELD(S) NAME", fields_data.inpstn_fields_name );
    create_h5_attr ( group_root, "FIELD(S) ID", fields_data.inpstn_flag_fields );
    create_h5_attr ( group_root, "NUMBER OF FIELDS", std::to_string(fields_data.n_flds) );
    create_h5_attr ( group_root, "REFERENCE FIELD", fields_data.inpstn_ref_field );
    create_h5_attr ( group02, "SNAPSHOTS FORMAT AND NUMBER", training_data.inpstn_snapshot_fmt+", "+std::to_string(training_data.n_snap) );
    create_h5_attr ( group02, "PARAMETERS", training_data.inpstn_name_parm );
    create_h5_attr ( group02, "FLOW TYPE, TIME FIELD ID (for unsteady)", training_data.inpstn_flow_type );
    create_h5_attr ( group02, "MESH FILE", training_data.inpstn_mesh_file );
//
    group02.close();
    group01.close();
    group_root.close();
    file.close();
//
}
//
// -------------------------------------------------------------------------------------------
void save_snap_matrix ( Eigen::MatrixXd &snap_matrix, flds_data &fields_data, int &field_id ) {
// -------------------------------------------------------------------------------------------
//
    std::cout << " Saving snapshots " << std::endl;
//  
//  Declare variables. Extended matrix --> fields are stored sequentially wrt rows
    H5File file;
    const std::string filename = { fields_data.database_name };
//
    int &fld_n = fields_data.n_flds;
    static Eigen::MatrixXd snap_flds;
//
//  Open file
    file.openFile ( filename, H5F_ACC_RDWR, FileAccPropList::DEFAULT );   
//
    if ( field_id == 0 ) {
        std::cout << "  Initialising extended matrix for snapshots " << std::endl;
        snap_flds.setZero( snap_matrix.rows() * fld_n, snap_matrix.cols() );
    }
//  Writing to temporary Eigen::MatrixX containing fields sequentially  
    int st = field_id * snap_matrix.rows();
    int end = ( field_id+1 ) * snap_matrix.rows();
    int i = 0;
    for ( st; st < end; st++ ) {
        snap_flds.row(st) = snap_matrix.row(i);
        i++;
    }
//
//  Saving to HDF5 database after all fields are stored
    if ( field_id == (fld_n-1) ) {
        EigenHDF5::save(file, "/TrainingSet/Snapshots", snap_flds);
        std::cout << "\n Snapshots written to database " << std::endl;
    }
//
    file.close();
//
}
//
// -------------------------------------------------------------------------------------------
Eigen::VectorXi save_mesh ( flds_data &fields_data, training_set_data &training_data ) { 
// -------------------------------------------------------------------------------------------
//
//  Checking the existence of a mesh file and read it
    CMesh mesh(fields_data);
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
//  Open H5 file, save mesh, close
    H5File file;
    const std::string filename = { fields_data.database_name };
//
    file.openFile ( filename, H5F_ACC_RDWR, FileAccPropList::DEFAULT ); 
//
    EigenHDF5::save(file, "TrainingSet/Coords", mesh.GetCoords());
    EigenHDF5::save(file, "TrainingSet/Connectivity", mesh.GetConnectivity());
    EigenHDF5::save(file, "TrainingSet/PointID", mesh.GetPointID());
//
    std::cout << " Mesh saved successfully" << std::endl;
//
    file.close();
    return mesh.GetPointID();
//
}

// Long comment *
//---------------------------------------------------------------------------------------------------------------------
//  Read snapshots and create matrix of snapshots. This operation is done within a loop over the fields to be
//  processed. This is to manage memory requirements if many large snapshot files need to be processed. The exception
//  is in the case of field types such as VECTOR or CONSERVATIVE or PRIMITIVE etc. where the field are processed 
//  all togehter in the attempt to ensure consistency in capturing the time dynamics for unsteady problems. This 
//  is implemented in sucha  way that when a VECTOR or CONSERVATIVE etc. type is specified the variable setting.n_flds
//  is set to 1 field only so that the components are treated in a single snapshot matrix   
//---------------------------------------------------------------------------------------------------------------------
