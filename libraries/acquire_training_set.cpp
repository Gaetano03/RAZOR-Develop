/* ------------------------------------------------------------------------------
\file acquire_training_set.cpp
* \brief Subroutines and functions to acquire and generate snapshots set.
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
#include "acquire_training_set.hpp"
//
// ----------------------------------------------------------------------------------------
Eigen::MatrixXd generate_snap_matrix( const int fld_id, const int Ni, std::ifstream 
    filestream[], flds_data fields_data, training_set_data training_data ) {
// ----------------------------------------------------------------------------------------
//
    if ( fields_data.field_type == "SCALAR" ) {
        std::cout << std::endl;
        std::cout << " Building snapshot matrix for SCALAR field ID " << fields_data.fields[fld_id] << std::endl;
        std::cout << " - estimated memory for the snapshots matrix is " 
        << Ni*training_data.n_snap*8.0/1000000.0 << "Mb" << std::endl;
    } else {
        std::cout << std::endl;
        std::cout << " Building snapshot matrix for VECTOR field with components ID (";
        for ( int i=0; i<fields_data.fields.size() ; i++ ) {
            std::cout << fields_data.fields[i]; }
            std::cout << " )" << std::endl;
        std::cout << " - estimated memory for the snapshots matrix is " 
        << fields_data.fields.size()*Ni*training_data.n_snap*8.0/1000000.0 << "Mb" << std::endl;
    }
//
    int Nj = training_data.n_snap; 
    std::string fld_type = fields_data.field_type;
    std::vector<int> fields = fields_data.fields;
    std::string filefmt = training_data.snap_fmt;
    std::string reffield_type = fields_data.ref_field;
//
    Eigen::MatrixXd snap_matrix;
    Eigen::MatrixXd snap_field;
//
    if ( fld_type == "SCALAR" || fld_type == "QCRITERION" ) {
        snap_matrix = Eigen::MatrixXd::Zero(Ni, Nj);
    } else {
        snap_matrix = Eigen::MatrixXd::Zero(Ni*fields.size(), Nj);
        snap_field = Eigen::MatrixXd::Zero(Ni,fields.size());        
    }
//
//------------------------------------------------------------------------------------------------------------------
//  Loop over the number of snapshots. Read the field from the i-th snapshots
//  field is preprocessed if necessary (e.g. PRIMITIVE, VELOCITY, QCRITERION, etc.)
//  after being read. Then fill the snapshots matrix
//------------------------------------------------------------------------------------------------------------------
//
    for ( int i = 0; i < Nj; i++ ) {
//        
        switch ( read_field_type(fld_type) ) {
//
            case SCALAR: {
//
//              Identify field to be read                
                std::vector<int> subfield;
                subfield.push_back(fields[fld_id]);
//                
//              Read scalar field from i-th snapshot into snap_matrix directly (size_of_snapshot, number_of_snapshots)
                snap_matrix.col(i) = read_field( i, Ni, subfield, filestream, filefmt );
                break; }
//
            case VECTOR: {
//
//              Identify field components to be read                
                std::vector<int> subfield;
                for ( int j = 0; j < fields.size(); j++ ) {
                    subfield.push_back(fields[j]);
                }
//                
//              Read vector field from i-th snapshot into matrix snap_field (size_of_snapshot, number_of_field_components)            
                snap_field = read_field( i, Ni, subfield, filestream, filefmt );
//
//              Put all fields components from i-th snapshot into a single column in the snapshot matrix
                int counter = 0;
                for ( int j = 0; j < fields.size(); j++ ) {
                    snap_matrix.block(counter,i,Ni,1) = snap_field.col(j); 
                    counter += Ni;
                }                
                break; }
//
            case CONSERVATIVE: { std::cout << "-> Error: CONSERVATIVE variables not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
            case PRIMITIVE: { std::cout << "-> Error: PRIMITIVE variables not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
            case VELOCITY: { std::cout << "-> Error: VEOCITY field not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
            case AEROLOADS: { std::cout << "-> Error: AEROLOADS variables not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
            case QCRITERION: { std::cout << "-> Error: QCRITERION not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
//            
            default: { std::cout << "-> Error: Unknown field type" << std::endl; exit (EXIT_FAILURE); break; }
//
        }
    }
//
//------------------------------------------------------------------------------------------------------------------    
//  Preprocessing of snapshot matrix. Subtract from the existing snapshot fields a reference field. If nothign is 
//  specified, no preprocessing will be applied
//------------------------------------------------------------------------------------------------------------------
//
    switch ( read_refsol_type(reffield_type) ) {
//
        case NONE: { break; }
//        
        case MEAN_FIELD: {
//
            Eigen::VectorXd mean = snap_matrix.rowwise().mean();
            for ( int i = 0; i < Nj; i++ ) {
                snap_matrix.col(i) -= mean; }
            break; }
//
        case INITIAL_FIELD: { std::cout << "-> Error: INITIAL_FIELD reference solution not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
        case CONSTANT_FIELD: { std::cout << "-> Error: CONSTANT_FIELD reference solution not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
        case READ_FROM_FILE: { std::cout << "-> Error: READ_FROM_FILE reference solution not supported yet" << std::endl; exit (EXIT_FAILURE); break; }
//
    }
    return snap_matrix;
}




