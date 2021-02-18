/* ------------------------------------------------------------------------------
\file low_dimensional_solution.cpp
* \brief Subroutines and functions to generate a low dimensional solution.
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
#include "low_dimensional_solution.hpp"
//
// -------------------------------------------------------------------------------------------
void CModalReconstruction::features_scaling() {

    std::cout << "Features Scaling... " << std::endl;
    m_mean_params.resize(m_training_data.n_parm);
    m_std_params.resize(m_training_data.n_parm);
//
    std::fill(m_mean_params.begin(), m_mean_params.end(), 0.);
    std::fill(m_std_params.begin(), m_std_params.end(), 0);
//
    m_train_pnts_scaled = std::vector<std::vector<double> > (m_training_data.n_snap, std::vector<double> (m_training_data.n_parm, 0.));
//
    for ( int i = 0; i < m_training_data.n_parm; i++ ) {
        for ( int j = 0; j < m_training_data.n_snap; j++ )
            m_mean_params[i] += m_training_data.snap_pnts[j][i];
        m_mean_params[i] /= m_training_data.n_snap;

        for ( int j = 0; j < m_training_data.n_snap; j++ )
            m_std_params[i] += (m_training_data.snap_pnts[j][i] - m_mean_params[i])*(m_training_data.snap_pnts[j][i] - m_mean_params[i]);
        m_std_params[i] = std::sqrt(m_std_params[i]/(double)(m_training_data.n_snap-1));

        if (m_std_params[i] < 1e-15) std::cout << "WARNING - constant parameter, feature scaling inappropriate ... " << std::endl;

    }

    for ( int i = 0; i < m_training_data.n_snap; i++ ){
        for ( int j = 0; j < m_training_data.n_snap; j++ ){
            m_train_pnts_scaled[i][j] = (m_training_data.snap_pnts[i][j] - m_mean_params[j])/m_std_params[j];
        }
    }

    // Scaling also target points accordingly
    m_targ_pnts_scaled = std::vector<std::vector<double> > (m_target_data.target_pnts.size(), std::vector<double> (m_target_data.target_pnts[0].size(), 0.));
    for ( int i = 0; i < m_target_data.target_pnts.size(); i++ ){
        for ( int j = 0; j < m_target_data.n_parm; j++ )
            m_targ_pnts_scaled[i][j] = (m_target_data.target_pnts[i][j] - m_mean_params[j])/m_std_params[j];
    }

//    std::cout << "Vector of means: " << std::endl;
//    for ( auto m : m_mean_params) std::cout << m << "\t";
//    std::cout << std::endl;
//
//    std::cout << "Vector of std: " << std::endl;
//    for ( auto m : m_std_params) std::cout << m << "\t";
//    std::cout << std::endl;
}
//
// -------------------------------------------------------------------------------------------
void CModalReconstruction::save_solutions(){
// -------------------------------------------------------------------------------------------
    //Defining Number of grid points
    int Np = m_Rec_field.rows()/m_flds_data.n_flds;
    for ( int targ_id = 0; targ_id < m_target_data.n_targ; targ_id++ ){
        std::ofstream flow_data;
        flow_data.open(m_target_data.target_list[targ_id]);

        //Writing row of headers
        for ( int field_id = 0; field_id < m_flds_data.n_flds; field_id++ )
            flow_data << "\"" << m_flds_data.flds_name[field_id] << "\",";
        flow_data << std::endl;

        //Writing fields
        for ( int iPoint = 0; iPoint < Np; iPoint++ ){
            for ( int field_id = 0; field_id < m_flds_data.n_flds; field_id++ )
                flow_data << std::setprecision(12) << std::scientific << m_Rec_field(field_id*Np + iPoint,targ_id) <<  ",";
        flow_data << std::endl;
        }
        // Close file
        flow_data.close();
    }
}
//
// -------------------------------------------------------------------------------------------
void CInterpolation::load_lowdim(const std::string filename) {
// -------------------------------------------------------------------------------------------
    std::cout << "Loading low dimensional solution ..." << std::endl;
    H5File file;
//
//  Open and load H5 database.
    file.openFile ( filename, H5F_ACC_RDWR, FileAccPropList::DEFAULT );
    EigenHDF5::load(file, "/LDmodels/PODbasis", m_Phi_flds);
    EigenHDF5::load(file, "/LDmodels/PODcoeff", m_Coefs_flds);
//
    file.close();

    m_Rec_field.setConstant(m_Phi_flds.rows(), m_target_data.n_targ, 0.);
}
//
// -------------------------------------------------------------------------------------------
void CInterpolation::compute_lowdim_surrogate(const int field_id) {
// -------------------------------------------------------------------------------------------
    std::cout << "Building low dimensional surrogate... " << std::endl;
    double avgDt = 0.0;  //!< \brief Need to define this (maybe with k-mean algorithm also used in ML)
    // Getting RBF interpolation method
    std::string delimiter = "_";
    size_t pos = 0;
    std::string rbf_met = m_target_data.coeff_fmla;
    pos = rbf_met.find(delimiter);
    rbf_met.erase(0, pos + delimiter.length());

    features_scaling();
//
    RBF_CONSTANTS rbf_const {avgDt, 0.0};
    m_Coefs = m_Coefs_flds.block(field_id*m_training_data.n_snap, 0, m_training_data.n_snap, m_target_data.n_modes);

    for ( int i = 0; i < m_target_data.n_modes; i++ ){

        std::vector<double> coefs ;
        for (int j = 0 ; j < m_training_data.n_snap ; j++)
            coefs.push_back(m_Coefs(j,i));

        surrogates[field_id].push_back( rbf(m_train_pnts_scaled, coefs, get_key_rbf(rbf_met), rbf_const) );
        surrogates[field_id][i].build();
    }
}
//
// -------------------------------------------------------------------------------------------
void CInterpolation::compute_lowdim_sol(const int field_id) {
// -------------------------------------------------------------------------------------------
    std::cout << "Computing low dimensional solution... " << std::endl;
    //Scaling features accordingly
    m_coefs_interp.setZero(m_target_data.n_modes,m_target_data.target_pnts.size());

    for ( int i = 0; i < m_target_data.n_modes; i++ ) {
        for ( int j = 0; j < m_target_data.target_pnts.size(); j++ )
            surrogates[field_id][i].evaluate(m_targ_pnts_scaled[j], m_coefs_interp(i, j));
    }

}
//
// -------------------------------------------------------------------------------------------
void CInterpolation::compute_highdim_sol(const int field_id) {
// -------------------------------------------------------------------------------------------
    std::cout << "Computing high dimensional solution... " << std::endl;
    m_Phi = m_Phi_flds.block(field_id*(m_Phi_flds.rows()/m_flds_data.n_flds), 0, m_Phi_flds.rows()/m_flds_data.n_flds, m_target_data.n_modes);
    m_Rec_field.block(field_id*m_Phi.rows(),0,m_Phi.rows(),m_target_data.n_targ) = m_Phi.real()*m_coefs_interp;
}
//
// -------------------------------------------------------------------------------------------
int compute_low_dimensional_solution ( const std::string filedata ) {
// -------------------------------------------------------------------------------------------
//
//  Variables declaration
    target_set_data target_data;
    training_set_data training_data;
    aero_data aerodynamic_data;
    flds_data fields_data;
    //
//  Reading configuration file
    read_soldata ( filedata, target_data, aerodynamic_data);        
    read_procdata( filedata, fields_data, training_data);
    //print_soldata ( lowdim_data, training_data, error_data,
    //    modal_data, manifold_data );

    CModalReconstruction *ModalRec;
//

    const std::string filename = { target_data.database_name };

    ModalRec = new CInterpolation(target_data, training_data,fields_data);
    ModalRec->load_lowdim(filename);

    for ( int field_id = 0; field_id < fields_data.n_flds; field_id++ ) {
        std::cout << "For Field ID " << field_id << std::endl;

        ModalRec->compute_lowdim_surrogate(field_id);
        ModalRec->compute_lowdim_sol(field_id);
        ModalRec->compute_highdim_sol(field_id);

        std::cout << std::endl;
    }
    ModalRec->save_solutions();

    return 0;
}
//
//
// -------------------------------------------------------------------------------------------
smartuq::surrogate::RBF_FUNCTION get_key_rbf ( const std::string &key_string ) {
// -------------------------------------------------------------------------------------------
//
    if ( key_string == "LINEAR" )
        return smartuq::surrogate::LINEAR;
    else if ( key_string == "CUBIC" )
        return smartuq::surrogate::CUBIC;
    else if ( key_string == "GAUSSIAN" )
        return smartuq::surrogate::GAUSSIAN;
    else if ( key_string == "THIN_PLATE" )
        return smartuq::surrogate::THIN_PLATE;
    else if ( key_string == "MULTIQUADRATICS" )
        return smartuq::surrogate::MULTIQUADRATICS;
//
    else {
//
        std::cout << " " << std::endl;
        std::cout << "-> Error: Unknown REDUCTION STRATEGY type" << std::endl;
        std::cout << " " << std::endl;
        exit (EXIT_FAILURE);
//     
    }
//
}
