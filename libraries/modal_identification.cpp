/* ------------------------------------------------------------------------------
\file modal_extraction.cpp
* \brief Subroutines and functions to generate a low dimensional model.
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
#include "modal_identification.hpp"
//
//------------------------------------------------------------------------------            
// Functions
//------------------------------------------------------------------------------            
//
// -------------------------------------------------------------------------------------------
void sort_eigen( Eigen::VectorXcd &eval, Eigen::MatrixXcd &evec ) {
// -------------------------------------------------------------------------------------------
//
    unsigned int swap_count = 1;
    std::complex<double> temp;
    Eigen::VectorXcd temp_vec(evec.rows());
//
    while ( swap_count > 0 ) {
//
        swap_count = 0;
        for( unsigned int index = 1; index < eval.size(); index++ ) {
//
            if ( eval.real()(index) > eval.real()(index-1) ) {
//
                temp = eval(index-1);
                eval(index-1) = eval(index);
                eval(index) = temp;
//
                temp_vec = evec.col(index-1);
                evec.col(index-1) = evec.col(index);
                evec.col(index) = temp_vec;
//
                swap_count++;
            }
        }
    }
}
//
// -------------------------------------------------------------------------------------------
void CModalIdentification_POD::compute_correlation_matrix( Eigen::MatrixXd &snap_matrix, 
	modal_identification_data &modal_data ) {
// -------------------------------------------------------------------------------------------
//
	ns = snap_matrix.cols();
	nr = snap_matrix.rows();
	xcmatrix.setZero(ns, ns);
	xcmatrix = snap_matrix.transpose()*snap_matrix;
//
}
//
// -------------------------------------------------------------------------------------------
void CModalIdentification_POD::compute_modes( Eigen::MatrixXd &snap_matrix ) {
// -------------------------------------------------------------------------------------------
//
    std::cout << " Computing modes " << std::endl;

//  Solving eignevalues-eigenvector problem. Sirovich method of snapshots
	Eigen::EigenSolver<Eigen::MatrixXd>  eigsol(xcmatrix);
//
	eigenVal = eigsol.eigenvalues().real();
    eigenVec = eigsol.eigenvectors().real();
    sort_eigen( eigenVal, eigenVec );
//
//  Computing energy content associated to each mode
    nrg_content.setZero(ns);
    double sum = 0;
    for ( int i = 0; i < ns; i++ ) {
        sum += eigenVal(i).real() / eigenVal.real().sum();
        nrg_content(i) = sum;
    }
//
//  Resizing the matrix of modes to exclude those modes with an energy content lower than tol
    double tol = eigenVal(0).real()*1e-12;
    nm = 0;
    while ( nm < eigenVal.size() && eigenVal(nm).real() > tol )
            nm++;
//
    Eigen::MatrixXcd phi_tmp;
    phi_tmp = snap_matrix*eigenVec;
    phi.setZero(snap_matrix.rows(),nm);
//
    for ( int i = 0; i < nm; i++ )
       phi.col(i) = phi_tmp.col(i) / sqrt(eigenVal(i).real());
//
}
//
// -------------------------------------------------------------------------------------------
void CModalIdentification_POD::get_modal_coefficients() {
// -------------------------------------------------------------------------------------------
//
    std::cout << " Computing modal coefficients " << std::endl;    
//    
	alpha.setZero(ns,nm);
    for ( int i = 0; i < nm; i++ )
        alpha.col(i) = eigenVec.col(i).real() * sqrt(eigenVal(i).real());
//
}
//

