/* ------------------------------------------------------------------------------
\file modal_identification.hpp
* \brief Subroutines and functions to generate a low dimensional model.
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
#ifndef modal_identification_hpp
#define modal_identification_hpp
//
#include "read_data.hpp"
//
//------------------------------------------------------------------------------            
//  Classes for low-dimensional methods (Modal identification, Manifold learning)
//------------------------------------------------------------------------------            
//
//  /*!< \brief Main class for defining the modal identification basis a child class for each particular identification method (xPOD for POD and SPOD, xDMD for DMD, multi-resolutionDMD and recursive DMD)
class CModalIdentification {
//
protected:
	int ns;
	int nm;
	int nr;
//
	Eigen::VectorXd nrg_content;
    Eigen::VectorXcd eigenVal;
//
	Eigen::MatrixXd alpha;
	Eigen::MatrixXd xcmatrix;
	Eigen::MatrixXcd phi;
    Eigen::MatrixXcd eigenVec;
//
public:
//
//  Constuctor
	CModalIdentification() {
		ns = 0;
		nm = 0;
		nr = 0;
	}
//
    static inline Eigen::MatrixXcd phi_flds;
    static inline Eigen::MatrixXd alpha_flds;
//
	virtual void compute_correlation_matrix( Eigen::MatrixXd &snap_matrix, modal_identification_data &modal_data ) = 0;
	virtual void compute_modes( Eigen::MatrixXd &snap_matrix ) = 0;
	virtual void get_modal_coefficients() = 0;
	virtual void save_modal_representation(ld_model_data &lowdim_data, int i) = 0;
//
};
//
class CModalIdentification_POD : public CModalIdentification {
protected:
//
public:
//
//  Constuctor
	CModalIdentification_POD() : CModalIdentification() {
		ns = 0;
		nm = 0;
		nr = 0;		
	}
//
	void compute_correlation_matrix( Eigen::MatrixXd &snap_matrix, modal_identification_data &modal_data );
	void compute_modes( Eigen::MatrixXd &snap_matrix );
	void get_modal_coefficients();
    void save_modal_representation(ld_model_data &lowdim_data, int i);
};
//
/*
class CModalExtraction_SPOD : public CModalIdentification_POD {
protected:
//
public:
	void compute_correlation_matrix(Eigen::MatrixXd &snap_matrix, modal_identification_data &modal_data);
}
//
class CModalExtraction_DMD : public CModalIdentification {
protected:
//
public:
	Eigen::MatrixXcd compute_modes(Eigen::MatrixXd &snap_matrix, modal_identification_data &modal_data);
	Eigen::VectorXcd compute_modal_coefficients();
}
//
class CModalExtraction_rDMD : public CModalIdentification_DMD {
protected:
//
public:
	Eigen::MatrixXcd compute_modes(Eigen::MatrixXd &snap_matrix, modal_identification_data &modal_data);
	Eigen::VectorXcd compute_modal_coefficients();
}
//
class CModalExtraction_mrDMD : public CModalIdentification_DMD {
protected:
//
public:
	Eigen::MatrixXcd compute_modes(Eigen::MatrixXd &snap_matrix, modal_identification_data &modal_data);
	Eigen::VectorXcd compute_modal_coefficients();
}
*/
//
//------------------------------------------------------------------------------            
// Functions
//------------------------------------------------------------------------------            
//
//!< \brief Sorting eigevalues and eigenvectors
void sort_eigen( Eigen::VectorXcd &eval, Eigen::MatrixXcd &evec );
//
#endif
