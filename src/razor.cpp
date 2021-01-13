/* ------------------------------------------------------------------------------
\file razor.cpp
* \brief Driver routine for RAZOR: Computational reduction and Flow identification
*
* \version 2020.10 "Wolf"
* Copyright 2016-2020, Aerospace Centre of Excellence University of Strathclyde
*
* RAZOR is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* MODES is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with MODES. If not, see <http://www.gnu.org/licenses/>.
* ------------------------------------------------------------------------------*/
//
#include "low_dimensional_model.hpp"
#include "low_dimensional_solution.hpp"
//
int main ( int argc, char *argv[] ) {
//
    std::cout << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;    
    std::cout << "  RAZOR v2020.10 - Wolf" << std::endl;    
    std::cout << "  Copyright 2016-2020, Aerospace Centre, University of Strathclyde" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
//
    if ( argc != 3 ) {    
        std::cout << std::endl;    
        std::cout << " -> Usage: razor -[gs] configuration_file" << std::endl; 
        std::cout << "    " << std::endl;
        std::cout << "    Options:" << std::endl;
        std::cout << "    -g     low dimensional model generation" << std::endl;
        std::cout << "    -s     computing low dimensional solution" << std::endl;
        std::cout << std::endl;
        exit (EXIT_FAILURE); 
    }
//    
    std::string execmode = argv[1];
    std::string filedata = argv[2];
//
//------------------------------------------------------------------------------------
//  Selection of execution mode from command line input:
//  GENERATION OF LOW DIMENSIONAL MODEL = acquire snapshots and generate low-
//                                        dimensional representation
//  COMPUTATION OF LOW DIMENSIONAL SOLUTION = reconstruct low-dimensional solutions  
//                                            at specified untried conditions
//---------------------------------------------------------------------------------
//
    if ( execmode == "-g" ) {        
        int status = generate_low_dimensional_model(filedata);
    } else if ( execmode == "-s" ) {
        int status = compute_low_dimensional_solution(filedata);
//        
    } else {
        std::cout << std::endl;
        std::cout << "-> Error: Unknown execution mode" << std::endl;        
        std::cout << std::endl;
        exit (EXIT_FAILURE);        
    }
//
    std::cout << std::endl;
    std::cout << " End of execution " << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
//
    return 0;
}
