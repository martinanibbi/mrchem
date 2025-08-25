/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include "GSDriver.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "tensor/RankZeroOperator.h"
#include "qmfunctions/orbital_utils.h"

namespace mrchem {

class FockBuiulder;
class MomentumOperator;
class RankZeroOperator;

 /** @brief Calculates and stores the one- and two-electron integrals for given orbitals
 *
 * @param Phi: Vector of orbitals
 *
 * Calculates the one- and two-electron integrals for the orbitals in Phi, and stores them
 * in the class members one_body_integrals and two_body_integrals.
 * 
 */

void GSDriver::set_integrals(OrbitalVector &Phi, FockBuilder &F){
    // operators
    KineticOperator K(F.momentum());
    NuclearOperator V = *(F.getNuclearOperator());
    // set the one- and two-body integrals
    GSDriver::set_one_body_integrals(Phi, K, V);
    //GSDriver::set_two_body_integrals(Phi, F);
}


// Private

void GSDriver::set_one_body_integrals(OrbitalVector &Phi, KineticOperator &K, NuclearOperator &V){
    OrbitalVector KPhi = K(Phi);
    OrbitalVector VPhi = V(Phi);
    *(this->one_body_integrals) = orbital::calc_overlap_matrix(Phi, KPhi)
                                + orbital::calc_overlap_matrix(Phi, VPhi);
}

void GSDriver::set_two_body_integrals(OrbitalVector &Phi, FockBuilder &F){
    int n_orb = Phi.size();
    this->two_body_integrals.resize(n_orb, std::vector<std::vector<std::vector<double>>>(n_orb, std::vector<std::vector<double>>(n_orb, std::vector<double>(n_orb, 0.0))));
    //CoulombOperator *J = F.getCoulombOperator();
    //ExchangeOperator *K = F.getExchangeOperator();
    for (int p = 0; p < n_orb; p++){
        for (int q = 0; q < n_orb; q++){
            for (int r = 0; r < n_orb; r++){
                for (int s = 0; s < n_orb; s++){
                    this->two_body_integrals[p][q][r][s] = 1; //J->calc(p, q, r, s) - F.exact_exchange * K->calc(p, q, r, s);
                }
            }
        }
    }
}

} // namespace mrchem