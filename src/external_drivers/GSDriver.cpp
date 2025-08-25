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

#include "MRCPP/MWOperators"

#include "GSDriver.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmfunctions/orbital_utils.h"

namespace mrchem {

class FockBuilder;
//class PoissonOperator;

 /** @brief Calculates and stores the one- and two-electron integrals for given orbitals
 *
 * @param Phi: Vector of orbitals
 *
 * Calculates the one- and two-electron integrals for the orbitals in Phi, and stores them
 * in the class members one_body_integrals and two_body_integrals.
 * 
 */

void GSDriver::set_integrals(OrbitalVector &Phi, FockBuilder &F,  mrcpp::PoissonOperator &P){
    // operators
    KineticOperator K(F.momentum());
    NuclearOperator V = *(F.getNuclearOperator());
    // set the one- and two-body integrals
    GSDriver::set_one_body_integrals(Phi, K, V);
    GSDriver::set_two_body_integrals(Phi, P);
}


// Private

void GSDriver::set_one_body_integrals(OrbitalVector &Phi, KineticOperator &K, NuclearOperator &V){
    OrbitalVector KPhi = K(Phi);
    OrbitalVector VPhi = V(Phi);
    *(this->one_body_integrals) = orbital::calc_overlap_matrix(Phi, KPhi)
                                + orbital::calc_overlap_matrix(Phi, VPhi);
}

void GSDriver::set_two_body_integrals(OrbitalVector &Phi,  mrcpp::PoissonOperator &P){
    int n_orb = Phi.size();
    *(this->two_body_integrals) = ComplexTensorR4(n_orb, n_orb, n_orb, n_orb);
    this->two_body_integrals->setZero();

    for (int j = 0; j < n_orb; j++) {
        for (int l = 0; l < n_orb; l++) {
            // calculate rho_jl = Phi_j^+ Phi_l
            Orbital rho_jl = Phi[j].paramCopy();
            mrcpp::cplxfunc::multiply(rho_jl, Phi[j].dagger(), Phi[l], this->prec, true, true);
            // calculate V_jl = P(rho_jl)
            Orbital Vjl = rho_jl.paramCopy();
            if (rho_jl.hasReal()) {
                Vjl.alloc(NUMBER::Real);
                mrcpp::apply(this->prec, Vjl.real(), P, rho_jl.real());
            }
            if (rho_jl.hasImag()) {
                Vjl.alloc(NUMBER::Imag);
                mrcpp::apply(this->prec, Vjl.imag(), P, rho_jl.imag());
            }
            rho_jl.release();
            for (int i = 0; i < n_orb; i++) {
                // calculate <Phi_i|V_jl|
                Orbital tmp_i = Phi[i].paramCopy();
                mrcpp::cplxfunc::multiply(tmp_i, Phi[i].dagger(), Vjl, this->prec, true, true);
                for (int k = 0; k < n_orb; k++) {
                    // calculate (ij|kl) = <Phi_i|V_jl|Phi_k>
                    // BUG: add factor 4pi??
                    (*this->two_body_integrals)(i,j,k,l) = orbital::dot(tmp_i, Phi[k]);
                }
            }
        }
    }
}

} // namespace mrchem