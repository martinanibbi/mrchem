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

#include "GenericTwoOrbitalsPotential.h"
#include "MRCPP/MWOperators"
#include "qmfunctions/orbital_utils.h"


namespace mrchem {

GenericTwoOrbitalsPotential::GenericTwoOrbitalsPotential(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<mrchem::OrbitalVector> Phi, bool mpi_share)
        : QMPotential(1, mpi_share)
        , orbitals(Phi)
        , poisson(P) {}


Orbital GenericTwoOrbitalsPotential::g_jl(int j, int l, double prec){
    OrbitalVector Phi = *(this->orbitals);
    // calculate rho_jl = Phi_j^+ Phi_l
    Orbital rho_jl = Phi[j].paramCopy();
    mrcpp::cplxfunc::multiply(rho_jl, Phi[j].dagger(), Phi[l], prec, true, true);
    // calculate g_jl = P(rho_jl)
    Orbital g_jl = rho_jl.paramCopy();
    if (rho_jl.hasReal()) {
        g_jl.alloc(NUMBER::Real);
        mrcpp::apply(prec, g_jl.real(), *(this->poisson), rho_jl.real());
    }
    if (rho_jl.hasImag()) {
        g_jl.alloc(NUMBER::Imag);
        mrcpp::apply(prec, g_jl.imag(), *(this->poisson), rho_jl.imag());
    }
    rho_jl.release();
    return g_jl;
}

}