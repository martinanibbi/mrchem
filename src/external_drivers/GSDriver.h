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

#pragma once

#include <vector>

#include "mrchem.h"
#include "qmfunctions/Orbital.h"
#include "qmoperators/two_electron/FockBuilder.h"

namespace mrchem {

class FockBuilder;

class GSDriver{
public:
    GSDriver() = default;
    virtual ~GSDriver() = default;

    void set_integrals(OrbitalVector &Phi, FockBuilder &F);
    virtual void optimize() = 0;
    virtual void get_rdms() = 0;

protected:
    std::shared_ptr<ComplexMatrix> one_body_integrals{};
    std::vector<std::vector<std::vector<std::vector<double>>>> two_body_integrals;

private:
    void set_one_body_integrals(OrbitalVector &Phi, KineticOperator &K, NuclearOperator &V);
    void set_two_body_integrals(OrbitalVector &Phi, FockBuilder &F);
};

} // namespace mrchem