/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Regulator

Description
    Regulator is a class which allows user to read value of a specified
    scalar field in some location inside mesh and returns a value between
    -1.0 and 1.0. Returned value dependent on the regulator settings
    obtained from "regulatorProperties" inside the "constant" directory.

    regulatorProperties dict:
    \table
        Property        | Description                 | Required |
        patchName       | name of measurement patch   | yes      |
        fieldName       | measured field name         | yes       |
        targetValue     | desired value               | yes      |
        Kp              | proportional gain           | yes      |
        Ti              | integral gain               | yes      |
        Td              | differential gain           | yes      |
    \endtable

\*---------------------------------------------------------------------------*/

#ifndef Regulator_H
#define Regulator_H

#include "fvCFD.H"

const scalar REG_SIGNAL_MIN = -1.0;
const scalar REG_SIGNAL_MAX = 1.0;

class Regulator
{
private:
    const fvMesh &mesh_;       // Reference to the system mesh
    dictionary regulatorDict_; // Dictionary with all regulator settings
    word regulatedFieldName_;  // Controlled process variable scalar field
    word targetPatchName_;     // Name of a patch from which the process value is read
    scalar targetValue_;       // Reference value used by the control loop
    scalar Kp_;                // Proportional control term
    scalar Ti_;                // Integral control term
    scalar Td_;                // Differential control term
    scalar error_;             // Measured error
    scalar errorIntegral_;     // Error integral w.r.t. time
    scalar oldError_;          // Old error
    label timeIndex_;          // Time index of the last update

public:
    // Initialize from mesh and dictionary
    Regulator(const fvMesh &mesh, const dictionary &dict);

    // Initialize from mesh and propertiesDict from "constant" directory
    Regulator(const fvMesh &mesh);

    // Construct a copy
    Regulator(const Regulator &reg);

    word fieldName() const { return regulatedFieldName_; }
    word targetPatchName() const { return targetPatchName_; }
    scalar targetValue() const { return targetValue_; }
    scalar Kp() const { return Kp_; }
    scalar Ti() const { return Ti_; }
    scalar Td() const { return Td_; }

    //- Read output signal from the regulator
    scalar read();

    //- Read value of the regulated field from regulation target patch
    scalar probeTargetPatch() const;

    //- Calculate face average of a value
    static scalar patchAverage(const word &fieldName, const fvPatch &patch);
};

#endif
