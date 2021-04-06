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


class Regulator
{
public:
    // Initialize from mesh and dictionary form "constant" directory
    Regulator(const fvMesh &mesh, const dictionary &dict);

    // Initialize from mesh and default dictionary
    Regulator(const fvMesh &mesh);

    // Construct a copy
    Regulator(const Regulator &reg);

    //- Operation mode enumeration
    enum operationMode
    {
        twoStep, // on/off control
        PID,
    };

    static const Enum<operationMode> operationModeNames;

    const word fieldName() const { return regulatedFieldName_; }
    const word targetPatchName() const { return targetPatchName_; }
    scalar targetValue() const { return targetValue_; }

    //- Read output signal from the regulator
    scalar read();

    //- Write to runtime dict
    void write(Ostream& os, const word dictName = "regulator") const;

    //- Read value of the regulated field from regulation target patch
    scalar probeTargetPatch() const;

    //- Calculate face average of a value
    static scalar patchAverage(const word &fieldName, const fvPatch &patch);

private:
    const fvMesh &mesh_;            // Reference to the system mesh
    const word regulatedFieldName_; // Controlled process variable scalar field
    const word targetPatchName_;    // Name of a patch from which the process value is read
    scalar targetValue_;            // Reference value used by the control loop
    operationMode mode_;            // Control type mode
    label timeIndex_;               // Time index of the last update
    scalar error_;                  // Measured error
    // on/off specific
    scalar outputSignal_;           // Current output signal (PID can be implemented without this member)
    // PID specific
    scalar Kp_;                     // Proportional control term
    scalar Ti_;                     // Integral control term
    scalar Td_;                     // Differential control term
    scalar outputMax_;              // Max output signal (saturation)
    scalar outputMin_;              // Min output signal
    scalar oldError_;               // Old error
    scalar errorIntegral_;          // Error integral w.r.t. time
};

#endif
