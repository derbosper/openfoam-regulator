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
        P               | proportional gain           | yes      |
        I               | integral gain               | yes      |
        D               | differential gain           | yes      |
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
    //- Reference to mesh
    const fvMesh& mesh_;

    //- Dictionary with all regulator settings
    dictionary regulatorDict_;

    //- Properties read from the dictionary
    word fieldName_;
    word patchName_;
    scalar targetValue_;
    scalar P_;
    scalar I_;
    scalar D_;

    //- Error
    scalar error_;

    //- Error integral w.r.t. time
    scalar errorIntegral_;

    //- Old error
    scalar oldError_;

    //- Time index of the last update
    label timeIndex_;

    //- Calculate face average of a value
    scalar patchAverage(const word& fieldName, const fvPatch& patch) const;

public:
    Regulator(const fvMesh &mesh);
    scalar read();
};

#endif
