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
    field using Sensor class and returns a value between outputMin and
    outputMax (default 0 and 1). Returned value is calculated by ControlMethod,
    where the regulation algorithm implementation is located.

    When modeling automatic regulation, the acutor is an OpenFOAM boundary
    condition. Usually a custom boundry condition needs to be implemented,
    inheriting or refering to any condition from regulatedPatchBase.

    The regulator is set up in the boundary condition as a subdict "regulator"
    with the following parameters:

    regulator dict:
    \table
        Property        | Description                       | Required | Default |
        targetValue     | desired value                     | yes      |         |
        field           | regulated field name              | yes      |         |
        mode            | regulation algorithm              | yes      |         |
        parameters      | mode-dependent dict with params   | depends  |         |
        sensor          | sensor dictionary                 | yes      |         |
    \endtable

\*---------------------------------------------------------------------------*/

#ifndef Regulator_H
#define Regulator_H

#include "fvCFD.H"
#include "sensor.H"
#include "controlMethod.H"


class Regulator
{
public:
    Regulator() = delete;
    ~Regulator() = default;

    // Initialize from mesh and dictionary form "constant" directory
    Regulator(const fvMesh &mesh, const dictionary &dict);

    // Initialize from mesh and default dictionary
    Regulator(const fvMesh &mesh);

    // Initialize as copy
    Regulator(const Regulator &);

    scalar targetValue() const
    {
        return targetValue_->value(mesh_.time().timeOutputValue());
    }

    // Name of the regulated field
    word fieldName() const
    {
        return sensor_->fieldName();
    }

    //- Read output signal from the regulator
    scalar read();

    //- Write to runtime dict
    void write(Ostream& os, const word dictName = "regulator") const;

private:
    // Reference to the system mesh
    const fvMesh &mesh_;

    // Sensor instance for reading field values
    const std::shared_ptr<Sensor> sensor_;

    // Control algorithm
    std::shared_ptr<ControlMethod> controlMethod_;

    // Reference value used by the control loop
    autoPtr<Function1<scalar>> targetValue_;

    // Time index of the last update
    label timeIndex_;
};

#endif
