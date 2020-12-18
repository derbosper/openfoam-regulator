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
    Foam::temperaturePIDControllerFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    This boundary condition tries to maintain a specified temperature of
    a downstream patch by regulating the temperature of itself.
    The regulation is done by implementing a PID control algorithm.

    \heading Patch usage

    \table
        Property    | Description                 | Required | Default value
        downstream  | downstream patch name       | yes      |
        targetT     | desired temperature         | yes      |
        value       | initial value of temperature| yes      |
        T           | temperature field name      | no       | T
        P           | proportional gain           | yes      |
        I           | integral gain               | yes      |
        D           | differential gain           | yes      |
        outputMin   | lower output signal limit   | yes      |
        outputMax   | upper output signal limit   | yes      |
    \endtable

    Example:

    \verbatim
    wall
    {
        type            temperaturePIDController;
        value           uniform 350;
        downstream      outlet;
        targetT         250;
        P               5;
        I               0.2;
        D               0.1;
        outputMin       -30;
        outputMax       30;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    temperaturePIDControllerFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef temperaturePIDControllerFvPatchScalarField_H
#define temperaturePIDControllerFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Switch.H"
#include "regulatorLibrary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
     Class temperaturePIDControllerFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class temperaturePIDControllerFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    Regulator regulator_;
    scalar m_; // temperature inertia momentum

public:

   //- Runtime type information
   TypeName("temperaturePIDController");


   // Constructors

        //- Construct from patch and internal field
        temperaturePIDControllerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        temperaturePIDControllerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchScalarField
        //  onto a new patch
        temperaturePIDControllerFvPatchScalarField
        (
            const temperaturePIDControllerFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        temperaturePIDControllerFvPatchScalarField
        (
            const temperaturePIDControllerFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new temperaturePIDControllerFvPatchScalarField ( *this )
            );
        }

        //- Construct as copy setting internal field reference
        temperaturePIDControllerFvPatchScalarField
        (
            const temperaturePIDControllerFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new temperaturePIDControllerFvPatchScalarField ( *this, iF )
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
