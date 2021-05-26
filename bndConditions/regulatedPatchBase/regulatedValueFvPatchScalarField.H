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
    Foam::regulatedValueFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a regulated value, which is dependend
    of the Regulator signal in the linear manner.

    The value is set accoring to formula:

        (maxValue - minValue) * regulator_signal + minValue

Usage
    \table
        Property     | Description                          | Required | Default
        maxValue     | reference value for signal = 1       | yes |
        minValue     | reference value for signal = 0       | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            regulatedValue;
        maxValue        10;
        minValue        1;
    }
    \endverbatim


SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    regulatedValueFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef regulatedValueFvPatchScalarField_H
#define regulatedValueFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Switch.H"
#include "regulator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
     Class regulatedValueFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class regulatedValueFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    Regulator regulator_;
    scalar minValue_;
    scalar maxValue_;

public:

   //- Runtime type information
   TypeName("regulatedValue");


   // Constructors

        //- Construct from patch and internal field
        regulatedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regulatedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchScalarField
        //  onto a new patch
        regulatedValueFvPatchScalarField
        (
            const regulatedValueFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        regulatedValueFvPatchScalarField
        (
            const regulatedValueFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new regulatedValueFvPatchScalarField ( *this )
            );
        }

        //- Construct as copy setting internal field reference
        regulatedValueFvPatchScalarField
        (
            const regulatedValueFvPatchScalarField&,
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
                new regulatedValueFvPatchScalarField ( *this, iF )
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
