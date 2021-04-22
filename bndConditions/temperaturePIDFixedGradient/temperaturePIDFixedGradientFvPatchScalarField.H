/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::temperaturePIDFixedGradientFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a total temperature condition.

Usage
    \table
        Property     | Description             | Required    | Default value
        Q            | Power [W]               | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            temperaturePIDFixedGradient;
        Q               1000;
    }
    \endverbatim

SourceFiles
    temperaturePIDFixedGradientFvPatchScalarField.C

See also
    Foam::fixedValueFvPatchField

\*---------------------------------------------------------------------------*/

#ifndef temperaturePIDFixedGradientFvPatchScalarField_H
#define temperaturePIDFixedGradientFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"
#include "regulator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class temperaturePIDFixedGradientFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class temperaturePIDFixedGradientFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    Regulator regulator_;
    scalar Q_;  // Heat power [W]
    scalar kappa_; // Heat conducivity of surrounding medium [W/m/K]

public:

    //- Runtime type information
    TypeName("temperaturePIDFixedGradient");


    // Constructors

        //- Construct from patch and internal field
        temperaturePIDFixedGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        temperaturePIDFixedGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given temperaturePIDFixedGradientFvPatchScalarField
        //  onto a new patch
        temperaturePIDFixedGradientFvPatchScalarField
        (
            const temperaturePIDFixedGradientFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        temperaturePIDFixedGradientFvPatchScalarField
        (
            const temperaturePIDFixedGradientFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new temperaturePIDFixedGradientFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        temperaturePIDFixedGradientFvPatchScalarField
        (
            const temperaturePIDFixedGradientFvPatchScalarField&,
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
                new temperaturePIDFixedGradientFvPatchScalarField(*this, iF)
            );
        }

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
