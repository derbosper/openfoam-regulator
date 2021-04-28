/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::regulatedVelocityFvPatchVectorField

Group
    grpGenericBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a surface-normal vector boundary condition
    by its magnitude.

Usage
    \table
        Property     | Description                          | Required | Default
        maxValue     | reference value                      | yes |
        ramp         | time-based ramping                   | no  |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            regulatedVelocity;
        maxValue        uniform -10;           // 10 INTO the domain
    }
    \endverbatim

Note
    Sign conventions:
    - the value is positive for outward-pointing vectors

See also
    Foam::fixedValueFvPatchField

SourceFiles
    regulatedVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef regulatedVelocityFvPatchVectorField_H
#define regulatedVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class regulatedVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class regulatedVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        scalarField maxValue_;
        scalarField minValue_;


public:

    //- Runtime type information
    TypeName("regulatedVelocity");


    // Constructors

        //- Construct from patch and internal field
        regulatedVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regulatedVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  regulatedVelocityFvPatchVectorField
        //  onto a new patch
        regulatedVelocityFvPatchVectorField
        (
            const regulatedVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        regulatedVelocityFvPatchVectorField
        (
            const regulatedVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new regulatedVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        regulatedVelocityFvPatchVectorField
        (
            const regulatedVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new regulatedVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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
