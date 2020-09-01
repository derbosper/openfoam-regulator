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
    Test boundary condition

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    temperaturePIDControllerFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef temperaturePIDControllerFvPatchScalarField_H
#define temperaturePIDControllerFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Switch.H"

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
    // Private data

        //- Name of the downstream face zone
        const word downstreamName_;

        //- Desired temperature difference at downstream
        const scalar targetT_;

        //- Name of the temperature field
        const word TName_;

        //- Name of the flux field
        const word phiName_;

        //- Proportional gain
        const scalar P_;

        //- Integral gain
        const scalar I_;

        //- Derivative gain
        const scalar D_;

        //- Controller output signal limits
        const scalar outputMax_;
        const scalar outputMin_;

        //- Temperature
        scalar T_;

        //- Error
        scalar error_;

        //- Error integral w.r.t. time
        scalar errorIntegral_;

        //- Old error
        scalar oldError_;

        //- Time index of the last update
        label timeIndex_;


    // Private member functions

        //- Calculate the average on a face patch
        const scalar patchAverage
        (
            const word& fieldName,
            const fvPatch& patch
        );


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
