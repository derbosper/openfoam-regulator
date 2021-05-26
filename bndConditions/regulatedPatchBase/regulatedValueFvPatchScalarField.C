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

\*---------------------------------------------------------------------------*/

#include "regulatedValueFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "sensor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regulatedValueFvPatchScalarField::
regulatedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    regulator_(p.boundaryMesh().mesh()),
    minValue_(0.),
    maxValue_(0.)
{}


Foam::regulatedValueFvPatchScalarField::
regulatedValueFvPatchScalarField
(
    const regulatedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    regulator_(ptf.regulator_),
    minValue_(ptf.minValue_),
    maxValue_(ptf.maxValue_)
{}


Foam::regulatedValueFvPatchScalarField::
regulatedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    regulator_(p.boundaryMesh().mesh(), dict.subDict("regulator")),
    minValue_(dict.get<scalar>("minValue")),
    maxValue_(dict.get<scalar>("maxValue"))
{}


Foam::regulatedValueFvPatchScalarField::
regulatedValueFvPatchScalarField
(
    const regulatedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    regulator_(ptf.regulator_)
{}


Foam::regulatedValueFvPatchScalarField::
regulatedValueFvPatchScalarField
(
    const regulatedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    regulator_(ptf.regulator_),
    minValue_(ptf.minValue_),
    maxValue_(ptf.maxValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regulatedValueFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator==((maxValue_ - minValue_) * regulator_.read() + minValue_);
    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::regulatedValueFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("minValue", minValue_);
    os.writeEntry("maxValue", maxValue_);
    regulator_.write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       regulatedValueFvPatchScalarField
   );
}


// ************************************************************************* //
