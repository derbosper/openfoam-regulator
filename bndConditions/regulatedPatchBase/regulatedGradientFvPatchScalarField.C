/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "regulatedGradientFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "sensor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regulatedGradientFvPatchScalarField::regulatedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    regulator_(p.boundaryMesh().mesh()),
    minValue_(0.0),
    maxValue_(0.0)
{}


Foam::regulatedGradientFvPatchScalarField::regulatedGradientFvPatchScalarField
(
    const regulatedGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    regulator_(ptf.regulator_),
    minValue_(ptf.minValue_),
    maxValue_(ptf.maxValue_)
{}


Foam::regulatedGradientFvPatchScalarField::regulatedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    regulator_(p.boundaryMesh().mesh(), dict.subDict("regulator")),
    minValue_(dict.get<scalar>("minValue")),
    maxValue_(dict.get<scalar>("maxValue"))
{
    // Initialize patch with internal field value
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}


Foam::regulatedGradientFvPatchScalarField::regulatedGradientFvPatchScalarField
(
    const regulatedGradientFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    regulator_(tppsf.regulator_),
    minValue_(tppsf.minValue_),
    maxValue_(tppsf.maxValue_)
{}


Foam::regulatedGradientFvPatchScalarField::regulatedGradientFvPatchScalarField
(
    const regulatedGradientFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    regulator_(tppsf.regulator_),
    minValue_(tppsf.minValue_),
    maxValue_(tppsf.maxValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regulatedGradientFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Calculate output signal from regulator
    gradient() = (maxValue_ - minValue_) * regulator_.read() + minValue_;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::regulatedGradientFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
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
        regulatedGradientFvPatchScalarField
    );
}

// ************************************************************************* //
