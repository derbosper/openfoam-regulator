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

\*---------------------------------------------------------------------------*/

#include "regulatedValueFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regulatedValueFvPatchVectorField::
regulatedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(p.size()),
    minValue_(p.size()),
    regulator_(p.boundaryMesh().mesh())
{}


Foam::regulatedValueFvPatchVectorField::
regulatedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    maxValue_("maxValue", dict, p.size()),
    minValue_("minValue", dict, p.size()),
    regulator_(p.boundaryMesh().mesh(), dict.subDict("regulator"))
{
    tmp<vectorField> tvalues(maxValue_*patch().nf());
    fvPatchVectorField::operator=(tvalues);
}


Foam::regulatedValueFvPatchVectorField::
regulatedValueFvPatchVectorField
(
    const regulatedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(ptf.maxValue_, mapper, pTraits<scalar>::zero),
    minValue_(ptf.minValue_, mapper, pTraits<scalar>::zero),
    regulator_(ptf.regulator_)
{
    tmp<vectorField> tvalues(maxValue_*patch().nf());
    fvPatchVectorField::operator=(tvalues);
}


Foam::regulatedValueFvPatchVectorField::
regulatedValueFvPatchVectorField
(
    const regulatedValueFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    maxValue_(ptf.maxValue_),
    minValue_(ptf.minValue_),
    regulator_(ptf.regulator_)
{}


Foam::regulatedValueFvPatchVectorField::
regulatedValueFvPatchVectorField
(
    const regulatedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    maxValue_(ptf.maxValue_),
    minValue_(ptf.minValue_),
    regulator_(ptf.regulator_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regulatedValueFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField outputValue = (maxValue_ - minValue_) * regulator_.read() + minValue_;
    tmp<vectorField> tvalues = outputValue*patch().nf();
    fvPatchVectorField::operator=(tvalues);
    fvPatchVectorField::updateCoeffs();
}


void Foam::regulatedValueFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    maxValue_.writeEntry("maxValue", os);
    minValue_.writeEntry("minValue", os);
    regulator_.write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        regulatedValueFvPatchVectorField
    );
}

// ************************************************************************* //
