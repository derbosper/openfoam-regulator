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

#include "convectiveHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::convectiveHeatFluxFvPatchScalarField::convectiveHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_("T"),
    Tsur_(0),
    kappa_(0),
    h_(0)
{}


Foam::convectiveHeatFluxFvPatchScalarField::convectiveHeatFluxFvPatchScalarField
(
    const convectiveHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    Tsur_(ptf.Tsur_),
    kappa_(ptf.kappa_),
    h_(ptf.h_)
{}


Foam::convectiveHeatFluxFvPatchScalarField::convectiveHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    Tsur_(dict.get<scalar>("Tsur")),
    kappa_(dict.get<scalar>("kappa")),
    h_(dict.get<scalar>("h"))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}


Foam::convectiveHeatFluxFvPatchScalarField::convectiveHeatFluxFvPatchScalarField
(
    const convectiveHeatFluxFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    TName_(tppsf.TName_),
    Tsur_(tppsf.Tsur_),
    kappa_(tppsf.kappa_),
    h_(tppsf.h_)

{}


Foam::convectiveHeatFluxFvPatchScalarField::convectiveHeatFluxFvPatchScalarField
(
    const convectiveHeatFluxFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    TName_(tppsf.TName_),
    Tsur_(tppsf.Tsur_),
    kappa_(tppsf.kappa_),
    h_(tppsf.h_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::convectiveHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& T =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    gradient() = - h_ / kappa_ * (T - Tsur_);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::convectiveHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
    os.writeEntry<scalar>("Tsur", Tsur_);
    os.writeEntry<scalar>("kappa", kappa_);
    os.writeEntry<scalar>("h", h_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        convectiveHeatFluxFvPatchScalarField
    );
}

// ************************************************************************* //
