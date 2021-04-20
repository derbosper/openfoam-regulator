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

#include "temperaturePIDFixedGradientFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "sensor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperaturePIDFixedGradientFvPatchScalarField::temperaturePIDFixedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    regulator_(p.boundaryMesh().mesh()),
    Q_(0.0)
{}


Foam::temperaturePIDFixedGradientFvPatchScalarField::temperaturePIDFixedGradientFvPatchScalarField
(
    const temperaturePIDFixedGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    regulator_(ptf.regulator_),
    Q_(ptf.Q_)
{}


Foam::temperaturePIDFixedGradientFvPatchScalarField::temperaturePIDFixedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    regulator_(p.boundaryMesh().mesh(), dict.subDict("regulator")),
    Q_(dict.get<scalar>("Q"))
{
    // Initialize patch with internal field value
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}


Foam::temperaturePIDFixedGradientFvPatchScalarField::temperaturePIDFixedGradientFvPatchScalarField
(
    const temperaturePIDFixedGradientFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    regulator_(tppsf.regulator_),
    Q_(tppsf.Q_)
{}


Foam::temperaturePIDFixedGradientFvPatchScalarField::temperaturePIDFixedGradientFvPatchScalarField
(
    const temperaturePIDFixedGradientFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    regulator_(tppsf.regulator_),
    Q_(tppsf.Q_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::temperaturePIDFixedGradientFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Read thermal conductivity from transportProperties
    const dictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
    dimensionedScalar k(transportProperties.get<dimensionedScalar>("k"));

    // Calculate output signal from regulator
    const scalar outputSignal = regulator_.read();

    if (outputSignal > 0) {
        gradient() = Q_ / k.value();
    } else {
        gradient() = 0;
    }


    Info << "Regulator: value at inlet = " << Sensor::patchAverage(regulator_.fieldName(), patch().boundaryMesh()["inlet"]) << endl;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::temperaturePIDFixedGradientFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    regulator_.write(os);
    os.writeEntry("Q", Q_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        temperaturePIDFixedGradientFvPatchScalarField
    );
}

// ************************************************************************* //
