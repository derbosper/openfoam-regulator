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

#include "temperaturePIDFixedValueFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "sensor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperaturePIDFixedValueFvPatchScalarField::
temperaturePIDFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    regulator_(p.boundaryMesh().mesh())
{}


Foam::temperaturePIDFixedValueFvPatchScalarField::
temperaturePIDFixedValueFvPatchScalarField
(
    const temperaturePIDFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    regulator_(ptf.regulator_)
{}


Foam::temperaturePIDFixedValueFvPatchScalarField::
temperaturePIDFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    regulator_(p.boundaryMesh().mesh(), dict.subDict("regulator"))
{}


Foam::temperaturePIDFixedValueFvPatchScalarField::
temperaturePIDFixedValueFvPatchScalarField
(
    const temperaturePIDFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    regulator_(ptf.regulator_)
{}


Foam::temperaturePIDFixedValueFvPatchScalarField::
temperaturePIDFixedValueFvPatchScalarField
(
    const temperaturePIDFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    regulator_(ptf.regulator_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::temperaturePIDFixedValueFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get temperature of this boundary
    const fvPatch& thisPatch = patch().boundaryMesh()[patch().name()];
    const scalar thisTemp = Sensor::patchAverage("T", thisPatch);

    // Calculate output signal
    const scalar outputSignal = regulator_.read();

    // Adjust the temperature of a wall based on output signal of the regulator
    // Assume that the temperature change is 50*sign(signal) for the testing
    // purposes, with no physical meaning
    // -- sign(0.) gives 1, so we need to substract small number to achieve no
    // -- no control signal when regulator output is 0.
    const scalar controlSignal = 50*sign(outputSignal - SMALL);
    operator==(thisTemp + controlSignal);

    Info << "Wall temperature: " << thisTemp << endl;
    Info << "Control signal: " << controlSignal << endl;

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::temperaturePIDFixedValueFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       temperaturePIDFixedValueFvPatchScalarField
   );
}


// ************************************************************************* //
