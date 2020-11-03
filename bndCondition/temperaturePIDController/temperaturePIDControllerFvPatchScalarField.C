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

#include "temperaturePIDControllerFvPatchScalarField.H"
#include "pressurePIDControlInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "linear.H"
#include "syncTools.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperaturePIDControllerFvPatchScalarField::
temperaturePIDControllerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    regulator_(p.boundaryMesh().mesh())
{}


Foam::temperaturePIDControllerFvPatchScalarField::
temperaturePIDControllerFvPatchScalarField
(
    const temperaturePIDControllerFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    regulator_(ptf.regulator_)
{}


Foam::temperaturePIDControllerFvPatchScalarField::
temperaturePIDControllerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    regulator_(p.boundaryMesh().mesh())
{}


Foam::temperaturePIDControllerFvPatchScalarField::
temperaturePIDControllerFvPatchScalarField
(
    const temperaturePIDControllerFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    regulator_(ptf.regulator_)
{}


Foam::temperaturePIDControllerFvPatchScalarField::
temperaturePIDControllerFvPatchScalarField
(
    const temperaturePIDControllerFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    regulator_(ptf.regulator_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::temperaturePIDControllerFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get this temperature field
    const fvPatch& thisPatch = patch().boundaryMesh()[patch().name()];
    const scalar thisTemp = Regulator::patchAverage(regulator_.fieldName(), thisPatch);

    // Calculate output signal
    const scalar outputSignal = regulator_.read();

    // Set patch temperature
    operator==(thisTemp + 10*outputSignal);

    Info << "Wall temperature: " << thisTemp << endl;
    Info << "Output signal: " << outputSignal << endl;

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::temperaturePIDControllerFvPatchScalarField::write
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
       temperaturePIDControllerFvPatchScalarField
   );
}


// ************************************************************************* //
