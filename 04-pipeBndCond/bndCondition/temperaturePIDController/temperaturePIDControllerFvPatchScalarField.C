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

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const Foam::scalar Foam::temperaturePIDControllerFvPatchScalarField::patchAverage
(
    const word& fieldName,
    const fvPatch& patch
)
{
    const fvPatchField<scalar>& field =
            patch.lookupPatchField<volScalarField, scalar>(fieldName);

    return gSum(field*patch.magSf())/gSum(patch.magSf());
}

const Foam::surfaceScalarField&
Foam::temperaturePIDControllerFvPatchScalarField::faceTemperature() const
{
    const word TfName(TName_ + "f");

    const volScalarField& T = db().lookupObject<volScalarField>(TName_);

    surfaceScalarField* TfPtr = db().getObjectPtr<surfaceScalarField>(TfName);

    if (!TfPtr)
    {
        TfPtr = new surfaceScalarField(TfName, linearInterpolate(T));
        TfPtr->store();
    }

    surfaceScalarField& Tf = *TfPtr;

    if (!Tf.upToDate(T))
    {
        Tf = linearInterpolate(T);
    }

    return Tf;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperaturePIDControllerFvPatchScalarField::
temperaturePIDControllerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    downstreamName_(word::null),
    targetT_(0),
    TName_("T"),
    phiName_("phi"),
    P_(0),
    I_(0),
    D_(0),
    T_(0),
    error_(0),
    errorIntegral_(0),
    oldT_(0),
    oldError_(0),
    oldErrorIntegral_(0),
    timeIndex_(db().time().timeIndex())
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
    downstreamName_(ptf.downstreamName_),
    targetT_(ptf.targetT_),
    TName_(ptf.TName_),
    phiName_(ptf.phiName_),
    P_(ptf.P_),
    I_(ptf.I_),
    D_(ptf.D_),
    T_(ptf.T_),
    error_(ptf.error_),
    errorIntegral_(ptf.errorIntegral_),
    oldT_(ptf.oldT_),
    oldError_(ptf.oldError_),
    oldErrorIntegral_(ptf.oldErrorIntegral_),
    timeIndex_(ptf.timeIndex_)
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
    downstreamName_(dict.lookup("downstream")),
    targetT_(dict.get<scalar>("targetT")),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    P_(dict.get<scalar>("P")),
    I_(dict.get<scalar>("I")),
    D_(dict.get<scalar>("D")),
    T_(0.),
    error_(dict.lookupOrDefault<scalar>("error", 0)),
    errorIntegral_(dict.lookupOrDefault<scalar>("errorIntegral", 0)),
    oldT_(0),
    oldError_(0),
    oldErrorIntegral_(0),
    timeIndex_(db().time().timeIndex())
{}


Foam::temperaturePIDControllerFvPatchScalarField::
temperaturePIDControllerFvPatchScalarField
(
    const temperaturePIDControllerFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    downstreamName_(ptf.downstreamName_),
    targetT_(ptf.targetT_),
    TName_(ptf.TName_),
    phiName_(ptf.phiName_),
    P_(ptf.P_),
    I_(ptf.I_),
    D_(ptf.D_),
    T_(ptf.T_),
    error_(ptf.error_),
    errorIntegral_(ptf.errorIntegral_),
    oldT_(ptf.oldT_),
    oldError_(ptf.oldError_),
    oldErrorIntegral_(ptf.oldErrorIntegral_),
    timeIndex_(ptf.timeIndex_)
{}


Foam::temperaturePIDControllerFvPatchScalarField::
temperaturePIDControllerFvPatchScalarField
(
    const temperaturePIDControllerFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    downstreamName_(ptf.downstreamName_),
    targetT_(ptf.targetT_),
    TName_(ptf.TName_),
    phiName_(ptf.phiName_),
    P_(ptf.P_),
    I_(ptf.I_),
    D_(ptf.D_),
    T_(ptf.T_),
    error_(ptf.error_),
    errorIntegral_(ptf.errorIntegral_),
    oldT_(ptf.oldT_),
    oldError_(ptf.oldError_),
    oldErrorIntegral_(ptf.oldErrorIntegral_),
    timeIndex_(ptf.timeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::temperaturePIDControllerFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the mesh
    const fvMesh& mesh(patch().boundaryMesh().mesh());

    // Get the time step
    const scalar deltaT(db().time().deltaTValue());

    // Get the downstream temperature field
    const fvPatch& downstreamPatch = patch().boundaryMesh()[downstreamName_];
    const scalar downstreamAvg = patchAverage(TName_, downstreamPatch);
    Info << "Measured temperature: " << downstreamAvg << endl;

    // Get this temperature field
    const fvPatch& thisPatch = patch().boundaryMesh()[patch().name()];
    const scalar thisAvg = patchAverage(TName_, thisPatch);
    Info << "Wall temperature " << thisAvg << endl;

    // This is wrong but assings value OK
    error_ = targetT_ - thisAvg;
    operator==(thisAvg + error_);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::temperaturePIDControllerFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);

    os.writeEntry("targetT", targetT_);
    os.writeEntry("downstream", downstreamName_);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
    os.writeEntry("P", P_);

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
