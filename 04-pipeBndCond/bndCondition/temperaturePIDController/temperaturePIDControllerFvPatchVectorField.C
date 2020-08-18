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

#include "temperaturePIDControllerFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "linear.H"
#include "syncTools.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const double Foam::temperaturePIDControllerFvPatchVectorField::patchAverage
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
Foam::temperaturePIDControllerFvPatchVectorField::faceTemperature() const
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

Foam::temperaturePIDControllerFvPatchVectorField::
temperaturePIDControllerFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
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


Foam::temperaturePIDControllerFvPatchVectorField::
temperaturePIDControllerFvPatchVectorField
(
    const temperaturePIDControllerFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
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


Foam::temperaturePIDControllerFvPatchVectorField::
temperaturePIDControllerFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    downstreamName_(dict.lookup("downstream")),
    targetT_(dict.get<scalar>("targetT")),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    P_(dict.get<scalar>("P")),
    I_(dict.get<scalar>("I")),
    D_(dict.get<scalar>("D")),
    T_(0),
    error_(dict.lookupOrDefault<scalar>("error", 0)),
    errorIntegral_(dict.lookupOrDefault<scalar>("errorIntegral", 0)),
    oldT_(0),
    oldError_(0),
    oldErrorIntegral_(0),
    timeIndex_(db().time().timeIndex())
{
    // calls the = operator to assign the value to the faces held by this BC
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    // Get the mesh
    const fvMesh& mesh(patch().boundaryMesh().mesh());
    const fvPatch& targetPatch = mesh.boundary()[downstreamName_];
    T_ = patchAverage(TName_, targetPatch);
    updateCoeffs();
}


Foam::temperaturePIDControllerFvPatchVectorField::
temperaturePIDControllerFvPatchVectorField
(
    const temperaturePIDControllerFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
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


Foam::temperaturePIDControllerFvPatchVectorField::
temperaturePIDControllerFvPatchVectorField
(
    const temperaturePIDControllerFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
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

void Foam::temperaturePIDControllerFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }


    if (db().foundObject<volScalarField>(TName_))
    {
        Info << "found T " << TName_ << endl;
    }
    else
    {
        WarningInFunction
            << "The temperature field name, \"TName\", is \"" << TName_ << "\", "
            << "but a field of that name was not found." 
            << endl << endl;
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::temperaturePIDControllerFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);

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
       fvPatchVectorField,
       temperaturePIDControllerFvPatchVectorField
   );
}


// ************************************************************************* //
