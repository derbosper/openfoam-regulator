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

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

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
    deltaT_(1),
    TName_("T"),
    phiName_("phi"),
    P_(0),
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
    deltaT_(ptf.deltaT_),
    TName_(ptf.TName_),
    phiName_(ptf.phiName_),
    P_(ptf.P_),
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
    deltaT_(dict.get<scalar>("deltaT")),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    P_(dict.get<scalar>("P")),
{}


Foam::temperaturePIDControllerFvPatchVectorField::
temperaturePIDControllerFvPatchVectorField
(
    const temperaturePIDControllerFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    downstreamName_(ptf.downstreamName_),
    deltaT_(ptf.deltaT_),
    TName_(ptf.TName_),
    phiName_(ptf.phiName_),
    P_(ptf.P_),
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
    deltaT_(ptf.deltaT_),
    TName_(ptf.TName_),
    phiName_(ptf.phiName_),
    P_(ptf.P_),
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::temperaturePIDControllerFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the mesh
    const fvMesh& mesh(patch().boundaryMesh().mesh());

    // Get the time step
    const scalar deltaT(db().time().deltaTValue());

    // Get the flux field
    const surfaceScalarField& phi
    (
        db().lookupObject<surfaceScalarField>(phiName_)
    );

    
    scalar deltaT = deltaT_;
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

    os.writeEntry("deltaT", deltaT_);
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
