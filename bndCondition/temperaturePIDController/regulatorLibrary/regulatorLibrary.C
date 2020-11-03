#include "regulatorLibrary.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

scalar Regulator::patchAverage(const word &fieldName, const fvPatch &patch)
{
    const fvPatchField<scalar> &field =
        patch.lookupPatchField<volScalarField, scalar>(fieldName);

    return gSum(field * patch.magSf()) / gSum(patch.magSf());
}

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Regulator::Regulator(const fvMesh &mesh)
    : mesh_(mesh),
      error_(0.),         // TODO change initial value to lookupOrDefault("error", 0.)
      errorIntegral_(0.), // TODO as above
      oldError_(0.), timeIndex_(mesh.time().timeIndex())
{
    // Get access to a custom dictionary
    const word dictName("regulatorProperties");

    // Create and input-output object - this holds the path to the dict and its name
    IOobject dictIO(dictName,               // name of the file
                    mesh.time().constant(), // path to where the file is
                    mesh,                   // reference to the mesh needed by the constructor
                    IOobject::MUST_READ     // indicate that reading this dictionary
                                            // is compulsory
    );

    // Check the if the dictionary is present and follows the OF format
    if (!dictIO.typeHeaderOk<dictionary>(true))
        FatalErrorIn("regulatorLibrary.C")
            << "Cannot open specified refinement dictionary " << dictName
            << exit(FatalError);
    else
        Info << "Dictionary OK" << endl;

    // Initialise the dictionary object
    regulatorDict_ = IOdictionary(dictIO);

    // Read various pieces of information from the main part of the dictionary
    fieldName_ = regulatorDict_.getWord("fieldName");
    patchName_ = regulatorDict_.getWord("patchName");
    targetValue_ = regulatorDict_.getScalar("targetValue");
    P_ = regulatorDict_.getScalar("P");
    I_ = regulatorDict_.getScalar("I");
    D_ = regulatorDict_.getScalar("D");
}

// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * *//

scalar Regulator::read()
{
    Info << "Time index: " << mesh_.time().timeIndex() << endl;
    Info << "fieldName: " << fieldName_ << endl;
    Info << "patchName: " << patchName_ << endl;
    Info << "targetValue: " << targetValue_ << endl;
    Info << "P: " << P_ << endl;
    Info << "I: " << I_ << endl;
    Info << "D: " << D_ << endl;

    // Get the time step
    const scalar deltaT(mesh_.time().deltaTValue());

    // Update the old-time quantities
    if (timeIndex_ != mesh_.time().timeIndex())
    {
        timeIndex_ = mesh_.time().timeIndex();
        oldError_ = error_;
    }

    // Get the target patch average field value
    const fvPatch &targetPatch = mesh_.boundary()[patchName_];
    const scalar targetPatchValue = patchAverage(fieldName_, targetPatch);

    // Calculate errors
    error_ = targetValue_ - targetPatchValue;
    errorIntegral_ += 0.5*(error_ + oldError_)*deltaT;
    const scalar errorDifferential = (error_ - oldError_) / deltaT;

    // Calculate output signal
    const scalar outputSignal = P_*error_ + I_*errorIntegral_ + D_*errorDifferential;

    // Return result within defined SIGNAL_MIN and SIGNAL_MAX
    const scalar result = max(min(outputSignal, REG_SIGNAL_MAX), REG_SIGNAL_MIN);

    return result;
}
