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
      oldError_(0.),
      timeIndex_(mesh.time().timeIndex())
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
    regulatedFieldName_ = regulatorDict_.getWord("fieldName");
    targetPatchName_ = regulatorDict_.getWord("patchName");
    targetValue_ = regulatorDict_.getScalar("targetValue");
    Kp_ = regulatorDict_.getScalar("Kp");
    Ti_ = regulatorDict_.getScalar("Ti");
    Td_ = regulatorDict_.getScalar("Td");
}

// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * *//

scalar Regulator::read()
{
    // Get the time step
    const scalar deltaT(mesh_.time().deltaTValue());

    // Update the old-time quantities
    if (timeIndex_ != mesh_.time().timeIndex())
    {
        timeIndex_ = mesh_.time().timeIndex();
        oldError_ = error_;
    }

    // Get the target patch average field value
    const fvPatch &targetPatch = mesh_.boundary()[targetPatchName_];
    const scalar currentRegulatedPatchValue = patchAverage(regulatedFieldName_, targetPatch);

    // Calculate errors
    error_ = targetValue_ - currentRegulatedPatchValue;
    errorIntegral_ += 0.5*(error_ + oldError_)*deltaT;
    const scalar errorDifferential = (error_ - oldError_) / deltaT;

    // Calculate output signal
    // A negliable value is added to Ti_ to prevent division by 0
    const scalar outputSignal = Kp_*(error_ + 1/(Ti_ + 1e-7)*errorIntegral_ + Td_*errorDifferential);

    // Return result within defined SIGNAL_MIN and SIGNAL_MAX
    const scalar result = max(min(outputSignal, REG_SIGNAL_MAX), REG_SIGNAL_MIN);

    Info << "Time index: " << mesh_.time().timeIndex() << endl;
    Info << "regulatedFieldName: " << regulatedFieldName_ << endl;
    Info << "targetPatchName: " << targetPatchName_ << endl;
    Info << "targetValue: " << targetValue_ << endl;
    Info << "Kp: " << Kp_ << endl;
    Info << "Ti: " << Ti_ << endl;
    Info << "Td: " << Td_ << endl;
    Info << "currentRegulatedPatchValue: " << currentRegulatedPatchValue << endl;
    Info << "Error: " << error_ << endl;

    return result;
}
