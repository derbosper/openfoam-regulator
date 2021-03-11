#include "regulatorLibrary.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

scalar Regulator::patchAverage(const word &fieldName, const fvPatch &patch)
{
    const fvPatchField<scalar> &field =
        patch.lookupPatchField<volScalarField, scalar>(fieldName);

    return gSum(field * patch.magSf()) / gSum(patch.magSf());
}

const Foam::Enum<Regulator::operationMode>
    Regulator::operationModeNames({
        {operationMode::twoStep, "twoStep"},
        {operationMode::PID, "PID"},
        {operationMode::PIDRestricted, "PIDRestricted"},
    });

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Regulator::Regulator(const fvMesh &mesh, const dictionary &dict)
    : mesh_(mesh),
      regulatedFieldName_(dict.getWord("fieldName")),
      targetPatchName_(dict.getWord("patchName")),
      targetValue_(dict.getScalar("targetValue")),
      mode_(operationModeNames.get("mode", dict)),
      Kp_(dict.getScalar("Kp")),
      Ti_(dict.getScalar("Ti")),
      Td_(dict.getScalar("Td")),
      error_(0.),         // TODO change initial value to lookupOrDefault("error", 0.)
      errorIntegral_(0.), // TODO as above
      oldError_(0.),
      timeIndex_(mesh.time().timeIndex())
{}

Regulator::Regulator(const fvMesh &mesh)
    : mesh_(mesh),
    regulatedFieldName_(word::null),
    targetPatchName_(word::null),
    targetValue_(0),
    mode_(PID),
    Kp_(0),
    Ti_(0),
    Td_(0),
    error_(0),
    errorIntegral_(0),
    oldError_(0),
    timeIndex_(mesh.time().timeIndex())
{}

Regulator::Regulator(const Regulator &reg)
    : mesh_(reg.mesh_),
      regulatedFieldName_(reg.regulatedFieldName_),
      targetPatchName_(reg.targetPatchName_),
      targetValue_(reg.targetValue_),
      mode_(reg.mode_),
      Kp_(reg.Kp_),
      Ti_(reg.Ti_),
      Td_(reg.Td_),
      error_(reg.error_),
      errorIntegral_(reg.errorIntegral_),
      oldError_(reg.oldError_),
      timeIndex_(reg.timeIndex_)
{}

// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * *//

scalar Regulator::probeTargetPatch() const
{
    // Get the target patch average field value
    const fvPatch &targetPatch = mesh_.boundary()[targetPatchName_];
    const scalar result = patchAverage(regulatedFieldName_, targetPatch);
    return result;
}

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
    const scalar currentRegulatedPatchValue = probeTargetPatch();

    // Calculate errors
    scalar result = 0.;
    error_ = targetValue_ - currentRegulatedPatchValue;
    if (mode_ == twoStep)
    {
        result = error_ <= 0 ? 0 : 1;
    }
    else
    {
        errorIntegral_ += 0.5*(error_ + oldError_)*deltaT;
        const scalar errorDifferential = (error_ - oldError_) / deltaT;

        // Calculate output signal
        // A negliable value is added to Ti_ to prevent division by 0
        const scalar outputSignal = Kp_*(error_ + 1/(Ti_ + 1e-7)*errorIntegral_ + Td_*errorDifferential);

        if (mode_ == PIDRestricted)
        {
            // Return result within defined SIGNAL_MIN and SIGNAL_MAX
            result = max(min(outputSignal, REG_SIGNAL_MAX), REG_SIGNAL_MIN);
        }
        else if (mode_ == PID)
        {
            // Return raw computed signal from PID algorithm
            result = outputSignal;
        }
        else
        {
            FatalIOError << "Unknown regulator mode" << endl;
            exit(FatalIOError);
        }
    }

    Info << "Time index: " << mesh_.time().timeIndex() << endl;
    Info << "regulatedFieldName: " << regulatedFieldName_ << endl;
    Info << "targetPatchName: " << targetPatchName_ << endl;
    Info << "targetValue: " << targetValue_ << endl;
    Info << "currentRegulatedPatchValue: " << currentRegulatedPatchValue << endl;
    Info << "Error: " << error_ << endl;

    return result;
}
