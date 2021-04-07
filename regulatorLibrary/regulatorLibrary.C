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
    });

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Regulator::Regulator(const fvMesh &mesh, const dictionary &dict)
    : mesh_(mesh),
      regulatedFieldName_(dict.getWord("fieldName")),
      targetPatchName_(dict.getWord("patchName")),
      targetValue_(dict.getScalar("targetValue")),
      mode_(operationModeNames.get("mode", dict)),
      timeIndex_(mesh.time().timeIndex()),
      error_(0.),
      outputSignal_(0.),
      h_(0.),
      Kp_(0.),
      Ti_(0.),
      Td_(0.),
      outputMax_(1.),
      outputMin_(0.),
      oldError_(0.),
      errorIntegral_(0.)
{
    switch (mode_)
    {
        // Two step regulator returns either 0 or 1
        case twoStep:
        {
            if (dict.found("outputMax") || dict.found("outputMin"))
            {
                FatalIOError << "outputMin and outputMax values cannot be set "
                << "in twoStep mode, as they are always 0. and 1."
                << exit(FatalIOError);
            }
            h_ = dict.getOrDefault<scalar>("h", 0.);
            break;
        }
        // PID returns a value between outputMin and outputMax, defaults to (-1, 1)
        case PID:
        {
            Kp_ = dict.getScalar("Kp");
            Ti_ = dict.getScalar("Ti");
            Td_ = dict.getScalar("Td");
            outputMax_ = dict.getOrDefault<scalar>("outputMax", 1.);
            outputMin_ = dict.getOrDefault<scalar>("outputMin", -1.);
            break;
        }
        default:
        {
            FatalIOError << "Unknown regulator mode: " << mode_ << endl;
            exit(FatalIOError);
        }
    }
}

Regulator::Regulator(const fvMesh &mesh)
    : mesh_(mesh),
    regulatedFieldName_(word::null),
    targetPatchName_(word::null),
    targetValue_(0),
    mode_(PID),
    timeIndex_(mesh.time().timeIndex()),
    error_(0),
    outputSignal_(0),
    h_(0),
    Kp_(0),
    Ti_(0),
    Td_(0),
    outputMax_(1),
    outputMin_(0),
    oldError_(0),
    errorIntegral_(0)
{}

Regulator::Regulator(const Regulator &reg)
    : mesh_(reg.mesh_),
      regulatedFieldName_(reg.regulatedFieldName_),
      targetPatchName_(reg.targetPatchName_),
      targetValue_(reg.targetValue_),
      mode_(reg.mode_),
      timeIndex_(reg.timeIndex_),
      error_(reg.error_),
      outputSignal_(reg.outputSignal_),
      h_(reg.h_),
      Kp_(reg.Kp_),
      Ti_(reg.Ti_),
      Td_(reg.Td_),
      outputMax_(reg.outputMax_),
      outputMin_(reg.outputMin_),
      oldError_(reg.oldError_),
      errorIntegral_(reg.errorIntegral_)
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
    }

    // Get the target patch average field value
    const scalar currentRegulatedPatchValue = probeTargetPatch();

    // Calculate errors
    error_ = targetValue_ - currentRegulatedPatchValue;

    switch (mode_)
    {
        case twoStep:
        {
            //- Przy istnieniu histerezy wprowadzana jest zmienna zastępczna y'_z taka,
            //- że y'_z = y + 0.5h, dla sygnału wyjściowego 1 albo y - 0.5h dla sygnału 0
            const scalar targetCorrected = outputSignal_ > 0. ? targetValue_ + 0.5 * h_ : targetValue_ - 0.5 * h_;
            const scalar errorCorrected = targetCorrected - currentRegulatedPatchValue;
            outputSignal_ = errorCorrected <= 0 ? 0. : 1.;
            break;
        }
        case PID:
        {
            oldError_ = error_;
            errorIntegral_ += error_ * deltaT;
            const scalar errorDifferential = (error_ - oldError_) / deltaT;

            // Calculate output signal
            // A negliable value is added to Ti_ to prevent division by 0
            const scalar outputSignal = Kp_*(error_ + 1/(Ti_ + SMALL)*errorIntegral_ + Td_*errorDifferential);

            // Return result within defined regulator saturation: outputMax_ and outputMin_
            outputSignal_ = max(min(outputSignal, outputMax_), outputMin_);
            break;
        }
    }

    Info << "Regulator: regulatedFieldName = " << regulatedFieldName_ << endl;
    Info << "Regulator: targetPatchName = " << targetPatchName_ << endl;
    Info << "Regulator: targetValue = " << targetValue_ << endl;
    Info << "Regulator: mode = " << mode_ << endl;
    Info << "Regulator: currentRegulatedPatchValue = " << currentRegulatedPatchValue << endl;
    Info << "Regulator: error = " << error_ << endl;
    Info << "Regulator: outputSignal = " << outputSignal_ << endl;

    return outputSignal_;
}

void Regulator::write(Ostream& os, const word dictName) const
{
    os.beginBlock(dictName);
    os.writeEntry("fieldName", regulatedFieldName_);
    os.writeEntry("patchName", targetPatchName_);
    os.writeEntry("targetValue", targetValue_);
    os.writeEntry("mode", operationModeNames.get(mode_));
    os.endBlock();
}