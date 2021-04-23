#include "controlMethod.H"


// * * * * * * * * * * * * Factory  * * * * * * * * * * * * //
const Foam::Enum<ControlMethod::controlType>
    ControlMethod::controlTypeNames({
        {controlType::twoStep, "twoStep"},
        {controlType::PID, "PID"},
    });

ControlMethod* ControlMethod::create(const dictionary &dict)
{
    const controlType type = controlTypeNames.get("mode", dict);

    switch (type)
    {
    case twoStep:
        return new TwoStepControl(dict);
    case PID:
        return new PIDControl(dict);
    default:
        FatalIOErrorInFunction(dict)
            << "    Unknown control method " << type
            << exit(FatalIOError);
        return NULL;
    }
}

// * * * * * * * * * * * * Two Step Control  * * * * * * * * * * * * //
TwoStepControl::TwoStepControl(const dictionary &dict)
  : h_(dict.getOrDefault<scalar>("h", 0.)),
    outputSignal_(0.)
{}

scalar TwoStepControl::calculate(scalar current, scalar target, scalar deltaT)
{
    const scalar targetCorrected = outputSignal_ > 0. ? target + 0.5 * h_ : target - 0.5 * h_;
    const scalar errorCorrected = targetCorrected - current;
    outputSignal_ = errorCorrected <= 0 ? 0. : 1.;
    return outputSignal_;
}

// * * * * * * * * * * * * PID Control  * * * * * * * * * * * * //
PIDControl::PIDControl(const dictionary &dict)
  : Kp_(dict.getScalar("Kp")),
    Ti_(dict.getScalar("Ti")),
    Td_(dict.getScalar("Td")),
    outputMax_(dict.getOrDefault<scalar>("outputMax", 1.)),
    outputMin_(dict.getOrDefault<scalar>("outputMin", -1.)),
    oldError_(0.),
    errorIntegral_(0.)
{}

scalar PIDControl::calculate(scalar current, scalar target, scalar deltaT)
{
    scalar error = target - current;
    errorIntegral_ += error * deltaT;
    const scalar errorDifferential = (error - oldError_) / deltaT;
    oldError_ = error;

    // Calculate output signal
    // A negliable value is added to Ti_ to prevent division by 0
    const scalar outputSignal = Kp_*(error + 1/(Ti_ + SMALL)*errorIntegral_ + Td_*errorDifferential);

    // Return result within defined regulator saturation: outputMax_ and outputMin_
    return max(min(outputSignal, outputMax_), outputMin_);
}
