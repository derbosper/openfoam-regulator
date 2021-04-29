#include "controlMethod.H"

// * * * * * * * * * * * * Utility function  * * * * * * * * * * * * //
const dictionary ControlMethod::parameters(const dictionary &dict)
{
    return dict.subOrEmptyDict("parameters");
};

// * * * * * * * * * * * * Factory  * * * * * * * * * * * * //
const Foam::Enum<ControlMethod::controlType>
    ControlMethod::controlTypeNames({
        {controlType::twoStep, "twoStep"},
        {controlType::PID, "PID"},
    });

std::shared_ptr<ControlMethod>
ControlMethod::create(const dictionary &dict)
{
    const controlType type = controlTypeNames.get("mode", dict);

    switch (type)
    {
    case twoStep:
        return std::make_shared<TwoStepControl>(dict);
    case PID:
        return std::make_shared<PIDControl>(dict);
    default:
        FatalIOErrorInFunction(dict)
            << "    Unknown control method " << type
            << exit(FatalIOError);
        return nullptr;
    }
}

// * * * * * * * * * * * * Constructor  * * * * * * * * * * * * //
ControlMethod::ControlMethod(const dictionary &dict)
  : type_(controlTypeNames.get("mode", dict))
{}

void ControlMethod::write(Ostream &os) const
{
    os.writeEntry("mode", controlTypeNames.get(type_));
}

// * * * * * * * * * * * * Two Step Control  * * * * * * * * * * * * //
TwoStepControl::TwoStepControl(const dictionary &dict)
  : ControlMethod(dict),
    h_(parameters(dict).getOrDefault<scalar>("h", 0.)),
    outputSignal_(0.)
{}

scalar TwoStepControl::calculate(scalar current, scalar target, scalar deltaT)
{
    const scalar targetCorrected = outputSignal_ > 0. ? target + 0.5 * h_ : target - 0.5 * h_;
    const scalar errorCorrected = targetCorrected - current;
    outputSignal_ = errorCorrected <= 0 ? 0. : 1.;
    return outputSignal_;
}

void TwoStepControl::write(Ostream &os) const
{
    ControlMethod::write(os);
    os.beginBlock("parameters");
    os.writeEntryIfDifferent("h", 0.,  h_);
    os.endBlock();
}

// * * * * * * * * * * * * PID Control  * * * * * * * * * * * * //
PIDControl::PIDControl(const dictionary &dict)
  : ControlMethod(dict),
    Kp_(parameters(dict).getScalar("Kp")),
    Ti_(parameters(dict).getScalar("Ti")),
    Td_(parameters(dict).getScalar("Td")),
    outputMax_(parameters(dict).getOrDefault<scalar>("outputMax", 1.)),
    outputMin_(parameters(dict).getOrDefault<scalar>("outputMin", 0.)),
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
    const scalar outputSignal = Kp_*(error + 1/(Ti_ + ROOTSMALL)*errorIntegral_ + Td_*errorDifferential);

    // Return result within defined regulator saturation: outputMax_ and outputMin_
    return max(min(outputSignal, outputMax_), outputMin_);
}

void PIDControl::write(Ostream &os) const
{
    ControlMethod::write(os);
    os.beginBlock("parameters");
    os.writeEntry("Kp", Kp_);
    os.writeEntry("Ti", Ti_);
    os.writeEntry("Td", Td_);
    os.writeEntryIfDifferent("outputMax", 1., outputMax_);
    os.writeEntryIfDifferent("outputMin", 0., outputMin_);
    os.endBlock();
}
