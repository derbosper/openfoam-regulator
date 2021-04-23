#ifndef ControlMethod_H
#define ControlMethod_H

#include "fvCFD.H"

class ControlMethod
{
public:
    virtual ~ControlMethod() = default;
    ControlMethod() = default;

    static ControlMethod* create(const dictionary &);

    enum controlType
    {
        twoStep, // on/off control
        PID,
    };
    static const Enum<controlType> controlTypeNames;

    // Calculate output signal from current process value, target value
    // and current time delta
    virtual scalar calculate(scalar current, scalar target, scalar deltaT) = 0;
};

class TwoStepControl : public ControlMethod
{
public:
    TwoStepControl() = default;
    TwoStepControl(const dictionary &);

    scalar calculate(scalar current, scalar target, scalar deltaT);

private:
    // Wartość histerezy
    //- Przy istnieniu histerezy wprowadzana jest zmienna zastępczna y'_z taka,
    //- że y'_z = y + 0.5h, dla sygnału wyjściowego 1 albo y - 0.5h dla sygnału 0
    scalar h_;

    // Preserve old output signal for histeresis correction
    scalar outputSignal_;

};

class PIDControl : public ControlMethod
{
public:
    PIDControl() = default;
    PIDControl(const dictionary &);

    scalar calculate(scalar current, scalar target, scalar deltaT);

private:
    scalar Kp_;                     // Proportional control term
    scalar Ti_;                     // Integral control term
    scalar Td_;                     // Differential control term
    scalar outputMax_;              // Max output signal (saturation)
    scalar outputMin_;              // Min output signal
    scalar oldError_;               // Old error
    scalar errorIntegral_;          // Error integral w.r.t. time
};

#endif
