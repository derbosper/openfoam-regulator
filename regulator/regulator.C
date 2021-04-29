#include "regulator.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Regulator::Regulator(const fvMesh &mesh, const dictionary &dict)
    : mesh_(mesh),
      sensor_(Sensor::create(mesh, dict.subDict("sensor"))),
      controlMethod_(ControlMethod::create(dict)),
      targetValue_(dict.getScalar("targetValue")),
      timeIndex_(mesh.time().timeIndex())
{}

Regulator::Regulator(const fvMesh &mesh)
    : mesh_(mesh),
    sensor_(nullptr),
    controlMethod_(nullptr),
    targetValue_(0),
    timeIndex_(mesh.time().timeIndex())
{}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * *//

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
    const scalar sensorValue = sensor_->read();

    const scalar outputSignal = controlMethod_->calculate(sensorValue, targetValue_, deltaT);

    Info << "Regulator: targetValue = " << targetValue_ << endl;
    Info << "Regulator: sensorValue = " << sensorValue << endl;
    Info << "Regulator: error = " << targetValue_ - sensorValue  << endl;
    Info << "Regulator: outputSignal = " << outputSignal << endl;

    return outputSignal;
}

void Regulator::write(Ostream& os, const word dictName) const
{
    os.beginBlock(dictName);
    os.writeEntry("targetValue", targetValue_);
    controlMethod_->write(os);

    os.beginBlock("sensor");
    sensor_->write(os);
    os.endBlock();

    os.endBlock();
}
