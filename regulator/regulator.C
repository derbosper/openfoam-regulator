#include "regulator.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Regulator::Regulator(const fvMesh &mesh, const dictionary &dict)
    : mesh_(mesh),
      sensor_(Sensor::create(mesh, dict)),
      controlMethod_(ControlMethod::create(dict)),
      targetValue_(Function1<scalar>::New("targetValue", dict)),
      timeIndex_(mesh.time().timeIndex())
{}

Regulator::Regulator(const fvMesh &mesh)
    : mesh_(mesh),
    sensor_(nullptr),
    controlMethod_(nullptr),
    targetValue_(nullptr),
    timeIndex_(mesh.time().timeIndex())
{}

Regulator::Regulator(const Regulator& reg)
  : mesh_(reg.mesh_),
    sensor_(reg.sensor_),
    controlMethod_(reg.controlMethod_),
    targetValue_(reg.targetValue_.clone()),
    timeIndex_(reg.timeIndex_)
{}

// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * *//

scalar Regulator::read()
{
    // Get time data
    const scalar deltaT = mesh_.time().deltaTValue();
    const scalar t = mesh_.time().timeOutputValue();

    // Update the old-time quantities
    if (timeIndex_ != mesh_.time().timeIndex())
    {
        timeIndex_ = mesh_.time().timeIndex();
    }

    // Get the target patch average field value
    const scalar sensorValue = sensor_->read();

    const scalar outputSignal = controlMethod_->calculate(sensorValue, targetValue_->value(t), deltaT);

    Info << "Regulator: targetValue = " << targetValue_->value(t) << endl;
    Info << "Regulator: sensorValue = " << sensorValue << endl;
    Info << "Regulator: error = " << targetValue_->value(t) - sensorValue  << endl;
    Info << "Regulator: outputSignal = " << outputSignal << endl;

    return outputSignal;
}

void Regulator::write(Ostream& os, const word dictName) const
{
    os.beginBlock(dictName);
    targetValue_->writeData(os);
    os.writeEntry("field", sensor_->fieldName());

    controlMethod_->write(os);

    os.beginBlock("sensor");
    sensor_->write(os);
    os.endBlock();

    os.endBlock();
}
