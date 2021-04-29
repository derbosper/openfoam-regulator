#include "sensor.H"

// * * * * * * * * * * * * Factory  * * * * * * * * * * * * //
const Foam::Enum<Sensor::sensorType>
    Sensor::sensorTypeNames({
        {sensorType::patch, "patch"},
        {sensorType::points, "points"},
        {sensorType::volume, "volume"},
    });

std::shared_ptr<Sensor>
Sensor::create(const fvMesh& mesh, const dictionary& dict)
{
    const sensorType type = sensorTypeNames.get("type", dict);

    switch (type)
    {
    case patch:
        return std::make_shared<PatchSensor>(mesh, dict);
    case points:
        return std::make_shared<PointSensor>(mesh, dict);
    case volume:
        return std::make_shared<VolumeSensor>(mesh, dict);
    default:
        FatalIOErrorInFunction(dict)
            << "    Unknown Sensor type " << type
            << exit(FatalIOError);
        return nullptr;
    }
}


// * * * * * * * * * * * * Base Sensor  * * * * * * * * * * * * //
Sensor::Sensor(const fvMesh &mesh, const dictionary &dict):
    mesh_(mesh),
    fieldName_(dict.getWord("field")),
    type_(sensorTypeNames.get("type", dict))
{}

word Sensor::fieldName() const
{
    return fieldName_;
}

void Sensor::write(Ostream &os) const
{
    os.writeEntry("field", fieldName_);
    os.writeEntry("type", sensorTypeNames.get(type_));
}

// * * * * * * * * * * * * PointSensor  * * * * * * * * * * * * //
PointSensor::PointSensor(const fvMesh &mesh, const dictionary &dict):
    Sensor(mesh, dict),
    points_(dict.get<pointField>("points"))
{}

scalar PointSensor::read() const
{
    const volScalarField& field = mesh_.lookupObject<volScalarField>(fieldName_);
    scalar fieldSum = 0.0;
    forAll(points_, pointi)
    {
        const vector& location = points_[pointi];
        const label celli = mesh_.findCell(location);
        fieldSum += field[celli];
    }
    // Sync across all processes
    reduce(fieldSum, sumOp<scalar>());
    return fieldSum / points_.size();
}

void PointSensor::write(Ostream &os) const
{
    Sensor::write(os);
    os.writeEntry("points", points_);
}

// * * * * * * * * * * * * PatchSensor  * * * * * * * * * * * * //
PatchSensor::PatchSensor(const fvMesh &mesh, const dictionary &dict):
    Sensor(mesh, dict),
    patchName_(dict.getWord("patchName"))
{}

scalar PatchSensor::read() const
{
    const fvPatch &targetPatch = mesh_.boundary()[patchName_];
    return patchAverage(fieldName_, targetPatch);
}

void PatchSensor::write(Ostream &os) const
{
    Sensor::write(os);
    os.writeEntry("patchName", patchName_);
}

// * * * * * * * * * * * * VolumeSensor  * * * * * * * * * * * * //
VolumeSensor::VolumeSensor(const fvMesh &mesh, const dictionary &dict):
    Sensor(mesh, dict)
{}

scalar VolumeSensor::read() const
{
    const volScalarField &field = mesh_.lookupObject<volScalarField>(fieldName_);
    return field.weightedAverage(mesh_.V()).value();
}

void VolumeSensor::write(Ostream &os) const
{
    Sensor::write(os);
}

// * * * * * * * * * * * * Helper Functions  * * * * * * * * * * * * //
scalar patchAverage(const word &fieldName, const fvPatch &patch)
{
    const fvPatchField<scalar> &field =
        patch.lookupPatchField<volScalarField, scalar>(fieldName);

    return gSum(field * patch.magSf()) / gSum(patch.magSf());
}
