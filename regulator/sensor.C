#include "sensor.H"

// * * * * * * * * * * * * Helper Functions  * * * * * * * * * * * * //
scalar patchAverage(const word &fieldName, const fvPatch &patch)
{
    const fvPatchField<scalar> &field =
        patch.lookupPatchField<volScalarField, scalar>(fieldName);

    return gSum(field * patch.magSf()) / gSum(patch.magSf());
}

// * * * * * * * * * * * * Constructor  * * * * * * * * * * * * //
Sensor::Sensor(const fvMesh &mesh, const dictionary &dict):
    mesh_(mesh),
    fieldName_(dict.getWord("field"))
{}

// * * * * * * * * * * * * Public member functions  * * * * * * * * * * * * //
word Sensor::fieldName() const
{
    return fieldName_;
}

Sensor* Sensor::create(const fvMesh& mesh, const dictionary& dict)
{
    const sensorType type = sensorTypeNames.get("type", dict);

    switch (type)
    {
    case patch:
        return new PatchSensor(mesh, dict);
    case points:
        return new PointSensor(mesh, dict);
    default:
        return NULL;
    }
}

const Foam::Enum<Sensor::sensorType>
    Sensor::sensorTypeNames({
        {sensorType::patch, "patch"},
        {sensorType::points, "points"},
    });

// * * * * * * * * * * * * Point sensor  * * * * * * * * * * * * //
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
    reduce(fieldSum, sumOp<scalar>());
    return fieldSum / points_.size();
}

// * * * * * * * * * * * * Patch sensor  * * * * * * * * * * * * //
PatchSensor::PatchSensor(const fvMesh &mesh, const dictionary &dict):
    Sensor(mesh, dict),
    patchName_(dict.getWord("patchName"))
{}

scalar PatchSensor::read() const
{
    const fvPatch &targetPatch = mesh_.boundary()[patchName_];
    return patchAverage(fieldName_, targetPatch);
}
