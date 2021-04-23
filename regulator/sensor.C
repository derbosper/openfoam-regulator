#include "sensor.H"

// * * * * * * * * * * * * Helper Functions  * * * * * * * * * * * * //
scalar patchAverage(const word &fieldName, const fvPatch &patch)
{
    const fvPatchField<scalar> &field =
        patch.lookupPatchField<volScalarField, scalar>(fieldName);

    return gSum(field * patch.magSf()) / gSum(patch.magSf());
}

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

const Foam::Enum<Sensor::sensorType>
    Sensor::sensorTypeNames({
        {sensorType::patch, "patch"},
        {sensorType::points, "points"},
    });

// * * * * * * * * * * * * Constructor  * * * * * * * * * * * * //

Sensor::Sensor(const fvMesh& mesh):
    mesh_(mesh),
    fieldName_(word::null),
    type_(sensorType::patch),
    patchName_(word::null),
    points_(pointField::null())
{}

Sensor::Sensor(const fvMesh &mesh, const dictionary &dict):
    mesh_(mesh),
    fieldName_(dict.getWord("field")),
    type_(sensorTypeNames.get("type", dict)),
    patchName_(type_ == sensorType::patch ? dict.getWord("patchName") : word::null),
    points_(type_ == sensorType::points ? dict.get<pointField>("points") : pointField::null())
{}

scalar Sensor::read() const
{
    switch (type_)
    {
    case sensorType::patch:
    {
        const fvPatch &targetPatch = mesh_.boundary()[patchName_];
        return patchAverage(fieldName_, targetPatch);
    }
    case sensorType::points:
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
    default:
    {
        FatalError << "Unknown sensor type" << endl;
        return 1;
    }
    }
}

word Sensor::fieldName() const
{
    return fieldName_;
}
