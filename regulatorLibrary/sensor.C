#include "sensor.H"

// * * * * * * * * * * * * Utility Functions  * * * * * * * * * * * * //
scalar Sensor::patchAverage(const word &fieldName, const fvPatch &patch)
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
    patchName_(word::null)
{}

Sensor::Sensor(const fvMesh &mesh, const dictionary &dict):
    mesh_(mesh),
    fieldName_(dict.getWord("field")),
    type_(sensorTypeNames.get("type", dict)),
    patchName_(type_ == sensorType::patch ? dict.getWord("patchName") : word::null)
{}

Sensor::Sensor(const Sensor &rhs):
    mesh_(rhs.mesh_),
    fieldName_(rhs.fieldName_),
    type_(rhs.type_),
    patchName_(rhs.patchName_)
{}

scalar Sensor::read() const
{
    scalar result = 0;
    switch (type_)
    {
    case sensorType::patch:
    {
        const fvPatch &targetPatch = mesh_.boundary()[patchName_];
        result = patchAverage(fieldName_, targetPatch);
        break;
    }
    case sensorType::points:
    {
        result = 0;
        break;
    }
    default:
        result = 0;
    }
    return result;
}

word Sensor::fieldName() const
{
    return fieldName_;
}
