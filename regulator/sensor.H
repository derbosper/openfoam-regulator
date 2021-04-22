/* 
sensor
{
    field       T;
    type        points;
    points
    (
        (0.02 0.005 0.0)
        (0.05 0.005 0.0)
        (0.07 0.005 0.0)
    );
}

or

sensor
{
    field       T;
    type        patch;
    patchName   outlet;
}

*/

#ifndef Sensor_H
#define Sensor_H

#include "fvCFD.H"
#include "pointField.H"

class Sensor
{
public:
    // Initilize from mesh and dictionary
    Sensor(const fvMesh &mesh, const dictionary &dict);

    // Initialize from mesh
    Sensor(const fvMesh &mesh);

    // Copy
    Sensor(const Sensor &rhs);

    enum sensorType
    {
        patch,  // reads from speficied patch
        points, // reads from specified points
    };
    static const Enum<sensorType> sensorTypeNames;

    // Read current field value
    scalar read() const;

    // Write to runtime dict
    void write(Ostream &os) const;

    word fieldName() const;

private:
    const fvMesh &mesh_;   // Reference to the system mesh
    const word fieldName_; // Controlled process variable scalar field name
    const sensorType type_;
    const word patchName_; // Name of a patch from which the process value is read
    const pointField points_;   // List of points for from which the average process value is read
};

// * * * * * * * * * * * * Helper Functions  * * * * * * * * * * * * //
scalar patchAverage(const word &fieldName, const fvPatch &patch);

#endif
