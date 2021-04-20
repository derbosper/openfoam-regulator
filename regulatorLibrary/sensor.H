#ifndef Sensor_H
#define Sensor_H

#include "fvCFD.H"

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

    // Calculate face average of a value
    static scalar patchAverage(const word &fieldName, const fvPatch &patch);

private:
    const fvMesh &mesh_;   // Reference to the system mesh
    const word fieldName_; // Name of the field being read
    const sensorType type_;
    const word patchName_;  // For patch type
};

#endif
