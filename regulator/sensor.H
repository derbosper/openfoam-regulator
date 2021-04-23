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
    virtual ~Sensor() = default;
    Sensor() = delete;

    // Initilize from mesh and dictionary
    Sensor(const fvMesh &mesh, const dictionary &dict);

    // Read current field value
    virtual scalar read() const = 0;

    // Write to runtime dict
    void write(Ostream &os) const;

    word fieldName() const;

    // Returns sensor-type-dependent implementation of Sensor
    static Sensor* create(const fvMesh& mesh, const dictionary& dict);
    
    enum sensorType
    {
        patch,  // reads from speficied patch
        points, // reads from specified points
    };

    static const Enum<sensorType> sensorTypeNames;

protected:
    const fvMesh &mesh_;   // Reference to the system mesh
    const word fieldName_; // Controlled process variable scalar field name

};

class PointSensor : public Sensor
{
public:
    PointSensor() = delete;

    PointSensor(const fvMesh &mesh, const dictionary &dict);

    scalar read() const override;

    void write(Ostream &os) const;

private:
    const pointField points_;   // List of points for from which the average process value is read
};

class PatchSensor : public Sensor
{
public:
    PatchSensor() = delete;

    PatchSensor(const fvMesh &mesh, const dictionary &dict);

    scalar read() const override;

    void write(Ostream &os) const;

private:
    // Name of a patch from which the process value is read
    const word patchName_;
};

// * * * * * * * * * * * * Helper Functions  * * * * * * * * * * * * //
scalar patchAverage(const word &fieldName, const fvPatch &patch);

#endif
