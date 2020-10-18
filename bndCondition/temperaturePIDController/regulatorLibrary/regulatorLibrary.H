#ifndef Regulator_H
#define Regulator_H

#include "fvCFD.H"

class Regulator
{
public:
    Regulator(const fvMesh &mesh);
    const fvMesh& mesh_;
};
#endif