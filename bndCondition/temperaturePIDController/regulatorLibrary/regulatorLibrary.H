#ifndef Regulator_H
#define Regulator_H

#include "fvCFD.H"

class Regulator
{
private:
    const fvMesh& mesh_;
    dictionary regulatorDict_;

public:
    Regulator(const fvMesh &mesh);
    void read() const;
};

#endif
