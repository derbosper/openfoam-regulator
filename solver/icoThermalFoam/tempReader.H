// Class reads temperature from mesh in the runtime and prints average
// temperature of a patch to console.
#include "fvCFD.H"

scalar getPatchAverage(const volScalarField& field, const fvPatch& patch);

class tempReader
{    
    public:
        tempReader();
        ~tempReader();

        static void printValues(const volScalarField& field, const fvPatch& patch);
};
