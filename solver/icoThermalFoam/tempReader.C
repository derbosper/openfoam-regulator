
#include "tempReader.H"

tempReader::tempReader()
{}

tempReader::~tempReader()
{}

void tempReader::printValues(const volScalarField &field, const fvPatch& patch)
{
    scalar avgTemp = getPatchAverage(field, patch);
    Info << "Temperature read: " << avgTemp << endl;
}

scalar getPatchAverage(const volScalarField &field, const fvPatch &patch)
{
    label cellId;
    scalar total(0.);
    for (label cellI = 0; cellI < patch.size(); cellI++)
    {
        cellId = patch.faceCells()[cellI];
        total += field[cellId] / patch.size();
    }
    return total;
}