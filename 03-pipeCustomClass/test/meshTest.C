#include "fvCFD.H"

scalar getPatchAverage(volScalarField field, const fvPatch& patch);

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    
    #include "createTime.H"
    #include "createMesh.H"

    // get outlet from mesh database
    const fvPatch& outletPatch = mesh.boundary()["outlet"];

    // load temperature
    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    scalar avgTemp = getPatchAverage(T, outletPatch);
    Info << avgTemp << endl;
    Info<< "Time = " << runTime.timeName() << endl;

    return 0;
}

scalar getPatchAverage(volScalarField field, const fvPatch& patch)
{
    label cellId;
    scalar total(0.);
    for (label cellI = 0; cellI < patch.size(); cellI++)
    {
        cellId = patch.faceCells()[cellI];
        total += field[cellId];
    }
    return total / patch.size();
}