#include "regulatorLibrary.H"

Regulator::Regulator(const fvMesh &mesh):
    mesh_(mesh),
    error_(0.), // TODO change initial value to lookupOrDefault("error", 0.)
    errorIntegral_(0.),  // TODO as above
    oldError_(0.),
    timeIndex_(mesh.time().timeIndex())
{
    // Get access to a custom dictionary
    const word dictName("regulatorProperties");

    // Create and input-output object - this holds the path to the dict and its name
    IOobject dictIO(
        dictName,               // name of the file
        mesh.time().constant(), // path to where the file is
        mesh,                   // reference to the mesh needed by the constructor
        IOobject::MUST_READ     // indicate that reading this dictionary is compulsory
    );

    // Check the if the dictionary is present and follows the OF format
    if (!dictIO.typeHeaderOk<dictionary>(true))
        FatalErrorIn("regulatorLibrary.C") << "Cannot open specified refinement dictionary "
                                           << dictName << exit(FatalError);
    else
        Info << "Dictionary OK" << endl;

    // Initialise the dictionary object
    regulatorDict_ = IOdictionary(dictIO);

    // Read various pieces of information from the main part of the dictionary
    fieldName_ = regulatorDict_.getWord("fieldName");
    patchName_ = regulatorDict_.getWord("patchName");
    targetValue_ = regulatorDict_.getScalar("targetValue");
    P_ = regulatorDict_.getScalar("P");
    I_ = regulatorDict_.getScalar("I");
    D_ = regulatorDict_.getScalar("D");
}

scalar Regulator::read()
{       
    Info << "Time index: " << mesh_.time().timeIndex() << endl;
    Info << "fieldName: " << fieldName_ << endl;
    Info << "patchName: " << patchName_ << endl;
    Info << "targetValue: " << targetValue_ << endl;
    Info << "P: " << P_ << endl;
    Info << "I: " << I_ << endl;
    Info << "D: " << D_ << endl;

    // Get the time step
    const scalar deltaT(mesh_.time().deltaTValue());

    // Update the old-time quantities
    if (timeIndex_ != mesh_.time().timeIndex())
    {
        timeIndex_ = mesh_.time().timeIndex();
        oldError_ = error_;
    }

    return 1.0;
}
