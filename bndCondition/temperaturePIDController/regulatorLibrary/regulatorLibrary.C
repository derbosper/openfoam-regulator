#include "regulatorLibrary.H"

Regulator::Regulator(const fvMesh &mesh) : mesh_(mesh)
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
        Info << "Dictionary OK";

    // Initialise the dictionary object
    regulatorDict_ = IOdictionary(dictIO);
}

void Regulator::read() const
{
    // Read various pieces of information from the main part of the dictionary

    // Lookup which does not need to be told what type of variable we're looking for and
    // uses the standard C++ stringstream syntax
    word someWord;
    regulatorDict_.lookup("someWord") >> someWord;
    Info << "someWord: " << someWord;
    // This template method needs to know the type of the variable and can provide
    // a default value if the entry is not found in the dictionary
    scalar someScalar(regulatorDict_.lookupOrDefault<scalar>("someScalar", 1.0));
    Info << "someScalar: " << someScalar;
}
