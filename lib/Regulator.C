#include "Regulator.H"
#include <stdio.h>

Regulator::Regulator(const fvMesh &mesh) : mesh_(mesh)
{
    std::cout << "Initialized!!!" << std::endl;
    /*     // Get access to a custom dictionary
    dictionary regulationDict;
    const word dictName("regulationProperties");

    // Create and input-output object - this holds the path to the dict and its name
    IOobject dictIO(
        dictName,               // name of the file
        mesh.time().constant(), // path to where the file is
        mesh,                   // reference to the mesh needed by the constructor
        IOobject::MUST_READ     // indicate that reading this dictionary is compulsory
    );

    // Check the if the dictionary is present and follows the OF format
    if (!dictIO.typeHeaderOk<dictionary>(true))
        FatalErrorIn(args.executable()) << "Cannot open specified refinement dictionary "
                                        << dictName << exit(FatalError);
    else
        Info << "Dictionary OK";

    // Initialise the dictionary object
    regulationDict = IOdictionary(dictIO);
    // Read various pieces of information from the main part of the dictionary

    // Lookup which does not need to be told what type of variable we're looking for and
    // uses the standard C++ stringstream syntax
    word someWord;
    regulationDict.lookup("someWord") >> someWord;
    Info << "someWord: " << someWord;
    // This template method needs to know the type of the variable and can provide
    // a default value if the entry is not found in the dictionary
    scalar someScalar(customDict.lookupOrDefault<scalar>("someScalar", 1.0));
    Info << "someScalar: " << someScalar; */
}