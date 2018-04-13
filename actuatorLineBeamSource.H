/****************************************************************************\
This program is based on the openFOAM, and is developed by MaZhe.
The goal of this program is to build an ActuatorLine-Beam coupled source.
\****************************************************************************/

#ifndef actuatorLineBeamSource_H
#define actuatorLineBeamSource_H

#include "actuatorLineTurbine.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "vector.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "List.H"
#include "cellSetOption.H"
#include "interpolationCellPoint.H"

/******************************************class declaration******************************************/

namespace Foam
{
namespace fv
{

class actuatorLineBeamSource
:
    public cellSetOption
{

private:

    //Epsilon
    scalar epsilon_;

    //flags
    ALFBM::flagBit flagBit_;

    //static parameters of wind turbine
    std::vector<turbineInfo> turbinesInfo_;

    std::vector<bladeInfo> bladesInfo_;

    std::vector<airfoilInfo> airfoilsInfo_;

	std::vector<std::shared_ptr<actuatorLineTurbine>> turbines_;
    
    //- Disallow default bitwise copy construct
    actuatorLineBeamSource(const actuatorLineBeamSource&);

    //- Disallow default bitwise assignment
    void operator=(const actuatorLineBeamSource&);

    //member functions for initialization
    bool read(const dictionary& dict);
    void turbinesInfoRead();
    void bladesInfoRead();
    void airfoilsInfoRead();
    Foam::fv::bladeInfo& findBladeInfo(const word& bladeName);

    //member functions for result output
    void writeResult();

    //member functions for previous result read
    void readPreviousData();


public:

	TypeName("actuatorLineBeamSource");

    // Selectors

    //- Return a reference to the selected fvOption model
    static autoPtr<actuatorLineBeamSource> New
    (
        const word& name,
        const dictionary& dict,
        const fvMesh& mesh
    );

	//Constructor
	actuatorLineBeamSource
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~actuatorLineBeamSource();

    //- Source term to momentum equation
    virtual void addSup
    (
        fvMatrix<vector>& eqn,
        const label fieldI
    );

    //- Source term to compressable scalars
    virtual void addSup
    (
        const volScalarField& rho,
        fvMatrix<vector>& eqn,
        const label fieldI
    );

};

} // End namespace fv
} // End namespace Foam

#endif

