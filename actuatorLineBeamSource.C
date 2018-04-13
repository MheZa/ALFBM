/****************************************************************************\
This program is based on the openFOAM, and is developed by MaZhe.
The goal of this program is to build an ActuatorLine-Beam coupled source.
\****************************************************************************/

#include "actuatorLineBeamSource.H"
#include "addToRunTimeSelectionTable.H"
#include "geometricOneField.H"
#include "fvMatrices.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "simpleMatrix.H"

/******************************************function definition******************************************/

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorLineBeamSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        actuatorLineBeamSource,
        dictionary
    );
}
}

//constructor
Foam::fv::actuatorLineBeamSource::actuatorLineBeamSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh)
{
    read(dict_);
    for(auto tprobe=turbinesInfo_.begin();tprobe!=turbinesInfo_.end();tprobe++)
    {
        turbines_.push_back(std::make_shared<actuatorLineTurbine>(mesh.time(),flagBit_,(*tprobe),findBladeInfo((*tprobe).bladeName()),airfoilsInfo_));
    }
    readPreviousData();
    Info<<"The initialization of ALFBM succeed!"<<endl;
}

//- Destructor
Foam::fv::actuatorLineBeamSource::~actuatorLineBeamSource()
{}

inline bool Foam::fv::actuatorLineBeamSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);
        
        //read epsilon
        coeffs_.lookup("epsilon") >> epsilon_;

        //read flags
        flagBit_.read(coeffs_);

        //read informations of turbines to turbinesInfo_
        turbinesInfoRead();

        //read informations of blades to bladesInfo_
        bladesInfoRead();

        //read informations of airfoils to airfoilsInfo_
        airfoilsInfoRead();

        return true;
    }
    else
    {
        return false;
    }
}

inline void Foam::fv::actuatorLineBeamSource::turbinesInfoRead()
{
    dictionary turbineDict;
    turbineDict=coeffs_.subDict("turbine");
    forAll(turbineDict.toc(),i)
    {
        turbineInfo t(turbineDict.subDict(turbineDict.toc()[i]));
        //!!!useless copying here!!!
        turbinesInfo_.push_back(t);
    }
}

inline void Foam::fv::actuatorLineBeamSource::bladesInfoRead()
{
    dictionary bladeDict;
    bladeDict=coeffs_.subDict("blade");
    forAll(bladeDict.toc(),i)
    {        
        bladeInfo b(bladeDict.subDict(bladeDict.toc()[i]));
        //!!!useless copying here!!!
        bladesInfo_.push_back(b);
    }
}

inline void Foam::fv::actuatorLineBeamSource::airfoilsInfoRead()
{
    dictionary airfoilDict;
    airfoilDict=coeffs_.subDict("airfoil");
    forAll(airfoilDict.toc(),i)
    {
        airfoilInfo a(airfoilDict.subDict(airfoilDict.toc()[i]));
        //!!!useless copying here!!!
        airfoilsInfo_.push_back(a);
    }
}

inline Foam::fv::bladeInfo& Foam::fv::actuatorLineBeamSource::findBladeInfo(const word& bladeName)
{
    for(auto bprobe=bladesInfo_.begin();bprobe!=bladesInfo_.end();bprobe++)
    {
        if((*bprobe).bladeName()==bladeName)
        {
            return (*bprobe);
        }
    }
    Info<<"Error: error occur in function findBladeInfo when trying to find bladeinfo of "
        <<bladeName<<". Please check the bladeInfo folder."<<endl;
    return bladesInfo_[0];
}

void Foam::fv::actuatorLineBeamSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
	volVectorField force
    (
        IOobject(name_+":actuatorLineBeamSource", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedVector("zero", eqn.dimensions()/dimVolume, vector::zero)
    );

    for(auto tprobe=turbines_.begin();tprobe!=turbines_.end();tprobe++)
    {
        //correct working conditions
        (*(*tprobe)).correct();

        (*(*tprobe)).turbineRotate();

        (*(*tprobe)).finiteElementResultToAero();

        //display the working conditions of wind turbines
        (*(*tprobe)).workingConditionPrint();

        //read previous velocity for blades
        const volVectorField& Uin(eqn.psi());
        interpolationCellPoint<Foam::vector> UInterp(Uin);
        forAll((*(*tprobe)).bladeElementPosition(),j)
        {
            forAll((*(*tprobe)).bladeElementPosition()[j],k)
            {
                label temp;
                point bEP;
                    bEP=flagBit_.alphaVP()*(*(*tprobe)).bladeElementPosition()[j][k] + (1.0-flagBit_.alphaVP())*(*(*tprobe)).bladeElementPositionLast()[j][k];
                if(Pstream::parRun())
                {
                    temp=mesh_.findCell(bEP);
                    if(temp<=0)
                    {
                        (*(*tprobe)).bladeElementVelocity()[j][k]=VGREAT*vector::one;
                    }
                    else
                    {
                        (*(*tprobe)).bladeElementVelocity()[j][k]=UInterp.interpolate(bEP,temp);
                    }
                    reduce((*(*tprobe)).bladeElementVelocity()[j][k],minOp<vector>());
                }
                else
                {
                    temp=mesh_.findCell(bEP);
                    if(temp<=0)
                    {
                        Info<<"Error: Can not find position!"<<endl;
                    }
                    else
                    {
                        (*(*tprobe)).bladeElementVelocity()[j][k]=UInterp.interpolate(bEP,temp);
                    }
                }
            }
        }

        //read previous velocity for tower
        forAll((*(*tprobe)).towerElementPosition(),j)
        {
            label temp;
            if(Pstream::parRun())
            {
                temp=mesh_.findCell((*(*tprobe)).towerElementPosition()[j]);
                if(temp<=0)
                {
                    (*(*tprobe)).towerElementVelocity()[j]=VGREAT*vector::one;
                }
                else
                {
                    (*(*tprobe)).towerElementVelocity()[j]=UInterp.interpolate((*(*tprobe)).towerElementPosition()[j],temp);
                }
                reduce((*(*tprobe)).towerElementVelocity()[j],minOp<vector>());
            }
            else
            {
                temp=mesh_.findCell((*(*tprobe)).towerElementPosition()[j]);
                if(temp<=0)
                {
                    Info<<"Error: Can not find position!"<<endl;
                }
                else
                {
                    (*(*tprobe)).towerElementVelocity()[j]=UInterp.interpolate((*(*tprobe)).towerElementPosition()[j],temp);
                }
            }
        }

        (*(*tprobe)).velocityUpdate();

        (*(*tprobe)).aeroForceCalculation();

        (*(*tprobe)).forceUpdate();

        (*(*tprobe)).aeroResultToFiniteElement();

        (*(*tprobe)).turbineDeform();

    }

    //write results
    writeResult();
    
    //apply force to CFD
    forAll(cells_,i)
    {
        for(auto tprobe=turbines_.begin();tprobe!=turbines_.end();tprobe++)
        {
            forAll((*(*tprobe)).bladeElementPosition(),j)
            {
                forAll((*(*tprobe)).bladeElementPosition()[j],k)
                {
                    scalar dis = mag(mesh_.C()[cells_[i]] - (*(*tprobe)).bladeElementPosition()[j][k]);
                    if(dis<7*epsilon_)
                    {
                        scalar factor = Foam::exp(-Foam::sqr(dis/epsilon_))
                            / (Foam::pow(epsilon_, 3)
                            * Foam::pow(Foam::constant::mathematical::pi, 1.5));
                        force[cells_[i]] -= (*(*tprobe)).bladeElementForce()[j][k]*factor; 
                    }
                }
            }
            forAll((*(*tprobe)).towerElementPosition(),j)
            {
                scalar dis = mag(mesh_.C()[cells_[i]] - (*(*tprobe)).towerElementPosition()[j]);
                if(dis>=7*epsilon_)
                {
                    continue;
                }
                else
                {
                    scalar factor = Foam::exp(-Foam::sqr(dis/epsilon_))
                        / (Foam::pow(epsilon_, 3)
                        * Foam::pow(Foam::constant::mathematical::pi, 1.5));
                    force[cells_[i]] -= (*(*tprobe)).towerElementForce()[j]*factor; 
                }
            }
        }
    }
    // Add source to rhs of eqn
    eqn += force;
}

void Foam::fv::actuatorLineBeamSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    this->addSup(eqn, fieldI);
}

void Foam::fv::actuatorLineBeamSource::readPreviousData()
{
    //try to read previous results for Results file
    std::string dir;
    std::ifstream turbinedata;
    std::ifstream controllerdata;
    std::ifstream structuredata;
    std::ifstream airfoildata;

    if (Pstream::parRun())
    {
        dir = mesh_.time().path() / "../Results"
            / mesh_.time().timeName();
    }
    else
    {
        dir = mesh_.time().path() / "Results"
            / mesh_.time().timeName();
    }

    for(auto tprobe=turbines_.begin();tprobe!=turbines_.end();tprobe++)
    {
        turbinedata.open(dir+"/"+(*(*tprobe)).turbineI().turbineName()+"/turbineCondition.csv");
        controllerdata.open(dir+"/"+(*(*tprobe)).turbineI().turbineName()+"/controllerCondition.csv");
        structuredata.open(dir+"/"+(*(*tprobe)).turbineI().turbineName()+"/structureCondition.csv");
        airfoildata.open(dir+"/"+(*(*tprobe)).turbineI().turbineName()+"/airfoilCondition.csv");

        if(controllerdata && structuredata && airfoildata)
        {
            (*(*tprobe)).readTurbineResults(turbinedata);

            (*(*tprobe)).readControllerResults(controllerdata);

            (*(*tprobe)).readStructureResults(structuredata);

            (*(*tprobe)).readAirfoilResults(airfoildata);
            
            controllerdata.close();
            structuredata.close();
            airfoildata.close();

            Info<<"Successfully read previous results."<<endl; 
        }
        else
        {
            Info<<"No previous data found for wind turbine" << (*(*tprobe)).turbineI().turbineName() << "." <<endl;
            continue;
        }
    }
}

void Foam::fv::actuatorLineBeamSource::writeResult()
{
    //try to read previous results for Results file
    string dir;
    std::ofstream turbinedata;
    std::ofstream controllerdata;
    std::ofstream structuredata;
    std::ofstream airfoildata;
    std::ofstream forcedata;
    std::ofstream deflectiondata;
    std::ofstream modaldata;

    if (Pstream::parRun())
    {
        dir = mesh_.time().path() / "../Results"
            / mesh_.time().timeName();
    }
    else
    {
        dir = mesh_.time().path() / "Results"
            / mesh_.time().timeName();
    }

    for(auto tprobe=turbines_.begin();tprobe!=turbines_.end();tprobe++)
    {
        mkDir(dir+"/"+(*(*tprobe)).turbineI().turbineName());
        turbinedata.open(dir+"/"+(*(*tprobe)).turbineI().turbineName()+"/turbineCondition.csv");
        controllerdata.open(dir+"/"+(*(*tprobe)).turbineI().turbineName()+"/controllerCondition.csv");
        structuredata.open(dir+"/"+(*(*tprobe)).turbineI().turbineName()+"/structureCondition.csv");
        airfoildata.open(dir+"/"+(*(*tprobe)).turbineI().turbineName()+"/airfoilCondition.csv");
        modaldata.open(dir+"/"+(*(*tprobe)).turbineI().turbineName()+"/modalCondition.csv");
        forcedata.open(dir+"/"+(*(*tprobe)).turbineI().turbineName()+"/forceCondition.csv");
        deflectiondata.open(dir+"/"+(*(*tprobe)).turbineI().turbineName()+"/deflectionCondition.csv");

        if(turbinedata && controllerdata && structuredata && airfoildata && forcedata)
        {
            (*(*tprobe)).writeTurbineResults(turbinedata);

            (*(*tprobe)).writeControllerResults(controllerdata);

            (*(*tprobe)).writeStructureResults(structuredata);
            
            (*(*tprobe)).writeForceResults(forcedata);

            (*(*tprobe)).writeDeflectionResults(deflectiondata);
            
            (*(*tprobe)).writeModalResults(modaldata);

            (*(*tprobe)).writeAirfoilResults(airfoildata);
            
            turbinedata.close();
            controllerdata.close();
            structuredata.close();
            airfoildata.close();
            modaldata.close();
            forcedata.close();
            deflectiondata.close();
            Info<<"Previous data is writen for turbine "<<(*(*tprobe)).turbineI().turbineName()<<". "<<endl;
        }
        else
        {
            Info<< "Result write out error!"<<endl;
            continue;
        }
    }
}