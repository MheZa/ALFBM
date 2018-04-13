/****************************************************************************\
This program is based on the openFOAM, and is developed by MaZhe.
The goal of this program is to build an actuatorLineBlade class .
\****************************************************************************/

//This is the local coordinate system position of blade
//Only represent the pre-deformation and deformation of the blade
//The local coordinate is:
//rotate axis for Z direction,also the flap direction
//main direction of blade is X
//Y axis is also the edge direction 

#ifndef actuatorLineBlade_H
#define actuatorLineBlade_H

#include "fvMesh.H"
#include "List.H"
#include "vector.H"
#include "point.H"
#include "syncTools.H"
#include "bladeInfo.H"
#include "turbineInfo.H"
#include "actuatorLineElement.H"

/******************************************class declaration******************************************/

namespace Foam
{
namespace fv
{

class actuatorLineBlade
{

public:

//Constructor
    actuatorLineBlade
    (
        const Time& time,
        int bN,
        ALFBM::flagBit & f,
        turbineInfo& tI,
        bladeInfo& bI,
        std::vector<airfoilInfo>& aI
    );

//- Destructor
    ~actuatorLineBlade(){}

//APIs
    const bladeInfo & bladeI() const {return bladeInfo_;}

    const List<point> & elementPosition() const {return elementPosition_;}

    const List<vector> & forceAtElement() const {return forceAtElement_;}

    const List<vector> & momentAtElement() const {return momentAtElement_;}

    const int & bladeNumber() const {return bladeNumber_;}

//velocity update
    void velocityUpdate();

//force calculation
    void forceCalculation(scalar pitchedAngle);

//access to velocityAtElement_
    List<vector> & velocityAtElement() {return velocityAtElement_;}

//access to actuatorLineElements
    //std::vector<actuatorLineElement> & elementsAccess() {return elements_;}
//blade deform
    void turbineBladeDeform(ALFBM::fETurbine & fETurbine_);
//read results
    void readAirfoilResults(std::ifstream & airfoildata);
//write results
    void writeAirfoilResults(std::ofstream & airfoildata);

private:

//time	
    const Time& time_;

//blade number
    int bladeNumber_;

//flag bit
    ALFBM::flagBit & flagBit_;

//reference of rotorInfo
    turbineInfo& turbineInfo_;

//reference of bladeInfo
	bladeInfo& bladeInfo_;

//airfoilsInfo
    std::vector<airfoilInfo>& airfoilsInfo_;

//actuatorLineElements
    std::vector<std::shared_ptr<actuatorLineElement>> elements_;

//aero element position
    List<point> elementPosition_;

//velocity in blade local coordinate system
	List<vector> velocityAtElement_;

//aero force 
	List<vector> forceAtElement_;// force at actuator line element position

//aero moment
    List<vector> momentAtElement_;//moment at actuator line element position

//private member funtions
    airfoilInfo& findAirfoilInfo(word airfoilName);

    void aeroDataCorrect();

};

}//end namespace fv
}//end namespace Foam

/******************************************function definition******************************************/

inline Foam::fv::airfoilInfo& Foam::fv::actuatorLineBlade::findAirfoilInfo(word airfoilName)
{
    for(auto aprobe=airfoilsInfo_.begin();aprobe!=airfoilsInfo_.end();aprobe++)
    {
        if((*aprobe).airfoilName()==airfoilName)
        {
            return (*aprobe);
        }
    }
    Info<<"Error: error occur in function findAirfoilInfo when trying to find airfoilinfo of "
        <<airfoilName<<". Please check the folder of airfoilInfo."<<endl;
    return airfoilsInfo_[0];
}

inline void Foam::fv::actuatorLineBlade::aeroDataCorrect()
{
    //correct data form elements
    int i=0;
    for(auto eprobe=elements_.begin();eprobe!=elements_.end();eprobe++)
    {
        forceAtElement_[i]=(*(*eprobe)).forceTN();
        momentAtElement_[i]= (*(*eprobe)).forceM();
        i+=1;
    }
}

Foam::fv::actuatorLineBlade::actuatorLineBlade
(
    const Time& time,
    int bN,
    ALFBM::flagBit & f,
    turbineInfo& tI,
    bladeInfo& bI,
    std::vector<airfoilInfo>& aI
):
    time_(time),
    bladeNumber_(bN),
    flagBit_(f),
    turbineInfo_(tI),
    bladeInfo_(bI),
    airfoilsInfo_(aI)
{
    forAll(bladeInfo_.bladeAEP(),i)
    {
        elements_.push_back(std::make_shared<actuatorLineElement>(time_,i,flagBit_,turbineInfo_,bladeInfo_,findAirfoilInfo(bladeInfo_.airfoilName()[i])));
    }
    elementPosition_=bladeInfo_.bladeAEP();
    velocityAtElement_.setSize(elementPosition_.size(),vector::zero);
    forceAtElement_.setSize(elementPosition_.size(),vector::zero);
    momentAtElement_.setSize(elementPosition_.size(),vector::zero);

}

void Foam::fv::actuatorLineBlade::velocityUpdate()
{
    int i=0;
    for(auto eprobe=elements_.begin();eprobe!=elements_.end();eprobe++)
    {
        (*(*eprobe)).velocityUpdate(velocityAtElement_[i]);
        i+=1;
    }
}

//calculate blade element force according to the velocityAtElement_
void Foam::fv::actuatorLineBlade::forceCalculation(scalar pitchedAngle)
{
    for(auto eprobe=elements_.begin();eprobe!=elements_.end();eprobe++)
    {
        (*(*eprobe)).forceCalculate(pitchedAngle);
    }
    aeroDataCorrect();
}

void Foam::fv::actuatorLineBlade::turbineBladeDeform(ALFBM::fETurbine & fETurbine_)
{
    for(auto eprobe=elements_.begin();eprobe!=elements_.end();eprobe++)
    {
        (*(*eprobe)).turbineElementDeform(fETurbine_,bladeNumber_);
    }
}

void Foam::fv::actuatorLineBlade::readAirfoilResults(std::ifstream & airfoildata)
{
    for(auto eprobe=elements_.begin();eprobe!=elements_.end();eprobe++)
    {
        (*(*eprobe)).readResults(airfoildata);
    }
}

void Foam::fv::actuatorLineBlade::writeAirfoilResults(std::ofstream & airfoildata)
{
    for(auto eprobe=elements_.begin();eprobe!=elements_.end();eprobe++)
    {
        (*(*eprobe)).writeResults(airfoildata);
    }
    aeroDataCorrect();
}

#endif