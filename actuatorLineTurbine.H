/****************************************************************************\
This program is based on the openFOAM, and is developed by MaZhe.
The goal of this program is to build an actuatorLineTurbine class .
\****************************************************************************/

#ifndef actuatorLineTurbine_H
#define actuatorLineTurbine_H

#include "actuatorLineBlade.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "List.H"
#include "vector.H"
#include "tensor.H"
#include <vector>

/******************************************class declaration******************************************/

namespace Foam
{
namespace fv
{

class actuatorLineTurbine
{


public:

//Constructor
    actuatorLineTurbine
    (
        const Time & time,
        ALFBM::flagBit & f,
        turbineInfo & t,
        bladeInfo & b,
        std::vector<airfoilInfo> & a
    );

//- Destructor
    ~actuatorLineTurbine(){}

//turbine information and turbine working condition API
    const turbineInfo& turbineI() const {return turbineInfo_;}

    const scalar & turbineTorque() const {return torque_;}

    const vector & turbineThrust() const {return thrust_;}

    const scalar & turbinePower() const {return power_;}

//access used to read velocity from CFD result
    List<List<vector>> & bladeElementVelocity() {return bladeElementVelocity_;}

    List<vector> & towerElementVelocity() {return towerElementVelocity_;}

//actuator line model APIs
    const List<List<point>> & bladeElementPosition() const {return bladeElementPosition_;}

    const List<List<point>> & bladeElementPositionLast() const {return bladeElementPositionLast_;}

    const List<List<vector>> & bladeElementForce() const {return bladeElementForce_;}

    const List<point> & towerElementPosition() const {return towerElementPosition_;}

    const List<vector> & towerElementForce() const {return towerElementForce_;}

//calculate turbine working conditions of present time step
    void correct();

//initialize the blade element position
    void bladeElementInitial();

//transform the finite element result to actuator line model
    void finiteElementResultToAero();

//print out working condition of the wind turbine
    void workingConditionPrint();

//transform the velocities from global coordinate system to local coordinate system 
//store the velocity of local coordinate system in actuatorLineElement 
    void velocityUpdate();

//calculate forces in local coordinate system
    void aeroForceCalculation();

//transform the forces from local coordinate system to global system
//store the force of global coordinate system in actuatorLineTurbine
    void forceUpdate();

//apply the forces to finite element turbine model
    void aeroResultToFiniteElement();

//calculate the turbine deformation
    void turbineDeform();

//turbine rotate
    void turbineRotate();

//read & write results
    void writeTurbineResults(std::ofstream & turbinedata);
    void readTurbineResults(std::ifstream & turbinedata);
    void readControllerResults(std::ifstream & controllerdata);
    void writeControllerResults(std::ofstream & controllerdata);
    void readStructureResults(std::ifstream & structuredata);
    void writeStructureResults(std::ofstream & structuredata);
    void writeForceResults(std::ofstream & forcedata);
    void writeDeflectionResults(std::ofstream & deflectiondata);
    void writeModalResults(std::ofstream & modaldata);
    void readAirfoilResults(std::ifstream & airfoildata);
    void writeAirfoilResults(std::ofstream & airfoildata);

private:

//time
	const Time& time_;

//last time
    scalar lastTime_;

//flag bit
    ALFBM::flagBit & flagBit_;

//turbine static informaiton
    turbineInfo & turbineInfo_;

//vector of blades static informations
    bladeInfo & bladeInfo_;

//vector of airfoils static informations
    std::vector<airfoilInfo>& airfoilsInfo_;
	
//all actuator line blades
	std::vector<std::shared_ptr<actuatorLineBlade>> blades_;

//aero actuator line tower
    //actuatorLineTower tower_;

////controller
    controller controller_;

//finite element turbine model
    ALFBM::fETurbine fETurbine_;

//force and position in global coordinate system
//data for blades
    List<List<point>> bladeElementPosition_;
    List<List<point>> bladeElementPositionLast_;
    List<List<vector>> bladeElementForce_;
    List<List<vector>> bladeElementMoment_;
    List<List<vector>> bladeElementVelocity_;
//data for tower
    List<point> towerElementPosition_;
    List<vector> towerElementForce_;
    List<vector> towerElementVelocity_;

//wind turbine working conditions
    scalar torque_=0.0;
    vector thrust_=vector::zero;
    scalar power_=0.0;

//private member functions
    void initialParameters();
    void initialTrans();
    void correctTrans(scalar angle);
    void timeCorrect(scalar & variable, const scalar & startTime, const scalar & initialValue, const scalar & endTime, const scalar & finalValue);
    void readCSVLine(std::istringstream & data, scalar & s);
};

}//end namespace fv
}//end namespace Foam


/******************************************************************************************************************************\
|                                                                                                                              |
|                                                    function definition                                                       |
|                                                                                                                              |
\******************************************************************************************************************************/

/******************************************private member functions******************************************/

inline void Foam::fv::actuatorLineTurbine::initialParameters()
{
    List<vector> bladeElementTemp;
    bladeElementTemp.setSize(bladeInfo_.bladeAEP().size(),vector::zero);
    bladeElementPosition_.setSize(turbineInfo_.bladeNumber(),bladeElementTemp);
    bladeElementPositionLast_.setSize(turbineInfo_.bladeNumber(),bladeElementTemp);
    bladeElementForce_.setSize(turbineInfo_.bladeNumber(),bladeElementTemp);
    bladeElementMoment_.setSize(turbineInfo_.bladeNumber(),bladeElementTemp);
    bladeElementVelocity_.setSize(turbineInfo_.bladeNumber(),bladeElementTemp);
}

inline void Foam::fv::actuatorLineTurbine::readCSVLine(std::istringstream & line,scalar & s)
{
    std::string temp;
    getline(line,temp,',');
    std::istringstream d(temp);
    d>>s;
}

/******************************************public functions******************************************/

Foam::fv::actuatorLineTurbine::actuatorLineTurbine
(
    const Time& time,
    ALFBM::flagBit & f,
    turbineInfo & t,
    bladeInfo & b,
    std::vector<airfoilInfo>& a
):
    time_(time),
    lastTime_(time_.value()),
    flagBit_(f),
    turbineInfo_(t),
    bladeInfo_(b),
    airfoilsInfo_(a),
    controller_(time_,flagBit_,turbineInfo_),
    fETurbine_(time_,flagBit_,turbineInfo_,bladeInfo_,controller_)
{
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        blades_.push_back(std::make_shared<actuatorLineBlade>(time_,i,flagBit_,turbineInfo_,bladeInfo_,airfoilsInfo_));           
    }

    initialParameters();

    finiteElementResultToAero();

    bladeElementInitial();
}

void Foam::fv::actuatorLineTurbine::correct()
{    
    controller_.correct(torque_);
}

void Foam::fv::actuatorLineTurbine::bladeElementInitial()
{
    fETurbine_.turbineDeform(bladeElementPosition_,towerElementPosition_);
    bladeElementPositionLast_=bladeElementPosition_;
}

void Foam::fv::actuatorLineTurbine::finiteElementResultToAero()
{
    bladeElementPositionLast_=bladeElementPosition_;
    fETurbine_.turbineDeform(bladeElementPosition_,towerElementPosition_);
    if(flagBit_.preciseDeform())
    {
        for(auto bprobe=blades_.begin();bprobe!=blades_.end();bprobe++)
        {
            (*(*bprobe)).turbineBladeDeform(fETurbine_);
        }
    }
}

void Foam::fv::actuatorLineTurbine::workingConditionPrint()
{
    controller_.workingConditionPrint();
}

void Foam::fv::actuatorLineTurbine::velocityUpdate()
{
    int b=0;
    vector rS(0,0,controller_.rotateSpeed());
    for(auto bprobe=blades_.begin();bprobe!=blades_.end();bprobe++)
    {
        //rotate trans
        //global coordinate to local coordinate
        controller_.bladeCorrect(2*Foam::constant::mathematical::pi*(*(*bprobe)).bladeNumber()/turbineInfo_.bladeNumber());
        
        //debug
        if(flagBit_.debug01())
            controller_.printAll();
        
        forAll((*(*bprobe)).velocityAtElement(),i)
        {
            vector bladeElementRotateSpeed;
            bladeElementRotateSpeed=(controller_.yawTensor().T()&(controller_.rotorTensor().T()&rS))^(bladeElementPosition_[b][i]-turbineInfo_.shaftOrigin());
            (*(*bprobe)).velocityAtElement()[i]=controller_.preconeTensor()&(controller_.bladeTensor()&(controller_.rotorTensor()&(controller_.yawTensor()&(bladeElementVelocity_[b][i]-bladeElementRotateSpeed))));
        }
        (*(*bprobe)).velocityUpdate();
        b+=1;
    }
}

void Foam::fv::actuatorLineTurbine::aeroForceCalculation()
{
    for(auto bprobe=blades_.begin();bprobe!=blades_.end();bprobe++)
    {
        (*(*bprobe)).forceCalculation(controller_.pitchedAngle());
    }
    Info<<"All Rotor Force Calculated!"<<endl;
}

void Foam::fv::actuatorLineTurbine::forceUpdate()
{
    int bN=0;
    torque_=0.0;
    thrust_=vector::zero;
    vector axis=vector::zero;
    axis.z()=1;
    for(auto bprobe=blades_.begin();bprobe!=blades_.end();bprobe++)
    {
        controller_.bladeCorrect(2*Foam::constant::mathematical::pi*(*(*bprobe)).bladeNumber()/turbineInfo_.bladeNumber());
        tensor temp;
        temp=controller_.yawTensor().T()&controller_.rotorTensor().T()&controller_.bladeTensor().T()&controller_.preconeTensor().T();
        forAll((*(*bprobe)).elementPosition(),i)
        {
            bladeElementForce_[bN][i]=temp&(*(*bprobe)).forceAtElement()[i];
            bladeElementMoment_[bN][i]=temp&(*(*bprobe)).momentAtElement()[i];
            thrust_ += flagBit_.airDensity() * bladeElementForce_[bN][i];
            torque_ += flagBit_.airDensity() * ((bladeElementPosition_[bN][i]-turbineInfo_.shaftOrigin()) ^ (bladeElementForce_[bN][i]))&(controller_.yawTensor().T()&(controller_.rotorTensor().T()&axis));
        }

        bN+=1;
    }
    power_=torque_*controller_.rotateSpeed();
    Info<<"Force updated!"<<endl;
}

void Foam::fv::actuatorLineTurbine::aeroResultToFiniteElement()
{
    fETurbine_.loadCalculation(bladeElementForce_, bladeElementMoment_);
    Info<<"Structural Load Calculated!"<<endl;
}

void Foam::fv::actuatorLineTurbine::turbineDeform()
{
    if(flagBit_.deform())
    {
        if(flagBit_.dynamic())
        {
            if(flagBit_.nonlinear())
            {
                if(time_.value()==time_.deltaT().value())
                    fETurbine_.initialSolve();
                else
                    fETurbine_.implicitIterationSolve();
            }
            else
                fETurbine_.linearSolve();
        }
        else
        {
            if(flagBit_.nonlinear())
            {
                if(time_.value()==time_.deltaT().value())
                    fETurbine_.initialSolve();
                else
                    fETurbine_.staticInterationSolve();
            }
            else
                fETurbine_.linearSolve();
        }
    }
    else
    {
        fETurbine_.noDeformSolve();
    }

    Info<<"Turbine deformed!"<<endl;
}

void Foam::fv::actuatorLineTurbine::turbineRotate()
{
    if(time_.value()!=lastTime_)
    {
        controller_.timeCorrect();
        fETurbine_.turbineRotate();
        lastTime_ = time_.value();
        Info<<"Turbine rotated!"<<endl;
    }
}

void Foam::fv::actuatorLineTurbine::writeTurbineResults(std::ofstream & turbinedata)
{
    turbinedata<<"thrustX,thrustY,thrustZ,torque,power,"<<std::endl;
    turbinedata<<thrust_.x()<<","<<thrust_.y()<<","<<thrust_.z()<<","<<torque_<<","<<power_<<","<<std::endl;
}

void Foam::fv::actuatorLineTurbine::readTurbineResults(std::ifstream & turbinedata)
{
    std::string line;
    getline(turbinedata,line);
    getline(turbinedata,line);
    std::istringstream tdata(line);
    readCSVLine(tdata,thrust_.x());
    readCSVLine(tdata,thrust_.y());
    readCSVLine(tdata,thrust_.z());
    readCSVLine(tdata,torque_);
    readCSVLine(tdata,power_);
    Info<<"turbinedata of turbine "<<turbineInfo_.turbineName()<<" has been read. "<<endl;
    if(flagBit_.debug01())
    {
        Info<<"turbine data are : "<<endl;
        Info<< thrust_.x()<<","<<thrust_.y()<<","<<thrust_.z()<<","<<torque_<<","<<power_<<","<<endl;
    }
}

void Foam::fv::actuatorLineTurbine::readControllerResults(std::ifstream & controllerdata)
{
    controller_.readControllerResults(controllerdata);
}

void Foam::fv::actuatorLineTurbine::writeControllerResults(std::ofstream & controllerdata)
{
    controller_.writeControllerResults(controllerdata);
}

void Foam::fv::actuatorLineTurbine::readStructureResults(std::ifstream & structuredata)
{
    std::string line;
    getline(structuredata,line);
    fETurbine_.readStructureResults(structuredata);
    Info<<"structuredata of turbine "<<turbineInfo_.turbineName()<<" has been read. "<<endl;
    bladeElementInitial();
    fETurbine_.nodeReadInitial();
}

void Foam::fv::actuatorLineTurbine::writeStructureResults(std::ofstream & structuredata)
{
    structuredata<<"X,Y,Z,thetaX,thetaY,thetaZ,"<<std::endl;
    fETurbine_.writeStructureResults(structuredata);
}

void Foam::fv::actuatorLineTurbine::writeForceResults(std::ofstream & forcedata)
{
    forcedata<<"forceX,forceY,forceZ,momentX,momentY,momentZ,"<<std::endl;
    fETurbine_.writeForceResults(forcedata);
}

void Foam::fv::actuatorLineTurbine::writeDeflectionResults(std::ofstream & deflectiondata)
{
    deflectiondata<<"X,Y,Z,thetaX,thetaY,thetaZ,"<<std::endl;
    fETurbine_.writeDeflectionResults(deflectiondata);
}

void Foam::fv::actuatorLineTurbine::writeModalResults(std::ofstream & modaldata)
{
    fETurbine_.writeModalResults(modaldata);
}

void Foam::fv::actuatorLineTurbine::readAirfoilResults(std::ifstream & airfoildata)
{
    for(auto bprobe=blades_.begin();bprobe!=blades_.end();bprobe++)
    {
        (*(*bprobe)).readAirfoilResults(airfoildata);
    }
    Info<<"Airfoiledata of turbine "<<turbineInfo_.turbineName()<<" has been read. "<<endl;
}

void Foam::fv::actuatorLineTurbine::writeAirfoilResults(std::ofstream & airfoildata)
{
    for(auto bprobe=blades_.begin();bprobe!=blades_.end();bprobe++)
    {
        (*(*bprobe)).writeAirfoilResults(airfoildata);
    }
}

#endif