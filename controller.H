/****************************************************************************\
This program is based on the openFOAM, and is developed by MaZhe.
The goal of this program is to build an controller class .
\****************************************************************************/

#ifndef controller_H
#define controller_H

#include "scalar.H"
#include "tensor.H"
#include "dictionary.H"
#include <Eigen/Dense>

namespace Foam
{
namespace fv
{

class controller
{

public:

//Constructor
    controller
    (
        const Time & time,
        ALFBM::flagBit & f,
        turbineInfo & tI
    ):
        time_(time),
        flagBit_(f),
        turbineInfo_(tI)
    {
        controllerInitial();
        paramatersInitial();
    }

//APIs
    const scalar & rotateSpeed() const {return rotateSpeed_;}

    const scalar & rotatedAngle() const {return rotatedAngle_;}

    const scalar & yawedAngle() const {return yawedAngle_;}

    const scalar & pitchedAngle() const {return pitchedAngle_;}

    const tensor & rotateTensor() const {return rotateTensor_;}

    const tensor & bladeTensor() const {return bladeTensor_;}

    const tensor & yawTensor() const {return yawTensor_;}

    const tensor & rotorTensor() const {return rotorTensor_;}

    const tensor & preconeTensor() const {return preconeTensor_;}

	const Eigen::Matrix<double,3,3> & rotateMatrix() const {return rotateMatrix_;}

	const Eigen::Matrix<double,3,3> & bladeMatrix() const {return bladeMatrix_;}

	const Eigen::Matrix<double,3,3> & yawMatrix() const {return yawMatrix_;}

	const Eigen::Matrix<double,3,3> & rotorMatrix() const {return rotorMatrix_;}

	const Eigen::Matrix<double,3,3> & preconeMatrix() const {return preconeMatrix_;}

//update all outputs
    void correct(const scalar & torque);
    void bladeCorrect(scalar angle);
    void timeCorrect();

//read & write results
    void readControllerResults(std::ifstream & controldata);
    void writeControllerResults(std::ofstream & controldata);

//print out working condition of the wind turbine
    void workingConditionPrint();

//print out all the mid terms, used for debug
    void printAll();

private:

//time
    const Time & time_;

//flagbit
    ALFBM::flagBit & flagBit_;

//turbine information
    turbineInfo & turbineInfo_;

//controller paramaters
    scalar rotateAccelerate_;

    scalar deltaPitch_;

    scalar rotateAccI_;

//controller scalar output
    scalar rotateSpeed_;

    scalar rotatedAngle_;

    scalar yawedAngle_;

    scalar pitchedAngle_;
//controller tensor output
    tensor rotateTensor_;

    tensor yawTensor_;

    tensor rotorTensor_;

    tensor preconeTensor_;

//controller matrix output
    Eigen::Matrix<double,3,3> rotateMatrix_;

    Eigen::Matrix<double,3,3> yawMatrix_;

    Eigen::Matrix<double,3,3> rotorMatrix_;

    Eigen::Matrix<double,3,3> preconeMatrix_;

//these tensor and matrix are only used for blades, and can only be used after bladeCorrect
    tensor bladeTensor_;

    Eigen::Matrix<double,3,3> bladeMatrix_;



//private member functions
    void controllerInitial();
    void paramatersInitial();
    void interp(scalar & variable, const scalar & startTime, const scalar & initialValue, const scalar & endTime, const scalar & finalValue);
    void scalarCorrect(const scalar & torque);
    void tensorCorrect();
    void readCSVLine(std::istringstream & line,scalar & s);

};

}//end namespace fv
}//end namespace Foam


inline void Foam::fv::controller::controllerInitial()
{
    rotorTensor_=tensor::I;
    rotorTensor_.yy()=cos(-turbineInfo_.shaftTiltAngle());
    rotorTensor_.yz()=sin(-turbineInfo_.shaftTiltAngle());
    rotorTensor_.zy()=-rotorTensor_.yz();
    rotorTensor_.zz()=rotorTensor_.yy();

    preconeTensor_=tensor::I;
    preconeTensor_.xx()=cos(-turbineInfo_.rotorPrecone());
    preconeTensor_.xz()=-sin(-turbineInfo_.rotorPrecone());
    preconeTensor_.zx()=-preconeTensor_.xz();
    preconeTensor_.zz()=preconeTensor_.xx();

    rotorMatrix_<<1,0,0,
                0, cos(-turbineInfo_.shaftTiltAngle()), sin(-turbineInfo_.shaftTiltAngle()),
                0,-sin(-turbineInfo_.shaftTiltAngle()), cos(-turbineInfo_.shaftTiltAngle());

    preconeMatrix_<<cos(-turbineInfo_.rotorPrecone()), 0, -sin(-turbineInfo_.rotorPrecone()),
                0,1,0,
                sin(-turbineInfo_.rotorPrecone()), 0, cos(-turbineInfo_.rotorPrecone());
}

inline void Foam::fv::controller::paramatersInitial()
{
    if(turbineInfo_.autoControl())
    {
    	rotateSpeed_=turbineInfo_.rotateInitialSpeed();
    	rotatedAngle_=0.0;
    	yawedAngle_=turbineInfo_.yawedInitialAngle();
    	pitchedAngle_=turbineInfo_.pitchedInitialAngle();
    	rotateAccelerate_=0.0;
    	rotateAccI_=0.0;
    }
    else
    {
    	interp(rotateSpeed_,
            turbineInfo_.rotateSpeedChangeStartTime(),turbineInfo_.rotateInitialSpeed(),turbineInfo_.rotateFinalSpeed(),turbineInfo_.rotateFinalSpeed());
        interp(yawedAngle_,
            turbineInfo_.yawedAngleChangeStartTime(),turbineInfo_.yawedInitialAngle(),turbineInfo_.yawedAngleChangeEndTime(),turbineInfo_.yawedFinalAngle());
        interp(pitchedAngle_,
            turbineInfo_.pitchedAngleChangeStartTime(),turbineInfo_.pitchedInitialAngle(),turbineInfo_.pitchedAngleChangeEndTime(),turbineInfo_.pitchedFinalAngle());
    }
    tensorCorrect();
}

inline void Foam::fv::controller::interp(scalar & variable, const scalar & startTime, const scalar & initialValue, const scalar & endTime, const scalar & finalValue)
{
    scalar t=time_.value();
    if(t <= startTime)
        variable = initialValue;
    else if(t >= endTime)
        variable = finalValue;
    else
    {
        variable = initialValue
            +(finalValue-initialValue)*(t-startTime)/(endTime-startTime);
    }
}

inline void Foam::fv::controller::scalarCorrect(const scalar & torque)
{
    if(turbineInfo_.autoControl())
    {
        if(rotateSpeed_<=turbineInfo_.rSpeed01())
        {
            rotateAccelerate_ = torque /turbineInfo_.iDrivetrain();
            deltaPitch_ = -pitchedAngle_;
            rotateAccI_=0.0;
        }
        else if(rotateSpeed_>turbineInfo_.rSpeed01() && rotateSpeed_<=turbineInfo_.rSpeed02())
        {
            scalar genTorque = turbineInfo_.rGenTorque01() + (turbineInfo_.rGenTorque02()-turbineInfo_.rGenTorque01()) * (rotateSpeed_-turbineInfo_.rSpeed01()) / (turbineInfo_.rSpeed02()-turbineInfo_.rSpeed01());
            rotateAccelerate_ = (torque - genTorque) / turbineInfo_.iDrivetrain();
            deltaPitch_ = -pitchedAngle_;
            rotateAccI_=0.0;
        }
        else if(rotateSpeed_>turbineInfo_.rSpeed02() && rotateSpeed_<=turbineInfo_.rSpeed03())
        {
            scalar velocity = 63*rotateSpeed_/7.55;
            scalar genTorque = 0.5*12469*1.225*velocity*velocity*velocity*0.482/rotateSpeed_;
            rotateAccelerate_ = (torque - genTorque) / turbineInfo_.iDrivetrain();
            deltaPitch_ = -pitchedAngle_;
            rotateAccI_=0.0;
        }
        else if(rotateSpeed_>turbineInfo_.rSpeed03() && rotateSpeed_<=turbineInfo_.rSpeed04())
        {
            scalar genTorque =  turbineInfo_.rGenTorque03() + (turbineInfo_.rGenTorque04()-turbineInfo_.rGenTorque03()) * (rotateSpeed_-turbineInfo_.rSpeed03()) / (turbineInfo_.rSpeed04()-turbineInfo_.rSpeed03());
            rotateAccelerate_ = (torque - genTorque) / turbineInfo_.iDrivetrain();
            deltaPitch_ = -pitchedAngle_;
            rotateAccI_=0.0;
        }
        else
        {
            rotateAccelerate_ = (torque - 5296610 / rotateSpeed_) / turbineInfo_.iDrivetrain();
            Foam::scalar k_P = -2*turbineInfo_.iDrivetrain()*1.267109*0.7*0.6/(turbineInfo_.alphaPToTheta()*(1+pitchedAngle_/turbineInfo_.thetaK()));
            Foam::scalar k_I = -turbineInfo_.iDrivetrain()*1.267109*0.6*0.6/(turbineInfo_.alphaPToTheta()*(1+pitchedAngle_/turbineInfo_.thetaK()));          
            deltaPitch_=k_P*rotateAccelerate_ + k_I * rotateAccI_;
            if(deltaPitch_<0.0 && pitchedAngle_ <=0.0)
                deltaPitch_ = -pitchedAngle_;
        }
        if(deltaPitch_>0.139626*time_.deltaT().value())
            deltaPitch_=0.139626*time_.deltaT().value();
        if(deltaPitch_<-0.139626*time_.deltaT().value())
            deltaPitch_=-0.139626*time_.deltaT().value();
    }
    else
    {
        interp(rotateSpeed_,
            turbineInfo_.rotateSpeedChangeStartTime(),turbineInfo_.rotateInitialSpeed(),turbineInfo_.rotateFinalSpeed(),turbineInfo_.rotateFinalSpeed());
        interp(yawedAngle_,
            turbineInfo_.yawedAngleChangeStartTime(),turbineInfo_.yawedInitialAngle(),turbineInfo_.yawedAngleChangeEndTime(),turbineInfo_.yawedFinalAngle());
        interp(pitchedAngle_,
            turbineInfo_.pitchedAngleChangeStartTime(),turbineInfo_.pitchedInitialAngle(),turbineInfo_.pitchedAngleChangeEndTime(),turbineInfo_.pitchedFinalAngle());
    }
}

inline void Foam::fv::controller::tensorCorrect()
{
    rotateTensor_=tensor::I;
    rotateTensor_.xx()=cos(rotatedAngle_);
    rotateTensor_.xy()=sin(rotatedAngle_);
    rotateTensor_.yx()=-rotateTensor_.xy();
    rotateTensor_.yy()=rotateTensor_.xx();
    
    yawTensor_=tensor::I;
    yawTensor_.xx()=cos(yawedAngle_);
    yawTensor_.xz()=-sin(yawedAngle_);
    yawTensor_.zx()=-yawTensor_.xz();
    yawTensor_.zz()=yawTensor_.xx();

    rotateMatrix_<<cos(rotatedAngle_),sin(rotatedAngle_),0,
                -sin(rotatedAngle_),cos(rotatedAngle_),0,
                0,0,1;

    yawMatrix_<<cos(yawedAngle_),0,-sin(yawedAngle_),
                0,1,0,
                sin(yawedAngle_),0,cos(yawedAngle_);
}

inline void Foam::fv::controller::readCSVLine(std::istringstream & line,scalar & s)
{
    std::string temp;
    getline(line,temp,',');
    std::istringstream d(temp);
    d>>s;
}

void Foam::fv::controller::correct(const scalar & torque)
{
    scalarCorrect(torque);
    tensorCorrect();
}

void Foam::fv::controller::bladeCorrect(scalar angle)
{
    bladeTensor_=tensor::I;
    bladeTensor_.xx()=cos(rotatedAngle_+angle);
    bladeTensor_.xy()=sin(rotatedAngle_+angle);
    bladeTensor_.yx()=-bladeTensor_.xy();
    bladeTensor_.yy()=bladeTensor_.xx();

    bladeMatrix_<<cos(rotatedAngle_+angle),sin(rotatedAngle_+angle),0,
                -sin(rotatedAngle_+angle),cos(rotatedAngle_+angle),0,
                0,0,1;
}

void Foam::fv::controller::timeCorrect()
{
    if(turbineInfo_.autoControl())
    {
        pitchedAngle_ += deltaPitch_*180.0/Foam::constant::mathematical::pi;
        if(pitchedAngle_ <0.0)
            pitchedAngle_ = 0.0;
        rotateSpeed_ += rotateAccelerate_*time_.deltaT().value();
        if(rotateSpeed_>turbineInfo_.rSpeed04())
        	rotateAccI_ += rotateAccelerate_*time_.deltaT().value();
    }
    rotatedAngle_+=rotateSpeed_*time_.deltaT().value();
}

void Foam::fv::controller::workingConditionPrint()
{
    Info<<"**********************************************************"<<endl;
    Info<<"Working condition of turbine "<<turbineInfo_.turbineName()<<" is: "<<endl;
    Info<<"    rotateSpeed: "<<rotateSpeed_<<endl;
    Info<<"    rotatedAngle: "<<rotatedAngle_<<endl;
    Info<<"    yawedAngle: "<<yawedAngle_<<endl;
    Info<<"    pitchedAngle: "<<pitchedAngle_<<endl;
    Info<<"**********************************************************"<<endl;
}

void Foam::fv::controller::printAll()
{
    Info<<"****All*mid*terms*in*controller***************************"<<endl;
    
    Info<<rotorTensor_.xx()<<","<<rotorTensor_.xy()<<","<<rotorTensor_.xz()<<","<<endl;
    Info<<rotorTensor_.yx()<<","<<rotorTensor_.yy()<<","<<rotorTensor_.yz()<<","<<endl;
    Info<<rotorTensor_.zx()<<","<<rotorTensor_.zy()<<","<<rotorTensor_.zz()<<","<<endl;
    Info<< " "<<endl;
    Info<<preconeTensor_.xx()<<","<<preconeTensor_.xy()<<","<<preconeTensor_.xz()<<","<<endl;
    Info<<preconeTensor_.yx()<<","<<preconeTensor_.yy()<<","<<preconeTensor_.yz()<<","<<endl;
    Info<<preconeTensor_.zx()<<","<<preconeTensor_.zy()<<","<<preconeTensor_.zz()<<","<<endl;
    Info<< " "<<endl;
    Info<<rotateTensor_.xx()<<","<<rotateTensor_.xy()<<","<<rotateTensor_.xz()<<","<<endl;
    Info<<rotateTensor_.yx()<<","<<rotateTensor_.yy()<<","<<rotateTensor_.yz()<<","<<endl;
    Info<<rotateTensor_.zx()<<","<<rotateTensor_.zy()<<","<<rotateTensor_.zz()<<","<<endl;
    Info<< " "<<endl;
    Info<<bladeTensor_.xx()<<","<<bladeTensor_.xy()<<","<<bladeTensor_.xz()<<","<<endl;
    Info<<bladeTensor_.yx()<<","<<bladeTensor_.yy()<<","<<bladeTensor_.yz()<<","<<endl;
    Info<<bladeTensor_.zx()<<","<<bladeTensor_.zy()<<","<<bladeTensor_.zz()<<","<<endl;
    Info<< " "<<endl;
    Info<<yawTensor_.xx()<<","<<yawTensor_.xy()<<","<<yawTensor_.xz()<<","<<endl;
    Info<<yawTensor_.yx()<<","<<yawTensor_.yy()<<","<<yawTensor_.yz()<<","<<endl;
    Info<<yawTensor_.zx()<<","<<yawTensor_.zy()<<","<<yawTensor_.zz()<<","<<endl;

    Info<<"**********************************************************"<<endl;
}

void Foam::fv::controller::readControllerResults(std::ifstream & controldata)
{
	std::string line;
    getline(controldata,line);
    getline(controldata,line);
    std::istringstream tdata(line);
    
    readCSVLine(tdata,rotateSpeed_);
    readCSVLine(tdata,rotatedAngle_);
    readCSVLine(tdata,yawedAngle_);
    readCSVLine(tdata,pitchedAngle_);
    
    Info<<"controllerdata of turbine "<<turbineInfo_.turbineName()<<" has been read. "<<endl;
    if(flagBit_.debug01())
    {
        Info<<"Working condition read for turbine "<<turbineInfo_.turbineName()<<" is: "<<endl;
        Info<<"    rotateSpeed: "<<rotateSpeed_<<endl;
        Info<<"    rotatedAngle: "<<rotatedAngle_<<endl;
        Info<<"    yawedAngle: "<<yawedAngle_<<endl;
        Info<<"    pitchedAngle: "<<pitchedAngle_<<endl;
    }
}

void Foam::fv::controller::writeControllerResults(std::ofstream & controldata)
{
	controldata<<"rotateSpeed,rotatedAngle,yawedAngle,pitchedAngle,"<<std::endl;
    controldata<<rotateSpeed_<<","<<rotatedAngle_<<","<<yawedAngle_<<","<<pitchedAngle_<<","<<std::endl;
}

#endif