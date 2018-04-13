/****************************************************************************\
This program is based on the openFOAM, and is developed by MaZhe.
The goal of this program is to build an actuatorLineElement class .
\****************************************************************************/

#ifndef actuatorLineElement_H
#define actuatorLineElement_H

#include "airfoilInfo.H"
#include <memory>
#include "dynamicStallModel.H"
#include "flagBit.H"
#include "OFstream.H"
#include "fETurbine.H"

/******************************************class declaration******************************************/

namespace Foam
{
namespace fv
{

class actuatorLineElement
{
public:
	//Constructors
	actuatorLineElement
	(
		const Time& time,
        int eN,
        ALFBM::flagBit & f,
        turbineInfo & tI,
        bladeInfo & bI,
        airfoilInfo & airfoil
	);

    ~actuatorLineElement(){}
//transformation
    const tensor & bladeET() const {return bladeET_;}
//velocity
	const vector & velocity() const {return velocity_;}
//lift and drag force
	const vector & forceLD() const {return forceLiftDrag_;}
//moment
	const vector & forceM() const {return forceMement_;}
//tangential and normal force
	const vector & forceTN() const {return forceTorqueNormal_;}
//phi angle and attack angle
    const scalar & phiAngle() const {return phiAngle_;}
    const scalar & atkAngle() const {return atkAngle_;}
//blade element velocity update
    void velocityUpdate(vector velocityB);
//calculation of blade element force
	void forceCalculate(scalar pitchedAngle);
//element deform
    void turbineElementDeform(ALFBM::fETurbine & fETurbine_, int & b);
//read results
    void readResults(std::ifstream & airfoildata);
//write results
    void writeResults(std::ofstream & airfoildata);
  
private:

    const Time& time_;

    int elementNumber_;

    ALFBM::flagBit & flagBit_;

    turbineInfo& turbineInfo_;

    bladeInfo& bladeInfo_;

    airfoilInfo& airfoilInfo_;
//transformation form blade coordinates to blade element coordinates
    tensor bladeET_;

	vector velocity_;

	vector forceLiftDrag_;

	vector forceTorqueNormal_;

	vector forceMement_;

	scalar phiAngle_;

	scalar atkAngle_;

	List<scalar> listClCdCm_;

	dynamicStallModel dynamicStallModel_;

//sub function for Cl Cd and Cm calculation	
    void findClCdCm(scalar angle);

    void readCSVLine(std::istringstream & line,scalar & s);
    void readResult(std::ifstream & airfoildata,vector & v);
    void readResult(std::ifstream & airfoildata,scalar & s);

    void writeResult(std::ofstream & airfoildata,vector & v);
    void writeResult(std::ofstream & airfoildata,scalar & s);

};

}//end namespace fv
}//end namespace Foam

/******************************************function definition******************************************/

inline void Foam::fv::actuatorLineElement::findClCdCm(scalar angle)
{
    if(angle<=airfoilInfo_.profileData()[0][0])
    {
        listClCdCm_=airfoilInfo_.profileData()[0];
    }
    else if(angle>=airfoilInfo_.profileData()[airfoilInfo_.profileData().size()-1][0])
    {
        listClCdCm_=airfoilInfo_.profileData()[airfoilInfo_.profileData().size()-1];
    }
    else
    {
        forAll(airfoilInfo_.profileData(),i)
        {
            if(angle<=airfoilInfo_.profileData()[i+1][0] && angle>airfoilInfo_.profileData()[i][0])
            {
                listClCdCm_=((angle-airfoilInfo_.profileData()[i][0])/(airfoilInfo_.profileData()[i+1][0]-airfoilInfo_.profileData()[i][0]))*airfoilInfo_.profileData()[i+1]
                    +((airfoilInfo_.profileData()[i+1][0]-angle)/(airfoilInfo_.profileData()[i+1][0]-airfoilInfo_.profileData()[i][0]))*airfoilInfo_.profileData()[i];
                break;
            }
        }
    }
}

inline void Foam::fv::actuatorLineElement::readCSVLine(std::istringstream & line,scalar & s)
{
    std::string temp;
    getline(line,temp,',');
    std::istringstream d(temp);
    d>>s;
}

inline void Foam::fv::actuatorLineElement::readResult(std::ifstream & airfoildata, vector & v)
{
    std::string line;
    getline(airfoildata,line);
    std::istringstream vdata(line);
    readCSVLine(vdata,v.x());
    readCSVLine(vdata,v.y());
    readCSVLine(vdata,v.z());
}

inline void Foam::fv::actuatorLineElement::readResult(std::ifstream & airfoildata, scalar & s)
{
    std::string line;
    getline(airfoildata,line);
    std::istringstream sdata(line);
    readCSVLine(sdata,s);
}

inline void Foam::fv::actuatorLineElement::writeResult(std::ofstream & airfoildata,vector & v)
{
    airfoildata<<v.x()<<","<<v.y()<<","<<v.z()<<","<<std::endl;
}

inline void Foam::fv::actuatorLineElement::writeResult(std::ofstream & airfoildata,scalar & s)
{
    airfoildata<< s <<","<<std::endl;
}

Foam::fv::actuatorLineElement::actuatorLineElement
(
    const Time& time,
    int eN,
    ALFBM::flagBit & f,
    turbineInfo& tI,
    bladeInfo& bI,
    airfoilInfo& airfoil
):
    time_(time),
    elementNumber_(eN),
    flagBit_(f), 
    turbineInfo_(tI),
    bladeInfo_(bI),
    airfoilInfo_(airfoil),
    velocity_(vector::zero),
    forceLiftDrag_(vector::zero),
    forceTorqueNormal_(vector::zero),
    forceMement_(vector::zero),
    dynamicStallModel_(time,airfoil)
{
    listClCdCm_.setSize(airfoilInfo_.profileData()[0].size(),0.0);
    bladeET_=tensor::I;
}

void Foam::fv::actuatorLineElement::velocityUpdate(vector velocityB)
{
    velocity_=bladeET_&velocityB;
    if(flagBit_.debug04())
    {
        Info<<"bladeET for turbine: "<<turbineInfo_.turbineName()<<" blade: "<<bladeInfo_.bladeName()<<" airfoil:"<<airfoilInfo_.airfoilName()<<"."<<endl;
        Info<<bladeET_.xx()<<","<<bladeET_.xy()<<","<<bladeET_.xz()<<","<<endl;
        Info<<bladeET_.yx()<<","<<bladeET_.yy()<<","<<bladeET_.yz()<<","<<endl;
        Info<<bladeET_.zx()<<","<<bladeET_.zy()<<","<<bladeET_.zz()<<","<<endl;
    }
}

void Foam::fv::actuatorLineElement::forceCalculate(scalar pitchedAngle)
{
    scalar pi=Foam::constant::mathematical::pi;
    scalar velocity_n=fabs(velocity_[2]);
    //scalar velocity_t=sqrt(pow(velocity_[0],2)+pow(velocity_[1],2));
    scalar velocity_t=fabs(velocity_[1]);
    phiAngle_=atan2(velocity_n,velocity_t);
    atkAngle_=phiAngle_*180/pi-bladeInfo_.bladeAEI()[elementNumber_][1]-pitchedAngle;//phi-theta-pitch
    if(flagBit_.dynamicStall())
    {
        dynamicStallModel_.correct(mag(velocity_),atkAngle_,listClCdCm_[1],listClCdCm_[2],listClCdCm_[3]);
    }
    else
    {
        findClCdCm(atkAngle_);
    }
    List<scalar> forceTemp;
    forceTemp.setSize(airfoilInfo_.profileData()[0].size(),0.0);
    forceTemp = 0.5*bladeInfo_.bladeAEI()[elementNumber_][0]*bladeInfo_.bladeAEI()[elementNumber_][2]*(pow(velocity_t,2)+pow(velocity_n,2))*listClCdCm_;
    //0.5*chordlength*elementlength*velocity^2*(Cl,Cd,Cm)
    //This is the force from blade element to fluid, which is minus of BEM result
    if(flagBit_.tipCorrection())
    {
        double tipCorr=(2/pi)*acos( exp( -(turbineInfo_.bladeNumber()/2.0)
            * ( (63-bladeInfo_.bladeAEP()[elementNumber_][0])/(bladeInfo_.bladeAEP()[elementNumber_][0]*sin(phiAngle_)) ) ) );
        forceTemp = tipCorr * forceTemp;
    }
    forceLiftDrag_.x()=0.0;
    forceLiftDrag_.y()=forceTemp[1];
    forceLiftDrag_.z()=forceTemp[2];
    forceMement_.x()=forceTemp[3];
    forceMement_.y()=0.0;
    forceMement_.z()=0.0;
    forceTorqueNormal_.x()=0.0;
    forceTorqueNormal_.y()=forceLiftDrag_.y()*sin(phiAngle_)-forceLiftDrag_.z()*cos(phiAngle_);
    forceTorqueNormal_.z()=-(forceLiftDrag_.y()*cos(phiAngle_)+forceLiftDrag_.z()*sin(phiAngle_));
    
    forceTorqueNormal_ = bladeET_.T() & forceTorqueNormal_;
    forceMement_ = bladeET_.T() & forceMement_;
    //Info<<"force calculated for "<<rotorInfo_.rotorName()<<","<<bladeInfo_.bladeName()<<" element "<<elementNumber_<<","<<airfoilInfo_.airfoilName()<<endl;
}

void Foam::fv::actuatorLineElement::turbineElementDeform(ALFBM::fETurbine & fETurbine_, int & bladeNumber)
{
    fETurbine_.turbineElementDeform(bladeET_,bladeNumber,elementNumber_);
}

void Foam::fv::actuatorLineElement::readResults(std::ifstream & airfoildata)
{
    readResult(airfoildata,velocity_);
    readResult(airfoildata,forceLiftDrag_);
    readResult(airfoildata,forceTorqueNormal_);
    readResult(airfoildata,forceMement_);
    readResult(airfoildata,phiAngle_);
    readResult(airfoildata,atkAngle_);
}

void Foam::fv::actuatorLineElement::writeResults(std::ofstream & airfoildata)
{
    writeResult(airfoildata,velocity_);
    writeResult(airfoildata,forceLiftDrag_);
    writeResult(airfoildata,forceTorqueNormal_);
    writeResult(airfoildata,forceMement_);
    writeResult(airfoildata,phiAngle_);
    writeResult(airfoildata,atkAngle_);
}

#endif