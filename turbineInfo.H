/****************************************************************************\
This class is developed by MaZhe.

This class is used to storage the informations of rotors.
Only static informations including:
rotorName, towerHeight, gapTowerHub, towerNumber, towerNames
are storaged in this class.

rotorCondition class storage the working conditions of rotors.
\****************************************************************************/

#ifndef turbineInfo_H
#define turbineInfo_H

#include "scalar.H"
#include "List.H"
#include <vector>

namespace Foam
{
namespace fv
{

class turbineInfo
{

public:

//Constructor
    turbineInfo
    (
        dictionary& turbineInfoPath
    ):
        turbineName_(turbineInfoPath.dictName())
    {
        read(turbineInfoPath);
        shaftOriginCalculation();
        interfaceCalulate();
    }

//turbine informations API
    const word & turbineName() const {return turbineName_;}

    const point & towerRootPosition() const {return towerRootPosition_;}

    const point & shaftOrigin() const {return shaftOrigin_;}

    const Eigen::Matrix<double,3,1> & shaftOriginm() const {return shaftOriginm_;}

    const scalar & shaftTiltAngle() const {return shaftTiltAngle_;}

//hub informations
    const scalar & shaftLength() const {return shaftLength_;}

    const scalar & hubInertia() const {return hubInertia_;}

    const scalar & hubMass() const {return hubMass_;}

    const scalar & shaftTorsionStiffness() const {return shaftTorsionStiffness_;}

//nacelle informations
    const scalar & nacelleMassCenter() const {return nacelleMassCenter_;}

    const scalar & nacelleMass() const {return nacelleMass_;}

    const scalar & nacelleHeight() const {return nacelleHeight_;}

    const scalar & nacelleInertia() const {return nacelleInertia_;}

    const scalar & yawTorsionStiffness() const {return yawTorsionStiffness_;}

//rotor informations
    const scalar & rotorPrecone() const {return rotorPrecone_;}

    const int & bladeNumber() const {return bladeNumber_;}

    const word & bladeName() const {return bladeName_;}
//turbine auto control flagbit
    const bool & autoControl() const {return autoControl_;}

//turbine rotate controler informations
    const scalar & rotateInitialSpeed() const {return rotateInitialSpeed_;}

    const scalar & rotateFinalSpeed() const {return rotateFinalSpeed_;}

    const scalar & rotateSpeedChangeStartTime() const {return rotateSpeedChangeStartTime_;}

    const scalar & rotateSpeedChangeEndTime() const {return rotateSpeedChangeEndTime_;}

//turbine yaw controler informations
    const scalar & yawedInitialAngle() const {return yawedInitialAngle_;}

    const scalar & yawedFinalAngle() const {return yawedFinalAngle_;}

    const scalar & yawedAngleChangeStartTime() const {return yawedAngleChangeStartTime_;}

    const scalar & yawedAngleChangeEndTime() const {return yawedAngleChangeEndTime_;}

//turbine pitch controler informations
    const scalar & pitchedInitialAngle() const {return pitchedInitialAngle_;}

    const scalar & pitchedFinalAngle() const {return pitchedFinalAngle_;}

    const scalar & pitchedAngleChangeStartTime() const {return pitchedAngleChangeStartTime_;}

    const scalar & pitchedAngleChangeEndTime() const {return pitchedAngleChangeEndTime_;}

//turbine controller parameters
    const scalar & rSpeed01() const {return rSpeed01_;}

    const scalar & rSpeed02() const {return rSpeed02_;}

    const scalar & rSpeed03() const {return rSpeed03_;}

    const scalar & rSpeed04() const {return rSpeed04_;}

    const scalar & rGenTorque01() const {return rGenTorque01_;}

    const scalar & rGenTorque02() const {return rGenTorque02_;}

    const scalar & rGenTorque03() const {return rGenTorque03_;}

    const scalar & rGenTorque04() const {return rGenTorque04_;}

//controller parameters
    const scalar & iDrivetrain() const {return iDrivetrain_;}

    const scalar & alphaPToTheta() const {return alphaPToTheta_;}

    const scalar & thetaK() const {return thetaK_;}

//turbine aerodynamic characteristic APIs
    const List<point> & towerAEP() const {return towerAEP_;}

    const List<word> & towerAS() const {return towerAS_;}

    const List<List<scalar>> & towerAEI() const {return towerAEI_;}

//turbine finite element characteristic APIs
    const List<point> & towerNP() const {return towerNP_;}

    const List<List<int>> & towerEI() const {return towerEI_;}

    const List<List<scalar>> & towerSI() const {return towerSI_;}

/*
    List<point> towerPositionToAero(const Eigen::MatrixXd & towerPosition);
    std::vector<Eigen::MatrixM<double,6,1>> aeroForceMomentToTower(ListM<vector> & aeroForce, List<vector> & aeroMoment);
*/

private:
//turbine name
    word turbineName_;

//informations must be found
//tower root postion
    point towerRootPosition_;

//initial and default rotate speed
    scalar rotateInitialSpeed_;

//number of blades
    int bladeNumber_;

//blade name
    word bladeName_;

//informations which have default values
//tower tip position
    point shaftOrigin_;
    Eigen::Matrix<double,3,1> shaftOriginm_;

//shaft tilt angle
    scalar shaftTiltAngle_;

//rotor precone angle
    scalar rotorPrecone_;

//shaft length to hub
    scalar shaftLength_;

//hub inertia
    scalar hubInertia_;

//hub mass
    scalar hubMass_;

//shaft torsion stiffness
    scalar shaftTorsionStiffness_;

//shaft length to nacelle mass center
    scalar nacelleMassCenter_;

//nacelle mass
    scalar nacelleMass_;

//nacelle height
    scalar nacelleHeight_;

//nacelle inertia
    scalar nacelleInertia_;

//yaw torsion stiffness
    scalar yawTorsionStiffness_;

//autoControl or constControl
    bool autoControl_;

//control parameters
    scalar rSpeed01_;
    scalar rSpeed02_;
    scalar rSpeed03_;
    scalar rSpeed04_;

    scalar rGenTorque01_;
    scalar rGenTorque02_;
    scalar rGenTorque03_;
    scalar rGenTorque04_;

    scalar iDrivetrain_;
    scalar alphaPToTheta_;
    scalar thetaK_;

//parameters of constControl
    scalar rotateFinalSpeed_;

    scalar rotateSpeedChangeStartTime_;

    scalar rotateSpeedChangeEndTime_;

    scalar yawedInitialAngle_;

    scalar yawedFinalAngle_;

    scalar yawedAngleChangeStartTime_;

    scalar yawedAngleChangeEndTime_;

    scalar pitchedInitialAngle_;

    scalar pitchedFinalAngle_;

    scalar pitchedAngleChangeStartTime_;

    scalar pitchedAngleChangeEndTime_;

//parameters of autoControl
    

//tower informations
//aero informations of tower
//tower aero element position
    List<point> towerAEP_;

//tower aero section
    List<word> towerAS_;

//tower aero element info 
    List<List<scalar>> towerAEI_;

//structure information of tower
//tower node position
    List<point> towerNP_;

//tower element info
    List<List<int>> towerEI_;

//tower element section info
    List<List<scalar>> towerSI_;

//tower aero and finite element interface
    List<List<int>> aeroInTowerNumber_;

    List<scalar> aeroInTowerK_;

//private member functions
    bool read(dictionary& turbineDict);

    bool check();

    void shaftOriginCalculation();
    
    void interfaceCalulate();
};

}//end namespace fv
}//end namespace Foam


/******************************************************************************************************************************\
|                                                                                                                              |
|                                                    function definition                                                       |
|                                                                                                                              |
\******************************************************************************************************************************/

/******************************************private member functions******************************************/

inline bool Foam::fv::turbineInfo::read(dictionary& turbineDict)
{
    //must find profileData!
    if(turbineDict.found("towerRootPosition"))
    {
    //information must be found
        turbineDict.lookup("towerRootPosition") >> towerRootPosition_;
        turbineDict.lookup("rotateInitialSpeed") >> rotateInitialSpeed_;
        turbineDict.lookup("bladeNumber") >> bladeNumber_;
        turbineDict.lookup("bladeName") >> bladeName_;
    //informations which have default values
        shaftTiltAngle_ = turbineDict.lookupOrDefault("shaftTiltAngle",0.0);
        rotorPrecone_ = turbineDict.lookupOrDefault("rotorPrecone",0.0);  
        
        autoControl_= turbineDict.lookupOrDefault<Foam::Switch>("autoControl", false);
        //rotate informations
        rotateFinalSpeed_ = turbineDict.lookupOrDefault("rotateFinalSpeed",rotateInitialSpeed_);
        rotateSpeedChangeStartTime_ = turbineDict.lookupOrDefault("rotateSpeedChangeStartTime",0.0);
        rotateSpeedChangeEndTime_ = turbineDict.lookupOrDefault("rotateSpeedChangeEndTime",1.0);
        //yaw informations
        yawedInitialAngle_ = turbineDict.lookupOrDefault("yawedInitialAngle",0.0);
        yawedFinalAngle_ = turbineDict.lookupOrDefault("yawedFinalAngle",0.0);
        yawedAngleChangeStartTime_ = turbineDict.lookupOrDefault("yawedAngleChangeStartTime",0.0);
        yawedAngleChangeEndTime_ = turbineDict.lookupOrDefault("yawedAngleChangeEndTime",1.0);
        //pitch informations
        pitchedInitialAngle_ = turbineDict.lookupOrDefault("pitchedInitialAngle",0.0);
        pitchedFinalAngle_ = turbineDict.lookupOrDefault("pitchedFinalAngle",0.0);
        pitchedAngleChangeStartTime_ = turbineDict.lookupOrDefault("pitchedAngleChangeStartTime",0.0);
        pitchedAngleChangeEndTime_ = turbineDict.lookupOrDefault("pitchedAngleChangeEndTime",1.0);

        if(autoControl_)
        {
        	turbineDict.lookup("rSpeed01") >> rSpeed01_;
        	turbineDict.lookup("rSpeed02") >> rSpeed02_;
        	turbineDict.lookup("rSpeed03") >> rSpeed03_;
        	turbineDict.lookup("rSpeed04") >> rSpeed04_;

        	turbineDict.lookup("rGenTorque01") >> rGenTorque01_;
        	turbineDict.lookup("rGenTorque02") >> rGenTorque02_;
        	turbineDict.lookup("rGenTorque03") >> rGenTorque03_;
        	turbineDict.lookup("rGenTorque04") >> rGenTorque04_;

        	turbineDict.lookup("iDrivetrain") >> iDrivetrain_;
		    turbineDict.lookup("alphaPToTheta") >> alphaPToTheta_;
		    turbineDict.lookup("thetaK") >> thetaK_;
        }

        if(turbineDict.found("towerAEP"))
        {
            turbineDict.lookup("towerAEP") >> towerAEP_;
            turbineDict.lookup("towerAS") >> towerAS_;
            turbineDict.lookup("towerAEI") >> towerAEI_;
            turbineDict.lookup("towerNP") >> towerNP_;
            turbineDict.lookup("towerEI") >> towerEI_;
            turbineDict.lookup("towerSI") >> towerSI_;
        //hub and shaft informations
            shaftLength_ = turbineDict.lookupOrDefault("shaftLength",1.0);   
            turbineDict.lookup("hubMass") >> hubMass_;
            turbineDict.lookup("hubInertia") >> hubInertia_;
            turbineDict.lookup("shaftTorsionStiffness") >> shaftTorsionStiffness_;
        //nacelle and yaw informations
            nacelleMassCenter_ = turbineDict.lookupOrDefault("nacelleMassCenter",1.0);
            turbineDict.lookup("nacelleMass") >> nacelleMass_;
            turbineDict.lookup("nacelleInertia") >>nacelleInertia_;
            turbineDict.lookup("yawTorsionStiffness") >>yawTorsionStiffness_;
            turbineDict.lookup("nacelleHeight") >> nacelleHeight_;
            if(check())
            {
                return true;
            }
            else
            {
                Info<<"Tower data check failed for "<<turbineName_<<". Please correct the turbine data."<<endl;
                return false;
            }
        }
        else
        {
            Info<<"No towerData for "<<turbineDict.name()<<endl;
            return false;
        }

        return true;
    }
    else
    {
        Info<<"No rotor for "<<turbineDict.name()<<endl;
        return false;
    }
}

inline bool Foam::fv::turbineInfo::check()
{
    if(towerAS_.size()!=towerAEP_.size())
    {
        Info<<"Error: error occured in towerAS data of "<<turbineName_<<endl;
        return false;
    }
    else
    {
        if(towerAEI_.size()!=towerAEP_.size())
        {
            Info<<"Error: error occured in tower aero element data of "<<turbineName_<<endl;
            return false;
        }
        else
        {
            if(towerSI_.size()!=towerNP_.size())
            {
                Info<<"Error: error occured in tower element section data of "<<turbineName_<<endl;
                return false;
            }
            else
            {
                return true;
            }
        }
    }
}

inline void Foam::fv::turbineInfo::shaftOriginCalculation()
{
    point temp;
    temp.x()=towerNP_[towerNP_.size()-1][0];
    temp.y()=towerNP_[towerNP_.size()-1][1]+nacelleHeight_;
    temp.z()=towerNP_[towerNP_.size()-1][2];
    shaftOrigin_=temp+towerRootPosition_;
    shaftOriginm_ << shaftOrigin_.x(),shaftOrigin_.y(),shaftOrigin_.z();
}

inline void Foam::fv::turbineInfo::interfaceCalulate()
{
    List<int> temp;
    temp.setSize(2,0);
    aeroInTowerNumber_.setSize(towerAEP_.size(),temp);
    aeroInTowerK_.setSize(towerAEP_.size(),0.0);
    int j=1;
    forAll(towerAEP_,i)
    {
        while(j!=towerNP_.size())
        {
            vector temp01=towerAEP_[i]-towerNP_[j-1];
            vector temp02=towerAEP_[i]-towerNP_[j];
            if((temp01&temp02)<=0.0)
            {
                aeroInTowerNumber_[i][0]=j-1;
                aeroInTowerNumber_[i][1]=j;
                aeroInTowerK_[i]=mag(temp01)/mag(towerNP_[j]-towerNP_[j-1]);
                break;
            }
            else
            {
                j+=1;
            }
        }
        if(j==towerNP_.size())
        {
            Info<<"Error: Error occured in interfaceCalulate function for aero element number "<<i<<" of turbine "<<turbineName_<<endl;
        }
    }
    Info<<"interface number of "<<turbineName_<<" is "<<aeroInTowerK_<<endl;
}

#endif
