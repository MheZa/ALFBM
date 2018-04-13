/*****************************************************\
|                       ALFBM                         |
|                 finiteElementTurbine                |
|                       MaZhe                         |
\*****************************************************/


#ifndef fETurbine_H
#define fETurbine_H
#include "fEElement.H"
#include <fstream>
#include <iostream>
#include "Eigen/IterativeLinearSolvers"
#include "flagBit.H"
#include "controller.H"

namespace ALFBM
{

enum solverType{noDeformSolver,linearSolver,initialSolver,staticSolver,implicitSolver};

class fETurbine
{
public:

/*******************\
|    constructor    |
\*******************/

	fETurbine
    (
        const Foam::Time & time,
        flagBit & f,
        Foam::fv::turbineInfo & t,
        Foam::fv::bladeInfo & b,
        Foam::fv::controller & c
    );

    ~fETurbine(){}

/*******************\
|  public functions |
\*******************/

//
    void nodeReadInitial();
    
    void loadCalculation(Foam::List< Foam::List< Foam::vector>> & aeroForce,Foam::List< Foam::List< Foam::vector>> & aeroMoment);

    void turbineDeform(Foam::List< Foam::List< Foam::point>> & blades, Foam::List<Foam::point> & tower);

    void turbineElementDeform(Foam::tensor & bladeET, const int & b, const unsigned int & i);

    void turbineRotate();
	
	void noDeformSolve();

    void linearSolve();

    void initialSolve();

    void staticInterationSolve();
	
	void implicitIterationSolve();

//read results
    void readStructureResults(std::ifstream & structuredata);
//write results
    void writeModalResults(std::ofstream & modaldata);
    void writeForceResults(std::ofstream & forcedata);
    void writeDeflectionResults(std::ofstream & deflectiondata);
    void writeStructureResults(std::ofstream & structuredata);

    const int & nodeNumber() const {return nodeNumber_;}

private:

/*******************\
| private variables |
\*******************/

//time
    const Foam::Time & time_;

//flag bit
    flagBit & flagBit_;

//reference of turbine information
    Foam::fv::turbineInfo & turbineInfo_;

//reference of blade informations
    Foam::fv::bladeInfo & bladeInfo_;

//references of turbine controller
    Foam::fv::controller & controller_;

//vector of blades nodes
    std::vector<std::vector<fENode>> bladeNodes_;

//vector of nacelle nodes
    std::vector<fENode> nacelleNodes_;

//vector of tower nodes
    std::vector<fENode> towerNodes_;

//number of all nodes
    int nodeNumber_;

//vector of blades elements
    std::vector<std::vector<fEElement>> bladeElements_;

//vector of nacelle elements
    std::vector<fEElement> nacelleElements_;

//vector of tower elements
    std::vector<fEElement> towerElements_;

//vector of node number of blade nacelle connection
    std::vector<int> bladeNacelleConnection_;

//node number of nacelle tower connection
    int nacelleTowerConnection_;

//global stiffness matrix
    Eigen::MatrixXd turbineStiffness_;

//global mass matrix
    Eigen::MatrixXd turbineMass_;

//global damp matrix
    Eigen::MatrixXd turbineDamp_;

//global coriolis damp matrix
    Eigen::MatrixXd coriolisDamp_;

//global spin softening matrix
    Eigen::MatrixXd spinSoften_;

//global stress stiffening matrix
    Eigen::MatrixXd stressStiffen_;

//equivalent stiffness matrix which is used to apply boundary condition and natural frequency calculation
    Eigen::MatrixXd equivalentK_;

//penalty coefficient
    double penaC_;

//time step
    double deltaT_;

//parameters for implicit method
    double c0_;
    double c1_;
    double c2_;
    double c3_;
    double c4_;
    double c5_;
    double c6_;
    double c7_;

//solver type , initial solver by default 
    solverType sT_=initialSolver;

//aerodynamic load
    Eigen::MatrixXd aeroP_;

//centifugal load
    Eigen::MatrixXd centP_;

//gravity load
    Eigen::MatrixXd gravP_;

//restore load 
    Eigen::MatrixXd restoreP_;

//total load
    Eigen::MatrixXd load_;

//last total load 
    Eigen::MatrixXd loadLast_;

//equivalent load
    Eigen::MatrixXd equivalentP_;

//rigid node position for plot
    Eigen::MatrixXd rigidNP_;

//time step results
    Eigen::MatrixXd nP_;
    Eigen::MatrixXd nD_;
    Eigen::MatrixXd nDNext_;
    Eigen::MatrixXd nDDot_;
    Eigen::MatrixXd nDDualDot_;

//modal results
    Eigen::MatrixXcd modalFrequence_;
    Eigen::MatrixXcd modalVector_;

//root force results
    //tower root, nacalle root, blade root in global, blade root in local
    std::vector<Eigen::Matrix<double,6,1>> rootForce_;
//tip deflection results
    //tower tip, nacalle tip, blade tip in global, blade tip in local
    std::vector<Eigen::Matrix<double,6,1>> tipDeflection_;

/*******************\
| private functions |
\*******************/

//node initialization
    void nodeInitial();
//np construction
    //main function
    void nPInitial();
    //sub function
    void nPUpdate(std::vector<fENode> & N);
//element initialization
    void elementInitial();
//parameters initialization
    void parametersInitial();
//root force result initialization
    void rootForceInitial();
//tip deflection result initialization
    void tipDeflectionInitial();
//turbine model construction
    //basic matrix assemble, stiffness matrix, mass matrix and damp matrix
    void turbineBMA();
    //rotate matrix assemble, spin softening matrix and coriolis damp matrix
    void turbineRMA();
    //stress stiffen matrix assemble 
    void turbineSMA();
    //sub function for matrix assemble
    void matrixAssemble(Eigen::MatrixXd & M,  const Eigen::Matrix<double,12,12> & m, int n0, int n1);
    //sub function for basic matrix assemble
    void bMA(std::vector<fEElement> & elements);
    //sub function for stress stiffening matrix assemble
    void sSMA(std::vector<fEElement> & elements);
    //sub function for rotate matrix assemble
    void rMA(std::vector<fEElement> & elements);
//load construction
    //centrifugal load assemble
    void centPA(std::vector<fEElement> & elements);
    //centrifugal load calculation
    void centPCal();
    //gravity load assembel
    void gravPA(std::vector<fEElement> & elements);
    //gravity load calculation
    void gravPCal();
    //sub function for restore load assemble
    void restorePA(std::vector<fEElement> & elements);
    //restore load calculation
    void restorePCal(); 
    //sub function for load assemble
    void loadAssemble(Eigen::MatrixXd & L, const Eigen::Matrix<double,12,1> & l, int n0, int n1);
//rotate function
    //position rotate
    void positionR(Eigen::Matrix<double,6,1> & PA, Eigen::Matrix<double,3,3> omega);
    //force rotate   
    void forceR(Eigen::Matrix<double,6,1> & FA, Eigen::Matrix<double,3,3> omega);
//matrix for solver
    //calculation of equivalent K
    void equivalentKCal();
    //boundary condition
    void boundaryApply();
    //calculation of equivalent P
    void equivalentPCal();
//solve the equation
    //stress stiffening solve
    void stressStiffenSolve();
    //deformation equation solve
    void deformationSolve();
    //modal equation solve
    void modalSolve();
//post process
    //node iteration
    void nodeUpdate(std::vector<fENode> & N);
    void nodeIteration();
    //result iteration
    void resultIteration();
    //root force calculation
    void rootForceCal();
    //tip deflection calculation
    void tipDeflectionCal();
    //sub function of force calculation
    void forceCal(Eigen::Matrix<double,6,1> & F , std::vector<fEElement> & elements);
    void forceG2L(Eigen::Matrix<double,6,1> & Fout , const Eigen::Matrix<double,6,1> & Fin, const int i);
    //sub function of deflection calculation
    void deflectionG2L(Eigen::Matrix<double,6,1> & Fout , const Eigen::Matrix<double,6,1> & Fin, const int i);
    //
//restart function
    //sub function
    void nodeRead(std::vector< fENode> & N);
//read & write function
    void readCSVLine(std::istringstream & line, double & d);

    void readStructureResult(std::ifstream & in, const std::string & dataTitle, int & dataNumber, Eigen::MatrixXd & data);

    void writeStructureResult(std::ofstream & out, const std::string & dataTitle, const int & dataNumber, const Eigen::MatrixXd & data);

    void writeForceResult(std::ofstream & out, const std::string & dataTitle, const int & number, const int & dataNumber, const Eigen::Matrix<double,6,1> & data);
};

}//end namespace ALFBM

/******************************************************************************************************************************\
|                                                                                                                              |
|                                                    function definition                                                       |
|                                                                                                                              |
\******************************************************************************************************************************/

inline void ALFBM::fETurbine::nodeInitial()
{
//node number
    int nodeN=0;

//key vectors
    Foam::point globalp(Foam::vector::zero),gap(Foam::vector::zero),massc(Foam::vector::zero);
    gap.z() = turbineInfo_.shaftLength();
    massc.z() = turbineInfo_.nacelleMassCenter();

//tower nodes
    forAll(turbineInfo_.towerNP(),i)
    {
        fENode nodetemp(nodeN);
        nodetemp.nPSet(turbineInfo_.towerNP()[i]);
        towerNodes_.push_back(nodetemp);
        nodeN+=1;
    }

//shaft origin node
    fENode shaftori(nodeN);
    shaftori.nPSet(turbineInfo_.shaftOrigin());
    nacelleNodes_.push_back(shaftori);
    nodeN+=1;

//blade nodes
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        double installangle=2*Foam::constant::mathematical::pi*i/turbineInfo_.bladeNumber();
        controller_.bladeCorrect(installangle);
        std::vector< fENode> blade;
        forAll(bladeInfo_.bladeNP(),j)
        {
            fENode nodetemp(nodeN);
            globalp = turbineInfo_.shaftOrigin() + (controller_.yawTensor().T()& (controller_.rotorTensor().T()& (gap + (controller_.bladeTensor().T()& (controller_.preconeTensor().T()& bladeInfo_.bladeNP()[j])))));
            nodetemp.nPSet(globalp);
            blade.push_back(nodetemp);
            nodeN+=1;
        }
        bladeNodes_.push_back(blade);
    }

//hub node
    fENode hub(nodeN);
    globalp= turbineInfo_.shaftOrigin() + (controller_.yawTensor().T()& (controller_.rotorTensor().T()& gap));
    hub.nPSet(globalp);
    nacelleNodes_.push_back(hub);
    nodeN+=1;

//nacelle mass node
    fENode mc(nodeN);
    globalp= turbineInfo_.shaftOrigin() + (controller_.yawTensor().T()& (controller_.rotorTensor().T()& (-massc)));
    mc.nPSet(globalp);
    nacelleNodes_.push_back(mc);

    nodeNumber_=nodeN+1;
}

inline void ALFBM::fETurbine::nPUpdate(std::vector<fENode> & N)
{
    for(auto probe=N.begin();probe!=N.end();probe++)
    {
        nP_.block((*probe).nN()*6,0,6,1)=(*probe).nP();
    }
}

inline void ALFBM::fETurbine::nPInitial()
{
    nPUpdate(towerNodes_);
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        nPUpdate(bladeNodes_[i]);
    }
    nPUpdate(nacelleNodes_);
}

inline void ALFBM::fETurbine::parametersInitial()
{
    deltaT_=time_.deltaT().value();

    double alpha=0.25;
    double beta=0.5;
    c0_=1/(alpha*pow(deltaT_,2));
    c1_=beta/(alpha*deltaT_);
    c2_=1/(alpha*deltaT_);
    c3_=(1/(2*alpha))-1;
    c4_=(beta/alpha)-1;
    c5_=(deltaT_/2.0)*((beta/alpha)-2);
    c6_=deltaT_*(1-beta);
    c7_=beta*deltaT_;
    
    aeroP_.setZero(6*nodeNumber_,1);
    centP_.setZero(6*nodeNumber_,1);
    gravP_.setZero(6*nodeNumber_,1);
    restoreP_.setZero(6*nodeNumber_,1);
    load_.setZero(6*nodeNumber_,1);
    loadLast_.setZero(6*nodeNumber_,1);
    equivalentP_.setZero(6*nodeNumber_,1);
    nP_.setZero(6*nodeNumber_,1);
    nD_.setZero(6*nodeNumber_,1);
    nDNext_.setZero(6*nodeNumber_,1);
    nDDot_.setZero(6*nodeNumber_,1);
    nDDualDot_.setZero(6*nodeNumber_,1);
}

inline void ALFBM::fETurbine::rootForceInitial()
{
    Eigen::Matrix<double,6,1> temp;
    temp=Eigen::MatrixXd::Zero(6,1);
    for(int i=0;i<2*turbineInfo_.bladeNumber()+2;++i)
    {
        rootForce_.push_back(temp);
    }
}

inline void ALFBM::fETurbine::tipDeflectionInitial()
{
    Eigen::Matrix<double,6,1> temp;
    temp=Eigen::MatrixXd::Zero(6,1);
    for(int i=0;i<2*turbineInfo_.bladeNumber()+2;++i)
    {
        tipDeflection_.push_back(temp);
    }
}

inline void ALFBM::fETurbine::elementInitial()
{

//element number    
    int elementN=0;

//node mark
    int tn=towerNodes_.size()-1;

//tower elements
    forAll(turbineInfo_.towerEI(),i)
    {
        fEElement towerelement(elementN, towerNodes_[turbineInfo_.towerEI()[i][0]], towerNodes_[turbineInfo_.towerEI()[i][1]],0.01,
                                turbineInfo_.towerSI()[turbineInfo_.towerEI()[i][0]],turbineInfo_.towerSI()[turbineInfo_.towerEI()[i][1]]);
        towerElements_.push_back(towerelement);
        elementN+=1;
    }

//nacelle element connect tower tip and nacelle
    fEElement yawelement(elementN, towerNodes_[tn], nacelleNodes_[0], 0.01,
                            turbineInfo_.towerSI()[tn],false);
    yawelement.eSILAccess()(3,3)+=turbineInfo_.yawTorsionStiffness();
    yawelement.eSILAccess()(3,9)+=-turbineInfo_.yawTorsionStiffness();
    yawelement.eSILAccess()(9,3)+=-turbineInfo_.yawTorsionStiffness();
    yawelement.eSILAccess()(9,9)+=turbineInfo_.yawTorsionStiffness();
    yawelement.eMILAccess()(9,9)+=turbineInfo_.nacelleInertia();
    nacelleElements_.push_back(yawelement);
    elementN+=1;
    

//blade elements
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        std::vector< fEElement> temp;
        forAll(bladeInfo_.bladeEI(),j)
        {
            fEElement bladeelement(elementN, bladeNodes_[i][bladeInfo_.bladeEI()[j][0]], bladeNodes_[i][bladeInfo_.bladeEI()[j][1]],0.00477465,
                                bladeInfo_.bladeSI()[bladeInfo_.bladeEI()[j][0]],bladeInfo_.bladeSI()[bladeInfo_.bladeEI()[j][1]]);
            temp.push_back(bladeelement);
            elementN+=1;
        }
        bladeElements_.push_back(temp);

    //hub blade connecting elements
        fEElement hubtoblade(elementN, nacelleNodes_[1], bladeNodes_[i][0],0.01,
                                bladeInfo_.bladeSI()[0],true);
        nacelleElements_.push_back(hubtoblade);
        elementN+=1;
    }

//shaft front element
    fEElement shaftfront(elementN, nacelleNodes_[0], nacelleNodes_[1], 0.01,
                        turbineInfo_.towerSI()[tn],false);
    //change the mass distribution
    shaftfront.eSILAccess()(3,3)+=turbineInfo_.shaftTorsionStiffness();
    shaftfront.eSILAccess()(3,9)+=-turbineInfo_.shaftTorsionStiffness();
    shaftfront.eSILAccess()(9,3)+=-turbineInfo_.shaftTorsionStiffness();
    shaftfront.eSILAccess()(9,9)+=turbineInfo_.shaftTorsionStiffness();
    yawelement.eMILAccess()(6,6)+=turbineInfo_.hubMass();
    yawelement.eMILAccess()(7,7)+=turbineInfo_.hubMass();
    yawelement.eMILAccess()(8,8)+=turbineInfo_.hubMass();
    yawelement.eMILAccess()(9,9)+=turbineInfo_.hubInertia();
    nacelleElements_.push_back(shaftfront);
    elementN+=1;

    fEElement shaftback(elementN, nacelleNodes_[2], nacelleNodes_[0], 0.01,
                        turbineInfo_.towerSI()[tn],true);
    shaftback.eMILAccess()(0,0)+=turbineInfo_.nacelleMass();
    shaftback.eMILAccess()(1,1)+=turbineInfo_.nacelleMass();
    shaftback.eMILAccess()(2,2)+=turbineInfo_.nacelleMass();
    nacelleElements_.push_back(shaftback);
}

inline void ALFBM::fETurbine::matrixAssemble(Eigen::MatrixXd & M, const Eigen::Matrix<double,12,12> & m, int n0, int n1)
{
    M.block(n0*6,n0*6,6,6) += m.topLeftCorner(6,6);
    M.block(n1*6,n0*6,6,6) += m.bottomLeftCorner(6,6);
    M.block(n0*6,n1*6,6,6) += m.topRightCorner(6,6);
    M.block(n1*6,n1*6,6,6) += m.bottomRightCorner(6,6);
}

inline void ALFBM::fETurbine::bMA(std::vector<fEElement> & elements)
{
    for(auto probe=elements.begin();probe!=elements.end();probe++)
    {
        (*probe).bMC(controller_.pitchedAngle());
        matrixAssemble(turbineStiffness_,(*probe).eSIG(),(*probe).node0().nN(),(*probe).node1().nN());
        matrixAssemble(turbineMass_,(*probe).eMIG(),(*probe).node0().nN(),(*probe).node1().nN());
        matrixAssemble(turbineDamp_,(*probe).eCIG(),(*probe).node0().nN(),(*probe).node1().nN());
    }
}

inline void ALFBM::fETurbine::sSMA(std::vector<fEElement> & elements)
{
    for(auto probe=elements.begin();probe!=elements.end();probe++)
    {
        Eigen::Matrix<double,12,1> nDis;
        nDis.block(0,0,6,1)=nDNext_.block(6*(*probe).node0().nN(),0,6,1);
        nDis.block(6,0,6,1)=nDNext_.block(6*(*probe).node1().nN(),0,6,1);
        (*probe).sSMC(controller_.pitchedAngle(),nDis);
        matrixAssemble(stressStiffen_,(*probe).eStressStiffenIG(),(*probe).node0().nN(),(*probe).node1().nN());
    }
}

inline void ALFBM::fETurbine::rMA(std::vector<fEElement> & elements)
{
    Eigen::Matrix<double,3,1> localRS;
    localRS << 0,0,controller_.rotateSpeed();
    Eigen::Matrix<double,3,1> globalRS;
    globalRS = controller_.yawMatrix().transpose()* (controller_.rotorMatrix().transpose()* localRS);
    Eigen::Matrix<double,3,1> center;
    center << nacelleNodes_[0].nP()(0,0),nacelleNodes_[0].nP()(1,0),nacelleNodes_[0].nP()(2,0);
    for(auto probe=elements.begin();probe!=elements.end();probe++)
    {
        //Foam::Info<<globalRS<<","<<center<<Foam::endl;
        (*probe).rMC(controller_.pitchedAngle(),globalRS,center);
        //Foam::Info<<(*probe).eSpinSoftenIG()<<Foam::endl;
        matrixAssemble(spinSoften_,(*probe).eSpinSoftenIG(),(*probe).node0().nN(),(*probe).node1().nN());
        matrixAssemble(coriolisDamp_,(*probe).eRCIG(),(*probe).node0().nN(),(*probe).node1().nN());
    }
}

inline void ALFBM::fETurbine::turbineBMA()
{
    turbineStiffness_.setZero(6*nodeNumber_,6*nodeNumber_);
    turbineMass_.setZero(6*nodeNumber_,6*nodeNumber_);
    turbineDamp_.setZero(6*nodeNumber_,6*nodeNumber_);
    bMA(towerElements_);
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        bMA(bladeElements_[i]);
    }
    bMA(nacelleElements_);
}

inline void ALFBM::fETurbine::turbineRMA()
{
    spinSoften_.setZero(6*nodeNumber_,6*nodeNumber_);
    coriolisDamp_.setZero(6*nodeNumber_,6*nodeNumber_);
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        rMA(bladeElements_[i]);
    }
}

inline void ALFBM::fETurbine::turbineSMA()
{
    stressStiffen_.setZero(6*nodeNumber_,6*nodeNumber_);
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        sSMA(bladeElements_[i]);
    }
}

inline void ALFBM::fETurbine::centPA(std::vector<fEElement> & elements)
{
    for(auto probe=elements.begin();probe!=elements.end();probe++)
    {
        loadAssemble(centP_,(*probe).eCCIG(),(*probe).node0().nN(),(*probe).node1().nN());
    }
}

inline void ALFBM::fETurbine::gravPA(std::vector<fEElement> & elements)
{
    for(auto probe=elements.begin();probe!=elements.end();probe++)
    {
        loadAssemble(gravP_,(*probe).eGIG(),(*probe).node0().nN(),(*probe).node1().nN());
    }
}
inline void ALFBM::fETurbine::restorePA(std::vector<fEElement> & elements)
{
    for(auto probe=elements.begin();probe!=elements.end();probe++)
    {
        Eigen::Matrix<double,12,1> nDis;
        nDis.block(0,0,6,1)=nDNext_.block(6*(*probe).node0().nN(),0,6,1);
        nDis.block(6,0,6,1)=nDNext_.block(6*(*probe).node1().nN(),0,6,1);
        (*probe).rLC(controller_.pitchedAngle(),nDis);
        loadAssemble(restoreP_,(*probe).eRIG(),(*probe).node0().nN(),(*probe).node1().nN());
    }
}

inline void ALFBM::fETurbine::loadAssemble(Eigen::MatrixXd & L, const Eigen::Matrix<double,12,1> & l, int n0, int n1)
{
    L.block(n0*6,0,6,1) += l.block(0,0,6,1);
    L.block(n1*6,0,6,1) += l.block(6,0,6,1);
}

inline void ALFBM::fETurbine::centPCal()
{
    centP_=-0.5*spinSoften_*nP_;
    centPA(towerElements_);
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        centPA(bladeElements_[i]);
    }
    centPA(nacelleElements_);
}

inline void ALFBM::fETurbine::gravPCal()
{
    gravP_.setZero(6*nodeNumber_,1);
    gravPA(towerElements_);
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        gravPA(bladeElements_[i]);
    }
    gravPA(nacelleElements_);
}

inline void ALFBM::fETurbine::restorePCal()
{
    restoreP_.setZero(6*nodeNumber_,1);
    restorePA(towerElements_);
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        restorePA(bladeElements_[i]);
    }
    restorePA(nacelleElements_);
}

inline void ALFBM::fETurbine::nodeRead(std::vector< fENode> & N)
{
    if(flagBit_.nonlinear())
        for(auto probe=N.begin(); probe!=N.end(); probe++)
        {
            (*probe).nPSet(nP_.block((*probe).nN()*6,0,6,1));
        }
    else
    {
        for(auto probe=N.begin(); probe!=N.end(); probe++)
        {
            (*probe).nPSet(rigidNP_.block((*probe).nN()*6,0,6,1));
        }
    }
}

inline void ALFBM::fETurbine::nodeReadInitial()
{
    nodeRead(towerNodes_);
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        nodeRead(bladeNodes_[i]);
    }
    nodeRead(nacelleNodes_);
}

inline void ALFBM::fETurbine::nodeUpdate(std::vector< fENode> & N)
{
    for(auto probe=N.begin(); probe!=N.end(); probe++)
    {
        (*probe).nPAdd(nDNext_.block((*probe).nN()*6,0,6,1));
    }
}

inline void ALFBM::fETurbine::nodeIteration()
{
    nodeUpdate(towerNodes_);
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        nodeUpdate(bladeNodes_[i]);
    }
    nodeUpdate(nacelleNodes_);
}

inline void ALFBM::fETurbine::positionR(Eigen::Matrix<double,6,1> & PA, Eigen::Matrix<double,3,3> omega)
{
    PA.block(0,0,3,1)=nacelleNodes_[1].nP().block(0,0,3,1)+(controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*
        omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*( PA.block(0,0,3,1)- nacelleNodes_[1].nP().block(0,0,3,1)));

    PA.block(3,0,3,1)=controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*
        omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*PA.block(3,0,3,1);
}

inline void ALFBM::fETurbine::forceR(Eigen::Matrix<double,6,1> & FA, Eigen::Matrix<double,3,3> omega)
{
    FA.block(0,0,3,1)=controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*
        omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()* FA.block(0,0,3,1);

    FA.block(3,0,3,1)=controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*
        omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*FA.block(3,0,3,1);
}

inline void ALFBM::fETurbine::equivalentKCal()
{
    if(sT_==noDeformSolver||sT_==linearSolver||sT_==staticSolver||sT_==initialSolver)
    {
        equivalentK_=turbineStiffness_;        
    }
    else
    {
        equivalentK_=turbineStiffness_+turbineMass_*c0_;
    }
    if(flagBit_.spinSoften())
        equivalentK_+=spinSoften_;
    if(flagBit_.stressStiffen())
        equivalentK_+=stressStiffen_;
    penaC_=50000*equivalentK_.maxCoeff();
}

inline void ALFBM::fETurbine::boundaryApply()
{
    for(int i=0;i<6;++i)
        equivalentK_(6*nacelleNodes_[0].nN()+i,6*nacelleNodes_[0].nN()+i) += penaC_;
}

inline void ALFBM::fETurbine::resultIteration()
{
    if(sT_==initialSolver)
    {
        loadLast_=load_;
    }
    else if (sT_==implicitSolver)
    {
        loadLast_=load_;
        load_.setZero(nodeNumber_*6,1);

        Eigen::MatrixXd nDDotNext;
        Eigen::MatrixXd nDDualDotNext;

        nDDualDotNext=c0_*(nDNext_-nD_)-c2_*nDDot_-c3_*nDDualDot_;
        nDDotNext=nDDot_+c6_*nDDualDot_+c7_*nDDualDotNext;

        nDDot_=nDDotNext;
        nDDualDot_=nDDualDotNext;
    }
    else if (sT_==staticSolver)
    {
        loadLast_=load_;
        //load_.setZero(nodeNumber_*6,1);
    }
    else
    {

    }
    //nDLast_=nD_;
    nD_=nDNext_;
    nDNext_.setZero(nodeNumber_*6,1);

}

inline void ALFBM::fETurbine::equivalentPCal()
{
    if(sT_==linearSolver)
    {
        equivalentP_=load_;
    }
    else if(sT_==initialSolver)
    {
        equivalentP_=load_/50;
    }
    else if (sT_==staticSolver)
    {
        equivalentP_=load_/50;
    }
    else if(sT_==implicitSolver)
    { 
        //equivalentP_=(load_-loadLast_)+turbineMass_*nD_/pow(deltaT_,2);
        equivalentP_=(load_-loadLast_)+turbineMass_*(c0_*nD_+c2_*nDDot_+c3_*nDDualDot_);
        /*
        if(flagBit_.spinSoften())
            equivalentP_=(load_-loadLast_)-(turbineStiffness_+2*turbineRhoNOmegaOmegaN_)*nD_+turbineMass_*(2*nD_-nDLast_)/pow(deltaT_,2);
        */
        //Foam::Info<<"turbineStiffness_*nD_ is:"<<Foam::endl;
        //Foam::Info<<turbineStiffness_*nD_<<Foam::endl;
    }
    else
    {

    }
}

inline void ALFBM::fETurbine::stressStiffenSolve()
{
    nDNext_=turbineStiffness_.ldlt().solve(centP_);
    turbineSMA();
}

inline void ALFBM::fETurbine::deformationSolve()
{
    if(flagBit_.cg())
    {
    	Eigen::ConjugateGradient<Eigen::MatrixXd> cg;
    	cg.compute(equivalentK_);
    	nDNext_=cg.solve(equivalentP_);
    }
    else
    {
    	nDNext_=equivalentK_.ldlt().solve(equivalentP_);
    }
    if(sT_==initialSolver || sT_==staticSolver || sT_==implicitSolver)
        nodeIteration();
}

inline void ALFBM::fETurbine::modalSolve()
{
    if(flagBit_.modalSolve())
    {
    	Eigen::MatrixXd M_K_=turbineStiffness_;
	    if(flagBit_.spinSoften())
	        M_K_+=spinSoften_;
	    if(flagBit_.stressStiffen())
	        M_K_+=stressStiffen_;
	    penaC_=100000*M_K_.maxCoeff();

	    for(int i=0;i<6;++i)
	        M_K_(i,i) += penaC_;
	/*
	    for(int i=6*nacelleNodes_[0].nN();i<6*nacelleNodes_[0].nN()+6;++i)
	        M_K_(i,i) += penaC_;
	*/
	    if(flagBit_.damp())
	    {
	    	Eigen::MatrixXd A,B;
	    	A.setZero(12*nodeNumber_,12*nodeNumber_);
	    	B.setZero(12*nodeNumber_,12*nodeNumber_);

	    	A.topLeftCorner(6*nodeNumber_,6*nodeNumber_)=turbineDamp_;
	    	A.topRightCorner(6*nodeNumber_,6*nodeNumber_)=M_K_;
	    	A.bottomLeftCorner(6*nodeNumber_,6*nodeNumber_)=-M_K_;

	    	B.topLeftCorner(6*nodeNumber_,6*nodeNumber_)=turbineMass_;
	    	B.bottomRightCorner(6*nodeNumber_,6*nodeNumber_)=M_K_;

	    	Eigen::LLT<Eigen::MatrixXd> lltOfA(B);
		    Eigen::MatrixXd L =lltOfA.matrixL();

		    Eigen::MatrixXd C;
		    C=L.inverse()*A*(L.transpose()).inverse();

		    Eigen::EigenSolver<Eigen::MatrixXd> es(C);

		    modalFrequence_=-es.eigenvalues();
		    modalVector_=es.eigenvectors();
	    }
	    else
	    {
		    M_K_=turbineMass_.inverse()*M_K_;

		    Eigen::EigenSolver<Eigen::MatrixXd> es(M_K_);

		    modalFrequence_=es.eigenvalues();
		    modalVector_=es.eigenvectors();
	    }
	    if(flagBit_.debug04())
	    {
	        std::cout<<"********************************"<<std::endl;
	        std::cout<<modalFrequence_<<std::endl;
	        std::cout<<"********************************"<<std::endl;
	    }
    }
}

inline void ALFBM::fETurbine::rootForceCal()
{
    rootForce_[0].setZero(6,1);
    rootForce_[1].setZero(6,1);
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        forceCal(rootForce_[i+2],bladeElements_[i]);
        rootForce_[1]+=rootForce_[i+2];
    }
    Eigen::Matrix<double,6,1> temp;
    forceCal(temp,towerElements_);
    rootForce_[0]=rootForce_[1]+temp;
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        forceG2L(rootForce_[i+2+turbineInfo_.bladeNumber()],rootForce_[i+2],i);
    }
}

inline void ALFBM::fETurbine::tipDeflectionCal()
{
    tipDeflection_[0]=nP_.block(6*nacelleNodes_[0].nN(),0,6,1)-rigidNP_.block(6*nacelleNodes_[0].nN(),0,6,1);
    tipDeflection_[1]=nP_.block(6*nacelleNodes_[1].nN(),0,6,1)-rigidNP_.block(6*nacelleNodes_[1].nN(),0,6,1);
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        tipDeflection_[2+i]=nP_.block(6*bladeNodes_[i][bladeNodes_[i].size()-1].nN(),0,6,1)-rigidNP_.block(6*bladeNodes_[i][bladeNodes_[i].size()-1].nN(),0,6,1);
    }
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        deflectionG2L(tipDeflection_[2+i+turbineInfo_.bladeNumber()],tipDeflection_[i+2],i);
    }
}

inline void ALFBM::fETurbine::forceCal(Eigen::Matrix<double,6,1> & F , std::vector<fEElement> & elements)
{
    F.setZero(6,1);
    for(auto probe=elements.begin();probe!=elements.end();probe++)
    {
        if(probe == (elements.end()-1))
        {
            F += load_.block(6*(*probe).node0().nN(),0,6,1);
            F += load_.block(6*(*probe).node1().nN(),0,6,1);
        }
        else
        {
            F += load_.block(6*(*probe).node0().nN(),0,6,1);
        }
    }
}

inline void ALFBM::fETurbine::forceG2L(Eigen::Matrix<double,6,1> & Fout , const Eigen::Matrix<double,6,1> & Fin, const int i)
{
    double installangle=2*Foam::constant::mathematical::pi*i/turbineInfo_.bladeNumber();
    controller_.bladeCorrect(installangle);
    Fout.block(0,0,3,1) = controller_.preconeMatrix() * ( controller_.bladeMatrix() * ( controller_.rotorMatrix() * ( controller_.yawMatrix() * Fin.block(0,0,3,1) ) ) );
    Fout.block(3,0,3,1) = controller_.preconeMatrix() * ( controller_.bladeMatrix() * ( controller_.rotorMatrix() * ( controller_.yawMatrix() * Fin.block(3,0,3,1) ) ) );
}

inline void ALFBM::fETurbine::deflectionG2L(Eigen::Matrix<double,6,1> & Fout , const Eigen::Matrix<double,6,1> & Fin, const int i)
{
    double installangle=2*Foam::constant::mathematical::pi*i/turbineInfo_.bladeNumber();
    controller_.bladeCorrect(installangle);
    Fout.block(0,0,3,1) =  controller_.bladeMatrix() * ( controller_.rotorMatrix() * ( controller_.yawMatrix() * Fin.block(0,0,3,1) ) ) ;
    Fout.block(3,0,3,1) =  controller_.bladeMatrix() * ( controller_.rotorMatrix() * ( controller_.yawMatrix() * Fin.block(3,0,3,1) ) ) ;
}

inline void ALFBM::fETurbine::readCSVLine(std::istringstream & line, double & d)
{
    std::string temp;
    getline(line,temp,',');
    std::istringstream tempd(temp);
    tempd>>d;
}

inline void ALFBM::fETurbine::readStructureResult(std::ifstream & in, const std::string & dataTitle, int & dataNumber, Eigen::MatrixXd & data)
{
    std::string line;
    getline(in,line);
    std::istringstream dataformfile(line);
    std::string title;
    dataformfile >> title;
    if(title==dataTitle)
    {
        for(int i=0; i<dataNumber; ++i)
        {
            getline(in,line);
            std::istringstream filedata(line);
            for(int j=0;j<6;++j)
            {
                readCSVLine(filedata,data(6*i+j,0));
            }
            //filedata >> data(6*i,0) >> data(6*i+1,0) >>data(6*i+2,0)
            //        >> data(6*i+3,0) >> data(6*i+4,0) >> data(6*i+5,0);
        }
    }
    else
    {
        Foam::Info<<"Data error for "<< dataTitle <<" of turbine "<<"."<<Foam::endl;
    }
}

inline void ALFBM::fETurbine::writeStructureResult(std::ofstream & out, const std::string & dataTitle, const int & dataNumber, const Eigen::MatrixXd & data)
{
    out << dataTitle << std::endl;
    for(int i=0; i<dataNumber; ++i)
    {
        out<<data(6*i,0)<<","<<data(6*i+1,0)<<","<<data(6*i+2,0)
                <<","<<data(6*i+3,0)<<","<<data(6*i+4,0)<<","<<data(6*i+5,0)<<","<<std::endl;
    }
}

inline void ALFBM::fETurbine::writeForceResult(std::ofstream & out, const std::string & dataTitle, const int & number, const int & dataNumber, const Eigen::Matrix<double,6,1> & data)
{
    out << dataTitle << ", "<<number << std::endl;
        out<<data(0,0)<<","<<data(1,0)<<","<<data(2,0)
                <<","<<data(3,0)<<","<<data(4,0)<<","<<data(5,0)<<","<<std::endl;
}

/******************************************public functions******************************************/
ALFBM::fETurbine::fETurbine
(
    const Foam::Time & time,
    flagBit & f,
    Foam::fv::turbineInfo & t,
    Foam::fv::bladeInfo & b,
    Foam::fv::controller & c
):
    time_(time),
    flagBit_(f),
    turbineInfo_(t),
    bladeInfo_(b),
    controller_(c)
{
    nodeInitial();
    parametersInitial();
    nPInitial();
    rigidNP_=nP_;
    elementInitial();
    rootForceInitial();
    tipDeflectionInitial();
}

void ALFBM::fETurbine::loadCalculation(Foam::List< Foam::List< Foam::vector>> & aeroForce,Foam::List< Foam::List< Foam::vector>> & aeroMoment)
{
    turbineRMA();
    centPCal();
    gravPCal();
    load_.setZero(6*nodeNumber_,1);
    aeroP_.setZero(6*nodeNumber_,1);
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        forAll(aeroForce[i],j)
        {
            Eigen::Matrix<double,6,1> F1;
            Eigen::Matrix<double,6,1> F2;
            Foam::vector forcetemp01;
            Foam::vector forcetemp02;
            Foam::vector momenttemp01;
            Foam::vector momenttemp02;
            forcetemp01=(1-bladeInfo_.aeroInBeamK()[j])*flagBit_.airDensity()*aeroForce[i][j];
            forcetemp02=bladeInfo_.aeroInBeamK()[j]*flagBit_.airDensity()*aeroForce[i][j];
            momenttemp01=(1-bladeInfo_.aeroInBeamK()[j])*flagBit_.airDensity()*aeroMoment[i][j];
            momenttemp02=bladeInfo_.aeroInBeamK()[j]*flagBit_.airDensity()*aeroMoment[i][j];
            F1<<forcetemp01[0],forcetemp01[1],forcetemp01[2],momenttemp01[0],momenttemp01[1],momenttemp01[2];
            F2<<forcetemp02[0],forcetemp02[1],forcetemp02[2],momenttemp02[0],momenttemp02[1],momenttemp02[2];
            aeroP_.block(6*bladeNodes_[i][bladeInfo_.aeroInBeamNumber()[j][0]].nN(),0,6,1) += F1;
            aeroP_.block(6*bladeNodes_[i][bladeInfo_.aeroInBeamNumber()[j][1]].nN(),0,6,1) += F2;
        }
    }
    
    load_=aeroP_+gravP_+centP_;
    if(flagBit_.nonlinear()&&flagBit_.dynamic())
    {
        load_-=restoreP_;
    }

    //Foam::Info<<aeroP_<<Foam::endl;
    //Foam::Info<<"********************"<<Foam::endl;
    //Foam::Info<<gravP_<<Foam::endl;
    //Foam::Info<<"********************"<<Foam::endl;
    //Foam::Info<<centP_<<Foam::endl;
}

void ALFBM::fETurbine::turbineRotate()
{
    double deltaAngle(controller_.rotateSpeed()*time_.deltaT().value());
    Eigen::Matrix<double,3,3> omega;
    omega<<cos(deltaAngle),sin(deltaAngle),0.0,
            -sin(deltaAngle),cos(deltaAngle),0.0,
            0.0,0.0,1.0;
    Eigen::Matrix<double,3,1> zero;
    zero<<0.0,0.0,0.0;
    //Foam::Info<<"check00"<<Foam::endl;
    for(unsigned int i=0;i<bladeNodes_.size();++i)
    {
        for(auto probe=bladeNodes_[i].begin();probe!=bladeNodes_[i].end();probe++)
        {
            Eigen::Matrix<double,6,1> temp;

            temp=(*probe).nP();
            positionR(temp,omega);
            (*probe).nPSet(temp);

            temp=nP_.block(6*(*probe).nN(),0,6,1);
            positionR(temp,omega);
            nP_.block(6*(*probe).nN(),0,6,1)=temp;

            temp=rigidNP_.block(6*(*probe).nN(),0,6,1);
            positionR(temp,omega);
            rigidNP_.block(6*(*probe).nN(),0,6,1)=temp;

            temp=loadLast_.block(6*(*probe).nN(),0,6,1);
            forceR(temp,omega);
            loadLast_.block(6*(*probe).nN(),0,6,1)=temp;
            
            temp=nD_.block(6*(*probe).nN(),0,6,1);
            forceR(temp,omega);
            nD_.block(6*(*probe).nN(),0,6,1)=temp;

            /*
            Eigen::Vector3d t;
            Eigen::Vector3d orit;

            Eigen::Matrix<double,6,1> temp((*probe).nP());
            t=(*probe).nP().block(0,0,3,1);
            orit=omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*(t- nacelleNodes_[1].nP().block(0,0,3,1));
            temp.block(0,0,3,1)=nacelleNodes_[1].nP().block(0,0,3,1)+(controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*orit);

            t=(*probe).nP().block(3,0,3,1);
            orit=omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*t;
            temp.block(3,0,3,1)=controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*orit;
            (*probe).nPSet(temp);
            
            t=nP_.block(6*(*probe).nN(),0,3,1);
            orit=omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*(t-nacelleNodes_[1].nP().block(0,0,3,1));
            nP_.block(6*(*probe).nN(),0,3,1)=nacelleNodes_[1].nP().block(0,0,3,1)+(controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*orit);

            t=nP_.block(6*(*probe).nN()+3,0,3,1);
            orit=omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*t;
            nP_.block(6*(*probe).nN()+3,0,3,1)=controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*orit;

            t=rigidNP_.block(6*(*probe).nN(),0,3,1);
            orit=omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*(t-nacelleNodes_[1].nP().block(0,0,3,1));
            rigidNP_.block(6*(*probe).nN(),0,3,1)=nacelleNodes_[1].nP().block(0,0,3,1)+(controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*orit);

            t=rigidNP_.block(6*(*probe).nN()+3,0,3,1);
            orit=omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*t;
            rigidNP_.block(6*(*probe).nN()+3,0,3,1)=controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*orit;

            t=loadLast_.block(6*(*probe).nN(),0,3,1);
            orit=omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*t;
            loadLast_.block(6*(*probe).nN(),0,3,1)=controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*orit;

            t=loadLast_.block(6*(*probe).nN()+3,0,3,1);
            orit=omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*t;
            loadLast_.block(6*(*probe).nN()+3,0,3,1)=controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*orit;

            t=nD_.block(6*(*probe).nN(),0,3,1);
            orit=omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*t;
            nD_.block(6*(*probe).nN(),0,3,1)=controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*orit;

            t=nD_.block(6*(*probe).nN()+3,0,3,1);
            orit=omega.transpose()*controller_.rotorMatrix()*controller_.yawMatrix()*t;
            nD_.block(6*(*probe).nN()+3,0,3,1)=controller_.yawMatrix().transpose()*controller_.rotorMatrix().transpose()*orit;
            */
        }
    }
    if (sT_==initialSolver || sT_==staticSolver || sT_==implicitSolver)
        nPInitial();
}

void ALFBM::fETurbine::turbineDeform(Foam::List< Foam::List< Foam::point>> & blades, Foam::List<Foam::point> & tower)
{
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        forAll(blades[i],j)
        {
            blades[i][j].x()=bladeInfo_.aeroInBeamK()[j]*nP_(6*bladeNodes_[i][bladeInfo_.aeroInBeamNumber()[j][0]].nN(),0)
                +(1-bladeInfo_.aeroInBeamK()[j])*nP_(6*bladeNodes_[i][bladeInfo_.aeroInBeamNumber()[j][1]].nN(),0);
            blades[i][j].y()=bladeInfo_.aeroInBeamK()[j]*nP_(6*bladeNodes_[i][bladeInfo_.aeroInBeamNumber()[j][0]].nN()+1,0)
                +(1-bladeInfo_.aeroInBeamK()[j])*nP_(6*bladeNodes_[i][bladeInfo_.aeroInBeamNumber()[j][1]].nN()+1,0);
            blades[i][j].z()=bladeInfo_.aeroInBeamK()[j]*nP_(6*bladeNodes_[i][bladeInfo_.aeroInBeamNumber()[j][0]].nN()+2,0)
                +(1-bladeInfo_.aeroInBeamK()[j])*nP_(6*bladeNodes_[i][bladeInfo_.aeroInBeamNumber()[j][1]].nN()+2,0);
        }
    }
}

void ALFBM::fETurbine::turbineElementDeform(Foam::tensor & bladeET, const int & b , const unsigned int & i)
{
    Eigen::Matrix<double,3,3> transX_;
    Eigen::Matrix<double,3,3> transY_;
    Eigen::Matrix<double,3,3> transZ_;

    Eigen::Matrix<double,3,1> Epos0;
    Epos0 = nP_.block(6*bladeNodes_[b][bladeInfo_.aeroInBeamNumber()[i][1]].nN(),0,3,1)-
            nP_.block(6*bladeNodes_[b][bladeInfo_.aeroInBeamNumber()[i][0]].nN(),0,3,1);

    double theta1=atan2(Epos0(1,0),Epos0(0,0));
    double theta2=atan2(Epos0(2,0),sqrt(pow(Epos0(0,0),2)+pow(Epos0(1,0),2)));

    transY_<<cos(theta2),0,-sin(-theta2),
            0,1.0,0,
            sin(-theta2),0,cos(theta2);
    transZ_<<cos(theta1),sin(theta1),0,
            -sin(theta1),cos(theta1),0,
            0,0,1.0;

    Eigen::Matrix<double,3,3> tempET0;
    tempET0=transY_*transZ_;

    Eigen::Matrix<double,3,1> Epos1;
    Epos1 = tempET0* 0.5* 
            (
                nP_.block(6*bladeNodes_[b][bladeInfo_.aeroInBeamNumber()[i][1]].nN()+3,0,3,1)+
                nP_.block(6*bladeNodes_[b][bladeInfo_.aeroInBeamNumber()[i][0]].nN()+3,0,3,1)
            );

    double theta3= Epos1(0,0);

    transX_<<1.0,0,0,
            0,cos(theta3),sin(theta3),
            0,-sin(theta3),cos(theta3);

    tempET0=transX_*transY_*transZ_;

    Epos0 = rigidNP_.block(6*bladeNodes_[b][bladeInfo_.aeroInBeamNumber()[i][1]].nN(),0,3,1)-
            rigidNP_.block(6*bladeNodes_[b][bladeInfo_.aeroInBeamNumber()[i][0]].nN(),0,3,1);

    theta1=atan2(Epos0(1,0),Epos0(0,0));
    theta2=atan2(Epos0(2,0),sqrt(pow(Epos0(0,0),2)+pow(Epos0(1,0),2)));

    transY_<<cos(theta2),0,-sin(-theta2),
            0,1.0,0,
            sin(-theta2),0,cos(theta2);
    transZ_<<cos(theta1),sin(theta1),0,
            -sin(theta1),cos(theta1),0,
            0,0,1.0;

    Eigen::Matrix<double,3,3> tempET1;
    tempET1=transY_*transZ_;

    Epos1 = tempET0* 0.5* 
            (
                rigidNP_.block(6*bladeNodes_[b][bladeInfo_.aeroInBeamNumber()[i][1]].nN()+3,0,3,1)+
                rigidNP_.block(6*bladeNodes_[b][bladeInfo_.aeroInBeamNumber()[i][0]].nN()+3,0,3,1)
            );

    theta3= Epos1(0,0);

    transX_<<1.0,0,0,
            0,cos(theta3),sin(theta3),
            0,-sin(theta3),cos(theta3);

    tempET1=transX_*transY_*transZ_;

    Eigen::Matrix<double,3,3> tempET;
    tempET=tempET0*tempET1.transpose();

    bladeET.xx() = tempET(0,0);
    bladeET.xy() = tempET(0,1);
    bladeET.xz() = tempET(0,2);
    bladeET.yx() = tempET(1,0);
    bladeET.yy() = tempET(1,1);
    bladeET.yz() = tempET(1,2);
    bladeET.zx() = tempET(2,0);
    bladeET.zy() = tempET(2,1);
    bladeET.zz() = tempET(2,2);
}

void ALFBM::fETurbine::noDeformSolve()
{
    sT_=noDeformSolver;
    turbineBMA();
    stressStiffenSolve();
    modalSolve();

    Foam::Info<<"No Deform Solver Used!"<<Foam::endl;
}

void ALFBM::fETurbine::linearSolve()
{
    sT_=linearSolver;
    nPInitial();
    turbineBMA();
    stressStiffenSolve();
    equivalentKCal();
    modalSolve();
    boundaryApply();
    equivalentPCal();
    deformationSolve();
    rootForceCal();
    nP_+=nDNext_;
    tipDeflectionCal();

    Foam::Info<<"Linear Solver Used!"<<Foam::endl;

}

void ALFBM::fETurbine::initialSolve()
{
	sT_=initialSolver;
    for(int i=0;i<50;++i)
    {
        nPInitial();
        turbineBMA();
        stressStiffenSolve();
        equivalentKCal();
        boundaryApply();
        equivalentPCal();
        deformationSolve();
        restorePCal();
        rootForceCal();
        resultIteration();
        tipDeflectionCal();
    }
    turbineBMA();
    stressStiffenSolve();
    modalSolve();

    Foam::Info<<"Initial Solver Used!"<<Foam::endl;
}

void ALFBM::fETurbine::staticInterationSolve()
{
    sT_=staticSolver;
    nP_=rigidNP_;
    nodeReadInitial();
    nPInitial();
    turbineBMA();
    stressStiffenSolve();
    for(int i=0;i<50;++i)
    {
        nPInitial();
        turbineBMA();
        equivalentKCal();
        boundaryApply();
        equivalentPCal();
        deformationSolve();
        rootForceCal();
        resultIteration();
        tipDeflectionCal();
    }
    turbineBMA();
    stressStiffenSolve();
    modalSolve();

    Foam::Info<<"Static Interation Solver Used!"<<Foam::endl;
}

void ALFBM::fETurbine::implicitIterationSolve()
{
	sT_=implicitSolver;
    nPInitial();
    turbineBMA();
    stressStiffenSolve();
    equivalentKCal();
    boundaryApply();
    equivalentPCal();
    deformationSolve();
    restorePCal();
    rootForceCal();
    resultIteration();
    tipDeflectionCal();
    modalSolve();

    Foam::Info<<"Explicit Interation Solver Used!"<<Foam::endl;
}

void ALFBM::fETurbine::readStructureResults(std::ifstream & structuredata)
{
    readStructureResult(structuredata, "RigidNodePosition", nodeNumber_, rigidNP_);
    readStructureResult(structuredata, "NodePosition", nodeNumber_, nP_);
    readStructureResult(structuredata, "LoadLast", nodeNumber_, loadLast_);
    readStructureResult(structuredata, "NodeDisplacementNext", nodeNumber_, nDNext_);
    readStructureResult(structuredata, "NodeDot", nodeNumber_, nDDot_);
    readStructureResult(structuredata, "NodeDualDot", nodeNumber_, nDDualDot_);
}

void ALFBM::fETurbine::writeStructureResults(std::ofstream & structuredata)
{
    writeStructureResult(structuredata, "RigidNodePosition", nodeNumber_, rigidNP_);
    writeStructureResult(structuredata, "NodePosition", nodeNumber_, nP_);
    writeStructureResult(structuredata, "LoadLast", nodeNumber_, loadLast_);
    writeStructureResult(structuredata, "NodeDisplacementNext", nodeNumber_, nDNext_);
    writeStructureResult(structuredata, "NodeDot", nodeNumber_, nDDot_);
    writeStructureResult(structuredata, "NodeDualDot", nodeNumber_, nDDualDot_);

    // !!! deflection !!! //
}

void ALFBM::fETurbine::writeModalResults(std::ofstream & modaldata)
{
    if(flagBit_.modalSolve())
    {
        for(int i=0;i<15;++i)
        {
            if(flagBit_.damp())
            {
                modaldata << "Frequence," << modalFrequence_(modalFrequence_.size()-1-i,0).imag()/(2*Foam::constant::mathematical::pi)<<std::endl;
                for(int j=0;j<nodeNumber_;++j)
                {
                    modaldata << modalVector_(6*j,modalFrequence_.size()-1-i).real()<<","
                                <<modalVector_(6*j+1,modalFrequence_.size()-1-i).real()<<","
                                <<modalVector_(6*j+2,modalFrequence_.size()-1-i).real()<<","
                                <<modalVector_(6*j+3,modalFrequence_.size()-1-i).real()<<","
                                <<modalVector_(6*j+4,modalFrequence_.size()-1-i).real()<<","
                                <<modalVector_(6*j+5,modalFrequence_.size()-1-i).real()<<","<<std::endl;
                }
            }
            else
            {
                modaldata << "Frequence," << sqrt(modalFrequence_(modalFrequence_.size()-1-i,0)).real()/(2*Foam::constant::mathematical::pi)<<std::endl;
                for(int j=0;j<nodeNumber_;++j)
                {
                    modaldata << modalVector_(6*j,modalFrequence_.size()-1-i).real()<<","
                                <<modalVector_(6*j+1,modalFrequence_.size()-1-i).real()<<","
                                <<modalVector_(6*j+2,modalFrequence_.size()-1-i).real()<<","
                                <<modalVector_(6*j+3,modalFrequence_.size()-1-i).real()<<","
                                <<modalVector_(6*j+4,modalFrequence_.size()-1-i).real()<<","
                                <<modalVector_(6*j+5,modalFrequence_.size()-1-i).real()<<","<<std::endl;
                }
            }
        }
    }
}

void ALFBM::fETurbine::writeForceResults(std::ofstream & forcedata)
{
    writeForceResult(forcedata, "towerRootForce", 0, 1, rootForce_[0]);
    writeForceResult(forcedata, "nacelleRootForce", 0, 1, rootForce_[1]);
    
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        writeForceResult(forcedata, "globalRootForce", i, 1, rootForce_[i+2]);
    }
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        writeForceResult(forcedata, "localRootForce", i, 1, rootForce_[i+2+turbineInfo_.bladeNumber()]);
    }
}

void ALFBM::fETurbine::writeDeflectionResults(std::ofstream & deflectiondata)
{
    writeForceResult(deflectiondata, "towerdeflection", 0, 1, tipDeflection_[0]);
    writeForceResult(deflectiondata, "nacelledeflection", 0, 1, tipDeflection_[1]);
    
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        writeForceResult(deflectiondata, "globaldeflection", i, 1, tipDeflection_[i+2]);
    }
    for(int i=0;i<turbineInfo_.bladeNumber();++i)
    {
        writeForceResult(deflectiondata, "localdeflection", i, 1, tipDeflection_[i+2+turbineInfo_.bladeNumber()]);
    }
}



#endif
