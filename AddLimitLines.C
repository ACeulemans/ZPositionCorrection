#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TAxis.h"
#include "TRandom.h"
#include "TString.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include <Math/RootFinderAlgorithms.h>
#include <TF1.h>
#include <TGraph.h>
#include <Math/RootFinder.h>
#include <Math/Functor.h>

using namespace ROOT::Math;


//Source/input/output options  (Don't forget to change name of output file as well if you change this)
TString DataFileName = "../IS587/R18-114_autocal_AllAlpha_newTresholds_sum_hists_13_02_23.root"; //name of input root file
TString TH2FileName = "EBISMode/Ex_vs_theta_ebis"; //name of Ex vs theta_CM th2 inside the datafile to open, for files from iss_sort
TString OutputFileName = "TestLimitLines.root"; //Name for main output file

//reaction information
//68Ni(d,p)
Double_t BeamEnergy = 408; //Beam kinetic energy in MeV
Double_t BeamMass = 67.93187; //Beam particle mass in amu
Double_t TargetMass = 2.01410; //Target mass in amu
Double_t EjectileMass = 1.00782; //Ejectile mass in amu
Double_t RecoilMass = 68.93561; //Recoil mass in amu (for ground state)
Int_t EjectileCharge = 1; //Charge of ejectile in multiples of electron charge
Int_t RecoilCharge = 28; //Charge of ejectile in multiples of electron charge
Double_t MagField = 2.05; //Magnetic field strength in Tesla

//Options for limit line calculations
//Lines are calculated for NrEnergiesForGraph energies of ejectile between MinEnergy and MaxEnergy
Int_t NrEnergiesForGraph = 1000; //Nr of points in limit lines TGraphs
Double_t MinEnergy = 0.;  //Minimum energy for ejectile limit calculations
Double_t MaxEnergy =  10.; //Maximum energy for ejectile limit calculations

//detector information
Double_t SiliconRadius = 28.75; //radius of si-array in mm (used for z-correction)
Double_t Edges[8] = {-61.9,-186.9,-187.4,-312.4,-312.9,-437.9,-438.4,-563.4}; //Edges of the silicon array wrt to the target in mm
Double_t DistToArray = 61.9; //Distance target-array in mm

//Drawing options 
Double_t MinCMAngleAxis = 0.; //Lower limit of com angle range (in degree) (used to set setrangeuser)
Double_t MaxCMAngleAxis = 60.; //Upper limit of com angle range (in degree) (used to set setrangeuser)
Double_t MinExEnAxis = -2.; //Lower limit of excitation energy range (in MeV) (used to set setrangeuser)
Double_t MaxExEnAxis = 10.; //Upper limit of excitation energy range (in MeV) (used to set setrangeuser)
Bool_t CreateOverlay = true; //If true draws the limit lines on top of the histogram, if false only plots the limit lines

//Add excitation energy list, the code will print the detector angular limits for these excitation energies
std::vector<Double_t> ExEnergiesList = {0, 0.300, 1.770, 2.500, 2.800, 3.100, 4.000, 4.350, 4.700}; //Excitation energies in MeV (for recoil particle)
std::vector<Double_t> LowerCoMbounds; //stores lowest CoM angle that can be detected for these excitation energies
std::vector<Double_t> UpperCoMbounds; //stores highest CoM angle that can be detected for these excitation energies


//Conversion factors (don't change)
Double_t amutoMeVperc2 = 931.5; // 1 amu = 931.5 MeV/c^2 
Double_t c = 299792458; //speed of light in m/s
Double_t eVperMeV = 1000000; //1MeV = 1 000 000 eV
Double_t nspers = 1000000000; //1s = 1 000 000 000 ns
Double_t mpermm = 0.001; // 1mm = 0.001 m





/*!  \brief Calculates velocity of com frame
*    
*    Returns the velocity of com frame as a fraction of the speed of light.
*/
Double_t GetCOMVelocity() {
    Double_t beamvel = sqrt( 2 * BeamEnergy / (BeamMass*amutoMeVperc2) ); //velocity of beam particle (in fraction of c)
    Double_t beammomentum = BeamMass * beamvel * amutoMeVperc2; //Beam momentum in MeV/c
    Double_t comvel = beammomentum / ((BeamMass+TargetMass)*amutoMeVperc2); //velocity of center of mass frame in fraction of c
    //std::cout << "The center of mass velocity is " << comvel << " c" << std::endl;
    return comvel;
}

/*!  \brief Calculates cyclotron period
*    
*    Returns the cyclotron period of the ejectile (in ns).
*/
Double_t GetCyclotronPeriod () {
    Double_t Cper = 2 * TMath::Pi() * EjectileMass * amutoMeVperc2 * eVperMeV *nspers / ( MagField * EjectileCharge * c * c); //cyclotron radius in ns
    //std::cout << "The cyclotron period is " << Cper << " ns" << std::endl;
    return Cper;
}

/*!  \brief Calculates ejectile lab velocity
*    
*    Returns ejectile velocity in fraction of c
*    @param LabEn lab energy of ejectile (in MeV)
*/
Double_t GetEjectileLabVelocity(Double_t LabEn) {

    //std::cout << "The lab energy of the ejectile is " << LabEn << " MeV" << std::endl;
    Double_t ejeclabvel = std::sqrt(2*LabEn/(EjectileMass*amutoMeVperc2)); //ejectile lab velocity in fraction of c
    //std::cout << "The lab velocity of the ejectile is " << ejeclabvel << " c" << std::endl;
    return ejeclabvel;
}

/*!  \brief Calculates ejectile lab angle
*    
*    Returns ejectile velocity in degree
*    @param LabVel lab velocity of ejectile (in fraction of c)
*    @param z z distance for a single loop (measured z on array corrected to x=y=0 axis) (in mm)
*    @param TCyc cyclotron period of the ejectile (in ns)
*/
Double_t GetEjectileLabAngle(Double_t LabVel, Double_t z, Double_t TCyc) {
    Double_t ejeclabangle = acos((z*mpermm)/(TCyc/nspers*LabVel*c));
    //std::cout << "The lab velocity of the ejectile is " << ejeclabvel << " c" << std::endl;
    return ejeclabangle;
}

/*!  \brief Calculates energy in the initial com frame
*    
*    Returns the energy in the initial com frame (in MeV)
*/
Double_t GetCOMEnergyInit() {
    Double_t ecomin = BeamEnergy * (TargetMass/(TargetMass+BeamMass)); //energy in initial com frame in MeV
    //std::cout << "The energy in the initial CoM frame is " << ecomin << " MeV" << std::endl;
    return ecomin;
}

/*!  \brief Calculates Q-value for reaction
*    
*    Returns the Q-value for reaction to ground state (in MeV)
*/
Double_t GetQValue() {
    Double_t qval = (BeamMass+TargetMass-EjectileMass-RecoilMass)*amutoMeVperc2; //q-value for reaction to GS in MeV
    //std::cout << "The Q-value for the reaction is " << qval << " MeV" << std::endl;
    return qval;
}

/*!  \brief Calculates ejectile com velocity
*    
*    Returns ejectile velocity in CM frame (in fraction of c)
*    @param vlab ejectile lab velocity (in fraction of c)
*    @param thetalab ejectile lab angle (in radian)
*    @param VCM center of mass velocity (in fraction of c)
*/
Double_t GetEjectileCMVelocity(Double_t vlab, Double_t thetalab, Double_t VCM) {

    //we calculate both components of the velocity in the cm frame
    Double_t V_cm_ejec_z = vlab * cos(thetalab) - VCM;
    Double_t V_cm_ejec_y = vlab * sin(thetalab);


    Double_t V_cm_ejec = sqrt(V_cm_ejec_z * V_cm_ejec_z + V_cm_ejec_y * V_cm_ejec_y);
    //std::cout << "The CM velocity of the ejectile is " << V_cm_ejec << " c" << std::endl;
    return V_cm_ejec;
}

/*!  \brief Calculates ejectile com angle
*    
*    Returns ejectile angle in CM frame (in radian)
*    @param vlab ejectile lab velocity (in fraction of c)
*    @param v_cm_ejec ejectile center of mass velocity (in fraction of c)
*    @param ThetaLab ejectile lab angle (in radian)
*/
Double_t GetEjectileCOMAngle(Double_t v_lab, Double_t v_cm_ejec, Double_t ThetaLab){
    
    //we calculate the y component of the lab velocity, which is equal to the y component of the cm veloctity of the ejectile
    Double_t V_cm_ejec_y = v_lab * sin(ThetaLab);


    Double_t ThetaCM = asin(V_cm_ejec_y/v_cm_ejec);
    //std::cout << "The CM angle of the ejectile is " << ThetaCM*180/TMath::Pi() << " degree" << std::endl;
    return ThetaCM;
}

/*!  \brief Calculates ejectile com energy
*    
*    Returns ejectile energy in CM frame (in MeV)
*    @param v_cm_ejec ejectile center of mass velocity (in fraction of c)
*/
Double_t GetEjectileCMEnergy(Double_t v_cm_ejec){

    Double_t Energy_ejec_cm = EjectileMass* amutoMeVperc2*  v_cm_ejec * v_cm_ejec / 2 ; 

    //std::cout << "The CM energy of the ejectile is " << Energy_ejec_cm << " MeV" << std::endl;
    return Energy_ejec_cm;
}

/*!  \brief Calculates recoil com energy
*    
*    Returns recoil energy in CM frame (in MeV)
*    @param v_cm_ejec ejectile center of mass velocity (in fraction of c)
*/
Double_t GetRecoilCMEnergy(Double_t v_cm_ejec){

    //First calculate the momentum of the ejectile
    Double_t Mom_ejec_cm = EjectileMass* amutoMeVperc2 * v_cm_ejec;


    Double_t Energy_rec_cm = Mom_ejec_cm * Mom_ejec_cm / (RecoilMass* amutoMeVperc2)/2; 
    //std::cout << "The CM energy of the recoil is " << Energy_rec_cm << " MeV" << std::endl;
    return Energy_rec_cm;
}

/*!  \brief Calculates z correction
*    
*    Returns reconstructed z intersection with x=y=0 in mm
*    @param z z distance for a single loop (measured z on array without correction) (in mm)
*    @param LabAngleEjec lab angle ejectile (in radian)
*    
*/
Double_t GetCorrectedZValue(Double_t z, Double_t LabAngleEjec) {
    if (std::isnan(LabAngleEjec)){
        return 0;
    }
    Double_t z_corrected = z + SiliconRadius/tan(LabAngleEjec);
    //std::cout << "The lab velocity of the ejectile is " << ejeclabvel << " c" << std::endl;
    return z_corrected;
}

struct my_f_params { double z_det; double h_det; double v_lab; double TCyc; };

my_f_params Current_Parameters = { 0., 0., 0., 0. };

double my_f (double x){

    //struct my_f_params * params = (struct my_f_params *)p;
    double Zd = (Current_Parameters.z_det);
    double Hd = (Current_Parameters.h_det);
    double Vl = (Current_Parameters.v_lab);
    double Tc = (Current_Parameters.TCyc);
    
    //std::cout << "Zd is: " << Zd << std::endl;
    //std::cout << "Hd is: " << Hd << std::endl;
    //std::cout << "Vl is: " << Vl << std::endl;
    //std::cout << "Tc is: " << Tc << std::endl;
    //std::cout << "Vl*c*Tc/nspers is: " << Vl*c*Tc/nspers << std::endl;
    //std::cout << "function at pi/2+0.001 is: " << cos(TMath::Pi()/2+0.001) - (1 / (Vl*c*Tc/nspers)) * (Zd + Hd / tan(TMath::Pi()/2+0.001))*mpermm << std::endl;
    //std::cout << "function at pi-0.1 is: " << cos(TMath::Pi()-0.1) - (1 / (Vl*c*Tc/nspers)) * (Zd + Hd / tan(TMath::Pi()-0.1))*mpermm << std::endl;

    return  cos(x) - (1 / (Vl*c*Tc/nspers)) * (Zd + Hd / tan(x))*mpermm ;
  }

/*!  \brief Checks if proton would have hit array side
*    
*    Some of the events in the Edep vs Z histogram do not correspond to a physical proton on array.
*    This function helps eliminate false events that should have hit the side of the array
*    or better those 'protons' with an energy too low to reach the detected z position.
*    Returns true if it should have hit side (and is thus a fake event)
*    @param LabAngleEjec lab angle ejectile (in radian)
*    
*/
Bool_t CheckIfHitsArraySide(Double_t LabAngleEjec){
    
    if(DistToArray*tan(TMath::Pi()-LabAngleEjec) <= SiliconRadius){
        //std::cout << "False proton event detected" << std::endl;
        return true;
    }
    else{
        return false;
    }

}

/*!  \brief Calculates corrected lab angle of ejectile
*    
*    Returns reconstructed lab angle of ejectile (in radian)
*    Check https://root-forum.cern.ch/t/solve-equation-with-rootfinder/42286/3
*    @param Z_det z distance for a single loop (measured z on array without correction) (in mm)
*    @param v_lab_ej ejectile lab velocity (in fraction of c)
*    @param TCyc cyclotron period of ejectile (in ns)
*    
*/
Double_t CalcCorrectedLabAngle(Double_t Z_det, Double_t v_lab_ej, Double_t TCyc){
    
    Double_t LabAngle_corr = -1;

    //update parameters
    Current_Parameters.z_det = Z_det;
    Current_Parameters.h_det = SiliconRadius;
    Current_Parameters.v_lab = v_lab_ej;
    Current_Parameters.TCyc = TCyc;

    //std::cout << "z_det is: " << Z_det << std::endl;
    //std::cout << "h_det is: " << SiliconRadius << std::endl;
    //std::cout << "v_lab is: " << v_lab_ej << std::endl;
    //std::cout << "TCyc is: " << TCyc << std::endl;

    //Set up the solver
    RootFinder *k = new RootFinder();
    k->SetMethod(RootFinder::kGSL_BRENT);
    ROOT::Math::Functor1D f(&my_f);

    //The function is always positive just after Pi/2, we need to find a negative endpoint with Pi/2<EndPoint<Pi
    double EndPoint = {TMath::Pi()};
    double TestIsNegative = {1};
    for (int i=0; i<10000; i++){
        TestIsNegative = my_f(TMath::Pi()/2 + TMath::Pi()/2*i/10000);
        if (TestIsNegative<0){
            EndPoint = TMath::Pi()/2 + TMath::Pi()/2*i/10000;
            break;
        }
    }
    //std::cout << "Used endpoint is "<< EndPoint <<endl;
    //std::cout << "Function value at endpoint is "<< my_f(EndPoint) <<endl;

    //Activate solver
    gSystem->RedirectOutput("/dev/null"); //Remove error printing (change if you are not on linux)
    k->SetFunction(f, TMath::Pi()/2+0.001, EndPoint);
    k->Solve();
    gSystem->RedirectOutput(0,0); //Restore normal output to terminal

 
    double c = -1;
    if(k->Status()==0){  //0 means no errors
        c = k->Root();
    }
    else{ //no good solution was found
        c = -1;
    }
    
    //std::cout << "Lab angle before correction is "<< v_lab_ej*180/TMath::Pi() <<endl;
    //std::cout << "New Lab angle is "<< c*180/TMath::Pi() <<endl;

    return c;
}

double GetCorrectedZValue(Double_t LabAngle, Double_t v_lab_ej, Double_t TCyc){
    return cos(LabAngle)*v_lab_ej*c*TCyc/nspers/mpermm;
}

/*!  \brief Calculates minimum and maximum com angle for each Ex
*    
*    These angles are maximal array hitting com angle (lower limit)
*    and the angle corresponding to the start of the array
*    @param SideArrGR TGraph for side hitting angle
*    @param StartArrGR TGraph for the array
*/
void GetAngularLimits(TGraph* SideArrGR, TGraph* StartArrGR, std::vector<Double_t> Energies){


    //Get number of excitation energies
    Int_t NrEnergies = Energies.size();
    
    for (Int_t i=0; i<NrEnergies; i++){
        Double_t CurrentExEnergyinKeV = Energies[i]*1000; //convert to keV
        LowerCoMbounds.push_back(-1);
        UpperCoMbounds.push_back(-1);
        //We find the closest matching energies
        Int_t ClosestIndexSide = 0;
        Int_t ClosestIndexStart = 0;
        Double_t MinEnergyDiffSide = 1000;  //large enough value so that we can find smallest
        Double_t MinEnergyDiffStart = 1000;
        Double_t CoMcurr = -1;
        Double_t ExEnCurr = -10;
        //Determine lower bound from side hitting graph
        for(Int_t j=0; j<SideArrGR->GetN(); j++){
            SideArrGR->GetPoint(j, CoMcurr, ExEnCurr);
            if(abs(ExEnCurr-CurrentExEnergyinKeV)<MinEnergyDiffSide){
                MinEnergyDiffSide = abs(ExEnCurr-CurrentExEnergyinKeV);
                LowerCoMbounds[i]=CoMcurr;
            }
        }
        //same for start array graph and upper bound
        for(Int_t j=0; j<StartArrGR->GetN(); j++){
            StartArrGR->GetPoint(j, CoMcurr, ExEnCurr);
            if(abs(ExEnCurr-CurrentExEnergyinKeV)<MinEnergyDiffStart){
                MinEnergyDiffStart = abs(ExEnCurr-CurrentExEnergyinKeV);
                UpperCoMbounds[i]=CoMcurr;
            }
        }
        std::cout << "For an excitation energy of " << Energies[i] << " MeV, the angular bounds are [" << LowerCoMbounds[i] << " ; " << UpperCoMbounds[i] << "]" << std::endl;
    }

}

/*!  \brief Calculates maximal com angle that hits side of array
*    
*    Returns maximal array hitting com angle
*    Used to determine limits of angular distribution
*    @param Elab Energy deposited in detector (in MeV)
*    @param ExEn Calculated excitation energy of recoil (in MeV), used only to pass output
*    @param SolutionFound Bool used to check wheter a valid solution is found, passed by reference
*/
Double_t GetMaxSideAnglevsEx(Double_t Elab, Double_t &ExEn, Bool_t &SolutionFound) {

//These four are only dependent on the reaction type
    Double_t TCyclotron = GetCyclotronPeriod();
    Double_t CMvelocity = GetCOMVelocity();
    Double_t CMEnergy_in = GetCOMEnergyInit();
    Double_t QValue = GetQValue();

    //The lab velocity of the proton can be calculated from the lab energy
    Double_t v_lab = GetEjectileLabVelocity(Elab);

    //Now we calculate the lab angle using the brent-decker algorithm
    Double_t ThetaLab = TMath::Pi()-atan2(SiliconRadius,DistToArray);


    //Define Theta_cm
    Double_t ThetaCM = -1;

    if(ThetaLab > 0){ //Theta lab solution was succesfully found

        SolutionFound = true;
    
        //The lab velocity combined with the angle can be converted to a center of mass angle
        Double_t v_cm_ejec = GetEjectileCMVelocity(v_lab, ThetaLab, CMvelocity);

        //Now the center of mass angle can be computed 
        ThetaCM = GetEjectileCOMAngle(v_lab, v_cm_ejec, ThetaLab);

        //Finally we also calculate the excitation energy, by first getting ejectile and recoil energies in CM frame
        Double_t EjectileCMEnergy = GetEjectileCMEnergy(v_cm_ejec);
        Double_t RecoilCMEnergy = GetRecoilCMEnergy(v_cm_ejec);

        Double_t ExcitationEnergy = CMEnergy_in - EjectileCMEnergy - RecoilCMEnergy + QValue;

        //Assign variables that are passed by reference
        ExEn = ExcitationEnergy;

    }

    else{ //No solution for theta lab was found
        SolutionFound = false;
        ExEn = -99;
        ThetaCM = -1;
    }

    return ThetaCM;

}


/*!  \brief Calculates ejectile com angle
*    
*    Returns ejectile com angle in degree
*    @param z_det measured z distance on array (in mm)
*    @param Elab Energy deposited in detector (in MeV)
*    @param Z_corr corrected z to axis (x=y=0)(in mm), used only to pass output
*    @param ExEn Calculated excitation energy of recoil (in MeV), used only to pass output
*    @param SolutionFound Bool used to check wheter a valid solution is found, passed by reference
*/
Double_t SolveReaction(Double_t z_det, Double_t Elab, Double_t &Z_corr, Double_t &ExEn, Bool_t &SolutionFound) {

    //These four are only dependent on the reaction type
    Double_t TCyclotron = GetCyclotronPeriod();
    Double_t CMvelocity = GetCOMVelocity();
    Double_t CMEnergy_in = GetCOMEnergyInit();
    Double_t QValue = GetQValue();

    //The lab velocity of the proton can be calculated from the lab energy
    Double_t v_lab = GetEjectileLabVelocity(Elab);

    //Now we calculate the lab angle using the brent-decker algorithm
    Double_t ThetaLab = CalcCorrectedLabAngle(z_det, v_lab, TCyclotron);

    //Check if this particle would hit the side of the array
    Bool_t HitsSide = CheckIfHitsArraySide(ThetaLab);

    //Define Theta_cm
    Double_t ThetaCM = -1;

    if((ThetaLab > 0) && (!HitsSide)){ //Theta lab solution was succesfully found

        SolutionFound = true;

        //Now we correct the z-value as it would be at the center line (x=y=0)
        Double_t z_corr = GetCorrectedZValue(z_det, ThetaLab);
    
        //The lab velocity combined with the angle can be converted to a center of mass angle
        Double_t v_cm_ejec = GetEjectileCMVelocity(v_lab, ThetaLab, CMvelocity);

        //Now the center of mass angle can be computed 
        ThetaCM = GetEjectileCOMAngle(v_lab, v_cm_ejec, ThetaLab);

        //Finally we also calculate the excitation energy, by first getting ejectile and recoil energies in CM frame
        Double_t EjectileCMEnergy = GetEjectileCMEnergy(v_cm_ejec);
        Double_t RecoilCMEnergy = GetRecoilCMEnergy(v_cm_ejec);

        Double_t ExcitationEnergy = CMEnergy_in - EjectileCMEnergy - RecoilCMEnergy + QValue;

        //Assign variables that are passed by reference
        Z_corr = z_corr;
        ExEn = ExcitationEnergy;

    }

    else{ //No solution for theta lab was found
        SolutionFound = false;
        Z_corr = 1000;
        ExEn = -99;
        ThetaCM = -1;
    }

    
    //std::cout << "The cyclotron period is " << TCyclotron << " ns" << std::endl;
    //std::cout << "The velocity of the CM frame is " << CMvelocity << " c" << std::endl;
    //std::cout << "The total energy in the initial CM frame is " << CMEnergy_in << " MeV" << std::endl;
    //std::cout << "The Q-value of the reaction is " << QValue << " MeV" << std::endl;
    //std::cout << "The total energy in the outgoing CM frame is " << CMEnergy_out << " MeV" << std::endl;
    //std::cout << "The ejectile CM velocity is " << v_cm_ejec << " c" << std::endl;
    //std::cout << "The detected z is " << z_det << " mm" << std::endl;
    //std::cout << "The CM angle is " << ThetaCM*180/TMath::Pi() << " degrees" << std::endl;
    //std::cout << "The lab angle is " << ThetaLab*180/TMath::Pi() << " degrees" << std::endl;
    //std::cout << "The corrected lab angle is " << ThetaLab_corr*180/TMath::Pi() << " degrees" << std::endl;
    //std::cout << "The corrected z is " << z_corr << " mm" << std::endl;
    //std::cout << "The single step corrected CM angle is " << ThetaCM_corr*180/TMath::Pi() << " degrees" << std::endl;
    //std::cout << "The brent step corrected z is " << z_corr_2 << " mm" << std::endl;
    //std::cout << "The brent corrected CM angle is " << ThetaCM_corr_2*180/TMath::Pi() << " degrees" << std::endl;
    //std::cout << "The ejectile com energy is " << EjectileCMEnergy << " MeV" << std::endl;
    //std::cout << "The recoil com energy is " << RecoilCMEnergy << " MeV" << std::endl;
    //std::cout << "The excitation energy is " << ExcitationEnergy << " MeV" << std::endl;

    return ThetaCM;
}


void AddLimitLines(){

  time_t t;
  time(&t);
  long double TimeInit = t;
  
  //Reading data file. 

  TString DataPath = DataFileName;

  //Check if input file is present and readable
  if (gSystem->AccessPathName(DataPath.Data())){
	std::string DataPathString(DataPath.Data());
    std::string ErrMsg = "Error with opening input file " + DataPathString + ", check if name and location are correct";
    throw std::runtime_error(ErrMsg);
  }
  //read input file
  TFile *InputFile = new TFile(DataPath,"READ");

  TH2F *Ex_thetaCM_Hist=(TH2F*)InputFile->Get(TH2FileName);
  if (!Ex_thetaCM_Hist){
      std::string TH2NameStd(TH2FileName); //convert TString to std::string
      std::string ErrMsg = "Histogram ( " + TH2NameStd + " ) not found inside file, check file content";
      throw std::runtime_error("Histogram not found inside file, check file content");
  }

  Int_t NrBinsX = Ex_thetaCM_Hist->GetNbinsX();
  Int_t NrBinsY = Ex_thetaCM_Hist->GetNbinsY();
  //std::cout << "Opened histogram contains " << NrBinsX << " x bins and " << NrBinsY << " y bins" << std::endl;

  
  //We create vectors to store the variables
  std::vector<Double_t> LeftLimitsCoMangle {};
  std::vector<Double_t> RightLimitsCoMangle {};
  std::vector<Double_t> StartArrayExEn {};
  std::vector<Double_t> StartArrayCoMangle {};
  std::vector<Double_t> EndArrayExEn {};
  std::vector<Double_t> EndArrayCoMangle {};
  std::vector<Double_t> SideArrayExEn {};
  std::vector<Double_t> SideArrayCoMangle {};


  //We also calculate the corresponding cm angles of the gaps in the si detector (for each energy) 
  for (Int_t n=0; n<=NrEnergiesForGraph; n++) {
       std::array<double,3> LeftLimits = {0.,0.,0.};
       std::array<double,3> RightLimits = {0.,0.,0.};
       std::array<double,3> ThisExEn = {0.,0.,0.};
       Double_t StartArrayCMa = 0;
       Double_t StartArrayEx = 0;
       Double_t EndArrayEx = 0;
       Double_t EndArrayCMa = 0;
       

       //Calculate max com angles that hit side
       Double_t CoMangleSide = -1;
       Double_t ExEnSide = -1;
       Bool_t SideSolFound = false;
       Double_t CurrentEnergy = MinEnergy + n*(MaxEnergy-MinEnergy)/NrEnergiesForGraph;
       CoMangleSide = GetMaxSideAnglevsEx(CurrentEnergy, ExEnSide, SideSolFound);  

       for(Int_t j=0; j<8; j++){
          
          Double_t ZCenter = Edges[j];
          Double_t ExEner = 0;
          Double_t ZCorrected = 0;
          Double_t CoMAnglecenter = 0;
          
          //std::cout << "n is "  << n << " j is " << j << std::endl;
          Bool_t GoodSolutionFound = false;
          CoMAnglecenter = SolveReaction(ZCenter, CurrentEnergy, ZCorrected, ExEner, GoodSolutionFound);  
          if(j == 1 || j == 3 || j==5){ //'left' sides of gaps
              int IndexInt = (int) (floor(j/2));
              LeftLimits[IndexInt] = CoMAnglecenter;
          }
          if(j == 2 || j == 4 || j==6){ //'right' sides of gaps
              int IndexInt = (j/2)-1;
              RightLimits[IndexInt] = CoMAnglecenter;
              ThisExEn[IndexInt] = ExEner;
          }
          if(j == 0){ //start of array
              StartArrayCMa = CoMAnglecenter;
              StartArrayEx = ExEner;
          }
          if(j == 7){ //end of array
              EndArrayCMa = CoMAnglecenter;
              EndArrayEx = ExEner;
          }
           
       }

       StartArrayExEn.push_back(StartArrayEx);
       StartArrayCoMangle.push_back(StartArrayCMa);
       EndArrayExEn.push_back(EndArrayEx);
       EndArrayCoMangle.push_back(EndArrayCMa);
       SideArrayExEn.push_back(ExEnSide);
       SideArrayCoMangle.push_back(CoMangleSide); 

  }


  //Create a tMultigraph to store limits of detector graphs 
  TMultiGraph* MG_StartEndArray = new TMultiGraph();
  MG_StartEndArray->SetTitle("Start and End of Array CoM angles for differen excitation energies; Center of Mass Angle [degree]; Excitation Energy [keV]");

  //TGraph for start of array
  TGraph* StartArrayGR = new TGraph();
  StartArrayGR->SetTitle("Start of array; Center of Mass Angle [degree]; Excitation Energy [keV]");
  StartArrayGR->SetMarkerStyle(20);
  StartArrayGR->SetMarkerColor(kBlue);

  //TGraph for end of array
  TGraph* EndArrayGR = new TGraph();
  EndArrayGR->SetTitle("End of array; Center of Mass Angle [degree]; Excitation Energy [keV]");
  EndArrayGR->SetMarkerStyle(20);
  EndArrayGR->SetMarkerColor(kRed);

  //TGraph for side of array
  TGraph* SideArrayGR = new TGraph();
  SideArrayGR->SetTitle("Side of array; Center of Mass Angle [degree]; Excitation Energy [keV]");
  SideArrayGR->SetMarkerStyle(20);
  SideArrayGR->SetMarkerColor(kGreen);

  //We fill the tgraphs
  for (Int_t j=0; j<StartArrayCoMangle.size(); j++) {
      StartArrayGR->SetPoint(j, StartArrayCoMangle[j]*180/TMath::Pi(), StartArrayExEn[j]*1000); //Convert energy to keV for drawing on hist
      EndArrayGR->SetPoint(j, EndArrayCoMangle[j]*180/TMath::Pi(), EndArrayExEn[j]*1000);
      SideArrayGR->SetPoint(j, SideArrayCoMangle[j]*180/TMath::Pi(), SideArrayExEn[j]*1000);
  }
  
  //Add graphs to multigraph
  MG_StartEndArray->Add(StartArrayGR,"px");
  MG_StartEndArray->Add(EndArrayGR,"px same");
  MG_StartEndArray->Add(SideArrayGR,"px same");

  //Create legend for multigraph
  TLegend *MG_Legend = new TLegend(0.6,0.7,0.9,0.9);	
  MG_Legend->SetTextSize(0.03);
  MG_Legend->AddEntry(StartArrayGR,"Start of Array","p");
  MG_Legend->AddEntry(EndArrayGR,"End of Array","p");
  MG_Legend->AddEntry(SideArrayGR,"Side of Array","p");

  //Create canvas for multigraph
  TCanvas *c_MG=new TCanvas("c_MG", "Canvas For Multigraph");

  //Draw graph over existing Ex vs theta_CM hist
  c_MG->cd();
  Ex_thetaCM_Hist->SetMinimum(0.);
  MG_StartEndArray->Draw("a");
  if(CreateOverlay) Ex_thetaCM_Hist->Draw("colz1 same");
  MG_StartEndArray->GetXaxis()->SetRangeUser(MinCMAngleAxis,MaxCMAngleAxis);
  MG_StartEndArray->GetYaxis()->SetRangeUser(MinExEnAxis*1000.,MaxExEnAxis*1000.);
  MG_Legend->Draw();

  //Get the angular limits for the given excitation energies
  GetAngularLimits(SideArrayGR, StartArrayGR, ExEnergiesList);

  //Save output histogram
  TFile *MyHistFile = new TFile(OutputFileName,"RECREATE");
  MyHistFile->cd();
  c_MG->Write();
  MyHistFile->Close();

  time(&t);
  long double TimeFin = t;
  std::cout << "Script was processed in " << (int)floor((TimeFin-TimeInit)/60) << " min and " << ((int)(TimeFin-TimeInit))%60 << " sec" << std::endl; 

}
