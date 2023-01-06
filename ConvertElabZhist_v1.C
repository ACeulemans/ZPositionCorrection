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
#include <Math/RootFinder.h>
#include <Math/Functor.h>

using namespace ROOT::Math;


//Source/input options  (Don't forget to change name of function/file as well if you change this)
//TString DataFileName = "R18-114_autocal_AllAlpha_sum_hists.root"; //file to read
TString DataFileName = "RGallium_resort_newcalib_notreshold_hists.root";
//TString DataFileName = "../Scripts/OUTPUT/R18-114_resortAllAlpha_sum_hists_RGallium_resort_newcalib_hists_E_vs_z_bgsub.root";

TString TH2FileName = "EBISMode/E_vs_z_ebis"; //name of th2 inside the datafile to open
//TString TH2FileName = "E_vs_Zdet_Ni_bgsub";

//reaction information
//68Ni(d,p)
Double_t BeamEnergy = 408; //Beam kinetic energy in MeV
Double_t BeamMass = 67.93187; //Beam particle mass in amu
Double_t TargetMass = 2.01410; //Target mass in amu
Double_t EjectileMass = 1.00782; //Ejectile mass in amu
Double_t RecoilMassGS = 68.93561; //Recoil mass in amu for GS
Double_t RecoilMass = RecoilMassGS; //Mass used will be overwritten to add excitation energy 
//input for ejectile radius calculation
Int_t EjectileCharge = 1; //Charge of ejectile in multiples of electron charge
Int_t RecoilCharge = 28; //Charge of ejectile in multiples of electron charge
std::string OutputFileName = "Ni68_d,p_Ni69_6MeV_2,05T"; 
Double_t MagField = 2.05; //Magnetic field strength in Tesla

//detector information
Double_t SiliconRadius = 28.75; //radius of si-array in mm (used for z-correction)
Double_t Edges[8] = {-61.9,-186.9,-190.4,-315.4,-318.9,-443.9,-447.4,-572.4}; //Edges of the silicon array

//Conversion factors
Double_t amutoMeVperc2 = 931.5; // 1 amu = 931.5 MeV/c^2 
Double_t c = 299792458; //speed of light in m/s
Double_t eVperMeV = 1000000; //1MeV = 1 000 000 eV
Double_t nspers = 1000000000; //1s = 1 000 000 000 ns
Double_t mpermm = 0.001; // 1mm = 0.001 m

//Output options
Bool_t SilenceSolver = false; //if true silences the solver errors generated when no solution is found

//Other declarations
TString DataFileStripped; //Definition for output file name


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
    /*
    const ROOT::Math::Roots::Brent* T = new ROOT::Math::Roots::Brent();
    ROOT::Math::BrentRootFinder* RF = new ROOT::Math::BrentRootFinder();

   TF1 F;
   struct my_f_params params = { Z_det, SiliconRadius, v_lab_ej, TCyc };

   F.function = &my_f;
   F.params = &params;
   */
    //struct my_f_params params = { Z_det, SiliconRadius, v_lab_ej, TCyc };

    //update parameters
    Current_Parameters.z_det = Z_det;
    Current_Parameters.h_det = SiliconRadius;
    Current_Parameters.v_lab = v_lab_ej;
    Current_Parameters.TCyc = TCyc;
    //std::cout << "z_det is: " << Z_det << std::endl;
    //std::cout << "h_det is: " << SiliconRadius << std::endl;
    //std::cout << "v_lab is: " << v_lab_ej << std::endl;
    //std::cout << "TCyc is: " << TCyc << std::endl;
    RootFinder *k = new RootFinder();
    k->SetMethod(RootFinder::kGSL_BRENT);
    ROOT::Math::Functor1D f(&my_f);
    //f.params = &params;
    double EndPoint = {TMath::Pi()};
    double TestIsNegative = {1};
    for (int i=0; i<10000; i++){
        TestIsNegative = my_f(TMath::Pi()/2 + TMath::Pi()/2*i/10000);
        if (TestIsNegative<0){
            if(Z_det < - 340){
            //std::cout << "Test angle is: " << (TMath::Pi()/2 + TMath::Pi()/2*i/10000)*180/TMath::Pi() << " and Function value is "<< TestIsNegative <<endl;
            }
            EndPoint = TMath::Pi()/2 + TMath::Pi()/2*i/10000;
            break;
        }
    }
    //std::cout << "Used endpoint is "<< EndPoint <<endl;
    //std::cout << "Function value at endpoint is "<< my_f(EndPoint) <<endl;
    k->SetFunction(f, TMath::Pi()/2+0.001, EndPoint);
    if(SilenceSolver){gROOT->ProcessLine( "gErrorIgnoreLevel = 6001;");} //Suppress all error messages
    k->Solve();
    if(SilenceSolver){gROOT->ProcessLine( "gErrorIgnoreLevel = 1;");} //back to normal

 
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

/*!  \brief Calculates ejectile com angle
*    
*    Returns ejectile com angle in degree
*    @param z_det measured z distance on array (in mm)
*    @param Elab Energy deposited in detector (in MeV)
*    @param Z_corr corrected z to axis (x=y=0)(in mm), used only to pass output
*    @param ExEn Calculated excitation energy of recoil (in MeV), used only to pass output
*/
Double_t SolveReaction(Double_t z_det, Double_t Elab, Double_t &Z_corr, Double_t &ExEn) {

    //These four are only dependent on the reaction type
    Double_t TCyclotron = GetCyclotronPeriod();
    Double_t CMvelocity = GetCOMVelocity();
    Double_t CMEnergy_in = GetCOMEnergyInit();
    Double_t QValue = GetQValue();

    //The lab velocity of the proton can be calculated from the lab energy
    Double_t v_lab = GetEjectileLabVelocity(Elab);

    //Now we calculate the lab angle using the brent-decker algorithm
    Double_t ThetaLab = CalcCorrectedLabAngle(z_det, v_lab, TCyclotron);

    //Define Theta_cm
    Double_t ThetaCM = -1;

    if(ThetaLab > 0){ //Theta lab solution was succesfully found

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
        Z_corr = 1000;
        ExEn = -1000;
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


void ConvertElabZhist_v1(){

  time_t t;
  time(&t);
  long double TimeInit = t;
  
  //1.---------Reading data file. 

  
  TString DataPath = "../IS587/"+ DataFileName;
  DataFileStripped = DataFileName.ReplaceAll(".root",""); //used for output file name composition in main function
  TFile *InputFile = new TFile(DataPath,"READ");

  if (!InputFile){
	  std::cout << "File containing histogram not found, check if file has correct name" << std::endl;
	  return;
  }
  TH2F *Edep_Z_Hist=(TH2F*)InputFile->Get(TH2FileName);
  if (!Edep_Z_Hist){
      std::cout << "Histogram not found inside file, check file content" << std::endl;
	  return;
  }

  Int_t NrBinsX = Edep_Z_Hist->GetNbinsX();
  Int_t NrBinsY = Edep_Z_Hist->GetNbinsY();
  std::cout << "Original histogram contains " << NrBinsX << " x bins and " << NrBinsY << " y bins" << std::endl;

  


  //We create vectors to store the variables
  std::vector<Double_t> CentersCoMangle {};
  std::vector<Double_t> CentersZCorr {};
  std::vector<Double_t> CentersExEn {};
  std::vector<Double_t> CentersEDep {};
  std::vector<Double_t> CountsVector {};
  std::vector<Double_t> LeftLimitsCoMangle {};
  std::vector<Double_t> RightLimitsCoMangle {};
  std::vector<std::array<double,3>> EmptyLeftLimitsCoMangle {}; //Saves the left and right limits for gaps of detector
  std::vector<std::array<double,3>> EmptyRightLimitsCoMangle {};
  std::vector<std::array<double,3>> GapsExEn {};

  for (Int_t i=1; i<=NrBinsY; i++) {
      for (Int_t j=1; j<=NrBinsX; j++) {
          Double_t ZCenter = ((TAxis*)Edep_Z_Hist->GetXaxis())->GetBinCenter(j);
          Double_t EnergyCenter = ((TAxis*)Edep_Z_Hist->GetYaxis())->GetBinCenter(i) / 1000;// divide by 1000 to convert from keV to MeV
          Double_t ExEner = 0;
          Double_t ZCorrected = 0;
          Double_t CoMAnglecenter = 0;
          CoMAnglecenter = SolveReaction(ZCenter, EnergyCenter, ZCorrected, ExEner);  
          CentersCoMangle.push_back(CoMAnglecenter); 
          CentersZCorr.push_back(ZCorrected); 
          CentersExEn.push_back(ExEner); 
          CentersEDep.push_back(EnergyCenter);
          CountsVector.push_back(Edep_Z_Hist->GetBinContent(j,i));
          //std::cout << "Current Z Center is: " << ZCenter << std::endl;
          //std::cout << "Current Energy Center is: " << EnergyCenter << std::endl;     
      }
  }


  //We also calculate the corresponding cm angles of the gaps in the si detector (for each energy) 
  for (Int_t n=1; n<=NrBinsY; n++) {
       std::array<double,3> LeftLimits = {0.,0.,0.};
       std::array<double,3> RightLimits = {0.,0.,0.};
       std::array<double,3> ThisExEn = {0.,0.,0.};
       for(Int_t j=0; j<8; j++){
          
          Double_t ZCenter = Edges[j];
          Double_t EnergyCenter = ((TAxis*)Edep_Z_Hist->GetYaxis())->GetBinCenter(n) / 1000;// divide by 1000 to convert from keV to MeV
          Double_t ExEner = 0;
          Double_t ZCorrected = 0;
          Double_t CoMAnglecenter = 0;
          
          //std::cout << "n is "  << n << " j is " << j << std::endl;
          
          CoMAnglecenter = SolveReaction(ZCenter, EnergyCenter, ZCorrected, ExEner);  
          if(j == 1 || j == 3 || j==5){ //'left' sides of gaps
              int IndexInt = (int) (floor(j/2));
              LeftLimits[IndexInt] = CoMAnglecenter;
          }
          if(j == 2 || j == 4 || j==6){ //'right' sides of gaps
              int IndexInt = (j/2)-1;
              LeftLimits[IndexInt] = CoMAnglecenter;
              ThisExEn[IndexInt] = ExEner;
          }
           
       }
       EmptyLeftLimitsCoMangle.push_back(LeftLimits);
       EmptyRightLimitsCoMangle.push_back(RightLimits);
       GapsExEn.push_back(ThisExEn);
  }


  //Create the histogram for angular distribution
  TH1D* AngularDistHist = new TH1D("Theta_CM", "Angular Distribution CM frame; CM angle [degree]; Counts", 60, 0, 60);

  //Create a histogram for the excitation energy
  TH1D* ExenHist = new TH1D("Ex", "Excitation Energy Spectrum ; Excitation Energy [MeV]; Counts", 1000, -5, 10);

  //Create a histogram for the corrected z position
  TH1D* ZcorrHist = new TH1D("Zcorr", "Corrected z Distance of Cyclotron Loop ; Corrected z distance [mm]; Counts", 1000, -700, 50);



  //Create a histogram for Ex vs Z (corrected Z)
  TH2D* ExvsZcorrHist = new TH2D("Ex_vs_Zcorr", "Ex vs Z_corr; Corrected z distance [mm]; Excitation Energy [MeV]", 1000, -700, 50,1000, -5, 10);

  //Create a histogram for Ex vs Theta_CM
  TH2D* ExvsThetaCMHist = new TH2D("Ex_vs_Theta_CM", "Ex vs Theta_CM; CM angle [degree]; Excitation Energy [MeV]", 60, 0, 60, 1000, -5, 10);

  //Create a histogram for Z vs Theta_CM
  TH2D* ZcorrvsThetaCMHist = new TH2D("ZCorr_vs_Theta_CM", "Zcorr vs Theta_CM; CM angle [degree]; Corrected z distance [mm]", 60, 0, 60, 1000, -700, 50);

  //We fill the histograms
  for (Int_t j=0; j<NrBinsX*NrBinsY; j++) {
      AngularDistHist->Fill(CentersCoMangle[j]*180/TMath::Pi(),CountsVector[j]);
      ExenHist->Fill(CentersExEn[j],CountsVector[j]);
      ZcorrHist->Fill(CentersZCorr[j],CountsVector[j]);
      ExvsZcorrHist->Fill(CentersZCorr[j], CentersExEn[j],CountsVector[j]);
      ExvsThetaCMHist->Fill(CentersCoMangle[j]*180/TMath::Pi(), CentersExEn[j],CountsVector[j]);
      ZcorrvsThetaCMHist->Fill(CentersCoMangle[j]*180/TMath::Pi(), CentersZCorr[j],CountsVector[j]);
  }


  TString HistFile = "OUTPUT/test.root";
  TFile *MyHistFile = new TFile(HistFile,"RECREATE");
  MyHistFile->cd();
  AngularDistHist->Write();
  ExenHist->Write();
  ZcorrHist->Write();
  ExvsZcorrHist->Write();
  ExvsThetaCMHist->Write();
  ZcorrvsThetaCMHist->Write();
  MyHistFile->Close();

  time(&t);
  long double TimeFin = t;
  //std::cout << "Events were processed in " << (int)floor((TimeFin-TimeInit)/60) << " min and " << ((int)(TimeFin-TimeInit))%60 << " sec" << std::endl; 


  InputFile->Close();

 
}
