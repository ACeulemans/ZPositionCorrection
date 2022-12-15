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
#include <Math/RootFinderAlgorithms.h>
#include <TF1.h>
#include <Math/RootFinder.h>
#include <Math/Functor.h>

using namespace ROOT::Math;


//Source/input options  (Don't forget to change name of function/file as well if you change this)
TString DataFileName = "R18-114_autocal_AllAlpha_sum_hists.root"; //file to read (must be located in ./ROOT_FILES)
TString sourcename = "Ba133_302keV"; //name of source
TString TH2FileName = "EBISMode/Ex_vs_z_ebis"; //name of th2 inside the datafile to open ("Scint Energy" for scintillators or "PadPlane Energy" for padplane)
Double_t ExcitationEnergy = 2500; //Excitation energy in keV (central value)
Double_t PeakFWHM = 100; //Bins with center in [ExcitationEnergy-PeakFWHM,ExcitationEnergy+PeakFWHM] are summed

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

//Conversion factors
Double_t amutoMeVperc2 = 931.5; // 1 amu = 931.5 MeV/c^2 
Double_t c = 299792458; //speed of light in m/s
Double_t eVperMeV = 1000000; //1MeV = 1 000 000 eV
Double_t nspers = 1000000000; //1s = 1 000 000 000 ns
Double_t mpermm = 0.001; // 1mm = 0.001 m

//Output options
Bool_t SaveHists = true; //If true saves all fitted histograms
Bool_t MakePlots = false; //If true histograms are drawn to canvas and displayed

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

/*!  \brief Calculates ejectile com angle
*    
*    Returns ejectile com angle in radian
*    @param TCyc cyclotron period of ejectile (in ns)
*    @param Vcm center of mass velocity (in fraction of c)
*    @param z z distance for a single loop (measured z on array corrected to x=y=0 axis) (in mm)
*    @param v_cm_ejec ejectile velocity in CM frame (in fraction of c)
*/
Double_t GetCoMAngle2(Double_t TCyc, Double_t Vcm, Double_t z, Double_t v_cm_ejec) {
    Double_t CoMAngle = acos((((z*mpermm)/(TCyc/nspers)/c)-Vcm)/v_cm_ejec);
    //std::cout << "The lab velocity of the ejectile is " << ejeclabvel << " c" << std::endl;
    return CoMAngle;
}

/*!  \brief Calculates ejectile com angle
*    
*    Returns ejectile com angle in radian
*    @param TCyc cyclotron period of ejectile (in ns)
*    @param Vcm center of mass velocity (in fraction of c)
*    @param z z distance for a single loop (measured z on array corrected to x=y=0 axis) (in mm)
*    @param v_cm_ejec ejectile velocity in CM frame (in fraction of c)
*/
Double_t GetCoMAngle3(Double_t TCyc, Double_t Vcm, Double_t z, Double_t v_cm_ejec, Double_t ThetaLab) {
    //Double_t CoMAngle = acos((((z*mpermm*tan(ThetaLab))/(TCyc/nspers))/(v_cm_ejec*c));
    Double_t CoMAngle = acos( -(z*mpermm*tan(ThetaLab)) / (TCyc/nspers*v_cm_ejec*c) );
    //std::cout << "The lab velocity of the ejectile is " << ejeclabvel << " c" << std::endl;
    return CoMAngle;
}

/*!  \brief Calculates ejectile com velocity
*    
*    Returns ejectile velocity in CM frame (in fraction of c)
*    @param CMEnergy_after total energy in CM frame after reaction (in MeV)
*/
Double_t GetEjectileCMVelocity(Double_t CMEnergy_after) {
    Double_t V_cm_ejec = sqrt(2*CMEnergy_after/(EjectileMass*amutoMeVperc2)*(1+(EjectileMass*EjectileMass)/(RecoilMass*RecoilMass)));
    //std::cout << "The lab velocity of the ejectile is " << ejeclabvel << " c" << std::endl;
    return V_cm_ejec;
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
*    @param LabAngle_in ejectile lab angle used for starting value (in radian)
*    @param Z_det z distance for a single loop (measured z on array without correction) (in mm)
*    @param v_lab_ej ejectile lab velocity (in fraction of c)
*    @param TCyc cyclotron period of ejectile (in ns)
*    
*/
Double_t CalcCorrectedLabAngle(Double_t LabAngle_in, Double_t Z_det, Double_t v_lab_ej, Double_t TCyc){
    
    Double_t LabAngle_corr = 0;
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
    std::cout << "v_lab is: " << v_lab_ej << std::endl;
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
            std::cout << "Test angle is: " << (TMath::Pi()/2 + TMath::Pi()/2*i/10000)*180/TMath::Pi() << " and Function value is "<< TestIsNegative <<endl;
            }
            EndPoint = TMath::Pi()/2 + TMath::Pi()/2*i/10000;
            break;
        }
    }
    //std::cout << "Used endpoint is "<< EndPoint <<endl;
    //std::cout << "Function value at endpoint is "<< my_f(EndPoint) <<endl;
    k->SetFunction(f, TMath::Pi()/2+0.001, EndPoint);
    k->Solve();

    double c = k->Root();
    //std::cout << "Lab angle before correction is "<< v_lab_ej*180/TMath::Pi() <<endl;
    std::cout << "New Lab angle is "<< c*180/TMath::Pi() <<endl;

    return c;
}

double GetCorrectedZValue(Double_t LabAngle, Double_t v_lab_ej, Double_t TCyc){
    return cos(LabAngle)*v_lab_ej*c*TCyc/nspers/mpermm;
}

/*!  \brief Calculates ejectile com angle
*    
*    Returns ejectile com angle in degree
*    @param z z distance for a single loop (measured z on array corrected to x=y=0 axis) (in mm)
*    @param Ex excitation energy of recoil (in MeV)
*/
Double_t GetCoMangle(Double_t z, Double_t Ex) {
    Double_t TCyclotron = GetCyclotronPeriod();
    Double_t CMvelocity = GetCOMVelocity();
    Double_t CMEnergy_in = GetCOMEnergyInit();
    Double_t QValue = GetQValue();
    Double_t CMEnergy_out = CMEnergy_in + QValue - Ex;
    Double_t Vcm_ejec = GetEjectileCMVelocity(CMEnergy_out);
    Double_t ThetaCM = TMath::Pi()-GetCoMAngle2(TCyclotron, CMvelocity, z, Vcm_ejec);

    //Double_t ThetaLab = GetEjectileLabAngle(Vcm_ejec-CMvelocity,z, TCyclotron);
    Double_t ThetaLab = GetEjectileLabAngle(std::sqrt(Vcm_ejec*Vcm_ejec*sin(ThetaCM)*sin(ThetaCM)+(Vcm_ejec*cos(ThetaCM)-CMvelocity)*(Vcm_ejec*cos(ThetaCM)-CMvelocity)),z, TCyclotron);
    Double_t ThetaLab_corr = CalcCorrectedLabAngle(ThetaLab, z, std::sqrt(Vcm_ejec*Vcm_ejec*sin(ThetaCM)*sin(ThetaCM)+(Vcm_ejec*cos(ThetaCM)-CMvelocity)*(Vcm_ejec*cos(ThetaCM)-CMvelocity)), TCyclotron); 
    Double_t z_corr = GetCorrectedZValue(z, ThetaLab);
    Double_t ThetaCM_corr = TMath::Pi()-GetCoMAngle2(TCyclotron, CMvelocity, z_corr, Vcm_ejec);
    //using brent dekker correction
    Double_t z_corr_2 = GetCorrectedZValue(ThetaLab_corr,std::sqrt(Vcm_ejec*Vcm_ejec*sin(ThetaCM_corr)*sin(ThetaCM_corr)+(Vcm_ejec*cos(ThetaCM_corr)-CMvelocity)*(Vcm_ejec*cos(ThetaCM_corr)-CMvelocity)),TCyclotron);
    Double_t ThetaCM_corr_2 = TMath::Pi()/2 - (TMath::Pi()-GetCoMAngle3(TCyclotron, CMvelocity, z_corr_2, Vcm_ejec, ThetaLab_corr));

    if (std::isnan(ThetaCM)){
        ThetaLab = std::numeric_limits<double>::quiet_NaN();
        z_corr = std::numeric_limits<double>::quiet_NaN();
        ThetaCM_corr = std::numeric_limits<double>::quiet_NaN();
    }
    
    //std::cout << "The cyclotron period is " << TCyclotron << " ns" << std::endl;
    //std::cout << "The velocity of the CM frame is " << CMvelocity << " c" << std::endl;
    //std::cout << "The total energy in the initial CM frame is " << CMEnergy_in << " MeV" << std::endl;
    //std::cout << "The Q-value of the reaction is " << QValue << " MeV" << std::endl;
    //std::cout << "The total energy in the outgoing CM frame is " << CMEnergy_out << " MeV" << std::endl;
    //std::cout << "The ejectile CM velocity is " << Vcm_ejec << " c" << std::endl;
    //std::cout << "The z is " << z << " mm" << std::endl;
    //std::cout << "The CM angle is " << ThetaCM*180/TMath::Pi() << " degrees" << std::endl;
    //std::cout << "The lab angle is " << ThetaLab*180/TMath::Pi() << " degrees" << std::endl;
    //std::cout << "The corrected lab angle is " << ThetaLab_corr*180/TMath::Pi() << " degrees" << std::endl;
    //std::cout << "The single step corrected z is " << z_corr << " mm" << std::endl;
    //std::cout << "The single step corrected CM angle is " << ThetaCM_corr*180/TMath::Pi() << " degrees" << std::endl;
    //std::cout << "The brent step corrected z is " << z_corr_2 << " mm" << std::endl;
    //std::cout << "The brent corrected CM angle is " << ThetaCM_corr_2*180/TMath::Pi() << " degrees" << std::endl;

    return ThetaCM_corr_2;
}


void ConvertExZhist(){

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
  TH2F *Ex_Z_Hist=(TH2F*)InputFile->Get(TH2FileName);
  if (!Ex_Z_Hist){
      std::cout << "Histogram not found inside file, check file content" << std::endl;
	  return;
  }

  Int_t NrBinsX = Ex_Z_Hist->GetNbinsX();
  Int_t NrBinsY = Ex_Z_Hist->GetNbinsY();
  //std::cout << "Original histogram contains " << NrBinsX << " x bins and " << NrBinsY << " y bins" << std::endl;

  //First we sum the requested bins

  
  //Create and initialize a vector for counts, left & right limit of bins 
  std::vector<Double_t> BinCenters {};
  std::vector<Double_t> Counts_sum_perZ {};
  std::vector<Double_t> LeftLimits {};
  std::vector<Double_t> RightLimits {};
  std::vector<Double_t> EmptyLeftLimits {}; //Saves the left and right limits for empty bins
  std::vector<Double_t> EmptyRightLimits {};

  for (Int_t j=1; j<=NrBinsX; j++) {
     Counts_sum_perZ.push_back(0);      
     LeftLimits.push_back(0);   
     RightLimits.push_back(0);   
     BinCenters.push_back(0); 
  }

  double BinCenterX {0};
  double BinCenterY {0};

  for (Int_t i=1; i<=NrBinsY; i++) {
      BinCenterY = ((TAxis*)Ex_Z_Hist->GetYaxis())->GetBinCenter(i);
      if((BinCenterY>=(ExcitationEnergy-PeakFWHM)) && (BinCenterY<=(ExcitationEnergy+PeakFWHM))){
          for (Int_t j=1; j<=NrBinsX; j++) {
              Counts_sum_perZ[j-1] += Ex_Z_Hist->GetBinContent(j,i);
              LeftLimits[j-1] = ((TAxis*)Ex_Z_Hist->GetXaxis())->GetBinLowEdge(j);
              RightLimits[j-1] = ((TAxis*)Ex_Z_Hist->GetXaxis())->GetBinUpEdge(j);
              BinCenters[j-1] = ((TAxis*)Ex_Z_Hist->GetXaxis())->GetBinCenter(j);
              if(Counts_sum_perZ[j]==0){
                  EmptyLeftLimits.push_back(LeftLimits[j-1]);
                  EmptyRightLimits.push_back(RightLimits[j-1]);
              }
       
          }
      }
  }

  //Now we need to calculate the center of mass velocity


  //Next we calculate the new histogram bin edges
  std::vector<Double_t> LeftLimitsCoMangle {};
  std::vector<Double_t> RightLimitsCoMangle {};
  std::vector<Double_t> CentersCoMangle {};
  std::vector<Double_t> EmptyLeftLimitsCoMangle {}; //Saves the left and right limits for empty bins
  std::vector<Double_t> EmptyRightLimitsCoMangle {};

  for (Int_t j=0; j<NrBinsX; j++) {
     //LeftLimitsCoMangle.push_back(GetCoMangle(LeftLimits[j],ExcitationEnergy/1000));      //division by 1000 to convert from keV to MeV
     //RightLimitsCoMangle.push_back(GetCoMangle(RightLimits[j],ExcitationEnergy/1000)); 
     CentersCoMangle.push_back(GetCoMangle(BinCenters[j],ExcitationEnergy/1000)); 
     std::cout << "CurrentBincenter is: " << BinCenters[j] << std::endl;
     std::cout << "CentersCoMangle is: " << CentersCoMangle[j]*180/TMath::Pi() << std::endl;     
  }

  for (Int_t j=0; j<EmptyLeftLimits.size(); j++) {
     //EmptyLeftLimitsCoMangle.push_back(GetCoMangle(EmptyLeftLimits[j],ExcitationEnergy/1000));      //division by 1000 to convert from keV to MeV
     //EmptyRightLimitsCoMangle.push_back(GetCoMangle(EmptyRightLimits[j],ExcitationEnergy/1000));   
  }

  //Create the histogram for angular distribution
  TH1D* AngularDistHist = new TH1D("Theta_CM", "Angular Distribution CM frame; CM angle [degree]; Counts", 60, 0, 60);

  //We fill the histogram
  for (Int_t j=0; j<NrBinsX; j++) {
      AngularDistHist->Fill(CentersCoMangle[j]*180/TMath::Pi(),Counts_sum_perZ[j]);
  }

  TString HistFile = "OUTPUT/test_brent.root";
  TFile *MyHistFile = new TFile(HistFile,"RECREATE");
  MyHistFile->cd();
  AngularDistHist->Write();
  MyHistFile->Close();





  

  time(&t);
  long double TimeFin = t;
  //std::cout << "Events were processed in " << (int)floor((TimeFin-TimeInit)/60) << " min and " << ((int)(TimeFin-TimeInit))%60 << " sec" << std::endl; 


  InputFile->Close();

 
}
