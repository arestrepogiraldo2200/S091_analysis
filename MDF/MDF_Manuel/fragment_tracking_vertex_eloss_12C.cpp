#include "libs.hh"
#include "mass.h"
#include <R3BMDFWrapper.h>
#include <R3BEventHeader.h>
#include <R3BRpcHitData.h>
#include <R3BTrackS522.h>
#include "R3BTrackS522.h"
#include "R3BLosHitData.h"
#include "R3BRpcHitData.h"
#include "R3BTofdHitData.h"
#include "R3BCalifaClusterData.h"
#include "R3BFootHitData.h"
#include "R3BMwpcHitData.h"
#include "R3BMwpcHitData.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BCoarseTimeStitch.h"

using namespace std;

struct vertex
{
    Double_t min_dist;
    TVector3 vertexPos;
};

struct Material{
	double ZoA;            // Atomic mass (in g/mol)
	double I;              // Mean excitation energy (in MeV)
	double rho;            // material density
	double weight_fraction; // Weight fraction in the mixture
};

//Storing indices of hits in TCA for potential track candidates
struct Track
{
	double vertex_x_foot;
	double vertex_y_foot;
	double vertex_z_foot;
	double f1_x;
	double f1_z;
	double f2_y;
	double f2_z;
	double f15_x;
	double f15_z;
	double f16_y;
	double f16_z;
	double f1_Q;
	double f2_Q;
	double f15_Q;
	double f16_Q;
	double tx_foot;
	double ty_foot;
	double f32_x;
	double f32_z;
	double f30_x;
	double f30_y;
	double f30_z;
	double last_x;
	double last_z;
};

double my_nan = numeric_limits<double>::quiet_NaN(); 

R3BMDFWrapper* MDF_PoQ_foot;
R3BMDFWrapper* MDF_TX0_foot;
R3BMDFWrapper* MDF_TY0_foot;
R3BMDFWrapper* MDF_FlightPath_foot;

R3BMDFWrapper* MDF_PoQ;
R3BMDFWrapper* MDF_TX0;
R3BMDFWrapper* MDF_TY0;
R3BMDFWrapper* MDF_FlightPath;

TRotation r;
TVector3 v3_localX;
TVector3 v3_localZ;
///measruments
const Double_t Angle_frag = -14* TMath::Pi() / 180.;//turning angle of GLAD to centr of RPC
const Double_t Target_to_GLAD_flange = 106.2;
const Double_t Target_to_TP = Target_to_GLAD_flange + 165.;
const double SPEED_OF_LIGHT = 29.9792458; // cm/ns
const double AMU = 0.9314940038;          // GeV/c2
const double FRS_BETA = 0.904252;

//some variable 
double mp = 0.93892; 
double me = 0.000511; 
double m_12C = Nuke_Mass_Tab[12][6]/1000. - me*6;

//MWPC alignment
double	offset_x_mwpc0 = 0.1543;
double	offset_y_mwpc0 = 0.6864;
double	offset_x_mwpc1 = -0.5384;
double	offset_y_mwpc1 = 2.7001;

////foor z offset
//double	offset_f1 = 20.305200;
//double	offset_f2 = 21.956;
//double	offset_f15 = 71.7636;
//double	offset_f16 = 75.2645;
//foor z offset
double	offset_f1 = 20.105200;
double	offset_f2 = 21.7256;
double	offset_f15 = 71.5636;
double	offset_f16 = 75.0645;

//FIB
const Double_t TP_to_F32 = 434.1;
const Double_t TP_to_F30 = 459.6;
const Double_t TP_to_F33 = 594.5;
const Double_t TP_to_F31 = 611.4;

const Double_t F31_offset = -25.6;
const Double_t F33_offset = 25;

const double fib32x_off = 3.1058277;
const double fib32z_off = 1.295688;
const double fib30y_off = - 0.580751;
const double fib30z_off = - 1.27848;
const double fib33x_off = 2.8319228;
const double fib33z_off = 0.0248933;
const double fib31x_off = 0.8904328;
const double fib31z_off = 5.3244533;

//const double fib32x_off = 3.99168;
////const double fib32z_off = 0;
////const double fib30y_off = - 0.45576;
////const double fib30z_off = 0;
////const double fib33x_off = 3.61499;
////const double fib33z_off = 0;
////const double fib31x_off = 3.10928;
////const double fib31z_off = 0;

TVector3 f32_position;
TVector3 f30_position;
TVector3 f33_position;
TVector3 f31_position;
TVector3 fib_angles;

TVector3 f32_point;
TVector3 f30_point;
TVector3 flast_point;

TVector3 f32_point_i;
TVector3 f30_point_i;
TVector3 flast_point_i;

std::vector<double > f1_hits, f2_hits, f15_hits, f16_hits;

std::vector<Track> tracks_out;//track candidates in FOOT
std::vector<Track> tracks_in;//track candidates in FOOT

inline TVector3 GetPositionFib(int fib_number){
	switch(fib_number){
		case 30:
			return f30_position;
			break;
		case 31:
			return f31_position;
			break;
		case 32:
			return f32_position;
			break;
		case 33:
			return f33_position;
			break;
	}
	return TVector3();
} // cm

inline TVector3 GetEulerAnglesFib()    { return fib_angles;    } // rad

void SetPositionFib(int fib_number, double x, double y, double z){
	switch(fib_number){
		case 30:
			f30_position.SetXYZ(x,y,z);
			break;
		case 31:
			f31_position.SetXYZ(x,y,z);
			break;
		case 32:
			f32_position.SetXYZ(x,y,z);
			break;
		case 33:
			f33_position.SetXYZ(x,y,z);
			break;
	}
}// cm
void SetEulerAnglesfib(double x, double y, double z)  { fib_angles.SetXYZ(x, y, z);  }// rad

//auxiliary functions
void TransformPoint(TVector3& point, TVector3* rotation, TVector3* translation);
bool IsGoodFootHit(R3BFootHitData* fhit);
bool MakeOutgoingTracks();
bool MakeIncomingTracks();
double CalculateEnergyLoss(double beta, double Z, double mass);
double EnergyLossBetheBloch(double beta, const std::vector<Material>& materials, double thickness, double charge, double mass);

vertex get_vertex(TVector3 slope1, TVector3 slope2, TVector3 intercept1, TVector3 intercept2);

TClonesArray *losCA = new TClonesArray("R3BLosHitData");
TClonesArray *tofdCA = new TClonesArray("R3BTofdHitData");
TClonesArray *footCA = new TClonesArray("R3BFootHitData");
TClonesArray *mwpc0CA = new TClonesArray("R3BMwpcHitData");
TClonesArray *mwpc1CA = new TClonesArray("R3BMwpcHitData");
TClonesArray *fib32CA = new TClonesArray("R3BFiberMAPMTHitData");
TClonesArray *fib33CA = new TClonesArray("R3BFiberMAPMTHitData");
TClonesArray *fib31CA = new TClonesArray("R3BFiberMAPMTHitData");
TClonesArray *fib30CA = new TClonesArray("R3BFiberMAPMTHitData");

R3BEventHeader *fEventHeader  = new R3BEventHeader();

R3BCoarseTimeStitch* fTimeStitch = new R3BCoarseTimeStitch();

TLorentzVector pbeam_tree;
TLorentzVector pfrag_mdf_foot_tree;
TLorentzVector pfrag_foot_tree;
TLorentzVector pfrag_mdf_vertex_tree;

TLorentzVector pbeam_lab;
TLorentzVector pfrag_mdf_foot_lab;
TLorentzVector pfrag_foot_lab;
TLorentzVector pfrag_mdf_vertex_lab;

bool good_incoming = false;
Double_t mwpc1_x = my_nan;
Double_t mwpc1_y = my_nan;
Double_t mwpc1_z = my_nan;
Double_t mwpc0_x = my_nan;
Double_t mwpc0_y = my_nan;
Double_t mwpc0_z = my_nan;

Double_t mul_frag_vertex = my_nan;
Double_t mul_frag_foot = my_nan;

Double_t Q_tree = my_nan;
Double_t AoZ_vertex = my_nan;
Double_t AoZ_foot = my_nan;

Double_t vertex_x = my_nan;
Double_t vertex_y = my_nan;
Double_t vertex_z = my_nan;

Double_t vertex_x_foot = my_nan;
Double_t vertex_y_foot = my_nan;
Double_t vertex_z_foot = my_nan;

Double_t tx_mwpc = my_nan;
Double_t ty_mwpc = my_nan;

Double_t tx_foot = my_nan;
Double_t ty_foot = my_nan;

Double_t tx_mdf = my_nan;
Double_t ty_mdf = my_nan;

TVector3 poq_vec_vertex;
TVector3 poq_vec_foot;
TVector3 rec_vertex;
Double_t min_dist = my_nan;
Double_t f1_Q = my_nan;

void reset_val(){
	good_incoming = false;
	mwpc1_x = my_nan;
	mwpc1_y = my_nan;
	mwpc1_z = my_nan;
	mwpc0_x = my_nan;
	mwpc0_y = my_nan;
	mwpc0_z = my_nan;

	mul_frag_vertex = my_nan;
	mul_frag_foot = my_nan;

	tx_foot = my_nan;
	ty_foot = my_nan;

	tx_mdf = my_nan;
	ty_mdf = my_nan;

	Q_tree = my_nan;
	AoZ_vertex = my_nan;
	AoZ_foot = my_nan;

	vertex_x = my_nan;
	vertex_y = my_nan;
	vertex_z = my_nan;

	vertex_x_foot = my_nan;
	vertex_y_foot = my_nan;
	vertex_z_foot = my_nan;

	tx_mwpc = my_nan;
	ty_mwpc = my_nan;

	poq_vec_vertex.SetXYZ(my_nan,my_nan,my_nan);
	poq_vec_vertex.SetMag(my_nan);

	poq_vec_foot.SetXYZ(my_nan,my_nan,my_nan);
	poq_vec_foot.SetMag(my_nan);

	rec_vertex.SetXYZ(my_nan,my_nan,my_nan);
	rec_vertex.SetMag(my_nan);

	min_dist = my_nan;

	pbeam_tree.SetXYZM(my_nan,my_nan,my_nan,my_nan);
	pfrag_mdf_vertex_tree.SetXYZM(my_nan,my_nan,my_nan,my_nan);
	pfrag_mdf_foot_tree.SetXYZM(my_nan,my_nan,my_nan,my_nan);
	pfrag_foot_tree.SetXYZM(my_nan,my_nan,my_nan,my_nan);

	pbeam_lab.SetXYZM(my_nan,my_nan,my_nan,my_nan);
       	pfrag_mdf_vertex_lab.SetXYZM(my_nan,my_nan,my_nan,my_nan);
       	pfrag_mdf_foot_lab.SetXYZM(my_nan,my_nan,my_nan,my_nan);
       	pfrag_foot_lab.SetXYZM(my_nan,my_nan,my_nan,my_nan);
	f1_Q = my_nan;
}

TH2F *tofd_Q1_vs_Q2 = new TH2F("tofd_Q1_vs_Q2","tofd_Q1_vs_Q2",800,0,8,800,0,8); 
TH1F *LOS_mul = new TH1F("LOS_mul","LOS_mul",10,0,10); 

TH1F *fiber32_cluster = new TH1F("fiber32_cluster","fiber32_cluster",20,0,20); 
TH1F *fiber30_cluster = new TH1F("fiber30_cluster","fiber30_cluster",20,0,20); 
TH1F *fiber31_cluster = new TH1F("fiber31_cluster","fiber31_cluster",20,0,20); 
TH1F *fiber33_cluster = new TH1F("fiber33_cluster","fiber33_cluster",20,0,20); 

TH1F *fiber32_time = new TH1F("fiber32_time","fiber32_time",500,-210,-100); 
TH1F *fiber30_time = new TH1F("fiber30_time","fiber30_time",500,-210,-100); 
TH1F *fiber31_time = new TH1F("fiber31_time","fiber31_time",500,-210,-100); 
TH1F *fiber33_time = new TH1F("fiber33_time","fiber33_time",500,-210,-100); 

TH2F *fiber32_time_vs_pos = new TH2F("fiber32_time_vs_pos","fiber32_time_vs_pos",500,-210,-100,180,-30,30); 
TH2F *fiber30_time_vs_pos = new TH2F("fiber30_time_vs_pos","fiber30_time_vs_pos",500,-210,-100,180,-30,30); 
TH2F *fiber31_time_vs_pos = new TH2F("fiber31_time_vs_pos","fiber31_time_vs_pos",500,-210,-100,180,-30,30); 
TH2F *fiber33_time_vs_pos = new TH2F("fiber33_time_vs_pos","fiber33_time_vs_pos",500,-210,-100,180,-30,30); 

std::map<std::string, TF1*> functionMap_Si;
std::map<std::string, TF1*> functionMap_H2;

void tree_merger()
{
	TFile *fit_file = new TFile("test_energy.root", "READ");
	if (!fit_file || fit_file->IsZombie()) {
		std::cerr << "Error opening file" << std::endl;
		return;
	}

	for (int mass = 1; mass <= 16; mass++) {
		for (int Q = 1; Q <= 6; Q++) {
			std::string title_fit_Si = "FIT_Si_Mass_" + std::to_string(mass) + "_Q_" + std::to_string(Q);
			std::string title_fit_H2 = "FIT_H2_Mass_" + std::to_string(mass) + "_Q_" + std::to_string(Q);
			TF1 *quadFit = nullptr;
			TF1 *linearFit = nullptr;
			fit_file->GetObject(title_fit_H2.c_str(), quadFit);
			fit_file->GetObject(title_fit_Si.c_str(), linearFit);
			if (quadFit) {
				functionMap_H2[title_fit_H2] = (TF1*)quadFit->Clone();
			}
			if (linearFit) {
				functionMap_Si[title_fit_Si] = (TF1*)linearFit->Clone();
			}
		}
	}

	SetPositionFib(
			32,
			TP_to_F32 * TMath::Sin(Angle_frag) + fib32x_off,
			0,
			Target_to_TP + TP_to_F32 * TMath::Cos(Angle_frag) + fib32z_off);

	SetPositionFib(
			30,
			TP_to_F30 * TMath::Sin(Angle_frag),
			+ fib30y_off,
			Target_to_TP + TP_to_F30 * TMath::Cos(Angle_frag) + fib32z_off);

	SetPositionFib(
			33,
			TP_to_F33 * TMath::Sin(Angle_frag) + F33_offset*TMath::Cos(14. * TMath::Pi() / 180.) + fib33x_off,
			0,
			Target_to_TP + TP_to_F33 * TMath::Cos(Angle_frag) + F33_offset*TMath::Sin(14. * TMath::Pi() / 180.) + fib33z_off);

	SetPositionFib(
			31,
			TP_to_F31 * TMath::Sin(Angle_frag) + F31_offset*TMath::Cos(14. * TMath::Pi() / 180.) + fib31x_off,
			0,
			Target_to_TP + TP_to_F31 * TMath::Cos(Angle_frag) + F31_offset*TMath::Sin(14. * TMath::Pi() / 180.) + fib31z_off);

	SetEulerAnglesfib(0, Angle_frag, 0);

	TString base_MDF = "/u/mxarepe/worky/PhD_analysis/MDF/MDF_VERTEX/output/";
	TString MDF_PoQ_file = base_MDF + "MDF_PoQ_2400A_20230627.txt";
	TString MDF_TX0_file = base_MDF + "MDF_TX0_2400A_20230627.txt";
	TString MDF_TY0_file = base_MDF + "MDF_TY0_2400A_20230627.txt";
	TString MDF_FlightPath_file = base_MDF + "MDF_FlightPath_2400A_20230627.txt";

	TString base_MDF_foot = "/u/mxarepe/worky/PhD_analysis/MDF/MDF_RK/output/";
	TString MDF_PoQ_foot_file = base_MDF_foot + "MDF_PoQ_2400A_20230711.txt";
	TString MDF_TX0_foot_file = base_MDF_foot + "MDF_TX0_2400A_20230711.txt";
	TString MDF_TY0_foot_file = base_MDF_foot + "MDF_TY0_2400A_20230711.txt";
	TString MDF_FlightPath_foot_file = base_MDF_foot + "MDF_FlightPath_2400A_20230711.txt";

	MDF_PoQ = new R3BMDFWrapper(MDF_PoQ_file);

	MDF_TX0 = new R3BMDFWrapper(MDF_TX0_file);

	MDF_TY0 = new R3BMDFWrapper(MDF_TY0_file);

	MDF_FlightPath = new R3BMDFWrapper(MDF_FlightPath_file);

	MDF_PoQ_foot = new R3BMDFWrapper(MDF_PoQ_foot_file);

	MDF_TX0_foot = new R3BMDFWrapper(MDF_TX0_foot_file);

	MDF_TY0_foot = new R3BMDFWrapper(MDF_TY0_foot_file);

	MDF_FlightPath_foot = new R3BMDFWrapper(MDF_FlightPath_foot_file);

	Double_t mdf_data[8]={0};   // data container for the MDF function
	Double_t mdf_data_foot[8]={0};   // data container for the MDF function

	TRandom rand(time(0));
	int limit =0;
	int counter0 = 0, counter1 = 0, counter2 = 0, counter3 = 0, counter4 = 0, counter5 = 0, counter6 = 0, counter7 = 0, counter8 = 0, counter9 = 0;
	for(int run = 140; run <=143; run++){
		if(run == 140){
			limit = 48;
			limit = 1;
		}
		if(run == 141){
			continue;
			limit = 4;
		}
		if(run == 142){
			continue;
			limit = 61;
		}
		if(run == 143){
			continue;
			limit = 253;
		}
		for(int nv =1; nv<=limit; nv++){

			TString filename_hit = TString::Format("/lustre/r3b/202205_s522/rpc/main%04d_%04d_TS.root",run, nv);
			TFile *inputFile1 = TFile::Open(filename_hit);
			TTree *tree1 = (TTree *)inputFile1->Get("evt");
			cout << "analysing " << filename_hit << endl;
			TString output_filename_hit = TString::Format("main%04d_%04d_bug_fixed_new_align.root",run, nv);
			//TString output_filename_hit = TString::Format("/lustre/land/mxarepe/MDF_tracked_files/Carbon_12/new_stitch/New_Frag/main%04d_%04d_bug_fixed_new_align.root",run, nv);
			TFile *_file2 = TFile::Open(output_filename_hit,"RECREATE");
			cout << "output file named : " << output_filename_hit << endl;
			TTree *outtree = new TTree("evt","General Tree");

			tree1->SetBranchStatus("*", 0);
			tree1->SetBranchStatus("EventHeader.*", 1);
			tree1->SetBranchStatus("EventHeader.*", 1);
			tree1->SetBranchStatus("TofdHit.*", 1);
			tree1->SetBranchStatus("LosHit.*", 1);
			tree1->SetBranchStatus("Mwpc0HitData.*", 1);
			tree1->SetBranchStatus("Mwpc1HitData.*", 1);
			tree1->SetBranchStatus("Fi30Hit.*", 1);
			tree1->SetBranchStatus("Fi31Hit.*", 1);
			tree1->SetBranchStatus("Fi32Hit.*", 1);
			tree1->SetBranchStatus("Fi33Hit.*", 1);
			tree1->SetBranchStatus("FootHitData.*", 1);

			tree1->SetBranchAddress("EventHeader.", &fEventHeader);

			TBranch  *loshitBranch = tree1->GetBranch("LosHit");
			loshitBranch->SetAddress(&losCA);

			TBranch  *tofdhitBranch = tree1->GetBranch("TofdHit");
			tofdhitBranch->SetAddress(&tofdCA);

			TBranch  *foothitBranch = tree1->GetBranch("FootHitData");
			foothitBranch->SetAddress(&footCA);

			TBranch  *mwpc0hitBranch = tree1->GetBranch("Mwpc0HitData");
			mwpc0hitBranch->SetAddress(&mwpc0CA);

			TBranch  *mwpc1hitBranch = tree1->GetBranch("Mwpc1HitData");
			mwpc1hitBranch->SetAddress(&mwpc1CA);

			TBranch  *fib30hitBranch = tree1->GetBranch("Fi30Hit");
			fib30hitBranch->SetAddress(&fib30CA);

			TBranch  *fib31hitBranch = tree1->GetBranch("Fi31Hit");
			fib31hitBranch->SetAddress(&fib31CA);

			TBranch  *fib32hitBranch = tree1->GetBranch("Fi32Hit");
			fib32hitBranch->SetAddress(&fib32CA);

			TBranch  *fib33hitBranch = tree1->GetBranch("Fi33Hit");
			fib33hitBranch->SetAddress(&fib33CA);

			outtree->Branch("good_incoming",&good_incoming);
			outtree->Branch("mul_frag_vertex",&mul_frag_vertex,"mul_frag_vertex/D");
			outtree->Branch("mul_frag_foot",&mul_frag_foot,"mul_frag_foot/D");

			outtree->Branch("pbeam", &pbeam_tree);
			outtree->Branch("pfrag_mdf_foot", &pfrag_mdf_foot_tree);
			outtree->Branch("pfrag_foot", &pfrag_foot_tree);
			outtree->Branch("pfrag_mdf_vertex", &pfrag_mdf_vertex_tree);

			outtree->Branch("pbeam_lab", &pbeam_lab);
			outtree->Branch("pfrag_mdf_foot_lab", &pfrag_mdf_foot_lab);
			outtree->Branch("pfrag_foot_lab", &pfrag_foot_lab);
			outtree->Branch("pfrag_mdf_vertex_lab", &pfrag_mdf_vertex_lab);

			outtree->Branch("Q",&Q_tree,"Q/D");

			outtree->Branch("AoZ_vertex",&AoZ_vertex,"AoZ_vertex/D");
			outtree->Branch("AoZ_foot",&AoZ_foot,"AoZ_foot/D");
			
			outtree->Branch("vertex_x_foot",&vertex_x_foot,"vertex_x_foot/D");
			outtree->Branch("vertex_y_foot",&vertex_y_foot,"vertex_y_foot/D");
			outtree->Branch("vertex_z_foot",&vertex_z_foot,"vertex_z_foot/D");

			outtree->Branch("tx_mwpc",&tx_mwpc,"tx_mwpc/D");
			outtree->Branch("ty_mwpc",&ty_mwpc,"ty_mwpc/D");

			outtree->Branch("tx_mdf",&tx_mdf,"tx_mdf/D");
			outtree->Branch("ty_mdf",&ty_mdf,"ty_mdf/D");

			outtree->Branch("tx_foot",&tx_foot,"tx_foot/D");
			outtree->Branch("ty_foot",&ty_foot,"ty_foot/D");

			outtree->Branch("vertex_x",&vertex_x,"vertex_x/D");
			outtree->Branch("vertex_y",&vertex_y,"vertex_y/D");
			outtree->Branch("vertex_z",&vertex_z,"vertex_z/D");

			outtree->Branch("rec_vertex",&rec_vertex);
			outtree->Branch("min_dist",&min_dist,"min_dist/D");
			outtree->Branch("f1_Q",&f1_Q,"f1_Q/D");

			Int_t nevt = tree1->GetEntries();
			Int_t nTofd=0, nLos=0, nFoot=0, nMwpc0=0, nMwpc1=0, nFrag = 0;
			//Create foot vectors
			vector <int> f1_hits, f2_hits, f15_hits, f16_hits;

			Double_t mdf_data[8]={0};   // data container for the MDF function
			for(Int_t evt = 0; evt < nevt; evt++)
			{
				reset_val();

				good_incoming=false;
				// Clear Hit level TCLones array
				losCA->Clear();
				tofdCA->Clear();
				footCA->Clear();
				mwpc1CA->Clear();
				mwpc0CA->Clear();
				//clear foot vectors
				f1_hits.clear(); f2_hits.clear(), f15_hits.clear(); f16_hits.clear(); 
				tree1->GetEntry(evt);
				counter0++;
				if(fEventHeader->GetTpat() & 0xf000 || fEventHeader->GetTpat() < 1)
				{
					outtree->Fill();
					continue;
				}
				counter1++;
				// hit level //
				// get event multiplicity for all dets
				nTofd = tofdCA->GetEntries();
				nLos = losCA->GetEntries();
				nFoot = footCA->GetEntries();
				nMwpc0 = mwpc0CA->GetEntries();
				nMwpc1 = mwpc1CA->GetEntries();
				////////////// GOOD incoming events ///////////////////////////////////////////////////////////////////////
				LOS_mul->Fill(nLos);
				//get mul1 events where los trigger time is not nan
				if(nLos != 1 || isnan(fEventHeader->GetTStart())){
					outtree->Fill();
					continue;
				}
				counter2++;
				// Los events with good charge
				bool los_bad_charge = false;
				for (Int_t j = 0; j<nLos; j++){
					auto los_hit = dynamic_cast<R3BLosHitData*>(losCA->At(j));
					if(los_hit->GetZ() < 5.5 || los_hit->GetZ() > 7.1){
						los_bad_charge =true;
						continue;
					}
				}
				if(los_bad_charge){
					outtree->Fill();
					continue;
				}
				counter3++;
				if(nMwpc0 == 0 || nMwpc1 ==0){
					outtree->Fill();
					continue;
				}
				counter4++;
				Double_t mw1_x = 0., mw1_y = 0., mw1_z = -223.8584 - 44.1;
				TVector3 mwpc1, mwpc0;
				auto *hit_m1 = (R3BMwpcHitData*)mwpc1CA->At(0);
				mw1_x = hit_m1->GetX();
				mw1_y = hit_m1->GetY();
				mw1_x = mw1_x*0.1 + offset_x_mwpc1;
				mw1_y = mw1_y*0.1 + offset_y_mwpc1;

				Double_t mw0_x = 0., mw0_y = 0., mw0_z = -104.7107 - 44.1;

				auto *hit_m0 = (R3BMwpcHitData*)mwpc0CA->At(0);
				mw0_x = hit_m0->GetX();
				mw0_y = hit_m0->GetY();
				mw0_x = mw0_x*0.1 + offset_x_mwpc0;
				mw0_y = mw0_y*0.1 + offset_y_mwpc0;

				if(isnan(mw0_x) || isnan(mw0_y) || isnan(mw1_x) || isnan(mw1_y)){
					outtree->Fill();
					continue;
				}
				counter5++;

				if(fabs(mw0_x) > 50 || fabs(mw0_y) > 50 || fabs(mw1_x) > 50 || fabs(mw1_y) > 50){
					outtree->Fill();
					continue;
				}
				counter6++;

				mwpc0.SetXYZ(mw0_x, mw0_y, mw0_z);
				mwpc1.SetXYZ(mw1_x, mw1_y, mw1_z);

				tx_mwpc = (mwpc0.X() - mwpc1.X())/(mwpc0.Z() - mwpc1.Z());
				ty_mwpc = (mwpc0.Y() - mwpc1.Y())/(mwpc0.Z() - mwpc1.Z());

				if(fabs(mwpc1.X() - tx_mwpc * mwpc1.Z() - 0.5) > 1 || fabs(mwpc1.Y() - ty_mwpc * mwpc1.Z()) > 0.7){
					outtree->Fill();
					continue;
				}

				vertex_x = mwpc1.X() - tx_mwpc * mwpc1.Z(); 
				vertex_y = mwpc1.Y() - ty_mwpc * mwpc1.Z();
				vertex_z = 0;

				counter7++;
				TLorentzVector PBeam, Ptg;
				TVector3 pbeam;
				TVector3 Z_lab(0,0,1);
				pbeam.SetXYZ(tx_mwpc,ty_mwpc,1);

				TVector3 Bz = pbeam.Unit();
				TVector3 By = Z_lab.Cross(pbeam.Unit()).Unit();
				TVector3 Bx = By.Cross(Bz);
				Ptg.SetXYZM(0.,0.,0.,mp);
				
				double beta_a_eloss = CalculateEnergyLoss(0.9042525, 6, 1000 * m_12C);
				double gamma_a_eloss = 1.0 / sqrt(1 - pow(beta_a_eloss, 2));

				pbeam.SetMag(m_12C*beta_a_eloss*gamma_a_eloss);
				pbeam_lab.SetVectM(pbeam,m_12C);

				//rotate pbeam
				TVector3 pbeam_r(pbeam.Dot(Bx),pbeam.Dot(By),pbeam.Dot(Bz));
				PBeam.SetVectM(pbeam_r,m_12C);
				pbeam_tree = PBeam;
				good_incoming = true;
				////////////////////////////////////////////////////////////////////////////////////////////////////////////
				
				//Tofd data
				double tofdq1=-99;
				double tofdq2=-101;
				int mul_tofd1=0;
				int mul_tofd2=0;
				bool is_good_tofd = false;
				bool is_good_tofd1 = false;
				bool is_good_tofd2 = false;
				for (auto i = 0; i < nTofd; ++i)
				{
					auto tofd_hit = static_cast<R3BTofdHitData*>(tofdCA->At(i));
					if (tofd_hit->GetDetId() == 1 && tofd_hit->GetTof() > 22 && tofd_hit->GetTof() < 32) // only hits from first plane, add Z later
					{
						if(tofd_hit->GetEloss() > 1.6){
							is_good_tofd1 = true;
							mul_tofd1++;
							tofdq1=tofd_hit->GetEloss();
						}
					}
					if (tofd_hit->GetDetId() == 2 && tofd_hit->GetTof() > 22 && tofd_hit->GetTof() < 32) // only hits from first plane, add Z later
					{
						if(tofd_hit->GetEloss() > 1.6){
							is_good_tofd2 = true;
							mul_tofd2++;
							tofdq2=tofd_hit->GetEloss();
						}
					}
				}
				if(is_good_tofd1 && is_good_tofd2 && mul_tofd1==1 && mul_tofd2==1 && fabs(tofdq1 - tofdq2)< 0.51){
					Q_tree = (tofdq1 + tofdq2)/2.;
					is_good_tofd = true;
				}
				if(is_good_tofd1 && !is_good_tofd2 && mul_tofd1==1){
					Q_tree = tofdq1;
					is_good_tofd = true;
				}
				if(!is_good_tofd1 && is_good_tofd2 && mul_tofd2==1){
					Q_tree = tofdq2;
					is_good_tofd = true;
				}
				if(!is_good_tofd){ 
					tofd_Q1_vs_Q2->Fill(tofdq1,tofdq2);
					outtree->Fill();
					continue;
				}
				tofd_Q1_vs_Q2->Fill(tofdq1,tofdq2);
				counter8++;
				if(!MakeOutgoingTracks()){ 
					outtree->Fill();
					continue;//at least one good track candidate in Fibers
				}
				counter9++;

				mul_frag_vertex=0;
				mul_frag_foot=0;

				TVector3 pfrag_mdf_vertex, pfrag_mdf_vertex_r;
				TVector3 pfrag_mdf_foot, pfrag_mdf_foot_r;
				TVector3 pfrag_foot, pfrag_foot_r;
				TLorentzVector PFrag_foot, PFrag_mdf_foot, PFrag_mdf_vertex;

				double mass = my_nan;

				double px_temp_foot = my_nan, py_temp_foot = my_nan, pz_temp_foot = my_nan, tx_temp = my_nan, ty_temp = my_nan;
				double px_temp_vertex = my_nan, py_temp_vertex = my_nan, pz_temp_vertex = my_nan;

				bool frag_exist = false;

				for(auto & tout : tracks_out){
					if(MakeIncomingTracks()){ 
						for(auto & tin : tracks_in){
							mdf_data_foot[0] = tin.f2_y;
							mdf_data_foot[1] = tin.f2_z;
							mdf_data_foot[2] = tin.f1_x;
							mdf_data_foot[3] = tin.f1_z;
							mdf_data_foot[4] = tout.f32_x;
							mdf_data_foot[5] = tout.f32_z;
							mdf_data_foot[6] = (tout.last_x - tout.f32_x)/(tout.last_z - tout.f32_z);
							mdf_data_foot[7] = (tout.f30_y  - vertex_y)/(tout.f30_z - vertex_z);
							// Calculate all required MDF values

							double flight_p = MDF_FlightPath_foot->MDF(mdf_data_foot);
							double poq = MDF_PoQ_foot->MDF(mdf_data_foot) * 2545./2400.;
							double tx0 = MDF_TX0_foot->MDF(mdf_data_foot);
							double ty0 = MDF_TY0_foot->MDF(mdf_data_foot);
							double tof = flight_p / beta_a_eloss / SPEED_OF_LIGHT;

							double aoz = poq / beta_a_eloss / gamma_a_eloss / AMU;
							if((int) ((aoz * (int) (Q_tree+0.5))+0.5) > 33 || (int) (Q_tree+0.5) > 12 || (int) (Q_tree+0.5) == 0){
								continue;
							}
							AoZ_foot = aoz;
							poq_vec_foot.SetXYZ(tx0,ty0,1);
							poq_vec_foot.SetMag(poq);

							px_temp_foot = poq_vec_foot.X()*((int) (Q_tree+0.5));
							py_temp_foot = poq_vec_foot.Y()*((int) (Q_tree+0.5));
							pz_temp_foot = poq_vec_foot.Z()*((int) (Q_tree+0.5));
							
							tx_temp = tin.tx_foot;
							ty_temp = tin.ty_foot;

							mass = Nuke_Mass_Tab[(int) ((AoZ_foot * (int) (Q_tree+0.5))+0.5)][(int) (Q_tree+0.5)]/1000. - me*( (int) (Q_tree+0.5));
							pfrag_mdf_foot.SetXYZ(px_temp_foot,py_temp_foot,pz_temp_foot);
							pfrag_mdf_foot.SetMag(sqrt(pow(px_temp_foot,2) + pow(py_temp_foot,2) + pow(pz_temp_foot,2)));
							pfrag_mdf_foot_r.SetXYZ(pfrag_mdf_foot.Dot(Bx),pfrag_mdf_foot.Dot(By),pfrag_mdf_foot.Dot(Bz));
							pfrag_mdf_foot_r.SetMag(sqrt(pow(px_temp_foot,2) + pow(py_temp_foot,2) + pow(pz_temp_foot,2)));

							pfrag_foot.SetXYZ(tx_temp,ty_temp,1);
							pfrag_foot.SetMag(sqrt(pow(px_temp_foot,2) + pow(py_temp_foot,2) + pow(pz_temp_foot,2)));
							pfrag_foot_r.SetXYZ(pfrag_foot.Dot(Bx),pfrag_foot.Dot(By),pfrag_foot.Dot(Bz));
							pfrag_foot_r.SetMag(sqrt(pow(px_temp_foot,2) + pow(py_temp_foot,2) + pow(pz_temp_foot,2)));
							/// TLorentz
							//using_angles foot
							PFrag_foot.SetVectM(pfrag_foot_r,mass);
							pfrag_foot_lab.SetVectM(pfrag_foot,mass);
							//using_angles mdf
							PFrag_mdf_foot.SetVectM(pfrag_mdf_foot_r,mass);
							pfrag_mdf_foot_lab.SetVectM(pfrag_mdf_foot,mass);
							//Boost
							PFrag_foot.Boost(-PBeam.BoostVector()); 
							pfrag_foot_tree= PFrag_foot;
							PFrag_mdf_foot.Boost(-PBeam.BoostVector()); 
							pfrag_mdf_foot_tree= PFrag_mdf_foot;
							tx_foot = tx_temp;
							ty_foot = ty_temp;

							mul_frag_foot++;

							TVector3 slope1(tx_temp,ty_temp,1);
							TVector3 slope2(tx_mwpc,ty_mwpc,1);
							TVector3 intercept2(vertex_x,vertex_y,vertex_z);
							TVector3 intercept1(tin.vertex_x_foot,tin.vertex_y_foot,tin.vertex_z_foot);
							vertex _vertex = get_vertex(slope2, slope1, intercept2, intercept1);
							rec_vertex = _vertex.vertexPos;
							vertex_x_foot = rec_vertex.X();
							vertex_y_foot = rec_vertex.Y();
							vertex_z_foot = rec_vertex.Z();
							min_dist = _vertex.min_dist;
							f1_Q=tin.f1_Q;
						}
					}
					mdf_data[0] = vertex_x;
					mdf_data[1] = vertex_y;
					mdf_data[2] = vertex_z;
					mdf_data[3] = tout.f32_x;
					mdf_data[4] = tout.f32_z;
					mdf_data[5] = (tout.last_x - tout.f32_x)/(tout.last_z - tout.f32_z);
					mdf_data[6] = (tout.f30_y  - vertex_y)/(tout.f30_z - vertex_z);
					// Calculate all required MDF values

					double flight_p = MDF_FlightPath->MDF(mdf_data);
					double poq = MDF_PoQ->MDF(mdf_data) * 2545./2400.;
					double tx0 = MDF_TX0->MDF(mdf_data);
					double ty0 = MDF_TY0->MDF(mdf_data);
					double tof = flight_p / beta_a_eloss / SPEED_OF_LIGHT;
					tx_mdf = tx0;
					ty_mdf = ty0;
					double aoz = poq / beta_a_eloss / gamma_a_eloss / AMU;

					if(tracks_out.size() >1){
					//	cout << "\n (Mass: " << (int) ((aoz * (int) (Q_tree+0.5))+0.5) << ", Charge: " << (int) (Q_tree+0.5) << ", AoZ: " << aoz <<", TY0: " << ty0 << ", TX0: " << tx0 <<")" << endl; 
					}

					AoZ_vertex = aoz;
					poq_vec_vertex.SetXYZ(tx0,ty0,1);
					poq_vec_vertex.SetMag(poq);
					
					if((int) ((aoz * (int) (Q_tree+0.5))+0.5) > 33 || (int) (Q_tree+0.5) > 12 || (int) (Q_tree+0.5) == 0){
						continue;
					}

					mass = Nuke_Mass_Tab[(int) ((aoz * (int) (Q_tree+0.5))+0.5)][(int) (Q_tree+0.5)]/1000. - me*( (int) (Q_tree+0.5));
					double energy = sqrt(pow(poq*((int) (Q_tree+0.5)),2) + pow(mass,2));
					double beta_out = poq*((int) (Q_tree+0.5))/energy;

					std::string title_fit_Si = "FIT_Si_Mass_" + std::to_string((int) ((aoz * (int) (Q_tree+0.5))+0.5)) + "_Q_" + std::to_string((int) (Q_tree+0.5));
					std::string title_fit_H2 = "FIT_H2_Mass_" + std::to_string((int) ((aoz * (int) (Q_tree+0.5))+0.5)) + "_Q_" + std::to_string((int) (Q_tree+0.5));
					auto it_Si = functionMap_Si.find(title_fit_Si);
					auto it_H2 = functionMap_H2.find(title_fit_H2);
					double beta_in = 0;
					if (it_Si != functionMap_Si.end() || it_H2 != functionMap_H2.end() ) {
						TF1 *linearFit = it_Si->second;
						TF1 *quadFit = it_H2->second;
						beta_in = linearFit->Eval(beta_out); 
						beta_in = linearFit->Eval(beta_in); 
						beta_in = linearFit->Eval(beta_in); 
						beta_in = linearFit->Eval(beta_in); 
						beta_in = quadFit->Eval(beta_in); 
					}
					else{
						mul_frag_vertex++;
						continue;
					}

					double gamma_in = 1.0 / sqrt(1-beta_in*beta_in);

					double mom_in = beta_in * gamma_in * mass;

					poq_vec_vertex.SetXYZ(tx0,ty0,1);
					poq_vec_vertex.SetMag(mom_in);

					px_temp_vertex = poq_vec_vertex.X();
					py_temp_vertex = poq_vec_vertex.Y();
					pz_temp_vertex = poq_vec_vertex.Z();


					pfrag_mdf_vertex.SetXYZ(px_temp_vertex,py_temp_vertex,pz_temp_vertex);
					pfrag_mdf_vertex.SetMag(sqrt(pow(px_temp_vertex,2) + pow(py_temp_vertex,2) + pow(pz_temp_vertex,2)));
					pfrag_mdf_vertex_r.SetXYZ(pfrag_mdf_vertex.Dot(Bx),pfrag_mdf_vertex.Dot(By),pfrag_mdf_vertex.Dot(Bz));
					pfrag_mdf_vertex_r.SetMag(sqrt(pow(px_temp_vertex,2) + pow(py_temp_vertex,2) + pow(pz_temp_vertex,2)));
					/// TLorentz
					PFrag_mdf_vertex.SetVectM(pfrag_mdf_vertex_r,mass);
					pfrag_mdf_vertex_lab.SetVectM(pfrag_mdf_vertex,mass);
					//Boost
					PFrag_mdf_vertex.Boost(-PBeam.BoostVector()); 
					pfrag_mdf_vertex_tree = PFrag_mdf_vertex;
					mul_frag_vertex++;
				}
				outtree->Fill();
			}
			inputFile1->Close();
			delete inputFile1;

			_file2->cd();
			outtree->Write();

			LOS_mul->Write(); 
			tofd_Q1_vs_Q2->Write(); 
			fiber31_time->Write();	
			fiber31_time_vs_pos->Write();	
			fiber30_time->Write();	
			fiber30_time_vs_pos->Write();	
			fiber32_time->Write();	
			fiber32_time_vs_pos->Write();	
			fiber33_time->Write();	
			fiber33_time_vs_pos->Write();	
			fiber31_cluster->Write();	
			fiber30_cluster->Write();	
			fiber32_cluster->Write();	
			fiber33_cluster->Write();	

			delete _file2;
			cout << "per file : " << counter0 << " " << counter1 << " " << counter2 << " " << counter3 << " " << counter4 << " " << counter5 << " " << counter6 << " " << counter7 << " " << counter8 << " " << counter9 << endl;
		}
		cout << " per run : " << counter0 << " " << counter1 << " " << counter2 << " " << counter3 << " " << counter4 << " " << counter5 << " " << counter6 << " " << counter7 << " " << counter8 << " " << counter9 << endl;
	}
	cout << "final : " << counter0 << " " << counter1 << " " << counter2 << " " << counter3 << " " << counter4 << " " << counter5 << " " << counter6 << " " << counter7 << " " << counter8 << " " << counter9 << endl;
}
void TransformPoint(TVector3& point, TVector3* rot, TVector3* trans)
{
	r.SetToIdentity();
	// First Euler rotation around Y axis
	r.RotateY(rot->Y());
	// get local X axis after first rotation
	v3_localX.SetMagThetaPhi(1, r.ThetaX(), r.PhiX());
	// Second Euler rotation around local X axis
	r.Rotate(rot->X(), v3_localX);
	// get local Z axis after second rotation
	v3_localZ.SetMagThetaPhi(1, r.ThetaZ(), r.PhiZ());
	// final rotation around local Z axis
	r.Rotate(rot->Z(), v3_localZ);
	point.Transform(r);
	point += (*trans);
	return;
}

bool IsGoodFootHit(R3BFootHitData* fhit)
{
	if(fhit->GetEnergy()>100 && fhit->GetEnergy()<1400)
		return true;
	else 
		return false;
}

bool SortFootData()
{
	f1_hits.clear(); f2_hits.clear(); f15_hits.clear(); f16_hits.clear();
	int mul_foot = footCA->GetEntriesFast();
	if(mul_foot==0) return false;
	for (auto f=0; f<mul_foot; ++f)
	{
		auto foot_hit = static_cast<R3BFootHitData*>(footCA->At(f));
		if(!IsGoodFootHit(foot_hit)) continue;
		switch(foot_hit->GetDetId())
		{
			case 1 ://Second FOOT1 for (Y)
				f1_hits.push_back(f);
				break;
			case 2 ://First FOOT2 (X)
				f2_hits.push_back(f);
				break;
			case 15 ://Last FOOT15 (X) 
				f15_hits.push_back(f);
				break;
			case 16 ://3rd FOOT16 (Y)
				f16_hits.push_back(f);
				break;
		}
	}
	if(f1_hits.empty() || f2_hits.empty() || f15_hits.empty() || f16_hits.empty())     
		return false; 
	return true;
}

bool MakeIncomingTracks()
{
	tracks_in.clear();
	if(!SortFootData()) return false;//at least 1 hit in every FOOT
	double tx0 = -999;
	double ty0 = -999; 
	double dx_vertex = -999;
	double dy_vertex = -999;
	double dx_angle = -999;
	double dy_angle = -999;	

	TVector3 f1_point; //cm
	TVector3 f2_point; //cm
	TVector3 f15_point; //cm
	TVector3 f16_point; //cm

	Track tr;
	//----- Make track candidates in FOOT and project them to the center of the target
	for (auto & f2 : f2_hits){
		auto f2_hit = static_cast<R3BFootHitData*>(footCA->At(f2));
		f2_point.SetXYZ(0, f2_hit->GetPosLab().Y() * 0.1, offset_f2); //cm
		for (auto & f1 : f1_hits){
			auto f1_hit = static_cast<R3BFootHitData*>(footCA->At(f1));
			f1_point.SetXYZ(f1_hit->GetPosLab().X() * 0.1 , 0, offset_f1); //cm
			for (auto & f16 : f16_hits){
				auto f16_hit = static_cast<R3BFootHitData*>(footCA->At(f16));
				f16_point.SetXYZ(0, f16_hit->GetPosLab().Y() * 0.1, offset_f16); //cm
				for (auto & f15 : f15_hits){
					auto f15_hit = static_cast<R3BFootHitData*>(footCA->At(f15));
					f15_point.SetXYZ(f15_hit->GetPosLab().X() * 0.1, 0, offset_f15); //cm
					//project foot track to the center of the target
					tx0 = (f15_point.X() - f1_point.X())/(f15_point.Z() - f1_point.Z());
					ty0 = (f16_point.Y() - f2_point.Y())/(f16_point.Z() - f2_point.Z());
					vertex_x_foot = f1_point.X() - tx0 * f1_point.Z();
					vertex_y_foot = f2_point.Y() - ty0 * f2_point.Z();
					vertex_z_foot = 0;
					//Condition on the vertex matching
					dx_vertex = vertex_x_foot - vertex_x;
					dy_vertex = vertex_y_foot - vertex_y;
					dx_angle = tx0 - tx_mwpc;
					dy_angle = ty0 - ty_mwpc;
					//s522
					if((dx_vertex>-0.6 && dx_vertex<0.6) && (dy_vertex>-0.7 && dy_vertex<0.7)){
						if(abs(dx_angle)<0.02 && abs(dy_angle)<0.04){
							tr.vertex_x_foot = vertex_x_foot;
							tr.vertex_y_foot = vertex_y_foot;
							tr.vertex_z_foot = vertex_z_foot;
							tr.f1_x   = f1_point.X();
							tr.f1_z   = f1_point.Z();
							tr.f2_y   = f2_point.Y();
							tr.f2_z   = f2_point.Z();
							tr.f15_x  = f15_point.X();
							tr.f15_z  = f15_point.Z();
							tr.f16_y  = f16_point.Y();
							tr.f16_z  = f16_point.Z();
							tr.tx_foot  = tx0;
							tr.ty_foot  = ty0;

							tr.f1_Q   = f1_hit->GetEnergy();
							tr.f2_Q   = f2_hit->GetEnergy();
							tr.f15_Q   = f15_hit->GetEnergy();
							tr.f16_Q   = f16_hit->GetEnergy();

							tracks_in.push_back(tr);
						}
						else{
							return false;
						}
					}
					else{
						return false;
					}
				}
			}
		}
	}
	return true;
}

vertex get_vertex(TVector3 slope1, TVector3 slope2, TVector3 intercept1, TVector3 intercept2)
{
	// a->slope , b->intercept
	// x-> XZ plane, y->YZ plane
	// right and left -> after (in beam foot) and before (mw) target

	double a1, a2, b1, b2, ap1, ap2, bp1, bp2;
	a1 = intercept2.X();
	a2 = intercept2.Y();
	b1 = slope2.X();
	b2 = slope2.Y();
	ap1 = intercept1.X();
	ap2 = intercept1.Y();
	bp1 = slope1.X();    
	bp2 = slope1.Y();    

	double alpha, beta, A, B, C;
	alpha = (bp1*(a1-ap1)+bp2*(a2-ap2))/(bp1*bp1 + bp2*bp2 + 1);
	beta = (bp1*b1+bp2*b2+1)/(bp1*bp1 + bp2*bp2 + 1);
	A = beta*(bp1*bp1 + bp2*bp2 + 1) - (bp1*b1 + bp2*b2 + 1);
	B = (b1*b1 + b2*b2 + 1) - beta*(bp1*b1+bp2*b2+1);
	C = beta*(bp1*(ap1-a1) + bp2*(ap2-a2)) - (b1*(ap1-a1) + b2*(ap2-a2));

	double sol1, solf1;
	double x1,y1,z1,xp,yp,zp;

	sol1 = -(A*alpha + C)/(A*beta + B);
	solf1 = alpha + beta* sol1;

	x1 = a1 + b1*sol1;
	y1 = a2 + b2*sol1;
	z1 = sol1;
	xp = ap1 + bp1*solf1;
	yp = ap2 + bp2*solf1;
	zp = solf1;

    	TVector3 Vertex_calc((x1+xp)/2.,(y1+yp)/2.,(z1+zp)/2.);

    	vertex recVertex;
	recVertex.vertexPos = Vertex_calc;
	recVertex.min_dist = TMath::Sqrt((x1-xp)*(x1-xp) + (y1-yp)*(y1-yp) + (z1-zp)*(z1-zp));

    	return recVertex;
}

bool IsGoodFiberHit(R3BFiberMAPMTHitData* fhit)
{
	if((fhit->GetEloss() > 0) && (fhit->GetEloss() < 1000) && 
			(fhit->GetTime() < 20000 && fhit->GetTime()>(-20000) )
	  )
		return true;
	else 
		return false;
}

bool MakeOutgoingTracks()
{
	tracks_out.clear();
	//tracks_out2.clear();
	Track tr;
	double angle_out, f30_slope, f30_offset, track_slope, track_offset;
	std::vector<double> f32x;
	std::vector<double> f30y;
	std::vector<double> flastx;
	std::vector<double> flast2x;
	std::vector<double> f32e;
	std::vector<double> f30e;
	std::vector<double> flaste;
	std::vector<double> flast2e;
	TVector3 f30_edge[2];//to extract z and x in f30
	for (auto i=0; i<fib32CA->GetEntriesFast(); ++i)
	{
		auto f32 = static_cast<R3BFiberMAPMTHitData*>(fib32CA->At(i));
		if(!IsGoodFiberHit(f32)) continue;
		double fitime = f32->GetTime_ns() - fEventHeader->GetTStart();
		//double fitime = fTimeStitch->GetTime(f32->GetTime_ns() - fEventHeader->GetTStart(), "clocktdc", "vftx");
		if((fitime > -190 && fitime < -170))
			//if((fitime > 13460 && fitime < 13490))
		{
			fiber32_time->Fill(fitime);	
			fiber32_time_vs_pos->Fill(fitime,f32->GetX());	
			f32x.push_back(f32->GetX());
			f32e.push_back(f32->GetEloss());
		}
	}
	for (auto i=0; i<fib30CA->GetEntriesFast(); ++i)
	{
		auto f30 = static_cast<R3BFiberMAPMTHitData*>(fib30CA->At(i));
		if(!IsGoodFiberHit(f30)) continue;
		double fitime = f30->GetTime_ns() - fEventHeader->GetTStart();
		//double fitime = fTimeStitch->GetTime(f30->GetTime_ns() - fEventHeader->GetTStart(), "clocktdc", "vftx");
		if((fitime > -185 && fitime < -165))
			//if((fitime > 13460 && fitime < 13490))
		{
			fiber30_time->Fill(fitime);	
			fiber30_time_vs_pos->Fill(fitime,f30->GetX());	
			f30y.push_back(f30->GetY());
			f30e.push_back(f30->GetEloss());
		}
	}
	for (auto i=0; i<fib33CA->GetEntriesFast(); ++i)
	{
		auto f33 = static_cast<R3BFiberMAPMTHitData*>(fib33CA->At(i));
		if(!IsGoodFiberHit(f33)) continue;

		double fitime = f33->GetTime_ns() - fEventHeader->GetTStart();
		//double fitime = fTimeStitch->GetTime(f33->GetTime_ns() - fEventHeader->GetTStart(), "clocktdc", "vftx");
		if((fitime > -200 && fitime < -175))
			//if((fitime > 13460 && fitime < 13490))
		{
			fiber33_time->Fill(fitime);	
			fiber33_time_vs_pos->Fill(fitime,f33->GetX());	
			flastx.push_back(f33->GetX());
			flaste.push_back(f33->GetEloss());
		}
	}
	for (auto i=0; i<fib31CA->GetEntriesFast(); ++i)
	{
		auto f31 = static_cast<R3BFiberMAPMTHitData*>(fib31CA->At(i));
		if(!IsGoodFiberHit(f31)) continue;
		double fitime = f31->GetTime_ns() - fEventHeader->GetTStart();
		//double fitime = fTimeStitch->GetTime(f31->GetTime_ns() - fEventHeader->GetTStart(), "clocktdc", "vftx");
		if((fitime > -200 && fitime < -175))
			//if((fitime > 13460 && fitime < 13490))
		{
			fiber31_time->Fill(fitime);	
			fiber31_time_vs_pos->Fill(fitime,f31->GetX());	
			flast2x.push_back(f31->GetX());
			flast2e.push_back(f31->GetEloss());
		}
	}
	if(f32x.size() < 1 || f30y.size() < 1 || (flastx.size() < 1 && flast2x.size() < 1)){
		return false;
	}
	Double_t cluster2[f32x.size()][f32x.size()][2];
	Bool_t set2[f32x.size()];
	Int_t num_hit2[f32x.size()];
	Int_t hit_clust_id2[f32x.size()];
	Int_t num_clust2 = 0;
	for(Int_t i = 0; i < f32x.size(); i++)
	{
		set2[i] = false;
		num_hit2[i] = 0;
		hit_clust_id2[i] = 0;
		for(Int_t j = 0; j < f32x.size(); j++)
		{
			cluster2[i][j][0] = 0./0.;
			cluster2[i][j][1] = 0./0.;
		}
	}

	for(Int_t i = 0; i < f32x.size(); i++)
	{
		if(!set2[i])
		{
			cluster2[num_clust2][num_hit2[num_clust2]][0] = f32x[i];
			cluster2[num_clust2][num_hit2[num_clust2]][1] = f32e[i];
			num_hit2[num_clust2]++;
			set2[i] = true;
			hit_clust_id2[i] = num_clust2;
			num_clust2++;
		}
		else{
			continue;
		}
		for(Int_t k = 0; k < num_hit2[hit_clust_id2[i]]; k++)
		{
			for(Int_t j = 0; j < f32x.size(); j++)
			{
				if(set2[j])
				{
					continue;
				}
				if(fabs(f32x[j]-cluster2[hit_clust_id2[i]][k][0]) < .25)
				{
					Int_t id = hit_clust_id2[i];
					cluster2[id][num_hit2[id]][0] = f32x[j];
					cluster2[id][num_hit2[id]][1] = f32e[j];
					num_hit2[id]++;
					set2[j] = true;
				}
			}
		}
	}
	Double_t cluster0[f30y.size()][f30y.size()][2];
	Bool_t set0[f30y.size()];
	Int_t num_hit0[f30y.size()];
	Int_t hit_clust_id0[f30y.size()];
	Int_t num_clust0 = 0;
	for(Int_t i = 0; i < f30y.size(); i++)
	{
		set0[i] = false;
		num_hit0[i] = 0;
		hit_clust_id0[i] = 0;
		for(Int_t j = 0; j < f30y.size(); j++)
		{
			cluster0[i][j][0] = 0./0.;
			cluster0[i][j][1] = 0./0.;
		}
	}

	for(Int_t i = 0; i < f30y.size(); i++)
	{
		if(!set0[i])
		{
			cluster0[num_clust0][num_hit0[num_clust0]][0] = f30y[i];
			cluster0[num_clust0][num_hit0[num_clust0]][1] = f30e[i];
			num_hit0[num_clust0]++;
			set0[i] = true;
			hit_clust_id0[i] = num_clust0;
			num_clust0++;
		}
		else
			continue;
		for(Int_t k = 0; k < num_hit0[hit_clust_id0[i]]; k++)
		{
			for(Int_t j = 0; j < f30y.size(); j++)
			{
				if(set0[j])
					continue;

				if(fabs(f30y[j]-cluster0[hit_clust_id0[i]][k][0]) < .25)
				{
					Int_t id = hit_clust_id0[i];
					cluster0[id][num_hit0[id]][0] = f30y[j];
					cluster0[id][num_hit0[id]][1] = f30e[j];
					num_hit0[id]++;
					set0[j] = true;
				}
			}
		}
	}
	Double_t cluster3[flastx.size()][flastx.size()][2];
	Bool_t set3[flastx.size()];
	Int_t num_hit3[flastx.size()];
	Int_t hit_clust_id3[flastx.size()];
	Int_t num_clust3 = 0;
	for(Int_t i = 0; i < flastx.size(); i++)
	{
		set3[i] = false;
		num_hit3[i] = 0;
		hit_clust_id3[i] = 0;
		for(Int_t j = 0; j < flastx.size(); j++)
		{
			cluster3[i][j][0] = 0./0.;
			cluster3[i][j][1] = 0./0.;
		}
	}

	for(Int_t i = 0; i < flastx.size(); i++)
	{
		if(!set3[i])
		{
			cluster3[num_clust3][num_hit3[num_clust3]][0] = flastx[i];
			cluster3[num_clust3][num_hit3[num_clust3]][1] = flaste[i];
			num_hit3[num_clust3]++;
			set3[i] = true;
			hit_clust_id3[i] = num_clust3;
			num_clust3++;
		}
		else
			continue;
		for(Int_t k = 0; k < num_hit3[hit_clust_id3[i]]; k++)
		{
			for(Int_t j = 0; j < flastx.size(); j++)
			{
				if(set3[j])
					continue;

				if(fabs(flastx[j]-cluster3[hit_clust_id3[i]][k][0]) < .25)
				{
					Int_t id = hit_clust_id3[i];
					cluster3[id][num_hit3[id]][0] = flastx[j];
					cluster3[id][num_hit3[id]][1] = flaste[j];
					num_hit3[id]++;
					set3[j] = true;
				}
			}
		}
	}
	Double_t cluster1[flast2x.size()][flast2x.size()][2];
	Bool_t set1[flast2x.size()];
	Int_t num_hit1[flast2x.size()];
	Int_t hit_clust_id1[flast2x.size()];
	Int_t num_clust1 = 0;
	for(Int_t i = 0; i < flast2x.size(); i++)
	{
		set1[i] = false;
		num_hit1[i] = 0;
		hit_clust_id1[i] = 0;
		for(Int_t j = 0; j < flast2x.size(); j++)
		{
			cluster1[i][j][0] = 0./0.;
			cluster1[i][j][1] = 0./0.;
		}
	}

	for(Int_t i = 0; i < flast2x.size(); i++)
	{
		if(!set1[i])
		{
			cluster1[num_clust1][num_hit1[num_clust1]][0] = flast2x[i];
			cluster1[num_clust1][num_hit1[num_clust1]][1] = flast2e[i];
			num_hit1[num_clust1]++;
			set1[i] = true;
			hit_clust_id1[i] = num_clust1;
			num_clust1++;
		}
		else
			continue;
		for(Int_t k = 0; k < num_hit1[hit_clust_id1[i]]; k++)
		{
			for(Int_t j = 0; j < flast2x.size(); j++)
			{
				if(set1[j])
					continue;

				if(fabs(flast2x[j]-cluster1[hit_clust_id1[i]][k][0]) < .25)
				{
					Int_t id = hit_clust_id1[i];
					cluster1[id][num_hit1[id]][0] = flast2x[j];
					cluster1[id][num_hit1[id]][1] = flast2e[j];
					num_hit1[id]++;
					set1[j] = true;
				}
			}
		}
	}
	Double_t fi32xcl[num_clust2];
	Double_t fi30ycl[num_clust0];
	Double_t filastxcl[num_clust3];
	Double_t filast2xcl[num_clust1];

	if(num_clust2 > 1 || num_clust0 > 1 || num_clust3 > 1 || num_clust1 > 1){

//		cout << "cluster 32 : " << num_clust2 << " with : ";
		for(Int_t i = 0; i < num_clust2; i++)
		{
//			cout << num_hit2[i] << " fibers in cluster : " << i << " ";
			fiber32_cluster->Fill(num_hit2[i]); 
		}
//		cout <<	" \n cluster 30 : " << num_clust0 << " with : ";
		for(Int_t i = 0; i < num_clust0; i++)
		{
//			cout << num_hit0[i] << " fibers in cluster : " << i << " ";
			fiber30_cluster->Fill(num_hit0[i]); 
		}
//		cout <<	" \n cluster 33 : " << num_clust3 << " with : "; 
		for(Int_t i = 0; i < num_clust3; i++)
		{
//			cout << num_hit3[i] << " fibers in cluster : " << i << " ";
			fiber33_cluster->Fill(num_hit3[i]); 
		}
//		cout <<	" \n cluster 31 : " << num_clust1 << " with : ";
		for(Int_t i = 0; i < num_clust1; i++)
		{
//			cout << num_hit1[i] << " fibers in cluster : " << i << " ";
			fiber31_cluster->Fill(num_hit1[i]); 
		}
	}
//	cout << endl;
	for(Int_t i = 0; i < num_clust2; i++)
	{
		fi32xcl[i] = 0.;
		Double_t sum_ener = 0.;
		for(Int_t j = 0; j < num_hit2[i]; j++)
		{
			fi32xcl[i] += cluster2[i][j][0]*cluster2[i][j][1];
			sum_ener +=cluster2[i][j][1];
		}
		if(sum_ener > 0.)
			fi32xcl[i] = (double)fi32xcl[i]/sum_ener;
	}
	for(Int_t i = 0; i < num_clust0; i++)
	{
		fi30ycl[i] = 0.;
		Double_t sum_ener = 0.;
		for(Int_t j = 0; j < num_hit0[i]; j++)
		{
			fi30ycl[i] += cluster0[i][j][0]*cluster0[i][j][1];
			sum_ener +=cluster0[i][j][1];
		}
		if(sum_ener > 0.)
			fi30ycl[i] = (double)fi30ycl[i]/sum_ener;
	}
	for(Int_t i = 0; i < num_clust3; i++)
	{
		filastxcl[i] = 0.;
		Double_t sum_ener = 0.;
		for(Int_t j = 0; j < num_hit3[i]; j++)
		{
			filastxcl[i] += cluster3[i][j][0]*cluster3[i][j][1];
			sum_ener +=cluster3[i][j][1];
		}
		if(sum_ener > 0.)
			filastxcl[i] = (double)filastxcl[i]/sum_ener;
	}
	for(Int_t i = 0; i < num_clust1; i++)
	{
		filast2xcl[i] = 0.;
		Double_t sum_ener = 0.;
		for(Int_t j = 0; j < num_hit1[i]; j++)
		{
			filast2xcl[i] += cluster1[i][j][0]*cluster1[i][j][1];
			sum_ener +=cluster1[i][j][1];
		}
		if(sum_ener > 0.)
			filast2xcl[i] = (double)filast2xcl[i]/sum_ener;
	}

	for(auto i = 0; i < num_clust2; i ++)
	{
		f32_point.SetXYZ(fi32xcl[i], 0, 0); //cm
		f32_point_i=f32_point;
		TransformPoint(f32_point, &fib_angles, &f32_position);
		tr.f32_x = f32_point.X();
		tr.f32_z = f32_point.Z();

		for (auto j=0; j<num_clust0; ++j)
		{
			f30_point.SetXYZ(0,fi30ycl[j],0); //cm
			f30_point_i=f30_point;
			TransformPoint(f30_point, &fib_angles, &f30_position);
			tr.f30_y = f30_point.Y();
			//make combination with every hit in fibers 33 and 31

			for (auto k = 0; k<num_clust3; ++k)
			{
				flast_point.SetXYZ(filastxcl[k], 0, 0); //cm
				flast_point_i=flast_point;
				TransformPoint(flast_point, &fib_angles, &f33_position);
				tr.last_x = flast_point.X();
				tr.last_z = flast_point.Z();
				angle_out = TMath::ATan((tr.last_x - tr.f32_x)/(tr.last_z - tr.f32_z)) * TMath::RadToDeg();
				//cout<<"a_33-> "<<angle_out<<endl;
				if(angle_out>(-10.) || angle_out<(-18.)) continue;
				// We need to extrapolate Z position in f30 because it was used for Y measurement
				// Define two (X,Z) points on the f30 plane:
				//Now track every combination of upstream and downstream tracks 
				f30_edge[0].SetXYZ(-1, 0, 0);
				f30_edge[1].SetXYZ(1, 0, 0);
				TransformPoint(f30_edge[0], &fib_angles, &f30_position);
				TransformPoint(f30_edge[1], &fib_angles, &f30_position);
				// Parameterize f30 plane
				f30_slope = (f30_edge[1].X() - f30_edge[0].X()) / (f30_edge[1].Z() - f30_edge[0].Z());
				f30_offset = f30_edge[0].X() - f30_slope * f30_edge[0].Z();
				track_slope  = (tr.last_x - tr.f32_x) / (tr.last_z - tr.f32_z);
				track_offset = (tr.last_x - track_slope * tr.last_z);
				// Extrapolate final X and Z position in f30
				tr.f30_z = ((track_offset - f30_offset) / (f30_slope - track_slope));// extrapolated
				tr.f30_x = (track_slope * tr.f30_z + track_offset);// extrapolated
				tracks_out.push_back(tr);
			}

			for (auto k = 0; k<num_clust1; ++k)
			{
				flast_point.SetXYZ(filast2xcl[k], 0, 0); //cm
				flast_point_i=flast_point;
				TransformPoint(flast_point, &fib_angles, &f31_position);
				tr.last_x = flast_point.X();
				tr.last_z = flast_point.Z();
				angle_out = TMath::ATan((tr.last_x - tr.f32_x)/(tr.last_z - tr.f32_z)) * TMath::RadToDeg();
				//cout<<"a_31-> "<<angle_out<<endl;
				if(angle_out>(-10.) || angle_out<(-18.)) continue;
				// We need to extrapolate Z position in f30 because it was used for Y measurement
				// Define two (X,Z) points on the f30 plane:
				//Now track every combination of upstream and downstream tracks 
				f30_edge[0].SetXYZ(-1, 0, 0);
				f30_edge[1].SetXYZ(1, 0, 0);
				TransformPoint(f30_edge[0], &fib_angles, &f30_position);
				TransformPoint(f30_edge[1], &fib_angles, &f30_position);
				// Parameterize f30 plane
				f30_slope = (f30_edge[1].X() - f30_edge[0].X()) / (f30_edge[1].Z() - f30_edge[0].Z());
				f30_offset = f30_edge[0].X() - f30_slope * f30_edge[0].Z();
				track_slope  = (tr.last_x - tr.f32_x) / (tr.last_z - tr.f32_z);
				track_offset = (tr.last_x - track_slope * tr.last_z);
				// Extrapolate final X and Z position in f30
				tr.f30_z = ((track_offset - f30_offset) / (f30_slope - track_slope));// extrapolated
				tr.f30_x = (track_slope * tr.f30_z + track_offset);// extrapolated
				tracks_out.push_back(tr);
			}
		}
	}
	if(tracks_out.empty()) return false;
	return true;
}

// Function to calculate energy loss for different materials
double CalculateEnergyLoss(double beta, double Z, double mass) {
	// Define the materials with their properties
	std::vector<Material> MWPC1Gas = {
		{18/39.95, 0.000188, 0.0017, 0.7841}, // Argon (78.41% by weight)
		{6/12.01, 0.000081, 0.0017, 0.0589},  // Carbon (5.89% by weight)
		{8/16.00, 0.000133, 0.0017, 0.1570}   // Oxygen (15.70% by weight)
	};
	std::vector<Material> Scint = {
		{1/1.0079, 0.0000192, 0.085292}, // H 
		{6/12.0107, 0.000081, 0.914708}   // C 
	};

	std::vector<Material> MusGas = { //P-75
		{1/1.0079, 0.0000192, 0.000957, 0.0739}, // H 
		{6/12.0107, 0.000081, 0.000957, 0.4435},   // C 
		{8/15.9994, 0.000095, 0.000957, 0.0039},   // O 
		{18/39.948, 0.000188, 0.000957, 0.0039}   // Ar 
	};

	std::vector<Material> materials = {
		{0.52037, 78.7e-6, 1.400, 1},// Mylar
		{1/1.008, 21.8e-6, 0.07080, 1}, // Liquid Hydrogen
		{22./47.867, 233.0e-6, 4.540, 1} // Ti
	};

	//std::cout << "Initial: " << beta << std::endl;
	// MWPC1 Window Mylar
	beta -= EnergyLossBetheBloch(beta, {materials[0]}, 0.01, Z, mass); 
	//std::cout << "MWPC1 Window Mylar: " << beta << std::endl;

	// MWPC1 Gas
	beta -= EnergyLossBetheBloch(beta, MWPC1Gas, 1, Z, mass); 
	//std::cout << "MWPC1 Gas: " << beta << std::endl;

	// MWPC1 Exit Window Mylar
	beta -= EnergyLossBetheBloch(beta, {materials[0]}, 0.02, Z, mass); 
	//std::cout << "MWPC1 Exit Window Mylar: " << beta << std::endl;

	// LOS Scintillator 3mm
	beta -= EnergyLossBetheBloch(beta, Scint, 0.3, Z, mass); 
	//std::cout << "LOS Scintillator 3mm: " << beta << std::endl;

	// MUSIC Window Mylar
	beta -= EnergyLossBetheBloch(beta, {materials[0]}, 0.01, Z, mass); 
	//std::cout << "MUSIC Window Mylar: " << beta << std::endl;

	// MUSIC Gas Freon (you need to add Freon to the materials list)
	beta -= EnergyLossBetheBloch(beta, MusGas, 550, Z, mass); 
	//std::cout << "MUSIC Gas Freon: " << beta << std::endl;

	// MUSIC Exit Window
	beta -= EnergyLossBetheBloch(beta, {materials[0]}, 0.01, Z, mass); 
	//std::cout << "MUSIC Exit Window: " << beta << std::endl;

	// MWPC0 Window Mylar
	beta -= EnergyLossBetheBloch(beta, {materials[0]}, 0.01, Z, mass); 
	//std::cout << "MWPC0 Window Mylar: " << beta << std::endl;

	// MWPC0 Gas (use the same as MWPC1Gas or define separately)
	beta -= EnergyLossBetheBloch(beta, MWPC1Gas, 1, Z, mass); 
	//std::cout << "MWPC0 Gas: " << beta << std::endl;

	// MWPC0 Exit Window Mylar
	beta -= EnergyLossBetheBloch(beta, {materials[0]}, 0.01, Z, mass); 
	//std::cout << "MWPC0 Exit Window Mylar: " << beta << std::endl;

	// ROLU Window Mylar
	beta -= EnergyLossBetheBloch(beta, {materials[0]}, 0.04, Z, mass); 
	//std::cout << "ROLU Window Mylar: " << beta << std::endl;

	// TC Window Mylar
	beta -= EnergyLossBetheBloch(beta, {materials[0]}, 0.0023, Z, mass); 
      	//std::cout << "TC Window Mylar: " << beta << std::endl;

	// Target Entrance
	beta -= EnergyLossBetheBloch(beta, {materials[0]}, 0.0072, Z, mass); 
	//std::cout << "Target Entrance: " << beta << std::endl;

	// TARGET H
	beta -= EnergyLossBetheBloch(beta, {materials[1]}, 2.5, Z, mass); 
	//std::cout << "TARGET H: " << beta << std::endl;

	return beta;
}

double EnergyLossBetheBloch(double beta, const std::vector<Material>& materials, double thickness, double charge, double mass) {
	// Constants
	double K = 0.307075;   // MeV mol^-1 cm^2 (constant factor)
	double me_MeV = 0.511;     // MeV/c^2 (mass of the electron)
	double gamma = 1.0 / sqrt(1 - pow(beta, 2));
	double beta_gamma = beta * gamma;
	double totalEnergyLoss = 0.0;
	// Loop through each component in the composite material
	for (const auto& component : materials) {
		// Maximum kinetic energy transferable to a free electron in a single collision (Tmax)
		double Tmax = (2 * me_MeV * pow(beta_gamma, 2)) /
			(1 + 2 * gamma * (me_MeV / (mass)) + pow(me_MeV / (mass), 2));

		// Bethe-Bloch formula for energy loss (dE/dx)
		double dEdx = K * component.rho * component.weight_fraction * (component.ZoA) * pow(charge / beta, 2) *
			0.5 * (log((2 * me_MeV * pow(beta_gamma, 2) * Tmax) / pow(component.I, 2)) - pow(beta, 2));

		// Add the energy loss for this component to the total energy loss
		totalEnergyLoss += dEdx * thickness;
	}
	double p_mag = (mass*beta*gamma);
	double etot = TMath::Sqrt(pow(p_mag,2) + mass * mass);
	double d_beta = 1. / beta * TMath::Power(mass / etot, 2) * totalEnergyLoss / etot;
	return d_beta;
}

int main(Int_t argc, Char_t* argv[]){
	TStopwatch t1; 
	t1.Start(); 
	tree_merger();
	cout << "\n-- Run time of tree_marger():\n";
	t1.Print();
	return 0;
}
