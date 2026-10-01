#include"libs.hh"
#include"definitions.hh"
#include"R3BMDFWrapper.h"
#include"R3BTrackingParticle.h"
#include"R3BGladFieldMap.h"
#include"R3BTPropagator.h"

using namespace std;


R3BGladFieldMap* magField = new R3BGladFieldMap("test_field","R");

int main(Int_t argc, Char_t* argv[])
{
	TStopwatch t1; 
	t1.Start(); 
	Bool_t is_PCA = kFALSE;
	Bool_t NeedHelp = kTRUE;

	gRandom = new TRandom3();
	gRandom->SetSeed(0);
	gStyle->SetPalette(kRainBow);

	magField->SetPosition(0.7871, 1.75-1.526, Target_to_GLAD_flange + 54.05 - 0.55580);//x,y,z in cm 
	magField->SetXAngle(-0.113); //deg
	magField->SetYAngle(-14.08); //deg
	magField->SetZAngle(-0.83); //deg
	magField->Init();

	cout << "\n-- Readout time of the field map file:\n";
	t1.Print();

	//magField->Print();
	R3BTPropagator * fPropagator = new R3BTPropagator(magField, kFALSE);

	TApplication* theApp = new TApplication("App", 0, 0);

	cout << "\n\n\t*************************************************" << endl;
	cout << "\t*                                               *" << endl;
	cout << "\t*  Multidimensional Fit for tracking in R3BRoot *" << endl;
	cout << "\t*                                               *" << endl;
	cout << "\t*************************************************" << endl;

	//--------------- Collecting data and perfroming MDF fitting ---------------/

	R3BMDFWrapper * PoQ = new R3BMDFWrapper("../MDF_VERTEX/output/MDF_PoQ_2400A_20230627.txt");
	R3BMDFWrapper * TX0 = new R3BMDFWrapper("../MDF_VERTEX/output/MDF_TX0_2400A_20230627.txt");
	R3BMDFWrapper * TY0 = new R3BMDFWrapper("../MDF_VERTEX/output/MDF_TY0_2400A_20230627.txt");
	//--------------- Performing MDF quality check ---------------/

	for(int nv = 0; nv < 6 ; nv++){
		Double_t GLAD_CURRENT = 2400 /*+ gRandom->Uniform(-600,800)*/; 
		switch(nv)
		{
			case 0 : 
				set_Z(6);
				set_A(12);
				set_EKIN(1.25);
				update_mass();
				GLAD_CURRENT = 2400; 
				break;
			case 1 :
				set_Z(6);
				set_A(16);
				set_EKIN(1.25);
				update_mass();
				GLAD_CURRENT = 2400; 
				break;
			case 2 :
				set_Z(5);
				set_A(11);
				set_EKIN(1.25);
				update_mass();
				GLAD_CURRENT = 2400; 
				break;
			case 3 : 
				set_Z(6);
				set_A(16);
				set_EKIN(0.560);
				update_mass();
				GLAD_CURRENT = 1850; 
				break;
			case 4 :
				set_Z(6);
				set_A(18);
				set_EKIN(0.560);
				update_mass();
				GLAD_CURRENT = 1907; 
				break;
			case 5 :
				set_Z(5);
				set_A(15);
				set_EKIN(0.560);
				update_mass();
				GLAD_CURRENT = 1907; 
				break;
		}
		string title_z = to_string(BEAM_Z);
		title_z.erase(title_z.find_last_not_of('0') +1, std::string::npos);
		title_z.erase(title_z.find_last_not_of('.') +1, std::string::npos);
 		string title_a = to_string(BEAM_A);
		title_a.erase(title_a.find_last_not_of('0') +1, std::string::npos);
		title_a.erase(title_a.find_last_not_of('.') +1, std::string::npos);
 		string title_e = to_string(BEAM_EKIN_GEV_U*1000);
		title_e.erase(title_e.find_last_not_of('0') +1, std::string::npos);
		title_e.erase(title_e.find_last_not_of('.') +1, std::string::npos);
 		string title_g = to_string(GLAD_CURRENT);
		title_g.erase(title_g.find_last_not_of('0') +1, std::string::npos);
		title_g.erase(title_g.find_last_not_of('.') +1, std::string::npos);

		TString title = "Tree_"+ title_g + "_Z_" + title_z + "_A_" + title_a + "_E_" + title_e;
		cout << title << endl;
		UInt_t counter=0;
		Double_t poq_value;
		Double_t tx0_value;
		Double_t ty0_value;
		std::vector<Data> data_array;
		Data data;
		cout << "\n-- Performing MDF quality check after fitting.... " << endl;
		TTree tree_out;
		Float_t poq_mdf={0};
		Float_t tx0_mdf={0};
		Float_t ty0_mdf={0};
		Float_t poq_sim={0};
		Float_t tx0_sim={0};
		Float_t ty0_sim={0};

		TFile *_file1 = TFile::Open(("output/"+ title +".root"),"RECREATE");
		tree_out.SetName(title);

		tree_out.Branch("poq_mdf", &poq_mdf , "poq_mdf/F");
		tree_out.Branch("tx0_mdf", &tx0_mdf , "tx0_mdf/F");
		tree_out.Branch("ty0_mdf", &ty0_mdf , "ty0_mdf/F");
		tree_out.Branch("poq_sim", &poq_sim , "poq_sim/F");
		tree_out.Branch("tx0_sim", &tx0_sim , "tx0_sim/F");
		tree_out.Branch("ty0_sim", &ty0_sim , "ty0_sim/F");

		magField->SetScale((-1)*GLAD_CURRENT/3583.81);
		while(counter<NEVENTS)
		{
			data.primary = Generate_Primary(); 

			if(!ExtractData(data, fPropagator))
			{
				cout << "\n\t-- WARNING! --> Failed to extract track data, trying different primary... \n";
			}
			else
			{
				poq_value = PoQ->MDF(data.edata);
				tx0_value = TX0->MDF(data.edata);
				ty0_value = TY0->MDF(data.edata);
			}
			poq_mdf=poq_value;
			tx0_mdf=tx0_value;
			ty0_mdf=ty0_value;
			poq_sim=data.value[0];
			tx0_sim=data.value[1];
			ty0_sim=data.value[2];
			tree_out.Fill();
			counter++;

		}
		cout << "\n-- FINISH --" << endl;
		tree_out.Write();
		_file1->Close();
	}
}

R3BTrackingParticle Generate_Primary()
{
	const Double_t  Ekin_Total = BEAM_EKIN_GEV_U * BEAM_MASS_GEV_C2/AMU;
	Double_t P, Theta, RCos, Phi, vx, vy, vz, Angle, Rand_R, beta;
	TVector3 PVector;
	P     = sqrt(Ekin_Total*(Ekin_Total+2*BEAM_MASS_GEV_C2)); //initial momentum
	auto Ptot  = P;
	P       = (Ptot+gRandom->Uniform((-1)*Ptot*SCALE_MIN, Ptot*SCALE_MAX));//min and max momentum spread
	while(1)
	{
		RCos    = gRandom->Uniform(TMath::Cos(BEAM_DTHETA_DEG*TMath::Pi()/180.),1.); 
		Theta   = TMath::ACos(RCos); 
		Phi     = gRandom->Uniform(0., 2.*TMath::Pi());
		const double Angle_R   = gRandom->Uniform(0., 2.*TMath::Pi());
		Rand_R  = gRandom->Uniform(0, BEAM_RADIUS_CM);
		vx	= Rand_R * cos(Angle_R);
		vy	= Rand_R * sin(Angle_R);
		vz      = gRandom->Uniform(-BEAM_THICK_CM/2.,BEAM_THICK_CM/2.);//somehow R3BTPropagator fails at z=0
		PVector.SetMagThetaPhi(P, Theta, Phi);
		//Do not use events which may potentially hit the wall of GLAD chamber
		if(fabs(vx + (Target_to_RP-vz) * PVector.X()/PVector.Z() )<13. && 
				fabs(vy + (Target_to_RP-vz) * PVector.Y()/PVector.Z() )<13.)
			break;
	}
	beta = 1. / TMath::Sqrt(1 + TMath::Power(BEAM_MASS_GEV_C2 / PVector.Mag(), 2));
	R3BTrackingParticle ion(BEAM_Z,
			vx,
			vy,
			vz,
			PVector.X(),
			PVector.Y(),
			PVector.Z(),
			beta,
			BEAM_MASS_GEV_C2);
	return ion;
}

//----------- Obtain track points as measured in the experiment ------------
bool ExtractData(Data &d, R3BTPropagator* fPropagator)
{
	double x,x1,x2,z,z1,z2,y,y1,y2,det_to_TP;//coordinates of the detector plane 
	TVector3 position = d.primary.GetStartPosition();
	TVector3 momentum = d.primary.GetStartMomentum();
	TVector3 v1, v2, v3;//points on the the plane for propagation
	R3BTrackingParticle * particle;
	TVector3 points[4];
	TVector3 mom_track[4];
	Double_t length[4];
	Double_t tof[4];
	int i = 0;
	particle = new R3BTrackingParticle(BEAM_Z,
			position.X(),
			position.Y(),
			position.Z(),
			momentum.X(),
			momentum.Y(),
			momentum.Z(),
			d.primary.GetBeta(),
			d.primary.GetMass());
	//Assuming 5 detectors where the trajectories are measured
	points[0] = particle->GetPosition();
	length[0] = particle->GetLength();
	mom_track[0] = particle->GetMomentum();
	for(i=1; i<=3; i++)
	{
		//Caluclate x, and z to define the detector plane
		switch(i)
		{
			case 1 ://First Fi32 (X) 
				det_to_TP = TP_to_F32;
				z1 = Target_to_TP + det_to_TP * TMath::Cos(Angle);
				x1 = 0. - (det_to_TP * TMath::Sin(Angle));
				z2 = z1 + 10 * TMath::Cos(Angle_det);
				x2 = x1 + 10 * TMath::Sin(Angle_det);
				break;

			case 2 ://Second Fi30 (Y)
				det_to_TP = TP_to_F30;
				z1 = Target_to_TP + det_to_TP * TMath::Cos(Angle);
				x1 = 0. - (det_to_TP * TMath::Sin(Angle));
				//2nd point is 10 cm from the central point along det surface
				z2 = z1 + 10 * TMath::Cos(Angle_det);
				x2 = x1 + 10 * TMath::Sin(Angle_det);
				break;

			case 3 ://Last two fibers
				det_to_TP = (TP_to_F31 + TP_to_F33)/2.;
				z1 = Target_to_TP + det_to_TP * TMath::Cos(Angle);
				x1 = 0. - (det_to_TP * TMath::Sin(Angle));
				//2nd point is 10 cm from the central point along det surface
				z2 = z1 + 10 * TMath::Cos(Angle_det);
				x2 = x1 + 10 * TMath::Sin(Angle_det);
				break;
		}

		v1.SetXYZ(x1,0,z1); //central det point 
		v2.SetXYZ(x1,Ymax,z1);//vertical det point above central
		v3.SetXYZ(x2,0,z2);//3rd point to define the plane
		double tof_value=0;

		if(!fPropagator->PropagateToPlaneRK_eloss(particle,v1,v2,v3,tof_value,0.1,true)) 
			//if(!fPropagator->PropagateToPlaneRK(particle,v1,v2,v3)) 
		{
			cout << "\n\t-- WARNING! --> Failed propagation...";
			return false;
		}
		tof[i] = tof_value;
		points[i] = particle->GetPosition();
		length[i] = particle->GetLength();
		mom_track[i] = particle->GetMomentum();
		if(isnan(points[i].X()) || isnan(points[i].Z()) || isnan(points[i].Y()))
		{
			cout << "\n\t-- WARNING! --> Propagated NaN values!";
			return false;
		}
	}
	i=0;
	double offset_target_x = 0;//gRandom->Uniform(-0.005,0.005);
	double offset_target_y = 0;//gRandom->Uniform(-0.005,0.005);
	double offset_target_z = 0;//gRandom->Uniform(-0.005,0.005);

	double pitch_fiber_X_32 =0;// gRandom->Uniform(-0.05,0.05);
	double pitch_fiber_Y_30 =0;// gRandom->Uniform(-0.05,0.05);
	double pitch_fiber_X_31 =0;// gRandom->Uniform(-0.05,0.05);

	d.edata[i]   =  points[0].X() + offset_target_x; 
	d.edata[++i] =  points[0].Y() + offset_target_y; 
	d.edata[++i] =  points[0].Z() + offset_target_z;

	d.edata[++i] =  points[1].X() + pitch_fiber_X_32; 
	d.edata[++i] =  points[1].Z();

	d.edata[++i] = ((points[3].X() + pitch_fiber_X_32) - (points[1].X() + pitch_fiber_X_31))/(points[3].Z() - points[1].Z());
	d.edata[++i] = ((points[2].Y() + pitch_fiber_Y_30) - (points[0].Y() + offset_target_y))/(points[2].Z() - (points[0].Z() + offset_target_z));
	Double_t FlightPath = length[1] - length[0];
	//cout <<  "\n\nFlight Path: " << FlightPath << endl;


	for(i=0; i<nVars; ++i)
	{
		if(TMath::IsNaN(d.edata[i]))
		{
			cout << "\nNaN values for tracking!!";
			return false;
		}
	}
	//PoQ
	d.value[0] = momentum.Mag()/BEAM_Z; 
        d.value[1] = mom_track[0].X()/mom_track[0].Z();
        d.value[2] = mom_track[0].Y()/mom_track[0].Z();
	delete particle;
	return true;
}


