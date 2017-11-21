#include <iostream>
#include <fstream>
#include <algorithm> //fill_n
#include <string>
#include <cstdlib>

#include "TMath.h"
#include "TCutG.h"
#include "TRandom2.h"

#include "ParticleTree.h"
#include "Event.h"
#include "Particle.h"
#include "dEdxCut.h"
#include "PPMCut.h"

using namespace std;

void RunPPMCut(TString inputfile, TString outputfile, TString system, TString energy)
{
	cout << "Running particle population matrix mode" << endl;

	TFile *input_rootfile = new TFile(inputfile);
	TTree* input_tree = (TTree*)input_rootfile->Get("events");

	ParticleTree output_tree(outputfile);

	Event *event = new Event();
	Particle *particle;
	input_tree->SetBranchAddress("event",&event);

	const Long64_t treeNentries = input_tree->GetEntries();
	Long64_t ev;
	UInt_t Npa;
	UInt_t part;
	UInt_t particles_in, particles_out;

	PPMCut partpopmatrix("PartPopMatrix.root",system, energy);

	float pt, p, angle;

	particles_in = particles_out = 0;

	for(ev=0; ev<treeNentries; ++ev)
	{
		if(!(ev%50))
			cout << "Event: " << ev << endl;

		input_tree->GetEntry(ev);
		Npa = event->GetNpa();
		output_tree.BeginEvent();

		particles_in += Npa;

		for(part=0; part<Npa; part++)
		{
			particle = event->GetParticle(part);
			pt = TMath::Sqrt(TMath::Power(particle->GetPx(),2)+TMath::Power(particle->GetPy(),2));
			p = TMath::Sqrt(TMath::Power(particle->GetPx(),2)+TMath::Power(particle->GetPy(),2)+TMath::Power(particle->GetPz(),2));
			angle = TMath::ATan2(particle->GetPy(), particle->GetPx());

			//CIECIE NA AKCEPTACJE
			if(partpopmatrix.PartPopMatrixCut(particle->GetCharge(),p,pt,angle))
			{
				output_tree.AddParticle(particle->GetCharge(),
						particle->GetPx(), particle->GetPy(), particle->GetPz(),
						particle->GetdEdx(), particle->GetdEdxVtpc1(), particle->GetdEdxVtpc2(), particle->GetdEdxMtpc(),
						particle->GetNdEdx(), particle->GetNdEdxVtpc1(), particle->GetNdEdxVtpc2(), particle->GetNdEdxMtpc());

				particles_out++;
			}
		}

		output_tree.EndEvent();
	}

	output_tree.Close();
	input_rootfile->Close();

	cout << "Particle population matrix cut summary\n------------" << endl
		<< "Particles before cut: " << particles_in << endl
		<< "Particles after cut: " << particles_out << endl
		<< "Cutted particles: " << (particles_in-particles_out) << endl
		<< "Ratio: " << ((Double_t)particles_out/particles_in) << endl;
}

Float_t choose_dedx(Particle *particle, TString system)
{
	static Int_t vtpc1_part;
	static Int_t vtpc2_part;
	static Int_t mtpc_part;

	if(!(system.CompareTo("pp")))
		return particle->GetdEdx();
	else if(!(system.CompareTo("PbPb")))
	{
		vtpc1_part = particle->GetNdEdxVtpc1();
		vtpc2_part = particle->GetNdEdxVtpc2();
		mtpc_part = particle->GetNdEdxMtpc();

		//std::cout << "dE/dx: VTPC1 part: " << vtpc1_part << "\tVTPC2 part: " << vtpc2_part << "\tMTPC part: " << mtpc_part << std::endl;
		if((vtpc1_part == 0) && (vtpc2_part == 0) && (mtpc_part == 0))
		{
			//std::cout << "WTF? Particle with no dE/dx information!" << std::endl;
			return 0;
		}
		else
		{
			if(mtpc_part > 0)
				return (particle->GetdEdxMtpc());
			else if(vtpc2_part >= vtpc1_part)
				return (particle->GetdEdxVtpc2());
			else
				return (particle->GetdEdxVtpc1());
		}
	}
	else
	{
		std::cout << "Unrecognized particle system" << std::endl;
		return 0;
	}
}

void RunDedxCut(TString inputfile, TString outputfile, TString system, Int_t energy)
{
	cout << "Running dE/dx mode with energy " << energy << endl;
	TCutG* cutg = initialise_dedx_cutg(system, energy);
	cout << "Graphcut intialized" << endl;
	
	TFile *input_rootfile = new TFile(inputfile);
	TTree* input_tree = (TTree*)input_rootfile->Get("events");

	ParticleTree output_tree(outputfile);

	Event *event = new Event();
	Particle *particle;
	input_tree->SetBranchAddress("event",&event);

	const Long64_t treeNentries = input_tree->GetEntries();
	Long64_t ev;
	Long_t particles_in = 0, particles_out = 0;
	UInt_t Npa;
	UInt_t part;

	Float_t local_dedx;
	Float_t dedx_uppercut = 3.;

	if(!(system.CompareTo("PbPb")))
	{	
		if((energy == 158 ) || (energy == 160))
			dedx_uppercut = 1.65;
		else if(energy == 20)
			dedx_uppercut = 1.6;
	}

	float p;

	cout << "Cut dE/dx > " << dedx_uppercut << " applied" << endl;

	for(ev=0; ev<treeNentries; ++ev)
	{
		if(!(ev%500))
			cout << "Event: " << ev << endl;

		input_tree->GetEntry(ev);
		Npa = event->GetNpa();
		output_tree.BeginEvent();

		for(part=0; part<Npa; part++)
		{
			particles_in++;
			particle = event->GetParticle(part);
			p = TMath::Sqrt(TMath::Power(particle->GetPx(),2)+TMath::Power(particle->GetPy(),2)+TMath::Power(particle->GetPz(),2));
			if(choose_dedx(particle, system) == 0)	//Using global dE/dx values instead of choosing local ones (i.e. by number of dE/dx points left in TPCs)
				continue;
		
			local_dedx = particle->GetdEdx();

			if(cutg->IsInside(p,local_dedx))
				continue;
			if(local_dedx > dedx_uppercut)
				continue;		//dodatkowy cut 3.02.2013

			output_tree.AddParticle(particle->GetCharge(),
					particle->GetPx(), particle->GetPy(), particle->GetPz(),
					particle->GetdEdx(), particle->GetdEdxVtpc1(), particle->GetdEdxVtpc2(), particle->GetdEdxMtpc(),
					particle->GetNdEdx(), particle->GetNdEdxVtpc1(), particle->GetNdEdxVtpc2(), particle->GetNdEdxMtpc());

			particles_out++;
		}
		output_tree.EndEvent();
	}

	output_tree.Close();
	input_rootfile->Close();

	cout << "dEdx cut summary\n------------" << endl
		<< "Particles before cut: " << particles_in << endl
		<< "Particles after cut: " << particles_out << endl
		<< "Cutted particles: " << (particles_in-particles_out) << endl
		<< "Ratio: " << ((Double_t)particles_out/particles_in) << endl;
}

Bool_t is_electron(double logP, double dEdx) {
	static const double a1 = 0.108696;
	static const double b1 = 1.358696;
	static const double a2 = -0.152174;
	//static const double b2 = 1.88;
	static const double b2 = 1.95;
	//static const double b2 = 2;

	if((logP >= 0.2) && (dEdx > 1.78))
		return kTRUE;
	else if((logP > -2) && (logP < 2))
	{
		if((dEdx > a1*logP + b1) && (dEdx < a2*logP + b2))
			return kTRUE;
		else
			return kFALSE;
	}
	else
		return kFALSE;
} 

void RunDedxCut2(TString inputfile, TString outputfile)
{
	TFile *input_rootfile = new TFile(inputfile);
	TTree* input_tree = (TTree*)input_rootfile->Get("events");

	ParticleTree output_tree(outputfile);

	Event *event = new Event();
	Particle *particle;
	input_tree->SetBranchAddress("event",&event);

	const Long64_t treeNentries = input_tree->GetEntries();
	Long_t particles_in = 0, particles_out = 0;
	Long64_t ev;
	UInt_t Npa;
	UInt_t part;

	//Float_t dedx_uppercut = 3.;

	float p;

	//cout << "Cut dE/dx > " << dedx_uppercut << " applied" << endl;

	for(ev=0; ev<treeNentries; ++ev)
	{
		if(!(ev%500))
			cout << "Event: " << ev << endl;

		input_tree->GetEntry(ev);
		Npa = event->GetNpa();
		output_tree.BeginEvent();

		particles_in += Npa;
		for(part=0; part<Npa; part++)
		{
			particle = event->GetParticle(part);
			p = TMath::Sqrt(TMath::Power(particle->GetPx(),2)+TMath::Power(particle->GetPy(),2)+TMath::Power(particle->GetPz(),2));

			if(is_electron(TMath::Log10(p),particle->GetdEdx()))
				continue;

			++particles_out;
			
			output_tree.AddParticle(particle->GetCharge(),
					particle->GetPx(), particle->GetPy(), particle->GetPz(),
					particle->GetdEdx(), particle->GetdEdxVtpc1(), particle->GetdEdxVtpc2(), particle->GetdEdxMtpc(),
					particle->GetNdEdx(), particle->GetNdEdxVtpc1(), particle->GetNdEdxVtpc2(), particle->GetNdEdxMtpc());
		}
		output_tree.EndEvent();
	}

	output_tree.Close();
	input_rootfile->Close();

	cout << "dEdx cut summary\n------------" << endl
		<< "Particles before cut: " << particles_in << endl
		<< "Particles after cut: " << particles_out << endl
		<< "Cutted particles: " << (particles_in-particles_out) << endl
		<< "Ratio: " << ((Double_t)particles_out/particles_in) << endl;
}

//pT cut. Default value is 1.5 GeV/c. Which means it takes everything between 0 and 1.5.
//For ptcut values lower than 0, the cut will take everything between abs(ptcut) and infinity.
void RunPtCut(TString inputfile, TString outputfile, const Float_t ptcut=1.5)
{
	if(ptcut > 0)
		cout << "Particles with pT > " << ptcut << " GeV/c will be rejected" << endl;
	else if(ptcut < 0)
		cout << "Particles with pT: 0 < " << ptcut << " GeV/c will be rejected" << endl;
	else
	{
		cout << "I don't see any sense in cutting all particles. Exiting" << endl;
		return;
	}

	TFile *input_rootfile = new TFile(inputfile);
	TTree* input_tree = (TTree*)input_rootfile->Get("events");

	ParticleTree output_tree(outputfile);

	Event *event = new Event();
	Particle *particle;
	input_tree->SetBranchAddress("event",&event);

	const Long64_t treeNentries = input_tree->GetEntries();
	Long_t particles_in = 0, particles_out = 0;
	Long64_t ev;
	UInt_t Npa;
	UInt_t part;

	Double_t pt;

	for(ev=0; ev<treeNentries; ++ev)
	{
		if(!(ev%500))
			cout << "Event: " << ev << endl;

		input_tree->GetEntry(ev);
		Npa = event->GetNpa();
		output_tree.BeginEvent();

		particles_in += Npa;
		for(part=0; part<Npa; part++)
		{
			particle = event->GetParticle(part);
			pt = TMath::Sqrt(TMath::Power(particle->GetPx(),2)+TMath::Power(particle->GetPy(),2));

			if(ptcut >= 0)
			{
				if(pt > ptcut)
					continue;
			}
			else if(ptcut < 0)
			{
				if(pt < TMath::Abs(ptcut))
					continue;
			}
			

			++particles_out;
			
			output_tree.AddParticle(*particle);
		}
		output_tree.EndEvent();
	}

	output_tree.Close();
	input_rootfile->Close();

	cout << "pT cut summary\n------------" << endl
		<< "Particles before cut: " << particles_in << endl
		<< "Particles after cut: " << particles_out << endl
		<< "Cutted particles: " << (particles_in-particles_out) << endl
		<< "Ratio: " << ((Double_t)particles_out/particles_in) << endl;
}
//pT cut

void RunYCut(TString inputfile, TString outputfile, const Double_t beam_momentum)
{
	TFile *input_rootfile = new TFile(inputfile);
	TTree* input_tree = (TTree*)input_rootfile->Get("events");

	ParticleTree output_tree(outputfile);

	Event *event = new Event();
	Particle *particle;
	input_tree->SetBranchAddress("event",&event);

	const Long64_t treeNentries = input_tree->GetEntries();
	Long_t particles_in = 0, particles_out = 0;
	Long64_t ev;
	UInt_t Npa;
	UInt_t part;

	const Double_t proton_mass = 0.938272013;
	const Double_t nucleon_mass = 0.9389186795;
	const Double_t beta = beam_momentum/(TMath::Sqrt(beam_momentum*beam_momentum+nucleon_mass*nucleon_mass)+nucleon_mass);
	const Double_t y_beam = 0.5*TMath::Log((1+beta)/(1-beta));

	cout << "y_target=0" << endl << "y_beam=" << y_beam << endl << "Everything outside 0.5 < y_prot < " << (y_beam - 0.5) << " will be rejected" << endl;

	Double_t E, p2, y_prot;

	for(ev=0; ev<treeNentries; ++ev)
	{
		if(!(ev%500))
			cout << "Event: " << ev << endl;

		input_tree->GetEntry(ev);
		Npa = event->GetNpa();
		output_tree.BeginEvent();

		particles_in += Npa;
		for(part=0; part<Npa; part++)
		{
			particle = event->GetParticle(part);
			p2 = TMath::Power(particle->GetPx(),2)+TMath::Power(particle->GetPy(),2)+TMath::Power(particle->GetPz(),2);
			E = TMath::Sqrt(proton_mass*proton_mass + p2);
			y_prot = 0.5*TMath::Log((E+particle->GetPz())/(E-particle->GetPz()));

			if((y_prot < 0.5) || (y_prot > y_beam-0.5))
				continue;

			++particles_out;
			
			output_tree.AddParticle(*particle);
		}
		output_tree.EndEvent();
	}

	output_tree.Close();
	input_rootfile->Close();

	cout << "y cut summary\n------------" << endl
		<< "Particles before cut: " << particles_in << endl
		<< "Particles after cut: " << particles_out << endl
		<< "Cutted particles: " << (particles_in-particles_out) << endl
		<< "Ratio: " << ((Double_t)particles_out/particles_in) << endl;

}

int main(int argc, char** argv)
{
	if(argc <= 1)
	{
		cout << "USAGE: cutter <inputfile> <outputfile> <cut_mode> [<energy> [<system>]]" << endl;
		return 0;
	}

	TString inputfile = argv[1];
	TString outputfile = argv[2];
	TString cut_mode = argv[3];
	TString energy = argv[4];
	TString system = argv[5];
	TString ttr_distance = argv[4];
	TString pt_cut_string, y_cut_string;

	cout << "cut mode:" << cut_mode << endl;
	if(!(cut_mode.CompareTo("PPM")))
	{
		if(argc != 6)
		{
			cout << "PPM cut requires additional arguments: 1.energy, 2.system" << endl;
			return 0;
		}
		RunPPMCut(inputfile, outputfile, system, energy);
	}
	else if(!(cut_mode.CompareTo("DEDX")))
	{
		if(argc != 6)
		{
			cout << "DEDX cut requires additional arguments: 1.energy, 2. system" << endl;
			return 0;
		}
		RunDedxCut(inputfile, outputfile, system, energy.Atoi());
	}
	else if(!(cut_mode.CompareTo("PT")))
	{
		if(argc == 5)
		{
			pt_cut_string = argv[4];
			RunPtCut(inputfile, outputfile, pt_cut_string.Atof());
		}
		else 
		{
			RunPtCut(inputfile, outputfile);	//default 1.5 GeV
		}
	}
	else if(!(cut_mode.CompareTo("Y")))
	{
		if(argc == 5)
		{
			y_cut_string = argv[4];
			cout << "Rapidity cut on spectators: y_target+0.5 < y < y_beam-0.5" << endl;
			RunYCut(inputfile, outputfile, y_cut_string.Atof());
		}
		else
			cout << "Y cut requires additional argument: beam momentum" << endl;

		return 0;
	}
}
