#include <iostream>
#include <algorithm> //fill_n
#include <string>
#include <cstdlib>

#include "TMath.h"
#include "TCutG.h"

#include "ParticleTree.h"
#include "Event.h"
#include "Particle.h"
#include "AccCut.h"
#include "dEdxCut.h"
#include "PPMCut.h"
#include "TTRCut.h"

using namespace std;

void RunAccCut(const int ener)
{
	cout << "Running acceptance mode" << endl;

	TFile *input_rootfile = new TFile("ParticleTree.root");
	TTree* input_tree = (TTree*)input_rootfile->Get("events");

	ParticleTree output_tree("ParticleTree_acc.root");

	Event *event = new Event();
	Particle *particle;
	input_tree->SetBranchAddress("event",&event);

	const Long64_t treeNentries = input_tree->GetEntries();
	Long64_t ev;
	UInt_t Npa;
	UInt_t part;

	AccCut acc_map("acceptance-map-medium.root",ener);

	float pt, E, p, y, angle;
	const float pion_mass = 0.13957018; //GeV/c^2
	const float nucleon_mass = 0.9389186795; //GeV/c^2
	const float beta = (ener/(ener+nucleon_mass));
	const float y_cms = 0.5*TMath::Log((1+beta)/(1-beta));

	for(ev=0; ev<treeNentries; ++ev)
	{
		if(!(ev%5000))
			cout << "Event: " << ev << endl;

			input_tree->GetEntry(ev);
			Npa = event->GetNpa();
			output_tree.BeginEvent();

			for(part=0; part<Npa; part++)
			{
				particle = event->GetParticle(part);
				pt = TMath::Sqrt(TMath::Power(particle->GetPx(),2)+TMath::Power(particle->GetPy(),2));
				p = TMath::Sqrt(TMath::Power(particle->GetPx(),2)+TMath::Power(particle->GetPy(),2)+TMath::Power(particle->GetPz(),2));
				E = TMath::Sqrt(pion_mass*pion_mass+p*p);
				y = 0.5*TMath::Log((E+particle->GetPz())/(E-particle->GetPz())) - y_cms;
				angle = TMath::ATan2(particle->GetPy(), particle->GetPx());

				//CIECIE NA AKCEPTACJE
				if(acc_map.acceptanceCut(particle->GetPx(),pt,particle->GetCharge(),y,angle))
					output_tree.AddParticle(particle->GetCharge(),
					particle->GetBx(), particle->GetBy(),
					particle->GetPx(), particle->GetPy(), particle->GetPz(),
					particle->GetdEdx(), particle->GetdEdxVtpc1(), particle->GetdEdxVtpc2(), particle->GetdEdxMtpc(),
					particle->GetNdEdx(), particle->GetNdEdxVtpc1(), particle->GetNdEdxVtpc2(), particle->GetNdEdxMtpc());
			}

			output_tree.EndEvent();
	}

	output_tree.Close();
	input_rootfile->Close();
}

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
					particle->GetBx(), particle->GetBy(),
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

void RunMultSplit(TString inputfile, TString outputfile, const TString mult_string)
{
	cout << "Running multiplicity splitter mode" << endl;
	const UInt_t multiplicity = atoi(mult_string);

	TFile *input_rootfile = new TFile(inputfile);
	TTree* input_tree = (TTree*)input_rootfile->Get("events");

	ParticleTree output_tree_all(outputfile + "_all.root");
	ParticleTree output_tree_pos(outputfile + "_pos.root");
	ParticleTree output_tree_neg(outputfile + "_neg.root");

	Event *event = new Event();
	Particle *particle;
	input_tree->SetBranchAddress("event",&event);

	const Long64_t treeNentries = input_tree->GetEntries();
	Long64_t ev;
	UInt_t Npa;
	UInt_t Npos;
	UInt_t Nneg;
	UInt_t part;
	UInt_t all_count = 0;
	UInt_t pos_count = 0;
	UInt_t neg_count = 0;
	
	//cout << "Event\t\tAll\tPos\tNeg" << endl;

	//DLA WSZYSTKICH
	for(ev=0; ev<treeNentries; ++ev)
	{
		if(!(ev%5000))
			cout << "Event: " << ev << endl;
		input_tree->GetEntry(ev);
		Npa = event->GetNpa();
		Npos = event->GetNpos();
		Nneg = event->GetNneg();

		if(Npa==multiplicity)
		{
			output_tree_all.BeginEvent();
			for(part=0; part<Npa; part++)
			{
				particle = event->GetParticle(part);
				//output_tree_all.AddParticle(particle->GetPid(), particle->GetCharge(), particle->GetBx(), particle->GetBy(), particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetDedx());
				output_tree_all.AddParticle(particle->GetCharge(),
						particle->GetBx(), particle->GetBy(),
						particle->GetPx(), particle->GetPy(), particle->GetPz(),
						particle->GetdEdx(), particle->GetdEdxVtpc1(), particle->GetdEdxVtpc2(), particle->GetdEdxMtpc(),
						particle->GetNdEdx(), particle->GetNdEdxVtpc1(), particle->GetNdEdxVtpc2(), particle->GetNdEdxMtpc());
			}
			++all_count;
			output_tree_all.EndEvent();
		}

		if(Npos==multiplicity)
		{
			output_tree_pos.BeginEvent();
			for(part=0; part<Npa; part++)
			{
				particle = event->GetParticle(part);
				if((particle->GetCharge())>0)
					//output_tree_pos.AddParticle(particle->GetPid(), particle->GetCharge(), particle->GetBx(), particle->GetBy(), particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetDedx());
				output_tree_pos.AddParticle(particle->GetCharge(),
						particle->GetBx(), particle->GetBy(),
						particle->GetPx(), particle->GetPy(), particle->GetPz(),
						particle->GetdEdx(), particle->GetdEdxVtpc1(), particle->GetdEdxVtpc2(), particle->GetdEdxMtpc(),
						particle->GetNdEdx(), particle->GetNdEdxVtpc1(), particle->GetNdEdxVtpc2(), particle->GetNdEdxMtpc());
			}
			++pos_count;
			output_tree_pos.EndEvent();
		}

		if(Nneg==multiplicity)
		{
			output_tree_neg.BeginEvent();
			for(part=0; part<Npa; part++)
			{
				particle = event->GetParticle(part);
				if((particle->GetCharge())<0)
					//output_tree_neg.AddParticle(particle->GetPid(), particle->GetCharge(), particle->GetBx(), particle->GetBy(), particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetDedx());
				output_tree_neg.AddParticle(particle->GetCharge(),
						particle->GetBx(), particle->GetBy(),
						particle->GetPx(), particle->GetPy(), particle->GetPz(),
						particle->GetdEdx(), particle->GetdEdxVtpc1(), particle->GetdEdxVtpc2(), particle->GetdEdxMtpc(),
						particle->GetNdEdx(), particle->GetNdEdxVtpc1(), particle->GetNdEdxVtpc2(), particle->GetNdEdxMtpc());
			}
			++neg_count;
			output_tree_neg.EndEvent();
		}
		//cout << "\r" << ev << "\t\t" << all_count << "\t" << pos_count << "\t" << neg_count;
	}
	
	cout << "\nEvents\t\tAll\tPos\tNeg" << endl << ev << "\t\t" << all_count << "\t" << pos_count << "\t" << neg_count << endl;

	output_tree_all.Close();
	output_tree_pos.Close();
	output_tree_neg.Close();
	input_rootfile->Close();
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
			std::cout << "WTF? Particle with no dE/dx information!" << std::endl;
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
			particle = event->GetParticle(part);
			p = TMath::Sqrt(TMath::Power(particle->GetPx(),2)+TMath::Power(particle->GetPy(),2)+TMath::Power(particle->GetPz(),2));
			local_dedx = choose_dedx(particle, system);
			if(cutg->IsInside(p,local_dedx))
				continue;

			if(local_dedx > dedx_uppercut)
				continue;		//dodatkowy cut 3.02.2013

			output_tree.AddParticle(particle->GetCharge(),
					particle->GetBx(), particle->GetBy(),
					particle->GetPx(), particle->GetPy(), particle->GetPz(),
					particle->GetdEdx(), particle->GetdEdxVtpc1(), particle->GetdEdxVtpc2(), particle->GetdEdxMtpc(),
					particle->GetNdEdx(), particle->GetNdEdxVtpc1(), particle->GetNdEdxVtpc2(), particle->GetNdEdxMtpc());
		}
		output_tree.EndEvent();
	}

	output_tree.Close();
	input_rootfile->Close();
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
					particle->GetBx(), particle->GetBy(),
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

void RunElasticCut(TString inputfile, TString outputfile, Int_t energy)
{
	TFile *input_rootfile = new TFile(inputfile);
	TTree* input_tree = (TTree*)input_rootfile->Get("events");

	ParticleTree output_tree(outputfile);

	Event *event = new Event();
	Particle *particle;
	input_tree->SetBranchAddress("event",&event);

	const Long64_t treeNentries = input_tree->GetEntries();
	Long64_t ev;
	Long_t particles_in = 0, particles_out = 0, events_out = 0;
	UInt_t part, Npa;

	Float_t p;
	Bool_t event_cancelled;

	for(ev=0; ev<treeNentries; ++ev)
	{
		if(!(ev%500))
			cout << "Event: " << ev << endl;

		input_tree->GetEntry(ev);

		Npa = event->GetNpa();

		event_cancelled = false;
		for(part=0; part<Npa; ++part)
		{
			particle = event->GetParticle(part);
			p = TMath::Sqrt(TMath::Power(particle->GetPx(),2)+TMath::Power(particle->GetPy(),2)+TMath::Power(particle->GetPz(),2));
			if((particle->isPositive()) && (p > energy - 3))
			{
				event_cancelled = true;
				break;
			}
		}

		particles_in+=Npa;

		if(event_cancelled)
			continue;

		events_out++;
		
		output_tree.BeginEvent();

		for(part=0; part<Npa; part++)
		{
			particle = event->GetParticle(part);
			
			output_tree.AddParticle(particle->GetCharge(),
					particle->GetBx(), particle->GetBy(),
					particle->GetPx(), particle->GetPy(), particle->GetPz(),
					particle->GetdEdx(), particle->GetdEdxVtpc1(), particle->GetdEdxVtpc2(), particle->GetdEdxMtpc(),
					particle->GetNdEdx(), particle->GetNdEdxVtpc1(), particle->GetNdEdxVtpc2(), particle->GetNdEdxMtpc());
		}
		output_tree.EndEvent();
		particles_out+=Npa;
	}

	output_tree.Close();
	input_rootfile->Close();

	cout << "Elastic cut summary\n------------" << endl
		<< "Events before cut: " << treeNentries  << endl
		<< "Events after cut: " << events_out << endl
		<< "Cutted events: " << (treeNentries-events_out) << endl
		<< "Ratio: " << ((Double_t)treeNentries/events_out) << "\n------------" << endl
		<< "Particles before cut: " << particles_in << endl
		<< "Particles after cut: " << particles_out << endl
		<< "Cutted particles: " << (particles_in-particles_out) << endl
		<< "Ratio: " << ((Double_t)particles_out/particles_in) << endl;
}

int RunTTRCut(TString inputfile, TString outputfile, double distance=1.6)
{
	TTRCut ttrcut;
	TFile *input_rootfile = new TFile(inputfile);
	TTree* input_tree = (TTree*)input_rootfile->Get("events");

	ParticleTree output_tree(outputfile);

	Event *event = new Event();
	Particle *particleA, *particleB;
	input_tree->SetBranchAddress("event",&event);

	const Long64_t treeNentries = input_tree->GetEntries();
	Long64_t ev;
	Long_t particles_in = 0, particles_out = 0;
	UInt_t partA, partB, Npa;

	Float_t distance_av;
	Bool_t track_ok;
	cout << "distance in function: " << distance << endl;

	for(ev=0; ev<treeNentries; ++ev)
	{
		if(!(ev%500))
			cout << "Event: " << ev << " " << event->GetNpa() << " particles" << endl;

		input_tree->GetEntry(ev);
		Npa = event->GetNpa();

		Bool_t ttr_flags[Npa];
		fill_n(ttr_flags,Npa,true);
		
		output_tree.BeginEvent();

		for(partA=0; partA<Npa; partA++)
		{
			++particles_in;
			if(!ttr_flags[partA])
			{
				//cout << "Ev: " << ev << " part: " << partA << " skipped" << endl;
				continue;
			}

			track_ok = true;
			particleA = event->GetParticle(partA);

			for(partB=partA+1; partB<Npa; ++partB)
			{
				if(!ttr_flags[partB])
				{
					//cout << "Ev: " << ev << " part: " << partB << " skipped" << endl;
					continue;
				}

				particleB = event->GetParticle(partB);

				///////////////////
				distance_av = ttrcut.calcAvDistance(particleA,particleB);
				//////////////////
				
				if(distance_av < distance)
				{
					track_ok = false;
					ttr_flags[partB] = false;
					//cout << "Ev: " << ev << " particles " << partA << " and " << partB << " will be cut" << endl;
					//cerr << "Average: " << distance_av << endl;
					break;
				}
			}

			if(track_ok)
			{
				//trackA accepted
				++particles_out;
				output_tree.AddParticle(*particleA);
			}
		}
		
		//cerr << "Event: " << ev << " | Cutted " << (Npa-output_tree.Check()) << " particles" << endl;
		output_tree.EndEvent();
	}

	input_rootfile->Close();
	output_tree.Close();

	cout << "TTR cut summary\n------------" << endl
		<< "Events: " << treeNentries << endl
		<< "Particles before cut: " << particles_in << endl
		<< "Particles after cut: " << particles_out << endl
		<< "Cutted particles: " << particles_in-particles_out << endl
		<< "Ratio: " << ((Double_t)particles_out/particles_in) << endl;

	return particles_out;
}

int main(int argc, char** argv)
{
	if(argc <= 1)
	{
		cout << "USAGE: cutter <inputfile> <outputfile> <cut_mode> [<energy>/<multsplit> [<system>]]" << endl;
		return 0;
	}

	TString inputfile = argv[1];
	TString outputfile = argv[2];
	TString cut_mode = argv[3];
	TString energy = argv[4];
	TString system = argv[5];
	TString ttr_distance = argv[4];
	TString mult_string;

	cout << "cut mode:" << cut_mode << endl;
	if(!(cut_mode.CompareTo("NA61ACC")))
	{
		cout << "WARNING: Energy fixed to 158 GeV!" << endl;
		RunAccCut(158);
	}
	else if(!(cut_mode.CompareTo("PPM")))
	{
		if(argc != 6)
		{
			cout << "PPM cut requires additional arguments: 1.energy, 2.system" << endl;
			return 0;
		}
		RunPPMCut(inputfile, outputfile, system, energy);
	}
	else if(!(cut_mode.CompareTo("MULTSPLIT")))
	{
		mult_string = argv[4];
		RunMultSplit(inputfile, outputfile, mult_string);
	}
	else if(!(cut_mode.CompareTo("DEDX")))
	{
		if(argc != 6)
		{
			cout << "DEDX cut requires additional arguments: 1.energy, 2. system" << endl;
			return 0;
		}
		//RunDedxCut(inputfile, outputfile, system, energy.Atoi());
		RunDedxCut2(inputfile, outputfile);
	}
	else if(!(cut_mode.CompareTo("ELASTIC")))
	{
		if(argc != 5)
		{
			cout << "ELASTIC cut requires additional argument: energy" << endl;
			return 0;
		}
		cout << "Elastic cut mode" << endl;
		RunElasticCut(inputfile,outputfile, energy.Atoi());
	}
	else if(!(cut_mode.CompareTo("TTR")))
	{
		cout << "Two-track resolution mode" << endl;
		if(argc == 5)
		{
			cout << "\tRunning with cut distance " << ttr_distance << " cm" << endl;
			RunTTRCut(inputfile,outputfile,ttr_distance.Atoi());
		}
		else
		{
			cout << "\tRunning with default cut distance" << endl;
			RunTTRCut(inputfile,outputfile);
		}
	}
}