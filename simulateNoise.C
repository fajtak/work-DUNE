#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TMath.h"

const double gNoiseRatePerMM2 = 0.03;
const double gSiPMAreaMM2 = 36;
const double gNoiseRatePerSipm = gNoiseRatePerMM2*gSiPMAreaMM2; // 10mHz/mm2 -> 360 mHz/SiPM (6x6 mm2)
const int gNSiPMsPerTile = 160;
const double gSinglePERatio = 0.90;
// const double gDoublePERatio = 0.10;
const double gTriplePERatio = 0.01;
const double gPulseDuration = 10e-6;
const double gTimeMeasInS = 36000;
const double gSPEWidth = 0.1;

TRandom3 r_randNumGen3;

TH1D* h_noiseCharge = new TH1D("h_noiseCharge","Charge distribution of noise hits; Q [p.e.]; f [Hz/0.1 p.e.]",100,0,10);
TH1D* h_noiseChargePerPE = new TH1D("h_noiseChargePerPE","Charge distribution of noise hits; Q [p.e.]; f [Hz/1 p.e.]",10,0,10);
TH1D* h_noiseChargeNonMerged = new TH1D("h_noiseChargeNonMerged","Charge distribution of (non-merged) noise hits; Q [p.e.]; f [Hz/1 p.e.]",10,0,10);

int DrawResults()
{
	TCanvas* c_noiseCharge = new TCanvas("c_noiseCharge","NoiseCharge",800,600);
	h_noiseCharge->Scale(1/gTimeMeasInS);
	h_noiseCharge->Draw();

	TCanvas* c_noiseChargePerPE = new TCanvas("c_noiseChargePerPE","NoiseChargePerPE",800,600);
	h_noiseChargePerPE->Scale(1/gTimeMeasInS);
	h_noiseChargePerPE->Draw();

	TCanvas* c_noiseChargeNonMerged = new TCanvas("c_noiseChargeNonMerged","NoiseChargeNonMerged",800,600);
	h_noiseChargeNonMerged->Scale(1/gTimeMeasInS);
	h_noiseChargeNonMerged->Draw();
	return 0;
}

int SaveResults()
{
	TString outputFileName = Form("rootResults/noiseSim_nr%2.2f_spe%2.2f_tpe%2.2f_pd%2.1f_spew%2.1f.root",gNoiseRatePerMM2,gSinglePERatio,gTriplePERatio,gPulseDuration*1e6,gSPEWidth);
	TFile* outputFile = new TFile(outputFileName,"RECREATE");

	h_noiseCharge->Write();
	h_noiseChargePerPE->Write();
	h_noiseChargeNonMerged->Write();

	outputFile->Close();
	return 0;
}

struct Hit
{
	double time;
	double nPE;
};

bool OrderByPE(Hit a, Hit b)
{
	return (a.nPE > b.nPE);
}

bool OrderByTime(Hit a, Hit b)
{
	return (a.time < b.time);
}

double GenerateNoiseQ(int nPE)
{
	return r_randNumGen3.Gaus(nPE,TMath::Sqrt(nPE)*gSPEWidth);
}

int CreateHits(int nHits, vector<Hit> &noiseHits)
{
	double nPE = 0;
	double time = 0;
	for (int i = 0; i < nHits; ++i)
	{
		double rand = r_randNumGen3.Uniform(1);
		if (rand < gSinglePERatio)
			nPE = 1;
		else if (rand > 1-gTriplePERatio)
			nPE = 3;
		else
			nPE = 2;
		time = r_randNumGen3.Uniform(gTimeMeasInS);
		noiseHits.push_back(Hit{time,GenerateNoiseQ(nPE)});
		h_noiseChargeNonMerged->Fill(nPE);
	}

	return 0;
}

int CreateNoiseSpectrum(vector<Hit> &noiseHits)
{
	for (int i = 0; i < noiseHits.size(); ++i)
	{
		int nHits = 1;
		double nPE = noiseHits[i].nPE;
		while(i+nHits < noiseHits.size() && noiseHits[i+nHits].time - noiseHits[i+nHits-1].time < gPulseDuration) {
		    nPE += noiseHits[i+nHits].nPE;
		    nHits++;
		}
		h_noiseCharge->Fill(nPE);
		h_noiseChargePerPE->Fill(nPE);
		i += nHits-1;
	}

	return 0;
}

int simulateNoise()
{
	int nHits = gNSiPMsPerTile*gNoiseRatePerSipm*gTimeMeasInS;
	cout << nHits << endl;

	vector<Hit> v_noiseHits;

	CreateHits(nHits,v_noiseHits);
	std::sort(v_noiseHits.begin(),v_noiseHits.end(),OrderByTime);

	CreateNoiseSpectrum(v_noiseHits);

	DrawResults();
	SaveResults();
	return 0;
}