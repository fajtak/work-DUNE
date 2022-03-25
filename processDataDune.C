#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TRandom2.h"
#include "TRandom1.h"
#include "TRandom.h"
#include "TFile.h"
#include "THStack.h"
#include "TVector3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TLatex.h"

const double gArapucaEfficiency = 0.035;
const bool gTilesHorizontally = false;
const bool gTilesVertically = false;
const bool gTilesOpposite = false;
const bool gTilesSmall = false;
const bool gReadAnodeTiles = false;
const bool gReadCathodeTiles = false;
const bool gReadCryoTiles = true;
const bool gMaskCentralTiles = true;
const bool gMaskMiddleTiles = false;
const bool gMaskNegativeTiles = false;

const bool gUseCathodeTileTransparency = true;
TRandom3 r_randNumGen3;

const int gPEThreshold = 5;
const int gMulti = 4;

const int gMultiThreshold2pe = 13;
const int gMultiThreshold3pe = 5;
const int gMultiThreshold4pe = 5;

struct Hit
{
	unsigned short int cellID;
	unsigned short int nPE;
};

struct  Event
{
	int eventID;
	float energy;
	float x;
	float y;
	float z;
	int nHits;
	vector<Hit> v_hits;
	unsigned int nPEWalls;
	unsigned int nPEAnode;
	unsigned int nPECathode;
};

const int nCubesX = 27;
const int nCubesY = 26;
const int nCubesZ = 30;
unsigned short int DetectedEvents[nCubesX][nCubesY][nCubesZ];
unsigned short int AllEvents[nCubesX][nCubesY][nCubesZ];

void PrintEvent(Event &event)
{
	cout << "EventID: " << event.eventID << " ( " << event.x << " , " << event.y << " , " << event.z << " ) nHits: " << event.nHits << " nPEWalls: " << event.nPEWalls << " nPEAnode: " << event.nPEAnode << " nPECathode: " << event.nPECathode << endl;
	for (unsigned int i = 0; i < event.v_hits.size(); ++i)
	{
		cout << "\t" << i << " " << event.v_hits[i].cellID << " " << event.v_hits[i].nPE << endl;
	}
}

bool OrderByPE (Hit a, Hit b)
{
	return (a.nPE > b.nPE);
}

TH1D* h_nPETotal = new TH1D("h_nPETotal","Npe_{Total}; N_{p.e.} [#]; NoE [#]",10000,0,10000);
TH1D* h_nPECryo = new TH1D("h_nPECryo","Npe_{Cryo}; N_{p.e.} [#]; NoE [#]",10000,0,10000);
TH1D* h_nPECathode = new TH1D("h_nPECathode","Npe_{Cathode}; N_{p.e.} [#]; NoE [#]",10000,0,10000);
TH1D* h_nPECryoPlus = new TH1D("h_nPECryoPlus","Npe_{CryoMinus}; N_{p.e.} [#]; NoE [#]",10000,0,10000);
TH1D* h_nPECryoMinus = new TH1D("h_nPECryoMinus","Npe_{CryoPlus}; N_{p.e.} [#]; NoE [#]",10000,0,10000);
TH1D* h_nPEAnode = new TH1D("h_nPEAnode","Npe_{Anode}; N_{p.e.} [#]; NoE [#]",10000,0,10000);
TH1D* h_PECryoRatio = new TH1D("h_PECryoRatio","nPECryo/nPETotal;Ratio [%]; NoE [#]",110,0,110);
TH1D* h_peRatio = new TH1D("h_peRatio",";Ratio [%]; NoE [#]",1100,0,110);

TH1D* h_nPEMaxMegacell = new TH1D("h_nPEMaxMegacell","Number of detected p.e. in max Megacell; N_{p.e.} [#]; NoE [#]",1000,0,1000);
TH2D* h_nPEMaxMegacellX	= new TH2D("h_nPEMaxMegacellX","Number of detected p.e. in max Megacell vs. X position; X [cm]; N_{p.e.} [#]",80,-100,700,1000,0,1000);
TH2D* h_nPEMaxMegacellY	= new TH2D("h_nPEMaxMegacellY","Number of detected p.e. in max Megacell vs. Y position; Y [cm]; N_{p.e.} [#]",70,0,700,1000,0,1000);
TH2D* h_nPEMaxMegacellZ	= new TH2D("h_nPEMaxMegacellZ","Number of detected p.e. in max Megacell vs. Z position; Z [cm]; N_{p.e.} [#]",300,0,3000,1000,0,1000);

TH1D* h_multiplicity = new TH1D("h_multiplicity","Hit multiplicity; N_{hits} [#]; NoE [#]",300,0,300);
TH2D* h_multiplicityVsTotalNPE = new TH2D("h_multiplicityVsTotalNPE",";Total Q [p.e.]; N_{hits} [#]; NoE [#]",1000,0,1000,300,0,300);
TH1D* h_multiplicityFiltered = new TH1D("h_multiplicityFiltered","Hit multiplicity Filtered; N_{hits} [#]; NoE [#]",300,0,300);
TH1D* h_multiplicityAbove5pe = new TH1D("h_multiplicityAbove5pe","Hit multiplicity above 5 p.e.; N_{hits} [#]; NoE [#]",300,0,300);
TH1D* h_multiplicityAbove4pe = new TH1D("h_multiplicityAbove4pe","Hit multiplicity above 4 p.e.; N_{hits} [#]; NoE [#]",300,0,300);
TH1D* h_multiplicityAbove3pe = new TH1D("h_multiplicityAbove3pe","Hit multiplicity above 3 p.e.; N_{hits} [#]; NoE [#]",300,0,300);
TH1D* h_multiplicityAbove2pe = new TH1D("h_multiplicityAbove2pe","Hit multiplicity above 2 p.e.; N_{hits} [#]; NoE [#]",300,0,300);
TH2D* h_multiplicityX	= new TH2D("h_multiplicityX","Hit multiplicity vs. X position; X [cm]; N_{hits} [#]",80,-100,700,1000,0,1000);
TH2D* h_multiplicityXAboveNpe = new TH2D("h_multiplicityXAboveNpe",Form("Hit multiplicity above %d p.e. vs. X position; X [cm]; N_{hits} [#]",gPEThreshold),80,-100,700,1000,0,1000);
TH2D* h_multiplicityYAboveNpe = new TH2D("h_multiplicityYAboveNpe",Form("Hit multiplicity above %d p.e. vs. Y position; Y [cm]; N_{hits} [#]",gPEThreshold),70,0,700,1000,0,1000);

TH1D* h_energy = new TH1D("h_energy","Energy distribution;E [MeV]; NoE [#]",200,0,20);
TH1D* h_x = new TH1D("h_x","X distribution;X [cm]; NoE [#]",80,-0,800);
TH1D* h_xDet = new TH1D("h_xDet","X distribution;X [cm]; NoE [#]",80,0,800);
TH1D* h_y = new TH1D("h_y","Y distribution;Y [cm]; NoE [#]",70,0,700);
TH1D* h_z = new TH1D("h_z","Z distribution;Z [cm]; NoE [#]",300,0,3000);
TH1D* h_nHits = new TH1D("h_nHits","NHits distribution;N_{hits} [#]; NoE [#]",2000,0,2000);

TH1D* h_nHitsDiff = new TH1D("h_nHitsDiff","NHitsNonMerged - NHitsMerged; #Delta N_{hits} [#]; NoE [#]",100,0,100);

TH1D* h_cellID = new TH1D("h_cellID","Cell ID; CellID [#]; NoE [#]",50000,0,50000);

TH1D* h_detEffi5pe = new TH1D("h_detEffi5pe","M_{T} = 1, N_{pe} = 5;X [cm]; Detection Efficiency [%]",80,-100,700);
TH1D* h_detEffi10pe = new TH1D("h_detEffi10pe","10 p.e.;X [cm]; Detection Efficiency [%]",80,-100,700);
TH1D* h_detEffi15pe = new TH1D("h_detEffi15pe","15 p.e.;X [cm]; Detection Efficiency [%]",80,-100,700);

TH1D* h_detEffi5peY = new TH1D("h_detEffi5peY","5 p.e.;Y [cm]; Detection Efficiency [%]",70,0,700);
TH1D* h_detEffi5peZ = new TH1D("h_detEffi5peZ","5 p.e.;Z [cm]; Detection Efficiency [%]",300,0,3000);

TH1D* h_detEffiNpeMulti1 = new TH1D("h_detEffiNpeMulti1",Form("M_{T} = 1, N_{pe} = %d;X [cm]; Detection Efficiency [%%]",gPEThreshold),80,-100,700);
TH1D* h_detEffiNpeMulti2 = new TH1D("h_detEffiNpeMulti2",Form("M_{T} = 2, N_{pe} = %d;X [cm]; Detection Efficiency [%%]",gPEThreshold),80,-100,700);
TH1D* h_detEffiNpeMulti3 = new TH1D("h_detEffiNpeMulti3",Form("M_{T} = 3, N_{pe} = %d;X [cm]; Detection Efficiency [%%]",gPEThreshold),80,-100,700);
TH1D* h_detEffiNpeMulti4 = new TH1D("h_detEffiNpeMulti4",Form("M_{T} = 4, N_{pe} = %d;X [cm]; Detection Efficiency [%%]",gPEThreshold),80,-100,700);
TH1D* h_detEffiNpeMulti5 = new TH1D("h_detEffiNpeMulti5",Form("M_{T} = 5, N_{pe} = %d;X [cm]; Detection Efficiency [%%]",gPEThreshold),80,-100,700);

TH1D* h_detEffiNpeMulti1Y = new TH1D("h_detEffiNpeMulti1Y",Form("M_{T} = 1, N_{pe} = %d;Y [cm]; Detection Efficiency [%%]",gPEThreshold),70,0,700);
TH1D* h_detEffiNpeMulti2Y = new TH1D("h_detEffiNpeMulti2Y",Form("M_{T} = 2, N_{pe} = %d;Y [cm]; Detection Efficiency [%%]",gPEThreshold),70,0,700);
TH1D* h_detEffiNpeMulti3Y = new TH1D("h_detEffiNpeMulti3Y",Form("M_{T} = 3, N_{pe} = %d;Y [cm]; Detection Efficiency [%%]",gPEThreshold),70,0,700);
TH1D* h_detEffiNpeMulti4Y = new TH1D("h_detEffiNpeMulti4Y",Form("M_{T} = 4, N_{pe} = %d;Y [cm]; Detection Efficiency [%%]",gPEThreshold),70,0,700);
TH1D* h_detEffiNpeMulti5Y = new TH1D("h_detEffiNpeMulti5Y",Form("M_{T} = 5, N_{pe} = %d;Y [cm]; Detection Efficiency [%%]",gPEThreshold),70,0,700);

TH2D* h_anodeEvents = new TH2D("h_anodeEvents","XY distribution of events with nPEAnode > 3; X [cm]; Y [cm]; NoE [#]",80,-100,700,70,0,700);

TH1D* h_detEffiFinal = new TH1D("h_detEffiFinal","M_{T} = 2, N_{pe} = 5, Multi > 30;X [cm]; Detection Efficiency [%]",80,0,800);
TH1D* h_detEffiFinalAll = new TH1D("h_detEffiFinalAll","M_{T} = 2, N_{pe} = 5, Multi > 30;X [cm]; Detection Efficiency [%]",80,0,800);
TH1D* h_detEffiFinalDet = new TH1D("h_detEffiFinalDet","M_{T} = 2, N_{pe} = 5, Multi > 30;X [cm]; Detection Efficiency [%]",80,0,800);

TH1D* h_mismatchX = new TH1D("h_mismatchX","Mismatch in X coordinate; #Delta X [cm]; NoE [#]",200,-1000,1000);
TH1D* h_mismatchY = new TH1D("h_mismatchY","Mismatch in Y coordinate; #Delta Y [cm]; NoE [#]",200,-1000,1000);
TH1D* h_mismatchZ = new TH1D("h_mismatchZ","Mismatch in Z coordinate; #Delta Z [cm]; NoE [#]",200,-1000,1000);

TH2D* h_mismatchXvsX = new TH2D("h_mismatchXvsX", "Mismatch in X coordinate vs x coordinate; X [cm]; #DeltaX [cm]; NoE [#]",80,-100,700,200,-1000,1000);
TH2D* h_mismatchYvsY = new TH2D("h_mismatchYvsY", "Mismatch in Y coordinate vs Y coordinate; Y [cm]; #DeltaY [cm]; NoE [#]",70,-0,700,200,-1000,1000);
TH2D* h_mismatchZvsZ = new TH2D("h_mismatchZvsZ", "Mismatch in Z coordinate vs Z coordinate; Z [cm]; #DeltaZ [cm]; NoE [#]",300,0,3000,200,-1000,1000);

TH1D* h_posXNotReco = new TH1D("h_posXNotReco","X position of not reconstructed events; X [cm]; Ratio of Untagged Events [%]",80,0,800);
TH1D* h_posYNotReco = new TH1D("h_posYNotReco","Y position of not reconstructed events; Y [cm]; Ratio of Untagged Events [%]",70,0,700);
TH1D* h_posZNotReco = new TH1D("h_posZNotReco","Z position of not reconstructed events; Z [cm]; Ratio of Untagged Events [%]",300,0,3000);

TH3D* h_averageLY = new TH3D("h_averageLY","Average LY per MeV; x [cm]; y [cm]; LY [p.e./MeV]",54,-675,675,26,0,650,1000,-0.5,999.5);

TGraph2D* h_posNotReco = new TGraph2D();

void DrawResults()
{
	THStack* s_nPE = new THStack("s_nPE","Number of detected p.e.; N_{p.e.} [#]; NoE [#]");
	s_nPE->Add(h_nPETotal);
	s_nPE->Add(h_nPECryo);
	h_nPECryo->SetLineColor(kGreen);
	s_nPE->Add(h_nPEAnode);
	h_nPEAnode->SetLineColor(kBlue);

	TCanvas* c_nPE = new TCanvas("c_nPE","NPE",800,600);
	gPad->SetGrid();
	s_nPE->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_nPE->Write();

	TCanvas* c_PECryoRatio = new TCanvas("c_PECryoRatio","PECryoRatio",800,600);
	h_PECryoRatio->Draw();

	TCanvas* c_basic = new TCanvas("c_basic","Basic",800,600);
	c_basic->Divide(2,2);
	c_basic->cd(1);
	h_energy->Draw();
	c_basic->cd(2);
	h_x->Draw();
	h_xDet->Draw("same");
	h_xDet->SetLineColor(kBlue);
	c_basic->cd(3);
	h_y->Draw();
	c_basic->cd(4);
	h_z->Draw();

	TCanvas* c_peRatio = new TCanvas("c_peRatio","PERatio",800,600);
	h_peRatio->Draw();

	TCanvas* c_nHitsDiff = new TCanvas("c_nHitsDiff","NHitsDifference",800,600);
	h_nHitsDiff->Draw();

	TCanvas* c_nHits = new TCanvas("c_nHits","NHits",800,600);
	h_nHits->Draw();

	TCanvas* c_nPECryoPlusMinus = new TCanvas("c_nPECryoPlusMinus","PECryoPlusMinus",800,600);
	h_nPECryoPlus->Draw();
	h_nPECryoMinus->Draw("same");
	h_nPECryoMinus->SetLineColor(kRed);

	TCanvas* c_cellID = new TCanvas("c_cellID","CellID",800,600);
	h_cellID->Draw();

	TCanvas* c_nPEMaxMegacell = new TCanvas("c_nPEMaxMegacell","NPEMegaCell",800,600);
	h_nPEMaxMegacell->Draw();

	TCanvas* c_nPEMaxMegacellX = new TCanvas("c_nPEMaxMegacellX","NPEMegaCellX",800,600);
	h_nPEMaxMegacellX->Draw("colz");

	TCanvas* c_nPEMaxMegacellY = new TCanvas("c_nPEMaxMegacellY","NPEMegaCellY",800,600);
	h_nPEMaxMegacellY->Draw("colz");

	TCanvas* c_nPEMaxMegacellZ = new TCanvas("c_nPEMaxMegacellZ","NPEMegaCellZ",800,600);
	h_nPEMaxMegacellZ->Draw("colz");

	TCanvas* c_multiplicity = new TCanvas("c_multiplicity","Multiplicity",800,600);
	h_multiplicity->Draw();
	h_multiplicityAbove3pe->Draw("same");
	h_multiplicityAbove3pe->SetLineColor(kRed);
	// h_multiplicityAbove5pe->Draw("same");
	// h_multiplicityAbove5pe->SetLineColor(kRed);
	// h_multiplicityFiltered->Draw("same");
	// h_multiplicityFiltered->SetLineColor(kBlue);
	TCanvas* c_multiplicityVsTotalNPE = new TCanvas("c_multiplicityVsTotalNPE","TotalNPEvsMultiplicity",800,600);
	h_multiplicityVsTotalNPE->Draw("colz");

	THStack* s_multi = new THStack("s_multi","Detection Multiplicity; N_{tiles} [#]; Number of Entries [#] ");
	s_multi->Add(h_multiplicityAbove2pe,"HIST");
	h_multiplicityAbove2pe->SetLineColor(kRed);
	h_multiplicityAbove2pe->SetLineWidth(4);
	s_multi->Add(h_multiplicityAbove3pe,"HIST");
	h_multiplicityAbove3pe->SetLineColor(kBlue);
	h_multiplicityAbove3pe->SetLineWidth(4);
	s_multi->Add(h_multiplicityAbove4pe,"HIST");
	h_multiplicityAbove4pe->SetLineColor(kGreen);
	h_multiplicityAbove4pe->SetLineWidth(4);
	s_multi->Add(h_multiplicityAbove5pe,"HIST");
	h_multiplicityAbove5pe->SetLineWidth(4);


	TCanvas* c_multi = new TCanvas("c_multi","MultiplicityComparison",800,600);
	gPad->SetGrid();
	s_multi->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_multi->Write();

	THStack* s_multiCumulative = new THStack("s_multiCumulative","Detection Multiplicity; N_{tiles} [#]; Tagging efficiency [1] ");
	TH1* tempHist = h_multiplicityAbove2pe->GetCumulative(false);
	tempHist->Scale(1/h_multiplicityAbove2pe->Integral());
	s_multiCumulative->Add(tempHist,"HIST");
	// h_multiplicityAbove2pe->SetLineColor(kRed);
	tempHist = h_multiplicityAbove3pe->GetCumulative(false);
	tempHist->Scale(1/h_multiplicityAbove3pe->Integral());
	s_multiCumulative->Add(tempHist,"HIST");
	// h_multiplicityAbove3pe->SetLineColor(kCyan);
	tempHist = h_multiplicityAbove4pe->GetCumulative(false);
	tempHist->Scale(1/h_multiplicityAbove4pe->Integral());
	s_multiCumulative->Add(tempHist,"HIST");
	// h_multiplicityAbove4pe->SetLineColor(kGreen);
	tempHist = h_multiplicityAbove5pe->GetCumulative(false);
	tempHist->Scale(1/h_multiplicityAbove5pe->Integral());
	s_multiCumulative->Add(tempHist,"HIST");


	TCanvas* c_multiCumulative = new TCanvas("c_multiCumulative","MultiplicityComparisonCumulative",800,600);
	gPad->SetGrid();
	s_multiCumulative->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_multiCumulative->Write();

	TCanvas* c_multiplicityX = new TCanvas("c_multiplicityX","MultiplicityX",800,600);
	h_multiplicityX->Draw("colz");

	TCanvas* c_multiplicityXAboveNpe = new TCanvas("c_multiplicityXAboveNpe",Form("MultiplicityXAbove%dpe",gPEThreshold),800,600);
	h_multiplicityXAboveNpe->Draw("colz");

	TCanvas* c_detEffi5pe = new TCanvas("c_detEffi5pe","DetectionEfficiency5pe",800,600);
	h_detEffi5pe->Draw("HIST");
	h_detEffi10pe->Draw("HISTSAME");
	h_detEffi10pe->SetLineColor(kRed);
	h_detEffi15pe->Draw("HISTSAME");
	h_detEffi15pe->SetLineColor(kBlue);

	THStack* s_detEffiX5peMulti = new THStack("s_detEffiX5peMulti","Detection efficiency vs. X position; X [cm]; Detection efficiency [%]");
	s_detEffiX5peMulti->Add(h_detEffiNpeMulti1,"HIST");
	s_detEffiX5peMulti->Add(h_detEffiNpeMulti2,"HIST");
	h_detEffiNpeMulti2->SetLineColor(kCyan);
	s_detEffiX5peMulti->Add(h_detEffiNpeMulti3,"HIST");
	h_detEffiNpeMulti3->SetLineColor(kRed);
	s_detEffiX5peMulti->Add(h_detEffiNpeMulti4,"HIST");
	h_detEffiNpeMulti4->SetLineColor(kBlue);
	s_detEffiX5peMulti->Add(h_detEffiNpeMulti5,"HIST");
	h_detEffiNpeMulti5->SetLineColor(kGreen);

	TCanvas* c_detEffiNpeMulti = new TCanvas("c_detEffiNpeMulti",Form("DetectionEfficiency%dpeMulti",gPEThreshold),800,600);
	gPad->SetGrid();
	s_detEffiX5peMulti->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_detEffiX5peMulti->Write();

	TCanvas* c_detEffiNpeMultiY = new TCanvas("c_detEffiNpeMultiY",Form("DetectionEfficiency%dpeMultiY",gPEThreshold),800,600);
	h_detEffiNpeMulti1Y->Draw("HIST");
	h_detEffiNpeMulti2Y->Draw("HISTSAME");
	h_detEffiNpeMulti2Y->SetLineColor(kCyan);
	h_detEffiNpeMulti3Y->Draw("HISTSAME");
	h_detEffiNpeMulti3Y->SetLineColor(kRed);
	h_detEffiNpeMulti4Y->Draw("HISTSAME");
	h_detEffiNpeMulti4Y->SetLineColor(kBlue);
	h_detEffiNpeMulti5Y->Draw("HISTSAME");
	h_detEffiNpeMulti5Y->SetLineColor(kGreen);

	TCanvas* c_detEffi5peY = new TCanvas("c_detEffi5peY","DetectionEfficiency5peY",800,600);
	h_detEffi5peY->Draw("HIST");

	TCanvas* c_detEffi5peZ = new TCanvas("c_detEffi5peZ","DetectionEfficiency5peZ",800,600);
	h_detEffi5peZ->Draw("HIST");

	TCanvas* c_anodeEvents = new TCanvas("c_anodeEvents","AnodeEvents",800,600);
	h_anodeEvents->Draw("COLZ");

	TCanvas* c_detEffiFinal = new TCanvas("c_detEffiFinal","DetectionEfficiencyFinal",800,600);
	h_detEffiFinal->Divide(h_detEffiFinalDet,h_detEffiFinalAll);
	h_detEffiFinal->Draw("");
	// h_detEffiNpeMulti4->Draw("SAME,HIST");

	TCanvas* c_positionMismatch = new TCanvas("c_positionMismatch","PositionMismatch",800,600);
	c_positionMismatch->Divide(2,2);
	c_positionMismatch->cd(1);
	h_mismatchX->Draw();
	c_positionMismatch->cd(2);
	h_mismatchY->Draw();
	c_positionMismatch->cd(3);
	h_mismatchZ->Draw();
	c_positionMismatch->cd(4);

	TCanvas* c_positionVsPositionMismatch = new TCanvas("c_positionVsPositionMismatch","PositionsVsPositionMismatch",800,600);
	c_positionVsPositionMismatch->Divide(2,2);
	c_positionVsPositionMismatch->cd(1);
	h_mismatchXvsX->Draw("colz");
	c_positionVsPositionMismatch->cd(2);
	h_mismatchYvsY->Draw("colz");
	c_positionVsPositionMismatch->cd(3);
	h_mismatchZvsZ->Draw("colz");
	c_positionVsPositionMismatch->cd(4);

	TCanvas* c_posNotReco = new TCanvas("c_posNotReco","PositionNotRecoEvents",800,600);
	c_posNotReco->Divide(2,2);
	c_posNotReco->cd(1);
	h_posXNotReco->Divide(h_x);
	h_posXNotReco->Scale(100);
	h_posXNotReco->Draw("HIST");
	c_posNotReco->cd(2);
	h_posYNotReco->Divide(h_y);
	h_posYNotReco->Scale(100);
	h_posYNotReco->Draw("HIST");
	c_posNotReco->cd(3);
	h_posZNotReco->Divide(h_z);
	h_posZNotReco->Scale(100);
	h_posZNotReco->Draw("HIST");
	c_posNotReco->cd(4);
	h_posNotReco->Draw("p0");
	h_posNotReco->SetTitle(";x [cm];y [cm];z [cm]");

	TCanvas* c_averageLY = new TCanvas("c_averageLY","AverageLY",800,600);
	c_averageLY->SetRightMargin(0.20);
	TProfile2D* p2d_averageLY_xy = h_averageLY->Project3DProfile("yx");
	TProfile2D* p2d_averageLY_xz = h_averageLY->Project3DProfile("zx");
	p2d_averageLY_xy->SetTitle(";x [cm]; y [cm]; <LY> [p.e./MeV]");
	p2d_averageLY_xy->Draw("colz");
	p2d_averageLY_xy->SetStats(false);
	p2d_averageLY_xy->SetMinimum(0);

	double avgMinLY = 1000000;
	double avgMaxLY = 0;
	double meanLY = 0;
	int nBins = 0;
	double minLY = 0;

	for (int i = 1; i <= p2d_averageLY_xy->GetXaxis()->GetNbins(); ++i)
	{
		for (int j = 1; j <= p2d_averageLY_xy->GetYaxis()->GetNbins(); ++j)
		{
			// cout << i << " " << j << " " << p2d_averageLY_xy->GetBinContent(i,j) << endl;
			if (p2d_averageLY_xy->GetBinContent(i,j) < avgMinLY)
			{
				avgMinLY = p2d_averageLY_xy->GetBinContent(i,j);
			}
			if (p2d_averageLY_xy->GetBinContent(i,j) > avgMaxLY)
			{
				avgMaxLY = p2d_averageLY_xy->GetBinContent(i,j);
			}
			meanLY += p2d_averageLY_xy->GetBinContent(i,j);
			nBins++;
			double recentMinLY = 0;
			for (int k = 1; k < h_averageLY->GetZaxis()->GetNbins(); ++k)
			{
				if (h_averageLY->GetBinContent(i,j,k) != 0)
				{
					recentMinLY = k-1;
					break;
				}
			}
			// cout << recentMinLY << endl;
			minLY += recentMinLY;
		}
	}
	// cout << nBins << endl;
	minLY /= nBins;
	meanLY /= nBins;

	cout << "<LY>_min = " << avgMinLY << endl;
	cout << "<LY>_max = " << avgMaxLY << endl;
	cout << "<LY> = " << meanLY << endl;
	cout << "<LY_min> = " << minLY << endl;

	TH1D* zProjection = h_averageLY->ProjectionZ();

	TLatex latex;
	latex.SetTextSize(0.035);
	latex.DrawLatex(50,600,Form("Eff_{PD} = %.01f %%",gArapucaEfficiency*100));
	latex.DrawLatex(50,550,Form("<LY> = %.02f p.e.",meanLY));
	latex.DrawLatex(50,500,Form("<LY>_{min} = %.02f p.e.",minLY));
	latex.DrawLatex(50,450,Form("N(LY < 0.5) = %.02f %%",zProjection->GetBinContent(1)/zProjection->Integral()*100));

	TCanvas* c_averageLYProj = new TCanvas("c_averageLYProj","AverageLYProjections",800,800);
	c_averageLYProj->Divide(2,2);
	c_averageLYProj->cd(1);
	p2d_averageLY_xy->ProfileX()->Draw();
	c_averageLYProj->cd(2);
	p2d_averageLY_xy->ProfileY()->Draw();
	c_averageLYProj->cd(3);
	p2d_averageLY_xz->ProfileY()->Draw();
	c_averageLYProj->cd(4);
	// h_averageLY->Project3D("z")->Draw();
	zProjection->Draw();

}

void SaveResults(int inputFile)
{
	TString outputFileName = Form("rootResults/effi%d_hor%d_vert%d_oppo%d_maCentral%d_inFile%d.root",(int)(gArapucaEfficiency*100),gTilesHorizontally,gTilesVertically,gTilesOpposite,gMaskCentralTiles,inputFile);
	TFile* outputFile = new TFile(outputFileName,"RECREATE");

	h_nPETotal->Write();
	h_nPECryo->Write();
	h_nPECathode->Write();
	h_nPECryoPlus->Write();
	h_nPECryoMinus->Write();
	h_PECryoRatio->Write();

	h_nPEMaxMegacell->Write();
	h_nPEMaxMegacellX->Write();
	h_nPEMaxMegacellY->Write();
	h_nPEMaxMegacellZ->Write();
	h_multiplicity->Write();
	h_multiplicityAbove5pe->Write();
	h_multiplicityAbove3pe->Write();
	h_multiplicityAbove2pe->Write();
	h_multiplicityFiltered->Write();
	h_multiplicityX->Write();
	h_multiplicityXAboveNpe->Write();

	h_energy->Write();
	h_x->Write();
	h_y->Write();
	h_z->Write();
	h_nHits->Write();

	h_cellID->Write();

	h_detEffi5pe->Write();
	h_detEffi10pe->Write();
	h_detEffi15pe->Write();

	h_detEffiNpeMulti2->Write();
	h_detEffiNpeMulti3->Write();
	h_detEffiNpeMulti4->Write();
	h_detEffiNpeMulti5->Write();

	h_detEffiNpeMulti2Y->Write();
	h_detEffiNpeMulti3Y->Write();
	h_detEffiNpeMulti4Y->Write();
	h_detEffiNpeMulti5Y->Write();

	h_detEffi5peY->Write();
	h_detEffi5peZ->Write();

	outputFile->Close();
}

double AnalyzePEDist(Event &event, unsigned short int &nPETotal, unsigned short int &nPECathode, unsigned short int &nPECryo)
{
	nPETotal = nPECathode = nPECryo = 0;
	unsigned short int nPECryoPlus = 0;
	unsigned short int nPECryoMinus = 0;
	for (unsigned int i = 0; i < event.v_hits.size(); ++i)
	{
		h_cellID->Fill(event.v_hits[i].cellID);
		nPETotal += event.v_hits[i].nPE;
		// if (event.v_hits[i].cellID < 10000)
			// nPECryo += event.v_hits[i].nPE;
		// else if (event.v_hits[i].cellID < 20000)
			// nPECathode += event.v_hits[i].nPE;
		if (event.v_hits[i].cellID <= 1200)
			nPECryoPlus += event.v_hits[i].nPE;
		if (event.v_hits[i].cellID > 1200 && event.v_hits[i].cellID <= 2400)
			nPECryoMinus += event.v_hits[i].nPE;
	}
	h_nPECryoPlus->Fill(nPECryoPlus);
	h_nPECryoMinus->Fill(nPECryoMinus);
	return nPETotal!=0?(double)nPECryo/nPETotal:-1;
}

unsigned short int GetNPEMegaCell(Event &event)
{
	unsigned short int maxPE = 0;
	for (unsigned int i = 0; i < event.v_hits.size(); ++i)
	{
		if (event.v_hits[i].cellID > 2400)
			continue;
		if (event.v_hits[i].nPE > maxPE)
		{
			maxPE = event.v_hits[i].nPE;
		}
	}
	return maxPE;
}

int CreateDetEffi(TH2D* h_nPEMaxMegacellX,TH2D* h_nPEMaxMegacellY,TH2D* h_nPEMaxMegacellZ)
{
	// cout << h_nPEMaxMegacellX->GetXaxis()->GetNbins() << endl;
	int yBins = h_nPEMaxMegacellX->GetYaxis()->GetNbins();
	for (int i = 0; i < h_nPEMaxMegacellX->GetXaxis()->GetNbins(); ++i)
	{
		int yIntegral = h_nPEMaxMegacellX->Integral(i+1,i+1,1,yBins);
		if (yIntegral != 0)
		{
			h_detEffi5pe->Fill(h_nPEMaxMegacellX->GetXaxis()->GetBinCenter(i+1),(double)h_nPEMaxMegacellX->Integral(i+1,i+1,6,yBins)/yIntegral);
			h_detEffi10pe->Fill(h_nPEMaxMegacellX->GetXaxis()->GetBinCenter(i+1),(double)h_nPEMaxMegacellX->Integral(i+1,i+1,11,yBins)/yIntegral);
			h_detEffi15pe->Fill(h_nPEMaxMegacellX->GetXaxis()->GetBinCenter(i+1),(double)h_nPEMaxMegacellX->Integral(i+1,i+1,16,yBins)/yIntegral);
		}
		// h_detEffi5pe->Fill(h_nPEMaxMegacellX->GetXaxis()->GetBinCenter(i+1),h_nPEMaxMegacellX->Integral(i+1,i+1,1,yBins));
	}

	int yBinsY = h_nPEMaxMegacellY->GetYaxis()->GetNbins();
	for (int i = 0; i < h_nPEMaxMegacellY->GetXaxis()->GetNbins(); ++i)
	{
		int yIntegral = h_nPEMaxMegacellY->Integral(i+1,i+1,1,yBinsY);
		if (yIntegral != 0)
		{
			h_detEffi5peY->Fill(h_nPEMaxMegacellY->GetXaxis()->GetBinCenter(i+1),(double)h_nPEMaxMegacellY->Integral(i+1,i+1,6,yBinsY)/yIntegral);
		}
		// h_detEffi5pe->Fill(h_nPEMaxMegacellX->GetXaxis()->GetBinCenter(i+1),h_nPEMaxMegacellX->Integral(i+1,i+1,1,yBins));
	}

	int yBinsZ = h_nPEMaxMegacellZ->GetYaxis()->GetNbins();
	for (int i = 0; i < h_nPEMaxMegacellZ->GetXaxis()->GetNbins(); ++i)
	{
		int yIntegral = h_nPEMaxMegacellZ->Integral(i+1,i+1,1,yBinsZ);
		if (yIntegral != 0)
		{
			h_detEffi5peZ->Fill(h_nPEMaxMegacellZ->GetXaxis()->GetBinCenter(i+1),(double)h_nPEMaxMegacellZ->Integral(i+1,i+1,6,yBinsZ)/yIntegral);
		}
		// h_detEffi5pe->Fill(h_nPEMaxMegacellX->GetXaxis()->GetBinCenter(i+1),h_nPEMaxMegacellX->Integral(i+1,i+1,1,yBins));
	}
	return 0;
}

int CreateDetEffiMulti(TH2D* h_multiplicityXAboveNpe,TH2D* h_multiplicityYAboveNpe)
{
	int yBins = h_multiplicityXAboveNpe->GetYaxis()->GetNbins();
	for (int i = 0; i < h_multiplicityXAboveNpe->GetXaxis()->GetNbins(); ++i)
	{
		int yIntegral = h_multiplicityXAboveNpe->Integral(i+1,i+1,1,yBins);
		if (yIntegral != 0)
		{
			h_detEffiNpeMulti1->Fill(h_multiplicityXAboveNpe->GetXaxis()->GetBinCenter(i+1),(double)h_multiplicityXAboveNpe->Integral(i+1,i+1,2,yBins)/yIntegral);
			h_detEffiNpeMulti2->Fill(h_multiplicityXAboveNpe->GetXaxis()->GetBinCenter(i+1),(double)h_multiplicityXAboveNpe->Integral(i+1,i+1,3,yBins)/yIntegral);
			h_detEffiNpeMulti3->Fill(h_multiplicityXAboveNpe->GetXaxis()->GetBinCenter(i+1),(double)h_multiplicityXAboveNpe->Integral(i+1,i+1,4,yBins)/yIntegral);
			h_detEffiNpeMulti4->Fill(h_multiplicityXAboveNpe->GetXaxis()->GetBinCenter(i+1),(double)h_multiplicityXAboveNpe->Integral(i+1,i+1,5,yBins)/yIntegral);
			h_detEffiNpeMulti5->Fill(h_multiplicityXAboveNpe->GetXaxis()->GetBinCenter(i+1),(double)h_multiplicityXAboveNpe->Integral(i+1,i+1,6,yBins)/yIntegral);
		}
		// h_detEffi5pe->Fill(h_nPEMaxMegacellX->GetXaxis()->GetBinCenter(i+1),h_nPEMaxMegacellX->Integral(i+1,i+1,1,yBins));
	}
	int yBinsY = h_multiplicityYAboveNpe->GetYaxis()->GetNbins();
	for (int i = 0; i < h_multiplicityYAboveNpe->GetXaxis()->GetNbins(); ++i)
	{
		int yIntegral = h_multiplicityYAboveNpe->Integral(i+1,i+1,1,yBinsY);
		if (yIntegral != 0)
		{
			h_detEffiNpeMulti1Y->Fill(h_multiplicityYAboveNpe->GetXaxis()->GetBinCenter(i+1),(double)h_multiplicityYAboveNpe->Integral(i+1,i+1,2,yBinsY)/yIntegral);
			h_detEffiNpeMulti2Y->Fill(h_multiplicityYAboveNpe->GetXaxis()->GetBinCenter(i+1),(double)h_multiplicityYAboveNpe->Integral(i+1,i+1,3,yBinsY)/yIntegral);
			h_detEffiNpeMulti3Y->Fill(h_multiplicityYAboveNpe->GetXaxis()->GetBinCenter(i+1),(double)h_multiplicityYAboveNpe->Integral(i+1,i+1,4,yBinsY)/yIntegral);
			h_detEffiNpeMulti4Y->Fill(h_multiplicityYAboveNpe->GetXaxis()->GetBinCenter(i+1),(double)h_multiplicityYAboveNpe->Integral(i+1,i+1,5,yBinsY)/yIntegral);
			h_detEffiNpeMulti5Y->Fill(h_multiplicityYAboveNpe->GetXaxis()->GetBinCenter(i+1),(double)h_multiplicityYAboveNpe->Integral(i+1,i+1,6,yBinsY)/yIntegral);
		}
		// h_detEffi5pe->Fill(h_nPEMaxMegacellX->GetXaxis()->GetBinCenter(i+1),h_nPEMaxMegacellX->Integral(i+1,i+1,1,yBinsY));
	}
	return 0;
}

int GetMultiplicity(Event &event)
{
	int nHits = 0;
	for (unsigned int i = 0; i < event.v_hits.size(); ++i)
	{
		// if (event.v_hits[i].cellID > 2400)
			// continue;
		nHits++;
	}
	return nHits;
}

int GetMultiplicityAbove(Event &event, int threshold)
{
	int nHits = 0;
	for (unsigned int i = 0; i < event.v_hits.size(); ++i)
	{
		if (event.v_hits[i].nPE < threshold)
			continue;
		nHits++;
	}
	return nHits;
}

bool SignalOpositteSides(Event &event,int threshold)
{
	bool positiveSide = false;
	bool negativeSide = false;
	for (unsigned int i = 0; i < event.v_hits.size(); ++i)
	{
		if (event.v_hits[i].cellID < 401 && event.v_hits[i].nPE >= threshold)
			negativeSide = true;
		if (event.v_hits[i].cellID > 400 && event.v_hits[i].cellID < 801 && event.v_hits[i].nPE >= threshold)
			positiveSide = true;
	}
	return (positiveSide && negativeSide);
}

double CalculateSensitiveVolume()
{
	int nCubes = 0;
	for (int i = 0; i < nCubesX; ++i)
	{
		for (int j = 0; j < nCubesY; ++j)
		{
			for (int k = 0; k < nCubesZ; ++k)
			{
				// cout << DetectedEvents[i][j][k] << "\t" << AllEvents[i][j][k] << "\t" << (double)DetectedEvents[i][j][k]/AllEvents[i][j][k] << "\t" << ((double)DetectedEvents[i][j][k]/AllEvents[i][j][k] > 0.99) << endl;
				if (AllEvents[i][j][k] != 0 && (double)DetectedEvents[i][j][k]/AllEvents[i][j][k] > 0.99)
					nCubes++;
			}
		}
	}
	cout << nCubes << endl;
	cout << "Efective volume with detection efficiency above 99%: " << nCubes*0.25*0.25*1 << " m3" << endl;
	return nCubes*0.25*0.25*1;
}

unsigned short int FindMaxCellID(Event &event)
{
	unsigned short int maxNPE = 0;
	unsigned short int maxCellID = -1;
	for (unsigned int i = 0; i < event.v_hits.size(); ++i)
	{
		if (event.v_hits[i].nPE > maxNPE)
		{
			maxNPE = event.v_hits[i].nPE;
			maxCellID = event.v_hits[i].cellID;
		}
	}
	return maxCellID;
}

TVector3 GetEventRecoPositionMean(Event &event, vector<TVector3> positions)
{
	int geometryCellID = -1;
	TVector3 recoPos = TVector3{0,0,0};
	double nPETotal = 0;
	double xPos = 0;
	double xNorm = 0;
	// cout << event.x << "\t" << event.y << "\t" << event.z << "\t" << event.nHits << "\t" << event.nPEWalls << endl;

	unsigned short int maxSignalCellID = FindMaxCellID(event);
	unsigned short int maxSignalColumn = ((maxSignalCellID-1)%400)/20;
	unsigned short int maxSignalRow = ((maxSignalCellID-1)%400)%20;
	unsigned short int maxDiffColumn = 19-maxSignalColumn;
	unsigned short int maxDiffRow = maxSignalRow;

	// cout << "Main: " << maxSignalCellID << "\t" << maxSignalColumn << "\t" << maxSignalRow << "\t" << maxDiffColumn << "\t" << maxDiffRow << endl;


	for (unsigned int i = 0; i < event.v_hits.size(); ++i)
	{
		geometryCellID = event.v_hits[i].cellID-1;
		// cout << "\t" << event.v_hits[i].cellID << "\t" << (short int)((geometryCellID%400)%20) << "\t" << (short int)maxSignalRow << "\t" << TMath::Abs(((short int)(geometryCellID%400)%20)-((short int)maxSignalRow)) << endl;
		// if (TMath::Abs((short int)((geometryCellID%400)/20)-((short int)maxSignalColumn)) > maxDiffColumn)
			// continue;
		// if (TMath::Abs((short int)((geometryCellID%400)%20)-((short int)maxSignalRow)) > maxDiffRow)
		// {
			// cout << "CONTINUE" << endl;
			// continue;
		// }
		// recoPos += positions[geometryCellID]*(TMath::Log(TMath::Sqrt(event.v_hits[i].nPE)*7.0)+1);
		xPos += positions[geometryCellID].X()*(TMath::Log(TMath::Sqrt(event.v_hits[i].nPE)*7.0)+1);
		// nPETotal += TMath::Log(TMath::Sqrt(event.v_hits[i].nPE)*21.0)+1;
		xNorm += TMath::Log(TMath::Sqrt(event.v_hits[i].nPE)*7.0)+1;
		recoPos += positions[geometryCellID]*event.v_hits[i].nPE;
		nPETotal += event.v_hits[i].nPE;
		// recoPos += positions[geometryCellID];
		// nPETotal += 1;
		// cout << "\t" << geometryCellID << "\t" << event.v_hits[i].nPE << endl;
		// positions[geometryCellID].Print();
	}
	recoPos = recoPos*(1.0/nPETotal);
	recoPos.SetX(xPos*(1.0/xNorm));
	// cout << "FINAL" << nPETotal << endl;
	// recoPos.Print();
	// if (nPETotal > 0 && event.y < 300 && event.z < 2000)
	return recoPos;
}

TVector3 GetEventRecoPositionMax(Event &event, vector<TVector3> positions)
{
	int geometryCellID = -1;
	TVector3 recoPos = TVector3{0,0,0};

	int chargeMax = 0;
	// cout << event.x << "\t" << event.y << "\t" << event.z << "\t" << event.nHits << "\t" << event.nPEWalls << endl;
	for (unsigned int i = 0; i < event.v_hits.size(); ++i)
	{
		if (event.v_hits[i].nPE > chargeMax)
		{
			recoPos = positions[event.v_hits[i].cellID-1];
			chargeMax = event.v_hits[i].nPE;
		}
		// cout << "\t" << geometryCellID << "\t" << event.v_hits[i].nPE << endl;
		// positions[geometryCellID].Print();
	}
	// cout << "FINAL" << nPETotal << endl;
	// recoPos.Print();
	return recoPos;
}

int AnalyzeData(vector<Event> &events, vector<TVector3> &positions, int fileID)
{
	unsigned short int nPETotal, nPECathode, nPECryo;
	double peCryoRatio;
	int nDetectedEvents = 0;
	int nDetectedEventsBelow4m = 0;
	int nEventsBelow4m = 0;

	for (unsigned int i = 0; i < events.size(); ++i)
	{
		// PrintEvent(events[i]);
		peCryoRatio = AnalyzePEDist(events[i],nPETotal,nPECathode,nPECryo);
		if (nPETotal != 0)
		{
			h_nPETotal->Fill(nPETotal);
			h_nPECathode->Fill(nPECathode);
			h_nPECryo->Fill(nPECryo);
			h_PECryoRatio->Fill(peCryoRatio*100);
			h_nPEAnode->Fill(events[i].nPEAnode);
			if (events[i].nPEAnode > 4)
				h_anodeEvents->Fill(events[i].x,events[i].y);
		}

		unsigned int short nPEMaxMegacell = GetNPEMegaCell(events[i]);
		h_nPEMaxMegacell->Fill(nPEMaxMegacell);
		h_nPEMaxMegacellX->Fill(675-events[i].x,nPEMaxMegacell);
		h_nPEMaxMegacellY->Fill(events[i].y,nPEMaxMegacell);
		h_nPEMaxMegacellZ->Fill(events[i].z,nPEMaxMegacell);
		int multi = GetMultiplicity(events[i]);
		int multiAbove5pe = GetMultiplicityAbove(events[i],5);
		int multiAbove3pe = GetMultiplicityAbove(events[i],3);
		int multiAbove2pe = GetMultiplicityAbove(events[i],2);
		int multiAboveNpe = GetMultiplicityAbove(events[i],gPEThreshold);
		h_multiplicity->Fill(multi);
		h_multiplicityAbove4pe->Fill(GetMultiplicityAbove(events[i],4));
		h_multiplicityAbove2pe->Fill(multiAbove2pe);
		h_multiplicityAbove3pe->Fill(multiAbove3pe);
		h_multiplicityAbove5pe->Fill(multiAbove5pe);
		if (multiAbove5pe > 0)
			h_multiplicityFiltered->Fill(multi);
		h_multiplicityX->Fill(675-events[i].x,multi);
		h_multiplicityXAboveNpe->Fill(675-events[i].x,multiAboveNpe);
		h_multiplicityYAboveNpe->Fill(events[i].y,multiAboveNpe);
		h_multiplicityVsTotalNPE->Fill(nPETotal,multiAbove2pe);

		if (multiAbove2pe < gMultiThreshold2pe)
		{
			h_posXNotReco->Fill(events[i].x);
			h_posYNotReco->Fill(events[i].y);
			h_posZNotReco->Fill(events[i].z);
			h_posNotReco->SetPoint(h_posNotReco->GetN(),events[i].x,events[i].y,events[i].z);
		}

		h_energy->Fill(events[i].energy);
		h_x->Fill(events[i].x);
		h_y->Fill(events[i].y);
		h_z->Fill(events[i].z);
		h_nHits->Fill(events[i].nHits);

		// AllEvents[(int)events[i].x/25][(int)events[i].y/25][(int)events[i].z/25] += 1;
		h_detEffiFinalAll->Fill(events[i].x);
		// if (multiAbove > 1 && multi > 30)
		// if (multiAbove > 1)
		// if (multiAboveNpe >= gMulti && multi > 30)
		bool bothSides = SignalOpositteSides(events[i],2);
		// if (multiAboveNpe >= gMulti && bothSides)
		if (multiAboveNpe >= gMulti)
		{
			nDetectedEvents++;
			h_detEffiFinalDet->Fill(events[i].x);
			if (675-events[i].x < 400)
			{
				nDetectedEventsBelow4m++;
				h_xDet->Fill(events[i].x);
			}
			// PrintEvent(events[i]);
			// DetectedEvents[(int)events[i].x/25][(int)events[i].y/25][(int)events[i].z/25] += 1;
		}

		if (675-events[i].x < 400)
		{
			nEventsBelow4m++;
		}

		// if (!gReadAnodeTiles && !gReadCathodeTiles)
		// {
		// 	TVector3 recoPos = GetEventRecoPositionMean(events[i],positions);
		// 	// TVector3 recoPos = GetEventRecoPositionMax(events[i],positions);
		// 	if (nPETotal > 0)
		// 	{
		// 		h_mismatchX->Fill(events[i].x-recoPos.x());
		// 		h_mismatchXvsX->Fill(events[i].x,events[i].x-recoPos.x());
		// 		h_mismatchY->Fill(events[i].y-recoPos.y());
		// 		h_mismatchYvsY->Fill(events[i].y,events[i].y-recoPos.y());
		// 		h_mismatchZ->Fill(events[i].z-recoPos.z());
		// 		h_mismatchZvsZ->Fill(events[i].z,events[i].z-recoPos.z());
		// 	}
		// }

		// if (events[i].y < 51 && events[i].y > 49)
		// {
			// cout << i << endl;
			// PrintEvent(events[i]);
		// }
		// cout << events[i].nPECathode << "\t" << events[i].nPEWalls << "\t" << events[i].energy << "\t" << (events[i].nPECathode+events[i].nPEWalls)/events[i].energy << endl;
		// if (events[i].z < 100 && events[i].z > 0)
		{
			h_averageLY->Fill(events[i].x,events[i].y,(events[i].nPECathode+events[i].nPEWalls)/events[i].energy);
			h_averageLY->Fill(-events[i].x,events[i].y,(events[i].nPECathode+events[i].nPEWalls)/events[i].energy);
		}
	}
	CreateDetEffi(h_nPEMaxMegacellX,h_nPEMaxMegacellY,h_nPEMaxMegacellZ);
	CreateDetEffiMulti(h_multiplicityXAboveNpe,h_multiplicityYAboveNpe);

	// CalculateSensitiveVolume();

	std::cout << "Detection Efficiency: " << (double)nDetectedEvents/events.size()*100.0 << " %" << std::endl;
	std::cout << "Detection Efficiency below 4 m: " << (double)nDetectedEventsBelow4m/nEventsBelow4m*100.0 << " %" << std::endl;

	return nDetectedEvents;
}

unsigned short int GetNPEDetected(int nPE, unsigned short int cellID)
{
	int nPEDetected = 0;
	for (int i = 0; i < nPE; ++i)
	{
		double randomNumber = r_randNumGen3.Uniform(1);
		if (gUseCathodeTileTransparency && cellID > 9999 && cellID < 19999)
		{
			if (randomNumber < gArapucaEfficiency*0.8)
				nPEDetected++;
		}else
		{
			if (randomNumber < gArapucaEfficiency)
				nPEDetected++;
		}
	}
	return nPEDetected;
}

int IsAlreadyInHits(vector <Hit> &hits, unsigned short int cellID)
{
	for (unsigned int i = 0; i < hits.size(); ++i)
	{
		if (hits[i].cellID == cellID)
			return i;
	}
	return -1;
}

int ReadData(TString filePath, vector<Event> &events, bool tilesHorizontally = false, bool tilesVertically = false, bool tilesOpposite = false)
{
	ifstream inputFile;
    inputFile.open(filePath);

    if (!inputFile)
    {
    	cerr << "Input file: " << filePath << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

  	if ((tilesVertically && tilesHorizontally) || (tilesVertically && tilesOpposite) || (tilesOpposite && tilesHorizontally))
  	{
  		cerr << "You can't merge megacells twice in the same time!" << endl;
  		return -2;
  	}

  	std::string line;
  	int eventID,nHits;
  	float energy,x,y,z;
  	unsigned short int cellID;
  	unsigned int nPE;

    while(std::getline(inputFile, line))
    {
    	// cout << line << endl;
    	// cout << line.at(0) << endl;
    	if (line.at(0) == ' ')
    		break;
    	// std::cout << line << std::endl;
        // Create a stringstream of the current line
        std::stringstream ss(line);
        ss.ignore();
        ss.ignore();
        ss >> eventID >> energy >> x >> y >> z >> nHits;
        events.push_back(Event{eventID,energy*1000,x,y,z,nHits,{},0,0,0});
        // std::cout << eventID << "\t" << energy << "\t" << x << "\t" << y << "\t" << z << "\t" << nHits << std::endl;
        int nHitsWithEfficiency = 0;
        unsigned int nPEGeneratedTotal = 0;
        unsigned int nPEDetectedTotal = 0;
        int nHitsWithoutMerging = 0;
        for (int i = 0; i < nHits; ++i)
        {
        	inputFile >> cellID >> nPE;
    		unsigned short int nPEDetected = GetNPEDetected(nPE,cellID);
        	// cout << cellID << "\t" << nPE << "\t" << nPEDetected << endl;
        	if (!gReadCathodeTiles && cellID > 9999 && cellID < 19999)
        		continue;
        	if (!gReadAnodeTiles && cellID > 19999)
        		continue;
        	if (!gReadCryoTiles && cellID < 9999)
        		continue;

    		nPEGeneratedTotal += nPE;
    		nPEDetectedTotal += nPEDetected;

    		if (nPEDetected > 0)
    		{
    			nHitsWithoutMerging++;
    			if (cellID > 9999)
    			{
					events.back().v_hits.push_back(Hit{cellID,nPEDetected});
					if (cellID > 19999)
						events.back().nPEAnode += nPEDetected;
					else
						events.back().nPECathode += nPEDetected;
					nHitsWithEfficiency++;
					continue;
    			}
    			if (tilesHorizontally)
    			{
    				int column = (cellID-1)/20;
    				unsigned short int searchedCellID = (column%3==0)?cellID:(cellID-(20*(column%3)));
    				int SNInHits = IsAlreadyInHits(events.back().v_hits,searchedCellID);
    				if (SNInHits != -1)
    					events.back().v_hits[SNInHits].nPE += nPEDetected;
    				else
    				{
    					nHitsWithEfficiency++;
    					events.back().v_hits.push_back(Hit{searchedCellID,nPEDetected});
    				}
    			} else if (tilesVertically)
    			{
    				int row = (cellID-1)%20;
    				unsigned short int searchedCellID = (row%3==0)?cellID:(cellID-((row%3)));
    				int SNInHits = IsAlreadyInHits(events.back().v_hits,searchedCellID);
    				if (SNInHits != -1)
    					events.back().v_hits[SNInHits].nPE += nPEDetected;
    				else
    				{
    					nHitsWithEfficiency++;
    					events.back().v_hits.push_back(Hit{searchedCellID,nPEDetected});
    				}
    			} else if (tilesOpposite)
    			{
    				unsigned short int searchedCellID = ((cellID-1)%1200)+1;
    				int SNInHits = IsAlreadyInHits(events.back().v_hits,searchedCellID);
    				if (SNInHits != -1)
    					events.back().v_hits[SNInHits].nPE += nPEDetected;
    				else
    				{
    					nHitsWithEfficiency++;
    					events.back().v_hits.push_back(Hit{searchedCellID,nPEDetected});
    				}
    			} else if (gTilesSmall)
    			{
    				int row = (cellID-1)%20;
    				unsigned short int searchedCellID = (row%5==0)?cellID:(cellID-((row%5)));
    				int SNInHits = IsAlreadyInHits(events.back().v_hits,searchedCellID);
    				if (SNInHits != -1)
    					events.back().v_hits[SNInHits].nPE += nPEDetected;
    				else
    				{
    					nHitsWithEfficiency++;
    					events.back().v_hits.push_back(Hit{searchedCellID,nPEDetected});
    				}
    			} else
    			{
    				int row = (cellID-1)%20;
    				if (gMaskCentralTiles && (row == 9 || row == 10))
    					continue;
    				if (gMaskMiddleTiles && (row > 3 && row < 16))
    					continue;
    				if (gMaskNegativeTiles && (row > 9))
    					continue;
        			nHitsWithEfficiency++;
		        	events.back().v_hits.push_back(Hit{cellID,nPEDetected});
    			}
    			events.back().nPEWalls += nPEDetected;
    		}
        }
        if (nHits % 10 == 0 && nHits != 0)
	        std::getline(inputFile, line);
        std::getline(inputFile, line);
    	events.back().nHits = nHitsWithEfficiency;
        // std::sort(events.back().v_hits,OrderByPE);
        // PrintEvent(events.back());
        if (nPEGeneratedTotal != 0)
        {
		    h_peRatio->Fill((double)nPEDetectedTotal/nPEGeneratedTotal*100.0);
        }
        h_nHitsDiff->Fill(nHitsWithoutMerging-nHitsWithEfficiency);
        if (events.back().y >= 650.0)
        {
        	events.pop_back();
        }
    }
	return 0;
}

int CreateMegaCellGeometry(vector<TVector3> &megaCellPos)
{
	int nMegaCells = 2400;
	double x,y,z;

	for (int i = 0; i < nMegaCells; ++i)
	{
		if (i < 1200)
			x = -675 - 60;
		else
			x = 675 + 60;
		int row = (i)%20;
		y = 618.5 - row*65;
		int vertTrip = (i%1200)/60;
		int column = (i%1200)/20;
		z = -2850 + 300*vertTrip;
		if (column%3 == 0)
			z -= 68;
		if (column%3 == 2)
			z += 68;

		if (gTilesVertically)
			y -= 65;
		if (gTilesHorizontally)
			z += 68;

		megaCellPos.push_back(TVector3(x,y,z));
	}

	return 0;
}

int CreateTileGeometry(vector<TVector3> &tilePos)
{
	int nMegaCells = 800;
	double x,y,z;

	for (int i = 0; i < nMegaCells; ++i)
	{
		if (i < 400)
			x = -675 - 60;
		else
			x = 675 + 60;
		int row = (i)%20;
		y = 618.5 - row*65;
		int column = (i%400)/20;
		z = -2850 + 300*column;

		tilePos.push_back(TVector3(x,y,z));
	}

	return 0;
}

int CreateGeometry(vector<TVector3> &megaCellPos, int fileID)
{
	if (fileID < 10)
		CreateMegaCellGeometry(megaCellPos);
	else
		CreateTileGeometry(megaCellPos);

	return 0;
}

void PrintGeometry(vector<TVector3> &megaCellPos)
{
	for (unsigned int i = 0; i < megaCellPos.size(); ++i)
	{
		cout << i << "\t" << megaCellPos[i].x() << "\t" << megaCellPos[i].y() << "\t" << megaCellPos[i].z() << endl;
	}
}

int processDataDune(int fileID = 0)
{
	gStyle->SetPalette(55);
	TString filePath = "/Data/duneData/";

	switch (fileID)
	{
		case 0:
			filePath +=	"dune_Ar39_cato_anode_fulla_tiles.data";
			break;
		case 1:
			filePath += "dune_cato_anode_5MeValla_tiles.data";
			break;
		case 2:
			filePath += "dune_cato_anode_10MeValla_tiles.data";
			break;
		case 3:
			filePath += "dune_cato_5MeV_alla_tiles.data";
			break;
		// Old MC with fixed bugs and anode tiles
		case 10:
			filePath += "dune_alltiles_Ar39_fulla_tiles.data";
			break;
		case 11:
			filePath += "dune_alltiles_5MeValla_tiles.data";
			break;
		case 12:
			filePath += "dune_alltiles_10MeValla_tiles.data";
			break;
		// New MC - zero cathode transparency, real FC, longer scattering Xe length
		case 20:
			filePath += "dune_fieldcage_Ar39_fulla_tiles.data";
			break;
		case 21:
			filePath += "dune_fieldcage_5MeValla_tiles.data";
			break;
		case 22:
			filePath += "dune_fieldcage_10MeValla_tiles.data";
			break;
		// back to uniform FC + 0 transparency of the cathode
		case 32:
			filePath += "dune_nofieldcage_10MeValla_tiles.data";
			break;
		// real FC + 80% cathode transparency
		case 42:
			filePath += "dune_fieldcage_ctrasp_10MeValla_tiles.data";
			break;
		// Small tiles with the new MC
		case 50:
			filePath += "dune_fieldcage_smalltiles_Ar39fulla_tiles.data";
			break;
		case 52:
			filePath += "dune_fieldcage_smalltiles_10MeValla_tiles.data";
			break;
	}

	vector<TVector3> v_megaCellPos;

	CreateGeometry(v_megaCellPos,fileID);
	// PrintGeometry(v_megaCellPos);

	vector<Event> v_myEvents;
	if (ReadData(filePath,v_myEvents,gTilesHorizontally,gTilesVertically,gTilesOpposite) < 0)
		return -1;

	// ReadData("/Data/duneData/dune_fieldcage_Ar39_fulla_tiles.data",v_myEvents,gTilesHorizontally,gTilesVertically,gTilesOpposite);
	// ReadData("/Data/duneData/dune_fieldcage_10MeValla_tiles.data",v_myEvents,gTilesHorizontally,gTilesVertically,gTilesOpposite);


	std::cout << "Number of read events: " << v_myEvents.size() << std::endl;

	int nDetectedEvents = AnalyzeData(v_myEvents,v_megaCellPos,fileID);



	DrawResults();
	SaveResults(fileID);

	return 0;
}