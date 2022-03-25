#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRatioPlot.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TLegend.h"

#include <iostream>

int compareDataDune()
{
	TString tilesVertically5 = "rootResults/effi5_hor0_vert0_oppo0_inFile11.root";
	TString tilesHorizontally5 = "rootResults/effi5_hor1_vert0_oppo0_inFile1.root";
	TString tilesOpposite5 = "rootResults/effi5_hor0_vert0_oppo1_inFile1.root";
	TString megaCells5 = "rootResults/effi5_hor0_vert0_oppo0_inFile1.root";
	TString tilesVertically3 = "rootResults/effi3_hor0_vert0_oppo0_inFile11.root";
	TString tilesHorizontally3 = "rootResults/effi3_hor1_vert0_oppo0_inFile1.root";
	TString tilesOpposite3 = "rootResults/effi3_hor0_vert0_oppo1_inFile1.root";
	TString megaCells3 = "rootResults/effi3_hor0_vert0_oppo0_inFile1.root";
	TString tilesVertically4 = "rootResults/effi4_hor0_vert0_oppo0_inFile11.root";
	TString tilesVertically6 = "rootResults/effi6_hor0_vert0_oppo0_inFile11.root";
	TString tilesVertically7 = "rootResults/effi7_hor0_vert0_oppo0_inFile11.root";

	TString tilesHorizontally5Bkg = "rootResults/effi5_hor1_vert0_oppo0_inFile0.root";
	TString tilesHorizontally3Bkg = "rootResults/effi3_hor1_vert0_oppo0_inFile0.root";
	TString tilesVertically3Bkg = "rootResults/effi3_hor0_vert1_oppo0_inFile0.root";

	TString tilesCryo3 = "rootResults/effi3_hor0_vert0_oppo0_inFile11.root";
	TString tilesCryo3_10Mev = "rootResults/effi3_hor0_vert0_oppo0_inFile12.root";
	TString tilesCryo3Bkg = "rootResults/effi3_hor0_vert0_oppo0_inFile10.root";

	TString histogramName = "h_detEffiNpeMulti4";

	TFile* f_tilesVertically5 = new TFile(tilesVertically5,"READ");
	TH1F* detEffi5peTV5  = (TH1F*)f_tilesVertically5->Get(histogramName);
	detEffi5peTV5->SetTitle("M_{T} = 4, N_{phe} = 5, Eff_{PD} = 5%");
	detEffi5peTV5->SetLineColor(kRed);

	TFile* f_tilesVertically3 = new TFile(tilesVertically3,"READ");
	TH1F* detEffi5peTV3  = (TH1F*)f_tilesVertically3->Get(histogramName);
	detEffi5peTV3->SetTitle("M_{T} = 4, N_{phe} = 5, Eff_{PD} = 3%");
	detEffi5peTV3->SetLineColor(kMagenta);

	TFile* f_tilesVertically4 = new TFile(tilesVertically4,"READ");
	TH1F* detEffi5peTV4  = (TH1F*)f_tilesVertically4->Get(histogramName);
	detEffi5peTV4->SetTitle("M_{T} = 4, N_{phe} = 5, Eff_{PD} = 4%");
	detEffi5peTV4->SetLineColor(kBlue);

	TFile* f_tilesVertically6 = new TFile(tilesVertically6,"READ");
	TH1F* detEffi5peTV6  = (TH1F*)f_tilesVertically6->Get(histogramName);
	detEffi5peTV6->SetTitle("M_{T} = 4, N_{phe} = 5, Eff_{PD} = 6%");
	detEffi5peTV6->SetLineColor(kCyan);

	TFile* f_tilesVertically7 = new TFile(tilesVertically7,"READ");
	TH1F* detEffi5peTV7  = (TH1F*)f_tilesVertically7->Get(histogramName);
	detEffi5peTV7->SetTitle("M_{T} = 4, N_{phe} = 5, Eff_{PD} = 7%");
	detEffi5peTV7->SetLineColor(kGreen);

	TFile* f_tilesHorizontally5 = new TFile(tilesHorizontally5,"READ");
	TH1F* detEffi5peTH5  = (TH1F*)f_tilesHorizontally5->Get(histogramName);
	detEffi5peTH5->SetTitle("Tiles Horizontally, 5%");
	detEffi5peTH5->SetLineColor(kCyan);
	TH1F* multiplicityTH5  = (TH1F*)f_tilesHorizontally5->Get("h_multiplicity");
	multiplicityTH5->SetTitle("Tiles Horizontally, 5%");
	multiplicityTH5->SetLineColor(kCyan);
	multiplicityTH5->SetMarkerColor(kCyan);
	TH1F* multiplicityFilteredTH5  = (TH1F*)f_tilesHorizontally5->Get("h_multiplicityFiltered");
	multiplicityFilteredTH5->SetTitle("Tiles Horizontally, 5%");
	multiplicityFilteredTH5->SetLineColor(kCyan);
	multiplicityFilteredTH5->SetMarkerColor(kCyan);

	TFile* f_tilesHorizontally3 = new TFile(tilesHorizontally3,"READ");
	TH1F* detEffi5peTH3  = (TH1F*)f_tilesHorizontally3->Get(histogramName);
	detEffi5peTH3->SetTitle("Tiles Horizontally, 3%");
	detEffi5peTH3->SetLineColor(kBlue);
	TH1F* multiplicityTH3  = (TH1F*)f_tilesHorizontally3->Get("h_multiplicity");
	multiplicityTH3->SetTitle("Tiles Horizontally, 3%");
	multiplicityTH3->SetLineColor(kBlue);
	multiplicityTH3->SetMarkerColor(kBlue);
	multiplicityTH3->SetMarkerStyle(kOpenCircle);
	multiplicityTH3->SetLineStyle(21);
	TH1F* multiplicityFilteredTH3  = (TH1F*)f_tilesHorizontally3->Get("h_multiplicityFiltered");
	multiplicityFilteredTH3->SetTitle("Tiles Horizontally, 3%");
	multiplicityFilteredTH3->SetLineColor(kBlue);
	multiplicityFilteredTH3->SetMarkerColor(kBlue);
	multiplicityFilteredTH3->SetMarkerStyle(kOpenCircle);
	multiplicityFilteredTH3->SetLineStyle(21);
	TH1F* multiplicityAbove3peTH3  = (TH1F*)f_tilesHorizontally3->Get("h_multiplicityAbove3pe");
	multiplicityAbove3peTH3->SetTitle("Tiles Horizontally, 3%");
	multiplicityAbove3peTH3->SetLineColor(kBlue);
	multiplicityAbove3peTH3->SetMarkerColor(kBlue);
	multiplicityAbove3peTH3->SetMarkerStyle(kOpenCircle);
	multiplicityAbove3peTH3->SetLineStyle(21);
	TH1F* multiplicityAbove2peTH3  = (TH1F*)f_tilesHorizontally3->Get("h_multiplicityAbove2pe");
	multiplicityAbove2peTH3->SetTitle("Tiles Horizontally, 3%");
	multiplicityAbove2peTH3->SetLineColor(kBlue);
	multiplicityAbove2peTH3->SetMarkerColor(kBlue);
	multiplicityAbove2peTH3->SetMarkerStyle(kOpenCircle);
	multiplicityAbove2peTH3->SetLineStyle(21);

	TFile* f_tilesHorizontally5Bkg = new TFile(tilesHorizontally5Bkg,"READ");
	TH1F* multiplicityTH5Bkg  = (TH1F*)f_tilesHorizontally5Bkg->Get("h_multiplicity");
	multiplicityTH5Bkg->SetTitle("Tiles Horizontally, 5%, 39Ar");
	multiplicityTH5Bkg->SetLineColor(kSpring);
	multiplicityTH5Bkg->SetMarkerColor(kSpring);
	TH1F* multiplicityFilteredTH5Bkg  = (TH1F*)f_tilesHorizontally5Bkg->Get("h_multiplicityFiltered");
	multiplicityFilteredTH5Bkg->SetTitle("Tiles Horizontally, 5%, 39Ar");
	multiplicityFilteredTH5Bkg->SetLineColor(kSpring);
	multiplicityFilteredTH5Bkg->SetMarkerColor(kSpring);

	TFile* f_tilesHorizontally3Bkg = new TFile(tilesHorizontally3Bkg,"READ");
	TH1F* detEffi5peTH3Bkg  = (TH1F*)f_tilesHorizontally3Bkg->Get(histogramName);
	detEffi5peTH3Bkg->SetTitle("Tiles Horizontally, 3%, 39Ar");
	detEffi5peTH3Bkg->SetLineColor(kSpring);
	TH1F* multiplicityTH3Bkg  = (TH1F*)f_tilesHorizontally3Bkg->Get("h_multiplicity");
	multiplicityTH3Bkg->SetTitle("Tiles Horizontally, 3%, 39Ar");
	multiplicityTH3Bkg->SetLineColor(kSpring);
	multiplicityTH3Bkg->SetMarkerColor(kSpring);
	multiplicityTH3Bkg->SetLineStyle(21);
	multiplicityTH3Bkg->SetMarkerStyle(kOpenCircle);
	TH1F* multiplicityAbove3peTH3Bkg  = (TH1F*)f_tilesHorizontally3Bkg->Get("h_multiplicityAbove3pe");
	multiplicityAbove3peTH3Bkg->SetTitle("Tiles Horizontally, 3%, 39Ar");
	multiplicityAbove3peTH3Bkg->SetLineColor(kSpring);
	multiplicityAbove3peTH3Bkg->SetMarkerColor(kSpring);
	multiplicityAbove3peTH3Bkg->SetLineStyle(21);
	multiplicityAbove3peTH3Bkg->SetMarkerStyle(kOpenCircle);
	TH1F* multiplicityAbove2peTH3Bkg  = (TH1F*)f_tilesHorizontally3Bkg->Get("h_multiplicityAbove2pe");
	multiplicityAbove2peTH3Bkg->SetTitle("Tiles Horizontally, 3%, 39Ar");
	multiplicityAbove2peTH3Bkg->SetLineColor(kSpring);
	multiplicityAbove2peTH3Bkg->SetMarkerColor(kSpring);
	multiplicityAbove2peTH3Bkg->SetLineStyle(21);
	multiplicityAbove2peTH3Bkg->SetMarkerStyle(kOpenCircle);
	TH1F* multiplicityFilteredTH3Bkg  = (TH1F*)f_tilesHorizontally3Bkg->Get("h_multiplicityFiltered");
	multiplicityFilteredTH3Bkg->SetTitle("Tiles Horizontally, 3%, 39Ar");
	multiplicityFilteredTH3Bkg->SetLineColor(kSpring);
	multiplicityFilteredTH3Bkg->SetMarkerColor(kSpring);
	multiplicityFilteredTH3Bkg->SetLineStyle(21);
	multiplicityFilteredTH3Bkg->SetMarkerStyle(kOpenCircle);

	TFile* f_tilesVertically3Bkg = new TFile(tilesVertically3Bkg,"READ");
	TH1F* detEffi5peTV3Bkg  = (TH1F*)f_tilesVertically3Bkg->Get(histogramName);
	detEffi5peTV3Bkg->SetTitle("Tiles Vertically, 3%, 39Ar");
	detEffi5peTV3Bkg->SetLineColor(kViolet);

	TFile* f_tilesOpposite5 = new TFile(tilesOpposite5,"READ");
	TH1F* detEffi5peTO5  = (TH1F*)f_tilesOpposite5->Get(histogramName);
	detEffi5peTO5->SetTitle("Tiles Opposite, 5%");
	detEffi5peTO5->SetLineColor(kBlue);

	TFile* f_megaCells5 = new TFile(megaCells5,"READ");
	TH1F* detEffi5peMC5  = (TH1F*)f_megaCells5->Get(histogramName);
	detEffi5peMC5->SetTitle("MegaCells, 5%");
	detEffi5peMC5->SetLineColor(kGreen);

	TFile* f_tilesCryo3 = new TFile(tilesCryo3,"READ");
	TH1F* detEffiCryo3  = (TH1F*)f_tilesCryo3->Get(histogramName);
	detEffiCryo3->SetTitle("5 MeV e-, 3%, M_{T} = 4, N_{phe} = 3");
	detEffiCryo3->SetLineColor(kBlue);

	TFile* f_tilesCryo3_10MeV = new TFile(tilesCryo3_10Mev,"READ");
	TH1F* detEffiCryo3_10MeV  = (TH1F*)f_tilesCryo3_10MeV->Get(histogramName);
	detEffiCryo3_10MeV->SetTitle("10 MeV e-, 3%, M_{T} = 4, N_{phe} = 3");
	detEffiCryo3_10MeV->SetLineColor(kRed);

	TFile* f_tilesCryo3Bkg = new TFile(tilesCryo3Bkg,"READ");
	TH1F* detEffiCryo3Bkg  = (TH1F*)f_tilesCryo3Bkg->Get(histogramName);
	detEffiCryo3Bkg->SetTitle("39Ar, 3%, M_{T} = 4, N_{phe} = 3");
	detEffiCryo3Bkg->SetLineColor(kGreen);



	THStack* s_detEffi5pe = new THStack("s_detEffi5pe","; X [cm]; Detection Efficiency [%]");
	s_detEffi5pe->Add(detEffi5peTV5,"HIST");
	s_detEffi5pe->Add(detEffi5peTH5,"HIST");
	// s_detEffi5pe->Add(detEffi5peTO5,"HIST");
	s_detEffi5pe->Add(detEffi5peMC5,"HIST");
	s_detEffi5pe->Add(detEffi5peTV3,"HIST");
	s_detEffi5pe->Add(detEffi5peTH3,"HIST");

	TCanvas* c_detEffi5pe = new TCanvas("c_detEffi5pe","DetEffi5pe",800,600);
	gPad->SetGrid();
	s_detEffi5pe->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_detEffi5pe->Write();

	THStack* s_detEffi5peVertical = new THStack("s_detEffi5peVertical","; X [cm]; Detection Efficiency [%]");
	s_detEffi5peVertical->Add(detEffi5peTV3,"HIST");
	s_detEffi5peVertical->Add(detEffi5peTV4,"HIST");
	s_detEffi5peVertical->Add(detEffi5peTV5,"HIST");
	s_detEffi5peVertical->Add(detEffi5peTV6,"HIST");
	s_detEffi5peVertical->Add(detEffi5peTV7,"HIST");

	TCanvas* c_detEffi5peVetical = new TCanvas("c_detEffi5peVetical","DetEffi5peVertical",800,600);
	gPad->SetGrid();
	s_detEffi5peVertical->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_detEffi5peVertical->Write();

	THStack* s_multiplicity = new THStack("s_multiplicity","; N_{hitTiles} [#]; NoE [#]");
	s_multiplicity->Add(multiplicityTH3,"E");
	s_multiplicity->Add(multiplicityTH5,"E");
	s_multiplicity->Add(multiplicityTH3Bkg,"E");
	s_multiplicity->Add(multiplicityTH5Bkg,"E");

	TCanvas* c_multiplicity = new TCanvas("c_multiplicity","Multiplicity",800,600);
	gPad->SetGrid();
	s_multiplicity->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_multiplicity->Write();

	THStack* s_multiplicityFiltered = new THStack("s_multiplicityFiltered","; N_{hitTiles} [#]; NoE [#]");
	s_multiplicityFiltered->Add(multiplicityFilteredTH3,"E");
	s_multiplicityFiltered->Add(multiplicityFilteredTH5,"E");
	s_multiplicityFiltered->Add(multiplicityFilteredTH3Bkg,"E");
	s_multiplicityFiltered->Add(multiplicityFilteredTH5Bkg,"E");

	TCanvas* c_multiplicityFiltered = new TCanvas("c_multiplicityFiltered","MultiplicityFiltered",800,600);
	gPad->SetGrid();
	s_multiplicityFiltered->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_multiplicityFiltered->Write();

	THStack* s_multiplicityAbove3pe = new THStack("s_multiplicityAbove3pe","; N_{hitTiles} [#]; NoE [#]");
	s_multiplicityAbove3pe->Add(multiplicityAbove3peTH3,"E");
	s_multiplicityAbove3pe->Add(multiplicityAbove3peTH3Bkg,"E");

	TCanvas* c_multiplicityAbove3pe = new TCanvas("c_multiplicityAbove3pe","MultiplicityAbove3pe",800,600);
	gPad->SetGrid();
	s_multiplicityAbove3pe->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_multiplicityAbove3pe->Write();

	THStack* s_multiplicityAbove2pe = new THStack("s_multiplicityAbove2pe","; N_{hitTiles} [#]; NoE [#]");
	s_multiplicityAbove2pe->Add(multiplicityAbove2peTH3,"E");
	s_multiplicityAbove2pe->Add(multiplicityAbove2peTH3Bkg,"E");

	TCanvas* c_multiplicityAbove2pe = new TCanvas("c_multiplicityAbove2pe","MultiplicityAbove2pe",800,600);
	gPad->SetGrid();
	s_multiplicityAbove2pe->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_multiplicityAbove2pe->Write();

	THStack* s_detEffi5peSB = new THStack("s_detEffi5peSB","; X [cm]; Detection Efficiency [%]");
	s_detEffi5peSB->Add(detEffi5peTH3,"HIST");
	s_detEffi5peSB->Add(detEffi5peTH3Bkg,"HIST");
	s_detEffi5peSB->Add(detEffi5peTV3Bkg,"HIST");

	TCanvas* c_detEffi5peSB = new TCanvas("c_detEffi5peSB","DetEffi5peSB",800,600);
	gPad->SetGrid();
	s_detEffi5peSB->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_detEffi5peSB->Write();

	THStack* s_detEffiSB = new THStack("s_detEffiSB","; X [cm]; Detection Efficiency [%]");
	s_detEffiSB->Add(detEffiCryo3,"HIST");
	s_detEffiSB->Add(detEffiCryo3_10MeV,"HIST");
	s_detEffiSB->Add(detEffiCryo3Bkg,"HIST");

	TCanvas* c_detEffiSB = new TCanvas("c_detEffiSB","DetEffiSB",800,600);
	gPad->SetGrid();
	s_detEffiSB->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_detEffi5peSB->Write();


	TFile* f_oldMC_noMask = new TFile("rootResults/effi3_hor0_vert0_oppo0_maCentral0_inFile12.root","READ");
	TH1F* mult_oldMC_noMask  = (TH1F*)f_oldMC_noMask->Get("h_multiplicityAbove2pe");
	mult_oldMC_noMask->SetTitle("Old MC");
	mult_oldMC_noMask->SetLineColor(kRed);
	mult_oldMC_noMask->Scale(1/mult_oldMC_noMask->Integral());

	TFile* f_oldMC_Mask = new TFile("rootResults/effi3_hor0_vert0_oppo0_maCentral1_inFile12.root","READ");
	TH1F* mult_oldMC_Mask  = (TH1F*)f_oldMC_Mask->Get("h_multiplicityAbove2pe");
	mult_oldMC_Mask->SetTitle("Old MC, masked two central tiles");
	mult_oldMC_Mask->SetLineColor(kRed);
	mult_oldMC_Mask->SetLineStyle(2);
	mult_oldMC_Mask->Scale(1/mult_oldMC_Mask->Integral());

	TFile* f_newMC_noMask = new TFile("rootResults/effi3_hor0_vert0_oppo0_maCentral0_inFile22.root","READ");
	TH1F* mult_newMC_noMask  = (TH1F*)f_newMC_noMask->Get("h_multiplicityAbove2pe");
	mult_newMC_noMask->SetTitle("New MC");
	mult_newMC_noMask->SetLineColor(kBlue);
	mult_newMC_noMask->Scale(1/mult_newMC_noMask->Integral());

	TFile* f_newMC_Mask = new TFile("rootResults/effi3_hor0_vert0_oppo0_maCentral1_inFile22.root","READ");
	TH1F* mult_newMC_Mask  = (TH1F*)f_newMC_Mask->Get("h_multiplicityAbove2pe");
	mult_newMC_Mask->SetTitle("New MC, masked two central tiles");
	mult_newMC_Mask->SetLineColor(kBlue);
	mult_newMC_Mask->SetLineStyle(2);
	mult_newMC_Mask->Scale(1/mult_newMC_Mask->Integral());

	TFile* f_uniformFC_noMask = new TFile("rootResults/effi3_hor0_vert0_oppo0_maCentral0_inFile32.root","READ");
	TH1F* mult_uniformFC_noMask  = (TH1F*)f_uniformFC_noMask->Get("h_multiplicityAbove2pe");
	mult_uniformFC_noMask->SetTitle("Ver3: Uniform FC, Cathode transparency 0%");
	mult_uniformFC_noMask->SetLineColor(kGreen);
	mult_uniformFC_noMask->Scale(1/mult_uniformFC_noMask->Integral());

	TFile* f_uniformFC_Mask = new TFile("rootResults/effi3_hor0_vert0_oppo0_maCentral1_inFile32.root","READ");
	TH1F* mult_uniformFC_Mask  = (TH1F*)f_uniformFC_Mask->Get("h_multiplicityAbove2pe");
	mult_uniformFC_Mask->SetTitle("Ver3: Uniform FC, Cathode transparency 0%, masked two central tiles");
	mult_uniformFC_Mask->SetLineColor(kGreen);
	mult_uniformFC_Mask->SetLineStyle(2);
	mult_uniformFC_Mask->Scale(1/mult_uniformFC_Mask->Integral());

	TFile* f_80cathTrans_noMask = new TFile("rootResults/effi3_hor0_vert0_oppo0_maCentral0_inFile42.root","READ");
	TH1F* mult_80cathTrans_noMask  = (TH1F*)f_80cathTrans_noMask->Get("h_multiplicityAbove2pe");
	mult_80cathTrans_noMask->SetTitle("Ver4: Real FC, Cathode transparency 80%");
	mult_80cathTrans_noMask->SetLineColor(kBlack);
	mult_80cathTrans_noMask->Scale(1/mult_80cathTrans_noMask->Integral());

	TFile* f_80cathTrans_Mask = new TFile("rootResults/effi3_hor0_vert0_oppo0_maCentral1_inFile42.root","READ");
	TH1F* mult_80cathTrans_Mask  = (TH1F*)f_80cathTrans_Mask->Get("h_multiplicityAbove2pe");
	mult_80cathTrans_Mask->SetTitle("Ver4: Real FC, Cathode transparency 80%, masked two central tiles");
	mult_80cathTrans_Mask->SetLineColor(kBlack);
	mult_80cathTrans_Mask->SetLineStyle(2);
	mult_80cathTrans_Mask->Scale(1/mult_80cathTrans_Mask->Integral());

	THStack* s_multComp = new THStack("s_multComp","; N_{tiles} [#]; NoE [#]");
	s_multComp->Add(mult_oldMC_noMask,"HIST");
	s_multComp->Add(mult_newMC_noMask,"HIST");
	s_multComp->Add(mult_uniformFC_noMask,"HIST");
	s_multComp->Add(mult_80cathTrans_noMask,"HIST");
	s_multComp->Add(mult_oldMC_Mask,"HIST");
	s_multComp->Add(mult_newMC_Mask,"HIST");
	s_multComp->Add(mult_uniformFC_Mask,"HIST");
	s_multComp->Add(mult_80cathTrans_Mask,"HIST");


	TCanvas* c_multComp = new TCanvas("c_multComp","MultiplicityComparison",800,600);
	gPad->SetGrid();
	s_multComp->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_detEffi5peSB->Write();


	return 0;
}