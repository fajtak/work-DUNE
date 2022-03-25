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

int compareNoiseDune()
{
	TString twoPE_1us_10mHz_w100m = "rootResults/noiseSim_nr0.01_spe0.90_tpe0.00_pd1.0_spew0.1.root";
	TString threePE_1us_10mHz_w100m = "rootResults/noiseSim_nr0.01_spe0.90_tpe0.01_pd1.0_spew0.1.root";
	TString threePE_10us_30mHz_w100m = "rootResults/noiseSim_nr0.03_spe0.90_tpe0.01_pd10.0_spew0.1.root";
	TString threePE_10us_10mHz_w100m = "rootResults/noiseSim_nr0.01_spe0.90_tpe0.01_pd10.0_spew0.1.root";


	TFile* f_twoPE_1us_10mHz = new TFile(twoPE_1us_10mHz_w100m,"READ");
	TH1F* h_nPEGen_twoPE_1us_10mHz  = (TH1F*)f_twoPE_1us_10mHz->Get("h_noiseChargeNonMerged");
	h_nPEGen_twoPE_1us_10mHz->SetTitle("90/10, 1us, 10mHz");
	h_nPEGen_twoPE_1us_10mHz->SetLineColor(kRed);
	h_nPEGen_twoPE_1us_10mHz->SetMarkerColor(kRed);
	TH1F* h_chargePerPE_twoPE_1us_10mHz  = (TH1F*)f_twoPE_1us_10mHz->Get("h_noiseChargePerPE");
	h_chargePerPE_twoPE_1us_10mHz->SetTitle("90/10, 1us, 10mHz");
	h_chargePerPE_twoPE_1us_10mHz->SetLineColor(kRed);
	h_chargePerPE_twoPE_1us_10mHz->SetMarkerColor(kRed);
	TH1F* h_charge_twoPE_1us_10mHz  = (TH1F*)f_twoPE_1us_10mHz->Get("h_noiseCharge");
	h_charge_twoPE_1us_10mHz->SetTitle("90/10, 1us, 10mHz");
	h_charge_twoPE_1us_10mHz->SetLineColor(kRed);
	h_charge_twoPE_1us_10mHz->SetMarkerColor(kRed);

	TFile* f_threePE_1us_10mHz = new TFile(threePE_1us_10mHz_w100m,"READ");
	TH1F* h_nPEGen_threePE_1us_10mHz  = (TH1F*)f_threePE_1us_10mHz->Get("h_noiseChargeNonMerged");
	h_nPEGen_threePE_1us_10mHz->SetTitle("90/9/1, 1us, 10mHz");
	h_nPEGen_threePE_1us_10mHz->SetLineColor(kBlue);
	h_nPEGen_threePE_1us_10mHz->SetMarkerColor(kBlue);
	TH1F* h_chargePerPE_threePE_1us_10mHz  = (TH1F*)f_threePE_1us_10mHz->Get("h_noiseChargePerPE");
	h_chargePerPE_threePE_1us_10mHz->SetTitle("90/9/1, 1us, 10mHz");
	h_chargePerPE_threePE_1us_10mHz->SetLineColor(kBlue);
	h_chargePerPE_threePE_1us_10mHz->SetMarkerColor(kBlue);
	TH1F* h_charge_threePE_1us_10mHz  = (TH1F*)f_threePE_1us_10mHz->Get("h_noiseCharge");
	h_charge_threePE_1us_10mHz->SetTitle("90/9/1, 1us, 10mHz");
	h_charge_threePE_1us_10mHz->SetLineColor(kBlue);
	h_charge_threePE_1us_10mHz->SetMarkerColor(kBlue);

	TFile* f_threePE_10us_30mHz = new TFile(threePE_10us_30mHz_w100m,"READ");
	TH1F* h_nPEGen_threePE_10us_30mHz  = (TH1F*)f_threePE_10us_30mHz->Get("h_noiseChargeNonMerged");
	h_nPEGen_threePE_10us_30mHz->SetTitle("90/9/1, 10us, 30mHz");
	h_nPEGen_threePE_10us_30mHz->SetLineColor(kGreen);
	h_nPEGen_threePE_10us_30mHz->SetMarkerColor(kGreen);
	TH1F* h_chargePerPE_threePE_10us_30mHz  = (TH1F*)f_threePE_10us_30mHz->Get("h_noiseChargePerPE");
	h_chargePerPE_threePE_10us_30mHz->SetTitle("90/9/1, 10us, 30mHz");
	h_chargePerPE_threePE_10us_30mHz->SetLineColor(kGreen);
	h_chargePerPE_threePE_10us_30mHz->SetMarkerColor(kGreen);
	TH1F* h_charge_threePE_10us_30mHz  = (TH1F*)f_threePE_10us_30mHz->Get("h_noiseCharge");
	h_charge_threePE_10us_30mHz->SetTitle("90/9/1, 10us, 30mHz");
	h_charge_threePE_10us_30mHz->SetLineColor(kGreen);
	h_charge_threePE_10us_30mHz->SetMarkerColor(kGreen);

	TFile* f_threePE_10us_10mHz = new TFile(threePE_10us_10mHz_w100m,"READ");
	TH1F* h_nPEGen_threePE_10us_10mHz  = (TH1F*)f_threePE_10us_10mHz->Get("h_noiseChargeNonMerged");
	h_nPEGen_threePE_10us_10mHz->SetTitle("90/9/1, 10us, 10mHz");
	h_nPEGen_threePE_10us_10mHz->SetLineColor(kViolet);
	h_nPEGen_threePE_10us_10mHz->SetMarkerColor(kViolet);
	TH1F* h_chargePerPE_threePE_10us_10mHz  = (TH1F*)f_threePE_10us_10mHz->Get("h_noiseChargePerPE");
	h_chargePerPE_threePE_10us_10mHz->SetTitle("90/9/1, 10us, 10mHz");
	h_chargePerPE_threePE_10us_10mHz->SetLineColor(kViolet);
	h_chargePerPE_threePE_10us_10mHz->SetMarkerColor(kViolet);
	TH1F* h_charge_threePE_10us_10mHz  = (TH1F*)f_threePE_10us_10mHz->Get("h_noiseCharge");
	h_charge_threePE_10us_10mHz->SetTitle("90/9/1, 10us, 10mHz");
	h_charge_threePE_10us_10mHz->SetLineColor(kViolet);
	h_charge_threePE_10us_10mHz->SetMarkerColor(kViolet);

	THStack* s_nPEGen_1us_10mHz = new THStack("s_nPEGen_1us_10mHz","; Q [p.e.]; f [Hz/1 p.e.]");
	s_nPEGen_1us_10mHz->Add(h_nPEGen_twoPE_1us_10mHz,"");
	s_nPEGen_1us_10mHz->Add(h_nPEGen_threePE_1us_10mHz,"");

	TCanvas* c_nPEGen_1us_10mHz = new TCanvas("c_nPEGen_1us_10mHz","NPEGen_1us_10mHz",800,600);
	gPad->SetGrid();
	s_nPEGen_1us_10mHz->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_nPEGen_1us_10mHz->Write();

	THStack* s_chargePerPE_1us_10mHz = new THStack("s_chargePerPE_1us_10mHz","; Q [p.e.]; f [Hz/1 p.e.]");
	s_chargePerPE_1us_10mHz->Add(h_chargePerPE_twoPE_1us_10mHz,"");
	s_chargePerPE_1us_10mHz->Add(h_chargePerPE_threePE_1us_10mHz,"");
	s_chargePerPE_1us_10mHz->Add(h_chargePerPE_threePE_10us_30mHz,"");
	s_chargePerPE_1us_10mHz->Add(h_chargePerPE_threePE_10us_10mHz,"");

	TCanvas* c_chargePerPE_1us_10mHz = new TCanvas("c_chargePerPE_1us_10mHz","ChargePerPE_1us_10mHz",800,600);
	gPad->SetGrid();
	s_chargePerPE_1us_10mHz->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_chargePerPE_1us_10mHz->Write();

	THStack* s_charge_1us_10mHz = new THStack("s_charge_1us_10mHz","; Q [p.e.]; f [Hz/0.1 p.e.]");
	s_charge_1us_10mHz->Add(h_charge_twoPE_1us_10mHz,"");
	s_charge_1us_10mHz->Add(h_charge_threePE_1us_10mHz,"");
	s_charge_1us_10mHz->Add(h_charge_threePE_10us_10mHz,"");

	TCanvas* c_charge_1us_10mHz = new TCanvas("c_charge_1us_10mHz","Charge_1us_10mHz",800,600);
	gPad->SetGrid();
	s_charge_1us_10mHz->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_charge_1us_10mHz->Write();

	return 0;
}