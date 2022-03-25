#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>

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

const double gArapucaEfficiency = 0.035;

const bool gReadAnodeTiles = false;
const bool gReadCathodeTiles = false;
const bool gReadCryoTiles = true;

const bool gMaskCentralTiles = true;
const bool gMaskMiddleTiles = false;
const bool gMaskNegativeTiles = false;

const bool gUseCathodeTileTransparency = true;
TRandom3 r_randNumGen3(0);

const double gBckgRate = 1.0*1e7;
const double gAcquisitionWindow = 1e-3;

TH2F* h_negativeTiles = new TH2F("h_negativeTiles","Distribution of hits per Event on the negative cryostat wall; columnID [#]; rowID [#]",20,0,20,20,0,20);
TH2F* h_positiveTiles = new TH2F("h_positiveTiles","Distribution of hits per Event on the positive cryostat wall; columnID [#]; rowID [#]",20,0,20,20,0,20);

TH2F* h_negativeTilesSignalWindow = new TH2F("h_negativeTilesSignalWindow","Distribution of hits per Event on the negative cryostat wall; columnID [#]; rowID [#]",20,0,20,20,0,20);
TH2F* h_positiveTilesSignalWindow = new TH2F("h_positiveTilesSignalWindow","Distribution of hits per Event on the positive cryostat wall; columnID [#]; rowID [#]",20,0,20,20,0,20);

TH2F* h_negativeTilesSignalOnly = new TH2F("h_negativeTilesSignalOnly","Distribution of hits per Event on the negative cryostat wall; columnID [#]; rowID [#]",20,0,20,20,0,20);
TH2F* h_positiveTilesSignalOnly = new TH2F("h_positiveTilesSignalOnly","Distribution of hits per Event on the positive cryostat wall; columnID [#]; rowID [#]",20,0,20,20,0,20);

TH1F* h_timelineQ = new TH1F("h_timelineQ","Timeline Q; T [#mu s]; Q [p.e.]",10000,-gAcquisitionWindow/2*1e6,gAcquisitionWindow/2*1e6);
TH1F* h_timelineNpe = new TH1F("h_timelineNpe","Timeline Q; T [#mu s]; N_{hits} > 2 p.e. [#]",10000,-gAcquisitionWindow/2*1e6,gAcquisitionWindow/2*1e6);

TH2F* h_timeVsClusterSize = new TH2F("h_timeVsClusterSize","Time vs cluster size; Time [#mu s]; Cluster size [#]; NoE [#]",10000,-gAcquisitionWindow/2*1e6,gAcquisitionWindow/2*1e6,40,0,40);

TH1F* h_deltaT = new TH1F("h_deltaT","DeltaT; #DeltaT [#mu s]; NoE [#]",10000,-gAcquisitionWindow/2*1e6,gAcquisitionWindow/2*1e6);
TH1F* h_deltaTCluster = new TH1F("h_deltaTCluster","DeltaT with clusters; #DeltaT [#mu s]; NoE [#]",10000,-gAcquisitionWindow/2*1e6,gAcquisitionWindow/2*1e6);
TH1F* h_maxNTiles = new TH1F("h_maxNTiles","Max N Tiles; max(N_{tiles}) [#]; NoE [#]",100,0,100);

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
	double time;
};

int DrawResults()
{
	TCanvas* c_eventSnapshotSignalOnly = new TCanvas("c_eventSnapshotSignalOnly","EventSnapshotSignalOnly",800,600);
	c_eventSnapshotSignalOnly->Divide(1,2);
	c_eventSnapshotSignalOnly->cd(1);
	gPad->SetRightMargin(0.1);
	h_negativeTilesSignalOnly->Draw("colz");
	h_negativeTilesSignalOnly->SetStats(false);
	c_eventSnapshotSignalOnly->cd(2);
	gPad->SetRightMargin(0.1);
	h_positiveTilesSignalOnly->Draw("colz");
	h_positiveTilesSignalOnly->SetStats(false);

	TCanvas* c_eventSnapshotSignalWindow = new TCanvas("c_eventSnapshotSignalWindow","EventSnapshotSignalWindow",800,600);
	c_eventSnapshotSignalWindow->Divide(1,2);
	c_eventSnapshotSignalWindow->cd(1);
	gPad->SetRightMargin(0.1);
	h_negativeTilesSignalWindow->Draw("colz");
	h_negativeTilesSignalWindow->SetStats(false);
	c_eventSnapshotSignalWindow->cd(2);
	gPad->SetRightMargin(0.1);
	h_positiveTilesSignalWindow->Draw("colz");
	h_positiveTilesSignalWindow->SetStats(false);

	TCanvas* c_eventSnapshot = new TCanvas("c_eventSnapshot","EventSnapshot",800,600);
	c_eventSnapshot->Divide(1,2);
	c_eventSnapshot->cd(1);
	gPad->SetRightMargin(0.1);
	gPad->SetRightMargin(0.18);
	h_negativeTiles->Draw("colz");
	h_negativeTiles->SetStats(false);
	c_eventSnapshot->cd(2);
	gPad->SetRightMargin(0.1);
	gPad->SetRightMargin(0.18);
	h_positiveTiles->Draw("colz");
	h_positiveTiles->SetStats(false);

	TCanvas* c_timelineQ = new TCanvas("c_timelineQ","TimelineQ",800,600);
	h_timelineQ->Draw("HIST");

	TCanvas* c_timelineNpe = new TCanvas("c_timelineNpe","TimelineNpe",800,600);
	h_timelineNpe->Draw("HIST");

	TCanvas* c_deltaT = new TCanvas("c_deltaT","DeltaT",800,600);
	h_deltaT->Draw("HIST");

	TCanvas* c_maxNTiles = new TCanvas("c_maxNTiles","MaxNTiles",800,600);
	h_maxNTiles->Draw("HIST");

	TCanvas* c_deltaTCluster = new TCanvas("c_deltaTCluster","DeltaTCluster",800,600);
	h_deltaTCluster->Draw("HIST");

	TCanvas* c_timeVsClusterSize = new TCanvas("c_timeVsClusterSize","TimeVsClusterSize",800,600);
	h_timeVsClusterSize->Draw("COLZ");

	return 0;
}

int DrawEventSnapshotSignalWindow()
{
	TCanvas* c_eventSnapshotSignalWindow = new TCanvas("c_eventSnapshotSignalWindow","EventSnapshotSignalWindow",800,600);
	c_eventSnapshotSignalWindow->Divide(1,2);
	c_eventSnapshotSignalWindow->cd();
	c_eventSnapshotSignalWindow->cd(1);
	gPad->SetRightMargin(0.1);
	h_negativeTilesSignalWindow->Draw("colz");
	gPad->Modified();
	h_negativeTilesSignalWindow->SetStats(false);
	c_eventSnapshotSignalWindow->cd(2);
	gPad->SetRightMargin(0.1);
	h_positiveTilesSignalWindow->Draw("colz");
	gPad->Modified();
	h_positiveTilesSignalWindow->SetStats(false);
	// c_eventSnapshotSignalWindow->Modified();
	c_eventSnapshotSignalWindow->Update();

	TCanvas* c_timelineNpe = new TCanvas("c_timelineNpe","TimelineNpe",800,600);
	c_timelineNpe->cd();
	h_timelineNpe->Draw("HIST");
	gPad->Modified();
	gPad->Update();
	// c_timelineNpe->Modified();
	c_timelineNpe->Update();

	int dummy;
	cin >> dummy;

	delete c_eventSnapshotSignalWindow;
	delete c_timelineNpe;

	return 0;
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

int ReadData(TString filePath, vector<Event> &events,unsigned int nReadEvents = std::numeric_limits<unsigned int>::max())
{
	ifstream inputFile;
    inputFile.open(filePath);

    if (!inputFile)
    {
    	cerr << "Input file: " << filePath << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

  	std::string line;
  	int eventID,nHits;
  	float energy,x,y,z;
  	unsigned short int cellID;
  	unsigned int nPE;

  	cout << events.size() << " " << nReadEvents << endl;

    while(std::getline(inputFile, line) && nReadEvents > events.size())
    // while(std::getline(inputFile, line))
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
        events.push_back(Event{eventID,energy*1000,x,y,z,nHits,{},0,0,0,0});
        // std::cout << eventID << "\t" << energy << "\t" << x << "\t" << y << "\t" << z << "\t" << nHits << std::endl;
        int nHitsWithEfficiency = 0;
        unsigned int nPEGeneratedTotal = 0;
        unsigned int nPEDetectedTotal = 0;
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

    		if (nPEDetected > 0)
    		{
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
				int row = (cellID-1)%20;
				if (gMaskCentralTiles && (row == 9 || row == 10))
					continue;
				if (gMaskMiddleTiles && (row > 3 && row < 16))
					continue;
				if (gMaskNegativeTiles && (row > 9))
					continue;
    			if (nPEDetected < 2)
    				continue;
    			nHitsWithEfficiency++;
	        	events.back().v_hits.push_back(Hit{cellID,nPEDetected});
    			events.back().nPEWalls += nPEDetected;
    		}
        }
        if (nHits % 10 == 0 && nHits != 0)
	        std::getline(inputFile, line);
        std::getline(inputFile, line);
    	events.back().nHits = nHitsWithEfficiency;
        // std::sort(events.back().v_hits,OrderByPE);
        // PrintEvent(events.back());
    }
	return 0;
}

int MirrorX(Event &event)
{
	for (unsigned int i = 0; i < event.v_hits.size(); ++i)
	{
		event.v_hits[i].cellID = event.v_hits[i].cellID < 401? event.v_hits[i].cellID+400:event.v_hits[i].cellID - 400;
	}
	return 0;
}

int MirrorY(Event &event)
{
	for (unsigned int i = 0; i < event.v_hits.size(); ++i)
	{
		int side = (event.v_hits[i].cellID-1)/400;
		int row = ((event.v_hits[i].cellID-1)%400)%20;
		int column = ((event.v_hits[i].cellID-1)%400)/20;

		event.v_hits[i].cellID = 400*side+column*20+(19-row)+1;
	}
	return 0;
}

int MirrorZ(Event &event)
{
	// cout << "Mirroring Z" << endl;
	for (unsigned int i = 0; i < event.v_hits.size(); ++i)
	{
		int side = (event.v_hits[i].cellID-1)/400;
		int row = ((event.v_hits[i].cellID-1)%400)%20;
		int column = ((event.v_hits[i].cellID-1)%400)/20;

		// cout << event.v_hits[i].cellID << " " << side << " " << row << " " << column << " " << 400*side+row+(19-column)*20+1 << endl;
		event.v_hits[i].cellID = 400*side+row+(19-column)*20+1;

	}
	return 0;
}

int MirrorBackgroundEvent(Event &backgroundEvent)
{
	double mirrorX = r_randNumGen3.Uniform(1);
	double mirrorY = r_randNumGen3.Uniform(1);
	double mirrorZ = r_randNumGen3.Uniform(1);

	if (mirrorX > 0.5)
		MirrorX(backgroundEvent);
	if (mirrorY > 0.5)
		MirrorY(backgroundEvent);
	if (mirrorZ > 0.5)
		MirrorZ(backgroundEvent);

	return 0;
}

int GenerateEvent(vector<Event> &v_signalEvents, vector<Event> &v_bckgEvents, vector<Event> &v_genEvent)
{
	v_genEvent.clear();
	int signalEventID = r_randNumGen3.Integer(v_signalEvents.size());
	v_genEvent.push_back(v_signalEvents[signalEventID]);
	v_genEvent.back().time = 0.00000001;

	int nBckgEvents = gBckgRate*gAcquisitionWindow;
	// cout << gBckgRate << " " << gAcquisitionWindow << " " << nBckgEvents << endl;

	for (int i = 0; i < nBckgEvents; ++i)
	{
		int bckgEventID = r_randNumGen3.Integer(v_bckgEvents.size());
		v_genEvent.push_back(v_bckgEvents[bckgEventID]);
		// v_genEvent.push_back(v_bckgEvents[i]);
		double bckgEventTime = r_randNumGen3.Uniform(gAcquisitionWindow)-gAcquisitionWindow/2;
		v_genEvent.back().time = bckgEventTime;
		MirrorBackgroundEvent(v_genEvent.back());
	}

	return 0;
}

int CreateSnapshot(vector<Event> v_genEvent)
{
	h_negativeTiles->Reset();
	h_positiveTiles->Reset();
	for (unsigned int i = 0; i < v_genEvent.size(); ++i)
	{
		for (unsigned int j = 0; j < v_genEvent[i].v_hits.size(); ++j)
		{
			if (v_genEvent[i].v_hits[j].cellID < 401)
			{
				if (v_genEvent[i].v_hits[j].nPE > 1)
					h_negativeTiles->Fill((v_genEvent[i].v_hits[j].cellID-1)/20,(v_genEvent[i].v_hits[j].cellID-1)%20,v_genEvent[i].v_hits[j].nPE);
			}
			else
			{
				if (v_genEvent[i].v_hits[j].nPE > 1)
					h_positiveTiles->Fill((v_genEvent[i].v_hits[j].cellID-401)/20,(v_genEvent[i].v_hits[j].cellID-401)%20,v_genEvent[i].v_hits[j].nPE);
			}
		}
	}
	return 0;
}

int CreateSnapshotSignalOnly(vector<Event> v_genEvent)
{
	h_negativeTilesSignalOnly->Reset();
	h_positiveTilesSignalOnly->Reset();
	for (unsigned int i = 0; i < 1; ++i)
	{
		for (unsigned int j = 0; j < v_genEvent[i].v_hits.size(); ++j)
		{
			if (v_genEvent[i].v_hits[j].cellID < 401)
				h_negativeTilesSignalOnly->Fill((v_genEvent[i].v_hits[j].cellID-1)/20,(v_genEvent[i].v_hits[j].cellID-1)%20,v_genEvent[i].v_hits[j].nPE);
			else
				h_positiveTilesSignalOnly->Fill((v_genEvent[i].v_hits[j].cellID-401)/20,(v_genEvent[i].v_hits[j].cellID-401)%20,v_genEvent[i].v_hits[j].nPE);
		}
	}
	return 0;
}

int CreateSnapshotSignalWindow(vector<Event> v_genEvent, double t0 = 0, double Twin = 0.2)
{
	h_negativeTilesSignalWindow->Reset();
	h_positiveTilesSignalWindow->Reset();
	for (unsigned int i = 0; i < v_genEvent.size(); ++i)
	{
		if (v_genEvent[i].time*1e6 - t0 <= -Twin || v_genEvent[i].time*1e6 - t0 >= Twin)
				continue;
		for (unsigned int j = 0; j < v_genEvent[i].v_hits.size(); ++j)
		{
			if (v_genEvent[i].v_hits[j].cellID < 401)
			{
				if (v_genEvent[i].v_hits[j].nPE > 1)
				{
					// cout << v_genEvent[i].v_hits[j].cellID << endl;
					h_negativeTilesSignalWindow->Fill((v_genEvent[i].v_hits[j].cellID-1)/20,(v_genEvent[i].v_hits[j].cellID-1)%20,v_genEvent[i].v_hits[j].nPE);
				}
			}
			else
			{
				if (v_genEvent[i].v_hits[j].nPE > 1)
				{
					// cout << v_genEvent[i].v_hits[j].cellID << endl;
					h_positiveTilesSignalWindow->Fill((v_genEvent[i].v_hits[j].cellID-401)/20,(v_genEvent[i].v_hits[j].cellID-401)%20,v_genEvent[i].v_hits[j].nPE);
				}
			}
		}
	}
	return 0;
}

struct XYPos
{
	int x;
	int y;
	double nPE;
};

bool NotInClusterAlready(vector<XYPos> &cluster, int x, int y)
{
	bool notInCluster = true;
	for (unsigned int i = 0; i < cluster.size(); ++i)
	{
		if (cluster[i].x == x && cluster[i].y == y)
			notInCluster = false;
	}
	return notInCluster;
}

int AddNewPixels(TH2F* hist2D,vector<XYPos> &cluster,int x, int y)
{
	int nNewPixels = 0;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			if (i == 1 && j == 1)
				continue;
			if (hist2D->GetBinContent(x-1+j,y-1+i) != 0 && NotInClusterAlready(cluster,x-1+j,y-1+i))
			{
				cluster.push_back(XYPos{x-1+j,y-1+i,hist2D->GetBinContent(x-1+j,y-1+i)});
				nNewPixels++;
			}
		}
	}
	return nNewPixels;
}

vector<XYPos> GetCluster(TH2F* hist2D,int x, int y)
{
	unsigned int vectorIDPointer = 0;
	vector<XYPos> cluster;
	cluster.push_back(XYPos{x,y,hist2D->GetBinContent(x,y)});

	while(vectorIDPointer != cluster.size())
	{
		AddNewPixels(hist2D,cluster,cluster[vectorIDPointer].x,cluster[vectorIDPointer].y);
		vectorIDPointer++;
	}

	return cluster;
}

double GetClusterZPosMean(vector<XYPos> cluster)
{
	double meanZPos = 0;
	for (unsigned int i = 0; i < cluster.size(); ++i)
	{
		meanZPos += -2850 +(cluster[i].x-1)*300;
	}
	return meanZPos /= cluster.size();
}

double GetClusterZPosMax(vector<XYPos> cluster)
{
	double maxValue = 0;
	double zPos;
	for (unsigned int i = 0; i < cluster.size(); ++i)
	{
		if (cluster[i].nPE > maxValue)
		{
			maxValue = cluster[i].nPE;
			zPos = -2850 +(cluster[i].x-1)*300;
		}
	}
	return zPos;
}

int FindLargestCluster(TH2F* hist2D,double &zPos)
{
	vector<XYPos> largestCluster;
	int nXBins = hist2D->GetXaxis()->GetNbins();
	int nYBins = hist2D->GetYaxis()->GetNbins();

	for (int i = 1; i <= nYBins; ++i)
	{
		for (int j = 1; j <= nXBins; ++j)
		{
			if (hist2D->GetBinContent(j,i) != 0)
			{
				if (!NotInClusterAlready(largestCluster,j,i))
					continue;
				vector<XYPos> cluster = GetCluster(hist2D,j,i);
				if (cluster.size() > largestCluster.size())
				{
					largestCluster = cluster;
					// cout << cluster.size();
				}
			}
		}
	}
	if (largestCluster.size() == 0)
		zPos = -99999;
	else
		zPos = GetClusterZPosMean(largestCluster);
	return largestCluster.size();
}

double CreateTimeLine(vector<Event> &v_genEvent)
{
	h_timelineQ->Reset();
	h_timelineNpe->Reset();
	for (unsigned int i = 0; i < v_genEvent.size(); ++i)
	{
		// ;
		h_timelineQ->Fill(v_genEvent[i].time*1000000,v_genEvent[i].nPEWalls);
		for (unsigned int j = 0; j < v_genEvent[i].v_hits.size(); ++j)
		{
			// if (v_genEvent[i].v_hits[j].nPE )
				h_timelineNpe->Fill(v_genEvent[i].time*1000000);
		}
	}
	h_maxNTiles->Fill(h_timelineNpe->GetMaximum());
	double T0 = h_timelineNpe->GetBinCenter(h_timelineNpe->GetMaximumBin());

	// if (TMath::Abs(T0) > 0.2)
	// {
	// 	cout << "T0: " << T0 << endl;
	// 	CreateSnapshotSignalWindow(v_genEvent,T0,0.05);
	// 	cout << FindLargestCluster(h_positiveTilesSignalWindow) << "\t" << FindLargestCluster(h_negativeTilesSignalWindow) << endl;
	// 	DrawEventSnapshotSignalWindow();
	// }

	return T0;
	// return 0;
}

int StudyTimeLine(vector<Event> &v_genEvent)
{
	int binWithLargestCluster = 0;
	int largestClusterInTWin = 0;
	int nBins = h_timelineNpe->GetXaxis()->GetNbins();
	for (int i = 1; i <= nBins; ++i)
	{
		if (h_timelineNpe->GetBinContent(i) > 13)
		{
			CreateSnapshotSignalWindow(v_genEvent,h_timelineNpe->GetBinCenter(i),0.05);
			double zPosNeg = 0;
			double zPosPos = 0;
			int largestClusterPos = FindLargestCluster(h_positiveTilesSignalWindow,zPosPos);
			int largestClusterNeg = FindLargestCluster(h_negativeTilesSignalWindow,zPosNeg);
			int largestCluster =  largestClusterPos > largestClusterNeg ? largestClusterPos : largestClusterNeg;
			double zPos =  largestClusterPos > largestClusterNeg ? zPosPos : zPosNeg;
			// cout << h_timelineNpe->GetBinCenter(i) << " " << v_genEvent[0].z << " " << largestClusterPos << " " << largestClusterNeg << " " << largestCluster << " " << zPosPos << " " << zPosNeg << " " << zPos << endl;
			// if (largestCluster > largestClusterInTWin && TMath::Abs(v_genEvent[0].z-zPos) < 300)
			if (largestCluster > largestClusterInTWin && TMath::Abs(v_genEvent[0].z-zPos) < 300)
			{
				largestClusterInTWin = largestCluster;
				binWithLargestCluster = i;
			}
			h_timeVsClusterSize->Fill(h_timelineNpe->GetBinCenter(i),largestCluster);
		}
	}
	h_deltaTCluster->Fill(h_timelineNpe->GetBinCenter(binWithLargestCluster));
	return 0;
}

int generateEvent(void)
{
	gStyle->SetPalette(55);

	vector<Event> v_signalEvents;
	ReadData("/Data/duneData/dune_fieldcage_10MeValla_tiles.data",v_signalEvents,10000);

	std::cout << "Number of signal events read: " << v_signalEvents.size() << std::endl;

	vector<Event> v_bckgEvents;
	ReadData("/Data/duneData/dune_fieldcage_Ar39_fulla_tiles.data",v_bckgEvents);

	std::cout << "Number of bckg events read: " << v_bckgEvents.size() << std::endl;

	vector<Event> v_genEvent;

	for (int i = 0; i < 1000; ++i)
	{
		if (i%10 == 0)
			cout << i << endl;
		GenerateEvent(v_signalEvents,v_bckgEvents,v_genEvent);
		h_deltaT->Fill(CreateTimeLine(v_genEvent));
		StudyTimeLine(v_genEvent);
	}
	// CreateSnapshot(v_genEvent);
	CreateSnapshotSignalOnly(v_genEvent);
	CreateSnapshotSignalWindow(v_genEvent);

	DrawResults();

	return 0;
}