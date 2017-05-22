/**********************************************************************************
 * In this code I want to fit the longitudinal shower shape 
 * to extract the alpha and beta parameters of gamma functions
 * this time using histograms prepeared by Lukas
 *
 * In this code I am going to save the alpha beta for different shower and energy in a root file.
 *   Reza Goldouzian, 11 Jan  2016
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TLine.h"
#include "TStopwatch.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include <cmath> 
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TPie.h"
#include "TF1.h"
#include "TVector3.h"
#include <TSystemDirectory.h>
#include <TSystemFile.h>



#define maxHDe 8 //energy points for hadrons
#define maxHDeta 1 //eta bins for hadrons
#define nGammaPar 2 //#pars for Gamma fn

using namespace TMath;

Double_t GammaProfile(Double_t *x,Double_t *par) {
double Gamma1=0;
Gamma1 = (par[1])*pow((x[0])*(par[1]),par[0]-1.0)*exp(-(x[0])*(par[1]))/tgamma(par[0]);
return Gamma1;
}



//function to fit Shower shape

void AverageHadronicProfile(){

   gROOT->ProcessLine("#include <vector>");
   using namespace std;
	//energy values
//	Double_t energies[] = {1., 2., 3., 5., 10, 50., 100.,150,200};
        Double_t energies[] = {10};
        const int N=1;

	std::vector<string> samples_;
	std::stringstream dname;
	for( int num = 0; num < 1; num++ )
	{
        dname << "/afs/cern.ch/work/r/rgoldouz/TASK2015/CMSSW_7_4_12/src/Lukastree/hadshowertuningdata_compressed_V1/"<<energies[num];
	samples_.push_back((dname.str()).c_str());
	dname.str(std::string());
	}


   std::vector<TH1F*> histsreal;

   for(unsigned int idx=0; idx<samples_.size(); ++idx){
	TSystemDirectory dir(samples_[idx].c_str(), samples_[idx].c_str());
	TList *files = dir.GetListOfFiles();
	TSystemFile *file;
	TString filename;
	TIter next(files);
	while ((file=(TSystemFile*)next())) {
        filename = file->GetName();
       if(file->IsDirectory()) continue;
	cout<<filename.Data()<<endl;
        std::stringstream name;
	name<<samples_[idx].c_str()<<"/"<<filename;

	int nbins = 320;
        double ecalradlen =0.89;
        double ecalintlen =19.5;
        double hcalradlen =1.49;
        double hcalintlen =16.42;
        TH1F *EMFshapeRL, *EMLshapeIL, *HADshapeIL;
        TH1F *h_H, *h_F, *h_L;

	//open file and tree
	TFile *input(0);	
	input = TFile::Open((name.str()).c_str());
	for (int Nshower = 0; Nshower < 500; ++Nshower){
        std::stringstream snameH;
        std::stringstream snameF;
        std::stringstream snameL;

	snameH<<"shower_"<<Nshower<<"/longEProf_fineBin_had";
        snameF<<"shower_"<<Nshower<<"/longEProf_fineBin_pi0_1";
        snameL<<"shower_"<<Nshower<<"/longEProf_fineBin_pi0_2";

	h_H = new TH1F("","" , nbins,120,280);
	h_H = (TH1F*)input->Get((snameH.str()).c_str());

        h_F = new TH1F("", "" , nbins,120,280);
        h_F = (TH1F*)input->Get((snameF.str()).c_str());

        h_L = new TH1F( "" , "" , nbins,120,280);
        h_L = (TH1F*)input->Get((snameL.str()).c_str());

        std::stringstream SS;
	SS<<"shower_"<<Nshower<<"/showerStart";
	TVector3 * showerStart = (TVector3*)input->Get((SS.str()).c_str());

	EMFshapeRL = new TH1F( (snameF.str()).c_str(),(snameF.str()).c_str() ,   50, 0, 50  );
	EMLshapeIL = new TH1F( (snameL.str()).c_str(),(snameL.str()).c_str() ,   20, 0  , 10  );
	HADshapeIL = new TH1F( (snameH.str()).c_str(),(snameH.str()).c_str() ,   50, 0  , 10    );


	int firstnonzerobin;
	firstnonzerobin = h_H->FindFirstBinAbove(0,1);

	if ((*showerStart)[2] < 177.7) continue;

        for (int b = 1; b < nbins; ++b){

//        if ( h_H->GetXaxis()->GetBinCenter(b) < (*showerStart)[2]) continue;

//	if (h_H->Integral()>0 && h_F->Integral()<h_H->Integral()/100 && h_L->Integral()<h_H->Integral()/100)
//	HADshapeIL->Fill((h_H->GetXaxis()->GetBinCenter(b) - h_H->GetXaxis()->GetBinCenter(firstnonzerobin))/hcalintlen, h_H->GetBinContent(b));
//        HADshapeIL->Fill((h_H->GetXaxis()->GetBinCenter(b) - firstnonzerobin)/hcalintlen, h_H->GetBinContent(b));

      HADshapeIL->Fill((h_H->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/hcalintlen, h_H->GetBinContent(b));

        if(h_F->GetXaxis()->GetBinCenter(firstnonzerobin)-(*showerStart)[2] < 2*hcalradlen)  EMFshapeRL->Fill((h_F->GetXaxis()->GetBinCenter(b) - (*showerStart)[2] )/hcalradlen, h_F->GetBinContent(b));
        else EMLshapeIL->Fill((h_F->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/hcalradlen, h_F->GetBinContent(b));
        EMLshapeIL->Fill((h_L->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/hcalintlen, h_L->GetBinContent(b));
}


	double totalE;
	double firstEME;
	double lateEME;
	double HADE;
	firstEME = EMFshapeRL->Integral();
	lateEME = EMLshapeIL->Integral();
	HADE = HADshapeIL->Integral();
	totalE = firstEME + lateEME + HADE;

        Double_t mH = HADshapeIL->GetMean();
        Double_t sH = pow(HADshapeIL->GetRMS(),2);

//if (energies[idx]==100 && h_H->Integral(firstnonzerobin,h_H->GetXaxis()->FindBin((*showerStart)[2]))/h_H->Integral()<0.05) histsreal.push_back(HADshapeIL);
//        if (HADshapeIL->Integral()>0) HADshapeIL->Scale(1/HADshapeIL->Integral());
//        if (EMFshapeRL->Integral()>0) EMFshapeRL->Scale(1/EMFshapeRL->Integral());
//        if (EMLshapeIL->Integral()>0) EMLshapeIL->Scale(1/EMLshapeIL->Integral());
        if (HADshapeIL->Integral()==0 || pow(mH,2)/sH>200 || mH/sH >200) continue;
//if (energies[idx]==100 && h_H->Integral(firstnonzerobin,h_H->GetXaxis()->FindBin((*showerStart)[2]))/h_H->Integral()<0.05) 
histsreal.push_back(HADshapeIL);
}
}
}

for (int i = 0; i < histsreal.size(); ++i){
histsreal[i]->Scale(1/histsreal[i]->Integral());
}

for (int i = 1; i < histsreal.size(); ++i){
histsreal[i]->Add(histsreal[i-1]);
}

   TCanvas* can;
   can = new TCanvas("ft","ft",700,500);
//   can->SetLogy(1);
   can->cd();
histsreal[histsreal.size()-1]->Scale(1/histsreal[histsreal.size()-1]->Integral());




/*
3 GeV
double meanalpha = -0.165;
double meanbeta = 0.119;
double varalpha = 0.8;
double varbeta = 0.32;
double corr = 0.848;

100 GeV
double meanalpha = 0.854;
double meanbeta = 0.460;
double varalpha = 0.600;
double varbeta = 0.266;
double corr = 0.745;
*/

double meanalpha =  0.390173;
double meanbeta = 0.434341;
double varalpha = 0.553233;
double varbeta = 0.50433;
double corr = 0.557659;

double x = sqrt((1+corr)/2);
double y = sqrt((1-corr)/2);


double alpha;
double beta;
int n=1000;

gRandom = new TRandom3();
TH1F *Eprofile;
std::vector<TH1F*> hists;

int nsudo=1000;
for (int i = 0; i < nsudo; ++i){
double z1=gRandom->Gaus(0, 1);
double z2=gRandom->Gaus(0, 1);

double rand;

alpha = meanalpha + varalpha * x * z1 + varalpha *y *z2;
beta = meanbeta + varbeta * x * z1 - varbeta *y *z2;


Eprofile = new TH1F("Eprofile","Eprofile",50,0,10);
funcH = new TF1("GammaH",GammaProfile,0  , 100,nGammaPar);
funcH->SetParameters(exp(alpha) ,exp(beta));
for (int j = 0; j < n; ++j){
rand = funcH->GetRandom();
Eprofile->Fill(rand);
}
Eprofile->Scale(1/Eprofile->Integral());
hists.push_back(Eprofile);
}

for (int i = 1; i < nsudo; ++i){
hists[i]->Add(hists[i-1]);
}

hists[nsudo-1]->Scale(1/hists[nsudo-1]->Integral());
hists[nsudo-1]->SetLineColor(2);

hists[nsudo-1]->Draw();
histsreal[histsreal.size()-1]->Draw("same");


} 
