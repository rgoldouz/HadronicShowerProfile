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

void Tree_alphabeta(){

   gROOT->ProcessLine("#include <vector>");
   using namespace std;
	//energy values
//	Double_t energies[] = {1., 2., 3., 5., 10, 50., 100.,150,200};
        Double_t energies[] = {10.};
        const int N=1;

	std::vector<string> samples_;
	std::stringstream dname;
	for( int num = 0; num < N; num++ )
	{
        dname << "/afs/cern.ch/work/r/rgoldouz/TASK2015/CMSSW_7_4_12/src/Lukastree/hadshowertuningdata_compressed_V1/"<<energies[num];
	samples_.push_back((dname.str()).c_str());
	dname.str(std::string());
	}


   std::vector<TH1F*> hists;

   vector<float>  E;
   vector<float> A_EM_F;
   vector<float> B_EM_F;
   vector<float> A_EM_L;
   vector<float> B_EM_L;
   vector<float> A_EM_HAD;
   vector<float> B_EM_HAD; 


   TFile f("HadronicShower_parameters.root","recreate");
   TTree *T     = new TTree("HadronicShower","HadronicShower");

   T->Branch("E",&E);
   T->Branch("A_EM_F",&A_EM_F);
   T->Branch("B_EM_F",&B_EM_F);
   T->Branch("A_EM_L",&A_EM_L);
   T->Branch("B_EM_L",&B_EM_L);
   T->Branch("A_EM_HAD",&A_EM_HAD);
   T->Branch("B_EM_HAD",&B_EM_HAD);
	double alphaFmean[samples_.size()];
        double alphaLmean[samples_.size()];
        double alphaHmean[samples_.size()];
        double betaFmean[samples_.size()];
        double betaLmean[samples_.size()];
        double betaHmean[samples_.size()];

   for(unsigned int idx=0; idx<samples_.size(); ++idx){
	TSystemDirectory dir(samples_[idx].c_str(), samples_[idx].c_str());
	TList *files = dir.GetListOfFiles();
	TSystemFile *file;
	TString filename;
	TIter next(files);
        TH1F *alpha_H, *alpha_F, *alpha_L;
        TH1F *beta_H, *beta_F, *beta_L;
        TH2F *abcorr_H, *abcorr_F, *abcorr_L;

//        alpha_H = new TH1F("alpha_H","alpha_H" ,   100, 0, 15  );
//        alpha_F = new TH1F("alpha_F","alpha_F" ,   100, 0, 50  );
//        alpha_L = new TH1F("alpha_L","alpha_L" ,   100, 0, 15  );

//	beta_H = new TH1F("beta_H","beta_H" ,   100, 0, 15  );
//        beta_F = new TH1F("beta_F","beta_F" ,   100, 0, 15  );
//        beta_L = new TH1F("beta_L","beta_L" ,   100, 0, 15  );

	abcorr_H = new TH2F("abcorr_H","abcorr_H",200,-10,10,200,-10,10); 
	abcorr_F = new TH2F("abcorr_F","abcorr_F",100,0,-10,100,0,-10);
	abcorr_L = new TH2F("abcorr_L","abcorr_L",200,-10,10,200,-10,10);


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
//        if (h_H->Integral(firstnonzerobin,h_H->GetXaxis()->FindBin((*showerStart)[2]))/h_H->Integral() > 0.05) continue;

        for (int b = 1; b < nbins; ++b){
//	if (h_H->Integral()>0 && h_F->Integral()<h_H->Integral()/100 && h_L->Integral()<h_H->Integral()/100)
//        if ( h_H->GetXaxis()->GetBinCenter(b) < (*showerStart)[2]) continue;
	HADshapeIL->Fill((h_H->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/hcalintlen, h_H->GetBinContent(b));

        if(h_F->GetXaxis()->GetBinCenter(firstnonzerobin)-(*showerStart)[2] < 2*hcalradlen)  EMFshapeRL->Fill((h_F->GetXaxis()->GetBinCenter(b) - (*showerStart)[2] )/hcalradlen, h_F->GetBinContent(b));
        else EMLshapeIL->Fill((h_F->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/hcalradlen, h_F->GetBinContent(b));
        EMLshapeIL->Fill((h_L->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/hcalintlen, h_L->GetBinContent(b));
/*
	if ( h_H->GetXaxis()->GetBinCenter(b) < (*showerStart)[2]) continue;

	if ( h_H->GetXaxis()->GetBinCenter(b) < 152 ){
	HADshapeIL->Fill((h_H->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/ecalintlen, h_H->GetBinContent(b));
	if (h_F->GetXaxis()->GetBinCenter(firstnonzerobin)-(*showerStart)[2] < 2*hcalradlen) EMFshapeRL->Fill((h_F->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/ecalradlen, h_F->GetBinContent(b));
	else EMLshapeIL->Fill((h_F->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/ecalintlen, h_L->GetBinContent(b));
        EMLshapeIL->Fill((h_L->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/ecalintlen, h_L->GetBinContent(b));
	}

        if ( h_H->GetXaxis()->GetBinCenter(b) > 177.7 ){
        HADshapeIL->Fill((h_H->GetXaxis()->GetBinCenter(b) - (*showerStart)[2] - 25.7)/hcalintlen, h_H->GetBinContent(b));
	if(h_F->GetXaxis()->GetBinCenter(firstnonzerobin)-(*showerStart)[2] < 2*hcalradlen)  EMFshapeRL->Fill((h_F->GetXaxis()->GetBinCenter(b) - (*showerStart)[2] - 25.7)/hcalradlen, h_F->GetBinContent(b));
	else EMLshapeIL->Fill((h_F->GetXaxis()->GetBinCenter(b) - (*showerStart)[2] - 25.7)/hcalradlen, h_F->GetBinContent(b));
        EMLshapeIL->Fill((h_L->GetXaxis()->GetBinCenter(b) - (*showerStart)[2] - 25.7)/hcalintlen, h_L->GetBinContent(b));
        }
*/	
}


	double totalE;
	double firstEME;
	double lateEME;
	double HADE;
	firstEME = EMFshapeRL->Integral();
	lateEME = EMLshapeIL->Integral();
	HADE = HADshapeIL->Integral();
	totalE = firstEME + lateEME + HADE;

        bool jj=false;
        if(HADshapeIL->Integral()>0) jj=true;


        if (HADshapeIL->Integral()>0) HADshapeIL->Scale(1/HADshapeIL->Integral());
        if (EMFshapeRL->Integral()>0) EMFshapeRL->Scale(1/EMFshapeRL->Integral());
        if (EMLshapeIL->Integral()>0) EMLshapeIL->Scale(1/EMLshapeIL->Integral());

	TF1 *funcH;
	TF1 *funcF;
	TF1 *funcL;


	TCanvas* can;
	TPad* pad;
	TLegend* leg;
	TPaveText* pave_par;

	//get values from histos
	Double_t mH = HADshapeIL->GetMean();
        Double_t mF = EMFshapeRL->GetMean();
        Double_t mL = EMLshapeIL->GetMean();


	Double_t sH = pow(HADshapeIL->GetRMS(),2);
        Double_t sF = pow(EMFshapeRL->GetRMS(),2);
        Double_t sL = pow(EMLshapeIL->GetRMS(),2);

	Int_t NH = HADshapeIL->GetEntries();
        Int_t NF = EMFshapeRL->GetEntries();
        Int_t NL = EMLshapeIL->GetEntries();
	//setup fitting function & do fit
	//


	E.push_back(energies[idx]);
	if (EMFshapeRL->Integral()>0 && pow(mF,2)/sF <200) A_EM_F.push_back(log(pow(mF,2)/sF));
	if (EMFshapeRL->Integral()>0 && mF/sF<200) B_EM_F.push_back(log(mF/sF));
	if (EMLshapeIL->Integral()>0 && pow(mL,2)/sL <200) A_EM_L.push_back(log(pow(mL,2)/sL));
	if (EMLshapeIL->Integral()>0 && mL/sL <200) B_EM_L.push_back(log(mL/sL));
	if (HADshapeIL->Integral()>0 && pow(mH,2)/sH<200) A_EM_HAD.push_back(log(pow(mH,2)/sH));
	if (HADshapeIL->Integral()>0 && mH/sH <200) B_EM_HAD.push_back(log(mH/sH));

//       if (HADshapeIL->Integral()>0 && pow(mH,2)/sH<200 && mH/sH <200) abcorr_H->Fill(pow(mH,2)/sH,mH/sH);
	if (HADshapeIL->Integral()>0 && pow(mH,2)/sH<200 && mH/sH <200 && jj) abcorr_H->Fill(log(pow(mH,2)/sH),log(mH/sH));

}
}
	alphaFmean [idx] = std::accumulate(std::begin(A_EM_F), std::end(A_EM_F), 0.0) / A_EM_F.size(); 
	alphaLmean [idx] = std::accumulate(std::begin(A_EM_L), std::end(A_EM_L), 0.0) / A_EM_L.size();
	alphaHmean [idx] = std::accumulate(std::begin(A_EM_HAD), std::end(A_EM_HAD), 0.0) / A_EM_HAD.size();
	betaFmean [idx] = std::accumulate(std::begin(B_EM_F), std::end(B_EM_F), 0.0) / A_EM_F.size();
	betaLmean [idx] = std::accumulate(std::begin(B_EM_L), std::end(B_EM_L), 0.0) / A_EM_L.size();
	betaHmean [idx] = std::accumulate(std::begin(B_EM_HAD), std::end(B_EM_HAD), 0.0) / A_EM_HAD.size();
 
        T->Fill();
cout<<"***********alpha mean****************    "<<abcorr_H->GetMean(1)<<endl;
cout<<"***********beta mean****************    "<<abcorr_H->GetMean(2)<<endl;
cout<<"***********alpha var****************    "<<abcorr_H->GetRMS(1)<<endl;
cout<<"***********beta var****************    "<<abcorr_H->GetRMS(2)<<endl;
cout<<"***********correlation****************    "<<abcorr_H->GetCorrelationFactor()<<endl;


if (energies[idx]==3){
 TCanvas *c1 = new TCanvas("c1","Root Canvas",900,20,540,550); // My Default Canvas 
abcorr_H->Draw();
c1->Print("alphabetaEdep.png","png");
}
}

for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<"alphaFmean [idx]"<<alphaFmean [idx]<<endl;
cout<<"alphaLmean [idx]"<<alphaLmean [idx]<<endl;
cout<<"alphaHmean [idx]"<<alphaHmean [idx]<<endl;
cout<<"betaFmean [idx]"<<betaFmean [idx]<<endl;
cout<<"betaLmean [idx]"<<betaLmean [idx]<<endl;
cout<<"betaHmean [idx]"<<betaHmean [idx]<<endl;
}


   TCanvas* can;
   can = new TCanvas("ft","ft",700,500);
   can->cd();
   can->SetLogy(1);
   can->SetLogx(1);

   TMultiGraph *mg = new TMultiGraph();

   TGraph *AFmean = new TGraph(N,energies,alphaFmean);
   AFmean->SetLineColor(2);
   AFmean->SetLineWidth(2);
   AFmean->SetMarkerColor(2);
   AFmean->SetMarkerStyle(20);

   TGraph *ALmean = new TGraph(N,energies,alphaLmean);
   ALmean->SetLineColor(3);
   ALmean->SetLineWidth(2);
   ALmean->SetMarkerColor(3);
   ALmean->SetMarkerStyle(20);

   TGraph *AHADmean = new TGraph(N,energies,alphaHmean);
   AHADmean->SetLineColor(4);
   AHADmean->SetLineWidth(2);
   AHADmean->SetMarkerColor(4);
   AHADmean->SetMarkerStyle(20);


   TGraph *BFmean = new TGraph(N,energies,betaFmean);
   BFmean->SetLineColor(5);
   BFmean->SetLineWidth(2);
   BFmean->SetMarkerColor(5);
   BFmean->SetMarkerStyle(20);

   TGraph *BLmean = new TGraph(N,energies,betaLmean);
   BLmean->SetLineColor(6);
   BLmean->SetLineWidth(2);
   BLmean->SetMarkerColor(6);
   BLmean->SetMarkerStyle(20);

   TGraph *BHADmean = new TGraph(N,energies,betaHmean);
   BHADmean->SetLineColor(7);
   BHADmean->SetLineWidth(2);
   BHADmean->SetMarkerColor(7);
   BHADmean->SetMarkerStyle(20);


   mg->Add(AFmean);
   mg->Add(ALmean);
   mg->Add(AHADmean);
   mg->Add(BFmean);
   mg->Add(BLmean);
   mg->Add(BHADmean);
//   mg->SetMinimum(0);
//   mg->SetMaximum(1);
   mg->Draw("APC");
   mg->GetXaxis()->SetTitle("E_{inc} [GeV]");
   mg->GetYaxis()->SetTitle("");

  TLegend* leg = new TLegend(0.9,0.4,1,0.9);
  leg->SetFillStyle ( 0);
  leg->SetFillColor ( 0);
  leg->SetBorderSize( 0);
  leg->AddEntry( AFmean, " #alpha_{First #pi^{0}} "                           , "L");
  leg->AddEntry(ALmean , " #alpha_{Late #pi^{0}}"            , "L");
  leg->AddEntry( AHADmean, " #alpha_{Hadronic}", "L");
  leg->AddEntry(BFmean, " #beta_{First #pi^{0}}", "L");
  leg->AddEntry(BLmean, " #beta_{Late #pi^{0}}", "L");
  leg->AddEntry(BHADmean, " #beta_{Hadronic}", "L");
  leg->Draw("same");

 f.Write();
 f.Close();



} 
