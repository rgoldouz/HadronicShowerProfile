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
#include "TRandom3.h"


#define maxHDe 8 //energy points for hadrons
#define maxHDeta 1 //eta bins for hadrons
#define nGammaPar 2 //#pars for Gamma fn

using namespace TMath;

Double_t GammaProfile(Double_t *x,Double_t *par) {
double Gamma1=0;
Gamma1 = (par[1])*pow((x[0])*(par[1]),par[0]-1.0)*exp(-(x[0])*(par[1]))/tgamma(par[0]);
return Gamma1;
}




//function to generate gamma distribution randomly
TH1F* RandomGammaProfile(double alpha, double beta, double range, int nbin, const char* name){
  int n=1000;
  double rand;
  TH1F *Eprofile;
  Float_t bins[] = {0,0.05,0.15,0.35,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.8,3.3,3.9,4.6,5.3,8};
  Eprofile = new TH1F(name,name,nbin,bins);
  funcH = new TF1("GammaH", GammaProfile, 0, range, nGammaPar);
  funcH->SetParameters(alpha ,beta);

  for (int j = 0; j < n; ++j){
    rand = funcH->GetRandom();
    Eprofile->Fill(rand);
  }
  Eprofile->Scale(1/Eprofile->Integral());
  return Eprofile;
}


//function to fit Shower shape

void Average_HadProfile(){
  gROOT->ProcessLine("#include <vector>");
  using namespace std;
  int nbins = 320;
  double ecalradlen =0.89;
  double ecalintlen =19.5;
  double hcalradlen =1.49;
  double hcalintlen =16.42;

  Float_t bins[] = {0,0.05,0.15,0.35,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.8,3.3,3.9,4.6,5.3,8};
  TF1 *funcH;
  TH1F *A, *B, *As, *Bs;
  TH1F *HADshapeIL;
  TH1F *h_H;
  TH2F *abcorr_H, *abcorrRand_H, *cABr, *cABg;
  TCanvas* can;
  TPad* pad;
  TH1F *Eprofile;
  gRandom = new TRandom3();
  double rand;
  A = new TH1F( "A", "A" ,   50, -5  , 5    );
  As = new TH1F( "As", "As" ,   50, -5  , 5    );
  B = new TH1F( "B", "B" ,   50, -5  , 5   );
  Bs = new TH1F( "Bs", "Bs" ,50,   -5  , 5    );
//energy values
  Double_t energies[] = {1., 2., 3., 5., 10, 50., 100.,150,200};
//  Double_t energies[] = {10, 50., 100.,150,200};
//  Double_t energies[] = {50};
  const int N=9;
  std::vector<string> samples_;
  std::stringstream dname;


//Add samples
  for( int num = 0; num < 9; num++ ){
    dname << "hadshowertuningdata_compressed_V3/"<<energies[num];
    samples_.push_back((dname.str()).c_str());
    dname.str(std::string());
  }

//Loop over samples
  for(unsigned int idx=0; idx<samples_.size(); ++idx){
 
    cout<<"INput pion energy= "<<energies[idx]<<endl;
    abcorr_H = new TH2F("abcorr_H","abcorr_H",200,-10,10,200,-10,10);
    abcorrRand_H = new TH2F("abcorrRand_H","abcorrRand_H",200,-10,10,200,-10,10);
    cABr = new TH2F("cABr","cABr",200,-10,10,200,-10,10);
    cABg = new TH2F("cABg","cABg",200,-10,10,200,-10,10);
    std::vector<TH1F*> histsreal;
    std::vector<TF1*> fu;
    std::vector<TH1F*> histsreal2;
    TSystemDirectory dir(samples_[idx].c_str(), samples_[idx].c_str());
    TList *files = dir.GetListOfFiles();
    TSystemFile *file;
    TString filename;
    TIter next(files);
//Loop over files
    while ((file=(TSystemFile*)next())) {
      filename = file->GetName();
      if(file->IsDirectory()) continue;
      std::stringstream name;
      name<<samples_[idx].c_str()<<"/"<<filename;
//open file and tree
      TFile *input(0);	
      input = TFile::Open((name.str()).c_str());
      for (int Nshower = 0; Nshower < 500; ++Nshower){
        std::stringstream snameH;
        snameH<<"shower_"<<Nshower<<"/longEProf_fineBin_had";
        h_H = new TH1F("","" , nbins,120,280);
//Get hadronic shower component only (Pi0 components are for the next steps ;))
        h_H = (TH1F*)input->Get((snameH.str()).c_str());
//Get Shower starting point for each shower        
        std::stringstream SS;
        SS<<"shower_"<<Nshower<<"/showerStart";
        TVector3 * showerStart = (TVector3*)input->Get((SS.str()).c_str());
        HADshapeIL = new TH1F( (snameH.str()).c_str(),(snameH.str()).c_str() ,   16,bins   );
//now focus on the showers started in HCAL[eta=0 and (*showerStart)[2] is the z component of the shower start 
        if ((*showerStart)[2] < 177.7) continue;
//loop over bins and save the profile in a histogram with new binning 
        for (int b = 1; b < nbins; ++b){
//do not consider hits before shower starting points 
          if ( h_H->GetXaxis()->GetBinCenter(b) < (*showerStart)[2]) continue;
          HADshapeIL->Fill((h_H->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/hcalintlen, h_H->GetBinContent(b));
        }
         if (HADshapeIL->Integral()< energies[idx] * 0.3) continue;
//condition to remove second peak in hadronic showers
        if (HADshapeIL->Integral(1,7) > HADshapeIL->Integral(7,16)) continue;
//find alpha and beta parameters from the mean and RMS of the histogram
        Double_t mH = HADshapeIL->GetMean();
        Double_t sH = pow(HADshapeIL->GetRMS(),2);
        abcorr_H->Fill(log(pow(mH,2)/sH),log(mH/sH));
//We found that log of alpha and beta are distributed like Gaussian
        A->Fill(log(pow(mH,2)/sH));
        B->Fill(log(mH/sH));
        cABr->Fill(pow(mH,2)/sH,mH/sH);
        histsreal.push_back(HADshapeIL);
//Now produce energy profile using estimated alpha and beta parameters
        std::stringstream name;
        name <<idx<<Nshower;
        funcH = new TF1("GammaH",GammaProfile,0  , 8,nGammaPar);
        funcH->SetParameters(pow(mH,2)/sH,mH/sH);
        Eprofile = new TH1F((name.str()).c_str(),(name.str()).c_str(),16,bins);
      for (int j = 0; j <1000 ; ++j){
        rand = funcH->GetRandom();
        Eprofile->Fill(rand);
       }

        histsreal2.push_back(Eprofile);
/* 
   can = new TCanvas("ft","ft",700,500);
   TPad *pad = new TPad("pad", "pad",0,0,1,1);
   pad->Draw();
   pad->cd();
   pad->Divide(2,1);
   std::stringstream name;
   name <<idx<<filename<<Nshower<<".png";
   pad->cd(1);
   HADshapeIL->Draw();
   pad->cd(2);
   h_H->Draw();    
   can->Print((name.str()).c_str(),"png");
   can = new TCanvas("ft","ft",700,500);
   std::stringstream name;
   name <<idx<<"_"<<Nshower<<".png";
   HADshapeIL->Draw();
   funcH = new TF1("GammaH",GammaProfile,0  , 5.5,nGammaPar);
   funcH->SetParameters(pow(mH,2)/sH,mH/sH);
   funcH->Draw("same");
   can->Print((name.str()).c_str(),"png");
*/
        }
    }
  
//Normalize shower profiles to 1 and add all histograms
    for (int i = 0; i < histsreal.size(); ++i){
      histsreal[i]->Scale(1/histsreal[i]->Integral());
      histsreal2[i]->Scale(1/histsreal2[i]->Integral());


//Draw each profile seperatly to see shower shape and fitted gamma function
/*
  TCanvas *can1;
  can1 = new TCanvas("ft","ft",700,500);
  can1->cd();
  histsreal2[i]->SetLineColor(2);
  histsreal[i]->Draw("E");
  histsreal2[i]->Draw("Csame");
  std::stringstream eout;
  eout <<i<<".png";
  can1->Print((eout.str()).c_str(),"png");
*/
  }
  
    for (int i = 1; i < histsreal.size(); ++i){
        histsreal[i]->Add(histsreal[i-1]);
        histsreal2[i]->Add(histsreal2[i-1]);
    }
  
//normalize histograms
    histsreal[histsreal.size()-1]->Scale(1/histsreal[histsreal.size()-1]->Integral());
    histsreal2[histsreal2.size()-1]->Scale(1/histsreal2[histsreal2.size()-1]->Integral());




// Now reproduce the shower profiles using the distribution of the alpha and beta parameters found in previous step from fullsim showers
    double meanalpha = abcorr_H->GetMean(1);
    double meanbeta = abcorr_H->GetMean(2);
    double RMSalpha = abcorr_H->GetRMS(1);
    double RMSbeta = abcorr_H->GetRMS(2);
    double corr = abcorr_H->GetCorrelationFactor();

    double x = sqrt((1+corr)/2);
    double y = sqrt((1-corr)/2);

    double alpha;
    double beta;
    int n=1000;
    int nsudo=1000;
    std::vector<TH1F*> hists;

    for (int i = 0; i < nsudo; ++i){
      double z1=gRandom->Gaus(0, 1);
      double z2=gRandom->Gaus(0, 1);
      alpha = meanalpha + RMSalpha * x * z1 + RMSalpha *y *z2;
      beta = meanbeta + RMSbeta * x * z1 - RMSbeta *y *z2;
      As->Fill(alpha);
      Bs->Fill(beta);   
      abcorrRand_H->Fill(alpha,beta);
      cABg->Fill(exp(alpha) ,exp(beta));
      std::stringstream name;
      name <<i;
      hists.push_back(RandomGammaProfile(exp(alpha) ,exp(beta), 8,16, (name.str()).c_str()));
    }

    for (int i = 1; i < hists.size(); ++i){
      hists[i]->Add(hists[i-1]);
    }

/*
cout<<"corr= "<< corr << "  found  "<< abcorrRand_H->GetCorrelationFactor()<<endl;
cout<<"meanalpha= "<< meanalpha<< " found  "<<abcorrRand_H->GetMean(1)<< endl;
cout<<"meanbetap= "<< meanbeta<< " found  "<<abcorrRand_H->GetMean(2)<< endl;
cout<<"RMSalpha= "<< RMSalpha<< " found  "<<abcorrRand_H->GetRMS(1)<< endl;
cout<<"RMSbetap= "<< RMSbeta<< " found  "<<abcorrRand_H->GetRMS(2)<< endl;


cout<<" +++++++++++++++++++++++++++++++++++  "<<endl;
cout<<"corr= "<<cABr->GetCorrelationFactor()<< "  found  "<<cABg->GetCorrelationFactor()<<endl;
cout<<"meanalpha= "<<cABr->GetMean(1)<< "  found  "<<cABg->GetMean(1)<< endl;
cout<<"meanbetap= "<<cABr->GetMean(2)<< "  found  "<<cABg->GetMean(2)<< endl;
cout<<"RMSalpha= "<<cABr->GetRMS(1)<< "  found  "<<cABg->GetRMS(1)<< endl;
cout<<"RMSbetap= "<<cABr->GetRMS(2)<< "  found  "<<cABg->GetRMS(2)<< endl;
*/

    hists[hists.size()-1]->Scale(1/hists[hists.size()-1]->Integral());
    hists[hists.size()-1]->SetLineColor(2);
    hists[hists.size()-1]->SetLineWidth(2);
    hists[hists.size()-1]->SetTitle("Average Hadronic energy profile");
    hists[hists.size()-1]->GetXaxis()->SetTitle("Interaction length");
    hists[hists.size()-1]->GetYaxis()->SetTitle("#frac{1}{E} #frac{dE}{dx}");
    hists[hists.size()-1]->SetStats(0);
    histsreal[histsreal.size()-1]->SetStats(0);
    histsreal2[histsreal2.size()-1]->SetStats(0);
    histsreal[histsreal.size()-1]->SetLineWidth(2);
    histsreal2[histsreal2.size()-1]->SetLineWidth(2);
    histsreal2[histsreal2.size()-1]->SetLineColor(8);

//Draw histograms
    can = new TCanvas("ft","ft",700,500);
    can->cd();
    hists[hists.size()-1]->Draw("C");
    histsreal[histsreal.size()-1]->Draw("Esame");
    histsreal2[histsreal2.size()-1]->Draw("Csame");

    TLegend* leg = new TLegend(0.5,0.6,0.9,0.9);
    leg->SetFillStyle ( 0);
    leg->SetFillColor ( 0);
    leg->SetBorderSize( 0);
    leg->AddEntry( histsreal[histsreal.size()-1], "Fullsim"                           , "lep");
    leg->AddEntry(histsreal2[histsreal2.size()-1] , "#Gamma(#alpha,#beta) from fullsim"            , "L");
    leg->AddEntry(hists[hists.size()-1], "#Gamma(#alpha,#beta) from #alpha & #beta distribution", "L");
    leg->Draw("same");

    std::stringstream textL;
    textL <<"Energy = "<<energies[idx];
    pave = new TPaveText(0.5,0.5,0.9,0.6,"NDC");
    pave->AddText((textL.str()).c_str());
    pave->Draw("same");
    std::stringstream nameout;
    nameout <<energies[idx]<<"_averageprofile.png";
    can->Print((nameout.str()).c_str(),"png");
    pave->Clear();
}
//abcorrRand_H->SetMarkerColor(kBlue);
//abcorrRand_H->Draw();
//abcorr_H->Draw();
//A->Scale(1/A->Integral());
//As->Scale(1/As->Integral());
//As->Draw("C");
//As->Draw("same");

//B->Scale(1/B->Integral());
//Bs->Scale(1/Bs->Integral());
//Bs->Draw("C");
//Bs->Draw("same");
}
