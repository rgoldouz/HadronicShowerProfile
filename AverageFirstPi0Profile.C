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

void AverageFirstPi0Profile(){
  gROOT->ProcessLine("#include <vector>");
  using namespace std;
  int nbins = 320;
  double ecalradlen =0.89;
  double ecalintlen =19.5;
  double hcalradlen =1.49;
  double hcalintlen =16.42;
  TF1 *funcH;
  TH1F *A, *B, *As, *Bs;
  TH1F *HADshapeIL;
  TH1F *h_H;
  TH2F *abcorr_H;
  TCanvas* can;
  TPad* pad;

  A = new TH1F( "A", "A" ,   50, -5  , 5    );
  As = new TH1F( "As", "As" ,   50, -5  , 5    );
  B = new TH1F( "B", "B" ,   50, -5  , 5   );
  Bs = new TH1F( "Bs", "Bs" ,50,   -5  , 5    );
//energy values
  Double_t energies[] = {1., 2., 3., 5., 10, 50., 100.,150,200};
//  Double_t energies[] = {3};
  const int N=9;
  std::vector<string> samples_;
  std::stringstream dname;
  for( int num = 0; num < 9; num++ ){
    dname << "/afs/cern.ch/work/r/rgoldouz/TASK2015/CMSSW_7_4_12/src/Lukastree/hadshowertuningdata_compressed_V3/"<<energies[num];
    samples_.push_back((dname.str()).c_str());
    dname.str(std::string());
  }
    abcorr_H = new TH2F("abcorr_H","abcorr_H",200,-10,10,200,-10,10);
    abcorrRand_H = new TH2F("abcorrRand_H","abcorrRand_H",200,-10,10,200,-10,10);
    std::vector<TH1F*> histsreal;
    std::vector<TF1*> fu;

  for(unsigned int idx=0; idx<samples_.size(); ++idx){

    cout<<energies[idx]<<endl;
    TSystemDirectory dir(samples_[idx].c_str(), samples_[idx].c_str());
    TList *files = dir.GetListOfFiles();
    TSystemFile *file;
    TString filename;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      filename = file->GetName();
      if(file->IsDirectory()) continue;
//      cout<<filename.Data()<<endl;
      std::stringstream name;
      name<<samples_[idx].c_str()<<"/"<<filename;
      //open file and tree
      TFile *input(0);	
      input = TFile::Open((name.str()).c_str());
      for (int Nshower = 0; Nshower < 500; ++Nshower){
        std::stringstream snameH;
        snameH<<"shower_"<<Nshower<<"/longEProf_fineBin_pi0_1";

        h_H = new TH1F("","" , nbins,120,280);
        h_H = (TH1F*)input->Get((snameH.str()).c_str());
        
        std::stringstream SS;
        SS<<"shower_"<<Nshower<<"/showerStart";
        TVector3 * showerStart = (TVector3*)input->Get((SS.str()).c_str());
        HADshapeIL = new TH1F( (snameH.str()).c_str(),(snameH.str()).c_str() ,  65, 0  , 65    );

        //now focus on the showers started in HCAL[eta=0 and (*showerStart)[2] is the z component of the shower start 
        if ((*showerStart)[2] < 177.7) continue;
        //Just keep fully hadronic showers
//        if (h_H->Integral() < 20 || h_H->Integral()>60) continue; 
        //loop over bins and save the profile in a histogram with new binning 
        for (int b = 1; b < nbins; ++b){
          //do not consider hits before shower starting points 
          if ( h_H->GetXaxis()->GetBinCenter(b) < (*showerStart)[2]) continue;
          HADshapeIL->Fill((h_H->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/ecalradlen, h_H->GetBinContent(b));
        }

        if (HADshapeIL->Integral() < 2 || HADshapeIL->Integral()>4) continue;
        Double_t mH = HADshapeIL->GetMean();
        Double_t sH = pow(HADshapeIL->GetRMS(),2);
        abcorr_H->Fill(log(pow(mH,2)/sH),log(mH/sH));
        A->Fill(log(pow(mH,2)/sH));
        B->Fill(log(mH/sH));
//cout<<HADshapeIL->Integral()<<endl;
//        abcorr_H->Fill(pow(mH,2)/sH,mH/sH);
        histsreal.push_back(HADshapeIL);

       can = new TCanvas("ft","ft",700,500);
        std::stringstream name;
        name <<i<<".png";
        histsreal[i]->Draw();
        Eprofile = new TH1F((name.str()).c_str(),(name.str()).c_str(),65,0,65);
        
        funcH->SetParameters(pow(mH,2)/sH,mH/sH);
        Eprofile->Draw("same")
        can->Print((name.str()).c_str(),"png");


/*        can = new TCanvas("ft","ft",700,500);
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
*/
        } 
    }
  }
cout<<histsreal.size()<<endl;
    for (int i = 0; i < histsreal.size(); ++i){
        histsreal[i]->Scale(1/histsreal[i]->Integral());
    }
/*    for (int i = 0; i < histsreal.size(); ++i){
        can = new TCanvas("ft","ft",700,500);
        std::stringstream name;
        name <<i<<".png";
        histsreal[i]->Draw();
        can->Print((name.str()).c_str(),"png");
    }
*/
    for (int i = 1; i < histsreal.size(); ++i){
        histsreal[i]->Add(histsreal[i-1]);
    }
    can = new TCanvas("ft","ft",700,500);
    can->cd();
    histsreal[histsreal.size()-1]->Scale(1/histsreal[histsreal.size()-1]->Integral());
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
    gRandom = new TRandom3();
    TH1F *Eprofile;
    std::vector<TH1F*> hists;

    for (int i = 0; i < nsudo; ++i){
      double z1=gRandom->Gaus(0, 1);
      double z2=gRandom->Gaus(0, 1);
      double rand;

      alpha = meanalpha + RMSalpha * x * z1 + RMSalpha *y *z2;
      beta = meanbeta + RMSbeta * x * z1 - RMSbeta *y *z2;
      As->Fill(alpha);
      Bs->Fill(beta);   
//      if (alpha<0 || beta<0) continue;
      abcorrRand_H->Fill(alpha,beta);
      std::stringstream name;
      name <<i;
      Eprofile = new TH1F((name.str()).c_str(),(name.str()).c_str(),65,0,65);
      funcH = new TF1("GammaH",GammaProfile,0  , 65,nGammaPar);
      funcH->SetParameters(exp(alpha) ,exp(beta));
//      funcH->SetParameters(alpha ,beta);
      for (int j = 0; j < n; ++j){
        rand = funcH->GetRandom();
        Eprofile->Fill(rand);
      }
//      Eprofile->Scale(1/Eprofile->Integral());
      hists.push_back(Eprofile);
    }

    for (int i = 1; i < hists.size(); ++i){
      hists[i]->Add(hists[i-1]);
    }



    hists[hists.size()-1]->Scale(1/hists[hists.size()-1]->Integral());
    hists[hists.size()-1]->SetLineColor(2);
    hists[hists.size()-1]->SetLineWidth(2);
    histsreal[histsreal.size()-1]->SetLineWidth(2);
//      hists[0]->Draw();
//    for (int i = 1; i < hists.size(); ++i){
//      hists[i]->SetLineColor(2);
//      hists[i]->Draw("same");
//    }

    hists[hists.size()-1]->Draw("C");
    histsreal[histsreal.size()-1]->Draw("Esame");

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
