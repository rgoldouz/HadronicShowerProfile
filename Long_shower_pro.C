/**********************************************************************************
 * In this code I want to fit the longitudinal shower shape 
 * to extract the alpha and beta parameters of gamma functions
 * this time using histograms prepeared by Lukas
 *
 *
 *   Reza Goldouzian, 10 Dec  2015
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
class Showershape {
	private:
		Double_t energy;
		Int_t imip;
		Int_t ieta;
		TF1 *fit;

	public:
	//constructors
	Showershape (Double_t en, Int_t ie) : energy(en),  ieta(ie) {}

	//set members
		void setFit(TF1* _fit) { fit = _fit; }
		void setEnergy(Double_t en) { energy = en; }
		void setEta(Int_t ie) { ieta = ie; }	

	//access members
	Double_t getEnergy() { return energy; }
		Int_t getEta() { return ieta; }
		TF1* getFit() { return fit; }
};


//Gamma  function
////parameters:
////alpha beta

Double_t GammaProfile(Double_t *x,Double_t *par) {
double Gamma1=0;
Gamma1 = (par[1])*pow((x[0])*(par[1]),par[0]-1.0)*exp(-(x[0])*(par[1]))/tgamma(par[0]);
return Gamma1;
}



//function to fit Shower shape

Showershape* get_shape(int num, int ieta, bool do_fit, bool do_show, bool do_print, bool do_batch=false){
	//energy values
	Double_t energies[] = {1., 2., 3., 5., 10, 50., 100.,150,200};
	if (num>=maxHDe || num<0) {
		std::cout << "num must be between 0 and " << maxHDe - 1 << std::endl;
		Showershape* theRes = new Showershape(0,0);
		return theRes;
	}

	if (ieta>maxHDeta || ieta<1) {
		std::cout << "ieta must be between 1 and " << maxHDeta << std::endl;
		Showershape* theRes = new Showershape(0,0);
		return theRes;
	}
	
	//make filenames	
	TCanvas* canab;
        TPaveText* pave;
	std::stringstream dname;
	dname << "/afs/cern.ch/work/r/rgoldouz/TASK2015/CMSSW_7_4_12/src/Lukastree/hadshowertuningdata_compressed_V1/"<<energies[num];
	TSystemDirectory dir((dname.str()).c_str(), (dname.str()).c_str());
	TList *files = dir.GetListOfFiles();
	TSystemFile *file;
	TString filename;
	TIter next(files);
        TH1F *alpha_H, *alpha_F, *alpha_L;
        TH1F *beta_H, *beta_F, *beta_L;
        TH2F *abcorr_H, *abcorr_F, *abcorr_L;

        alpha_H = new TH1F("alpha_H","alpha_H" ,   100, 0, 15  );
        alpha_F = new TH1F("alpha_F","alpha_F" ,   100, 0, 50  );
        alpha_L = new TH1F("alpha_L","alpha_L" ,   100, 0, 15  );

	beta_H = new TH1F("beta_H","beta_H" ,   100, 0, 15  );
        beta_F = new TH1F("beta_F","beta_F" ,   100, 0, 15  );
        beta_L = new TH1F("beta_L","beta_L" ,   100, 0, 15  );


	while ((file=(TSystemFile*)next())) {
        filename = file->GetName();
       if(file->IsDirectory()) continue;
	cout<<filename.Data()<<endl;
        std::stringstream name;
	name<<(dname.str()).c_str()<<"/"<<filename;

cout<<(name.str()).c_str()<<endl;

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

        h_F = new TH1F("","" , nbins,120,280);
        h_F = (TH1F*)input->Get((snameF.str()).c_str());

        h_L = new TH1F("","" , nbins,120,280);
        h_L = (TH1F*)input->Get((snameL.str()).c_str());

        std::stringstream SS;
	SS<<"shower_"<<Nshower<<"/showerStart";
	TVector3 * showerStart = (TVector3*)input->Get((SS.str()).c_str());

	
	EMFshapeRL = new TH1F( (snameF.str()).c_str(),(snameF.str()).c_str() ,   50, 0, 50  );
	EMLshapeIL = new TH1F( (snameL.str()).c_str(),(snameL.str()).c_str() ,   20, 0  , 10  );
	HADshapeIL = new TH1F( (snameH.str()).c_str(),(snameH.str()).c_str() ,   20, 0  , 10    );


	int firstnonzerobin;
	firstnonzerobin = h_F->FindFirstBinAbove(0,1);
	if ((*showerStart)[2] < 177) continue;
        for (int b = 1; b < nbins; ++b){

	if ( h_H->GetXaxis()->GetBinCenter(b) < (*showerStart)[2]) continue;

	if ( h_H->GetXaxis()->GetBinCenter(b) < 152 ){
	HADshapeIL->Fill((h_H->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/ecalintlen, h_H->GetBinContent(b));
	if (h_F->GetXaxis()->GetBinCenter(firstnonzerobin)-(*showerStart)[2] < 2*hcalradlen) EMFshapeRL->Fill((h_F->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/ecalradlen, h_F->GetBinContent(b));
	else EMLshapeIL->Fill((h_F->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/ecalintlen, h_L->GetBinContent(b));
        EMLshapeIL->Fill((h_L->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/ecalintlen, h_L->GetBinContent(b));
	}

        if ( h_H->GetXaxis()->GetBinCenter(b) > 177.7 ){
        HADshapeIL->Fill((h_H->GetXaxis()->GetBinCenter(b) - (*showerStart)[2] )/hcalintlen, h_H->GetBinContent(b));
	if(h_F->GetXaxis()->GetBinCenter(firstnonzerobin)-(*showerStart)[2] < 2*hcalradlen)  EMFshapeRL->Fill((h_F->GetXaxis()->GetBinCenter(b) - (*showerStart)[2])/hcalradlen, h_F->GetBinContent(b));
	else EMLshapeIL->Fill((h_F->GetXaxis()->GetBinCenter(b) - (*showerStart)[2] )/hcalradlen, h_F->GetBinContent(b));
        EMLshapeIL->Fill((h_L->GetXaxis()->GetBinCenter(b) - (*showerStart)[2] )/hcalintlen, h_L->GetBinContent(b));
        }
	}

	double totalE;
	double firstEME;
	double lateEME;
	double HADE;
	firstEME = EMFshapeRL->Integral();
	lateEME = EMLshapeIL->Integral();
	HADE = HADshapeIL->Integral();
	totalE = firstEME + lateEME + HADE;


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
	if (do_fit){
	funcH = new TF1("GammaH",GammaProfile,0  , 50,nGammaPar);	
        funcF = new TF1("GammaF",GammaProfile,0  , 50,nGammaPar);  
        funcL = new TF1("GammaL",GammaProfile,0  , 50,nGammaPar);  



	funcH->SetParameters(pow(mH,2)/sH,mH/sH);
        funcF->SetParameters(pow(mF,2)/sF,mF/sF);
        funcL->SetParameters(pow(mL,2)/sL,mL/sL);


	funcH->SetLineColor(kRed);
	funcH->SetMarkerColor(kRed);
	funcH->SetLineWidth(2);

        funcF->SetLineColor(kRed);
        funcF->SetMarkerColor(kRed);
        funcF->SetLineWidth(2);

        funcL->SetLineColor(kRed);
        funcL->SetMarkerColor(kRed);
        funcL->SetLineWidth(2);


   
 
	if (HADshapeIL->Integral()>0) alpha_H->Fill(log(pow(mH,2)/sH));
        if (EMFshapeRL->Integral()>0) alpha_F->Fill(log(pow(mF,2)/sF));
        if (EMLshapeIL->Integral()>0) alpha_L->Fill(log(pow(mL,2)/sL));

	if (HADshapeIL->Integral()>0) beta_H->Fill(log(mH/sH));
        if (EMFshapeRL->Integral()>0) beta_F->Fill(log(mF/sF));
        if (EMLshapeIL->Integral()>0) beta_L->Fill(log(mL/sL));

cout<<HADshapeIL->Integral()<<"     "<<funcH->Integral(0,50)<<endl;
	//fit
//	if (HADshapeIL->GetEntries()>100) HADshapeIL->Fit(funcH);
//	if (EMFshapeRL->GetEntries()>100) EMFshapeRL->Fit(funcF);
//        if (EMLshapeIL->GetEntries()>100) EMLshapeIL->Fit(funcL);
	}


	if (do_show){
        std::stringstream outfile;
	std::stringstream textH;
        std::stringstream textL;
        std::stringstream textF;

	textH <<"Hadronic profile [Energy = "<<energies[num]<<" , #alpha =  " << funcH->GetParameter(0)<< " , #beta =  " << funcH->GetParameter(1)<<"]  Energy fraction: "<<HADE/totalE;
        textF <<"First #pi^{0} profile Energy = "<<energies[num]<<" , #alpha =  " << funcF->GetParameter(0)<< " , #beta =  " << funcF->GetParameter(1)<<"]  Energy fraction: "<<firstEME/totalE;
        textL <<"Late #pi^{0} profile [Energy = "<<energies[num]<<" , #alpha =  " << funcL->GetParameter(0)<< " , #beta =  " << funcL->GetParameter(1)<<"]  Energy fraction: "<<lateEME/totalE;


        outfile <<filename<<"_E_"<<energies[num]<<"_event_"<<Nshower<<".png";

if (Nshower==144){
cout<<(*showerStart)[2]<<endl;
TCanvas *can1;
  can1 = new TCanvas("ft","ft",700,500);
   can1->cd();
h_H->Draw();
can1->Print("144.png","png");
}

   can = new TCanvas("ft","ft",700,500);
   can->cd();

   TPad *pad = new TPad("pad", "pad",0,0,1,1);
   pad->Draw();
   pad->cd();
   pad->Divide(1,3);

   pad->cd(1);
   EMFshapeRL->Draw("hist");
   funcF->Draw("same");
   EMFshapeRL->SetLineWidth(2);
   EMFshapeRL->SetLineColor(1);
   EMFshapeRL->GetXaxis()->SetTitle("z (radiation length)");
   EMFshapeRL->GetYaxis()->SetTitle("First #pi^{0} <dE/dz>");
   pave = new TPaveText(0.3,0.5,1,0.7,"NDC");
   pave->AddText((textF.str()).c_str());
   pave->SetTextSize(0.08);
   pave->Draw("same");
 
   pad->cd(2);   
   EMLshapeIL->Draw("hist");
   funcL->Draw("same");
   EMLshapeIL->SetLineWidth(2);
   EMLshapeIL->SetLineColor(1);
   EMLshapeIL->GetXaxis()->SetTitle("z (interaction length)");
   EMLshapeIL->GetYaxis()->SetTitle("Late #pi^{0} <dE/dz>");
   pave = new TPaveText(0.3,0.5,1,0.7,"NDC");
   pave->AddText((textL.str()).c_str());
   pave->SetTextSize(0.08);
   pave->Draw("same");

   pad->cd(3);
   HADshapeIL->Draw("hist");
   funcH->Draw("same");
   HADshapeIL->SetLineWidth(2);
   HADshapeIL->SetLineColor(1);
   HADshapeIL->GetXaxis()->SetTitle("z (Interaction length)");
   HADshapeIL->GetYaxis()->SetTitle("Hadronic component <dE/dz>");
   pave = new TPaveText(0.3,0.5,1,0.7,"NDC");
   pave->AddText((textH.str()).c_str());
   pave->Draw("same");
   pave->SetTextSize(0.08);
   can->Print((outfile.str()).c_str(),"png");

}

}
input->Close();
}

std::stringstream outfileab;
std::stringstream textab;


outfileab <<"alphabeta_"<<energies[num]<<".png";

canab = new TCanvas("ab","ab",700,500);
canab->cd();

TPad *padab = new TPad("padab", "padab",0,0,1,1);
padab->Draw();
padab->cd();
padab->Divide(3,2);

padab->cd(1);
alpha_H->Draw("hist");
alpha_H->SetLineWidth(2);
pave = new TPaveText(0.5,0.6,0.9,0.7,"NDC");
pave->AddText("#alpha of Hadronic #Gamma function");
pave->Draw("same");

padab->cd(2);
alpha_F->Draw("hist");
alpha_F->SetLineWidth(2);
pave = new TPaveText(0.5,0.6,0.9,0.7,"NDC");
pave->AddText("#alpha of First #pi^{0} #Gamma function");
pave->Draw("same");

padab->cd(3);
alpha_L->Draw("hist");
alpha_L->SetLineWidth(2);
pave = new TPaveText(0.5,0.6,0.9,0.7,"NDC");
pave->AddText("#alpha of Late #pi^{0} #Gamma function");
pave->Draw("same");

padab->cd(4);
beta_H->Draw("hist");
beta_H->SetLineWidth(2);
pave = new TPaveText(0.5,0.6,0.9,0.7,"NDC");
pave->AddText("#beta of Hadronic #Gamma function");
pave->Draw("same");


padab->cd(5);
beta_F->Draw("hist");
beta_F->SetLineWidth(2);
pave = new TPaveText(0.5,0.6,0.9,0.7,"NDC");
pave->AddText("#beta of First #pi^{0} #Gamma function");
pave->Draw("same");

padab->cd(6);
beta_L->Draw("hist");
beta_L->SetLineWidth(2);
pave = new TPaveText(0.5,0.6,0.9,0.7,"NDC");
pave->AddText("#beta of Late #pi^{0} #Gamma function");
pave->Draw("same");

canab->Print((outfileab.str()).c_str(),"png");
} 
