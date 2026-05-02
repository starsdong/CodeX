#include <fstream>
#include <iostream>

using namespace std;

void hf_review_style()
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.01);
  gStyle->SetTickLength(0.03,"X");
  gStyle->SetTickLength(0.03,"Y");
  gStyle->SetGridWidth(1);
  gStyle->SetGridStyle(2);
  gStyle->SetGridColor(kGray+1);

  gStyle->SetFillColor(10);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(2);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.18);

  gStyle->SetNdivisions(108,"X");
  gStyle->SetNdivisions(108,"Y");
  gStyle->SetLabelOffset(0.002,"X");
  gStyle->SetLabelOffset(0.01,"Y");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetTitleSize(0.07,"X");
  gStyle->SetTitleSize(0.07,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");

  gStyle->SetPalette(1);
  gStyle->SetMarkerSize(1.8);
  gStyle->SetMarkerStyle(20);
  gStyle->SetLegendFillColor(10);
}

void make_LcD0_ratio(const Int_t mPlotTh = 1)
{
//  gROOT->Reset();
  hf_review_style();
  gSystem->mkdir("fig", kTRUE);

  const TString kHFInputDir = "/Users/starsdong/Work/work/HF";
  const Int_t kPythiaColor = kGray+2;
  const Int_t kPythiaLineStyle = 2;
  const Int_t kPythiaCRLineStyle = 7;
  const Int_t kTsinghuaColor = kAzure+10;
  const Int_t kRappColor = kRed-4;
  const Int_t kCaoKoColor = kRed+1;

  const Int_t NConfig = 4;
  const Int_t NMax = 5;
  const Int_t N[NConfig] = {5, 4, 5, 3}; // ALICE_pp, ALICE_pPb, STAR_AuAu
  const Char_t *ExpName[NConfig] = {"ALICE","ALICE","ALICE", "STAR"};
  const Char_t *CollName[NConfig] = {"pp_7TeV","pPb_5p02TeV","PbPb_0_10_5p02TeV","AuAu_10_80_200GeV"};
  double x[NConfig][NMax], y[NConfig][NMax], ye[NConfig][NMax], yes_u[NConfig][NMax], yes_d[NConfig][NMax];
  TGraphErrors *gr[NConfig];

  ifstream inData;
  for(int i=0;i<NConfig;i++) {
    inData.open(Form("%s/dat/%s_LcD0_%s.txt", kHFInputDir.Data(), ExpName[i], CollName[i]));
    for(int j=0;j<N[i];j++) {
      inData >> x[i][j] >> y[i][j] >> ye[i][j] >> yes_u[i][j] >> yes_d[i][j];
    }
    inData.close();
    gr[i] = new TGraphErrors(N[i], x[i], y[i], 0, ye[i]);
    cout << Form("%s/dat/%s_LcD0_%s.txt", kHFInputDir.Data(), ExpName[i], CollName[i]) << endl;
    gr[i]->Print();
  }
  
  // Theory calculations
  // pythia
  TFile *fin = new TFile(Form("%s/root/pythia_D_ratio.root", kHFInputDir.Data()));
  TGraphErrors *gr_pythia_tmp = (TGraphErrors *)fin->Get("Ratio_Lc_D0");
  TGraph *gr_pythia = new TGraph(*gr_pythia_tmp);
  gr_pythia->SetLineWidth(2);
  gr_pythia->SetLineColor(kPythiaColor);
  gr_pythia->SetLineStyle(kPythiaLineStyle);
  fin->Close();

  fin = new TFile(Form("%s/root/PYTHIA8_D_ratio_200GeV_20191101_CR.root", kHFInputDir.Data()));
  TGraphErrors *gr_pythia_CR_tmp = (TGraphErrors *)fin->Get("gLc2D0");
  TGraph *gr_pythia_CR = new TGraph(*gr_pythia_CR_tmp);
  gr_pythia_CR->SetLineWidth(2);
  gr_pythia_CR->SetLineStyle(kPythiaCRLineStyle);
  gr_pythia_CR->SetLineColor(kPythiaColor);
  
  TGraph *gr_Tsinghua = new TGraph(Form("%s/dat/Tsinghua_LcD0_1_AuAu200GeV_10_80.txt", kHFInputDir.Data()), "%lg %lg");
  gr_Tsinghua->SetLineWidth(3);
  gr_Tsinghua->SetLineColor(kTsinghuaColor);
  gr_Tsinghua->SetLineStyle(1);

  TGraph *gr_CaoKo = new TGraph(Form("%s/dat/LcD0_Cao_1911.txt", kHFInputDir.Data()), "%lg %lg");
  gr_CaoKo->SetLineWidth(2);
  gr_CaoKo->SetLineColor(kCaoKoColor);
  gr_CaoKo->SetLineStyle(7);

  const Int_t N_Rapp = 17;
  Double_t pt_Rapp[N_Rapp], y1_Rapp[N_Rapp], y2_Rapp[N_Rapp], y_Rapp[N_Rapp], ye_Rapp[N_Rapp];
  inData.open(Form("%s/dat/LcD0_Rapp_1905.09216.txt", kHFInputDir.Data()));
  for(int i=0;i<N_Rapp;i++) {
    inData >> pt_Rapp[i] >> y1_Rapp[i] >> y2_Rapp[i];
    y_Rapp[i] = 0.5*(y1_Rapp[i]+y2_Rapp[i]);
    ye_Rapp[i] = 0.5*fabs(y1_Rapp[i]-y2_Rapp[i]);
  }
  inData.close();
  TGraphErrors *gr_Rapp = new TGraphErrors(N_Rapp, pt_Rapp, y_Rapp, 0, ye_Rapp);
  gr_Rapp->SetFillColor(kRed-9);
  gr_Rapp->SetFillStyle(3004);
  gr_Rapp->SetLineColor(kRappColor);
  gr_Rapp->SetLineWidth(2);
  
  

  TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  c1->SetFillColor(10);
  c1->SetFillStyle(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetFrameFillColor(10);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.16);
  c1->SetRightMargin(0.02);
  c1->SetTopMargin(0.02);
  c1->SetBottomMargin(0.16);
  // c1->SetTickx();
  // c1->SetTicky();
  c1->Draw();
  c1->cd();


  double x1 = 0.01;
  double x2 = 8.5;
  double y1 = 0.;
  double y2 = 1.9;
  TH1 *h0 = new TH1D("h0","",1,x1, x2);
  h0->SetMinimum(y1);
  h0->SetMaximum(y2);
  h0->GetXaxis()->SetNdivisions(208);
  h0->GetXaxis()->CenterTitle();
  h0->GetXaxis()->SetTitle("Transverse Momentum p_{T} (GeV/c)");
  h0->GetXaxis()->SetTitleOffset(1.1);
  h0->GetXaxis()->SetTitleSize(0.06);
  h0->GetXaxis()->SetLabelOffset(0.01);
  h0->GetXaxis()->SetLabelSize(0.05);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetTitleFont(42);
  h0->GetYaxis()->SetNdivisions(505);
  h0->GetYaxis()->SetTitle("#Lambda_{c}/D^{0}");
  h0->GetYaxis()->SetTitleOffset(1.0);
  h0->GetYaxis()->SetTitleSize(0.06);
  h0->GetYaxis()->SetLabelOffset(0.01);
  h0->GetYaxis()->SetLabelSize(0.05);
  h0->GetYaxis()->SetLabelFont(42);
  h0->GetYaxis()->SetTitleFont(42);
  h0->Draw("c");

  TLine *l1 = new TLine(x1,y1,x2,y1);
  l1->SetLineWidth(3);
  l1->Draw("same");
  TLine *l2 = new TLine(x1,y2,x2,y2);
  l2->SetLineWidth(3);
  l2->Draw("same");
  TLine *l3 = new TLine(x1,y1,x1,y2);
  l3->SetLineWidth(3);
  l3->Draw("same");
  TLine *l4 = new TLine(x2,y1,x2,y2);
  l4->SetLineWidth(3);
  l4->Draw("same");


  gr_pythia->Draw("c");
  gr_pythia_CR->Draw("c");

  
  if(mPlotTh) {
    gr_Rapp->Draw("e3");
    gr_Tsinghua->Draw("c");
    gr_CaoKo->Draw("c");
  }

  // plotting the data points
  const Int_t kStyle[NConfig] = {25, 26, 21, 20};
  const Int_t kColor[NConfig] = {1, 1, 1, 4};
  const Double_t kSize[NConfig] = {2.0, 2.0, 2.0, 1.8};
  const double xo = fabs(x2-x1)/100.;
  const double yo = fabs(y2-y1)/100.;
  for(int i=0;i<NConfig;i++) {

    if(i!=3) continue; // STAR Au+Au only
    for(int j=0;j<N[i];j++) {
      double x1 = x[i][j] - xo;
      double x2 = x[i][j] + xo;
      double y1 = y[i][j] - yes_d[i][j];
      double y2 = y[i][j] + yes_u[i][j];
      
      TLine *la = new TLine(x1, y1, x1, y1+yo);
      la->SetLineColor(kColor[i]);
      la->SetLineWidth(1);
      la->Draw("same");
      TLine *lb = new TLine(x2, y1, x2, y1+yo);
      lb->SetLineColor(kColor[i]);
      lb->SetLineWidth(1);
      lb->Draw("same");
      TLine *lc = new TLine(x1, y2, x1, y2-yo);
      lc->SetLineColor(kColor[i]);
      lc->SetLineWidth(1);
      lc->Draw("same");
      TLine *ld = new TLine(x2, y2, x2, y2-yo);
      ld->SetLineColor(kColor[i]);
      ld->SetLineWidth(1);
      ld->Draw("same");
      TLine *le = new TLine(x1, y1, x2, y1);
      le->SetLineColor(kColor[i]);
      le->SetLineWidth(2);
      le->Draw("same");
      TLine *lf = new TLine(x1, y2, x2, y2);
      lf->SetLineColor(kColor[i]);
      lf->SetLineWidth(2);
      lf->Draw("same");
    }
    gr[i]->SetMarkerStyle(kStyle[i]);
    gr[i]->SetMarkerColor(kColor[i]);
    gr[i]->SetMarkerSize(kSize[i]);
    gr[i]->SetLineColor(kColor[i]);
    gr[i]->SetLineWidth(2);
    gr[i]->Draw("p");
  }
  

   TLegend *leg = new TLegend(0.18, 0.88, 0.58, 0.95);
   leg->SetFillColor(10);
   leg->SetFillStyle(10);
   leg->SetLineStyle(4000);
   leg->SetLineColor(10);
   leg->SetLineWidth(0.);
   leg->SetTextFont(42);
   leg->SetTextSize(0.038);
   leg->SetMargin(0.20);
   leg->AddEntry(gr[3], "Au+Au@200GeV, 10-80%", "pl");
   leg->Draw();
  
   leg = new TLegend(0.69, 0.70, 0.99, 0.95);
   leg->SetFillColor(10);
   leg->SetFillStyle(10);
   leg->SetLineStyle(4000);
   leg->SetLineColor(10);
   leg->SetLineWidth(0.);
   leg->SetTextFont(42);
   leg->SetTextSize(0.032);
   leg->SetMargin(0.28);
   leg->AddEntry(gr_Tsinghua, "Tsinghua 10-80%", "l");
   leg->AddEntry(gr_Rapp, "He,Rapp 0-20%", "f");
   leg->AddEntry(gr_CaoKo, "Cao,Ko 10-80%", "l");
   leg->AddEntry(gr_pythia_CR, "PYTHIA8 CR", "l");
   leg->AddEntry(gr_pythia, "PYTHIA8", "l");
   leg->Draw();

   c1->Update();
  c1->cd();
  
  c1->SaveAs("fig/LcD0_ratio.pdf");
  c1->SaveAs("fig/LcD0_ratio.png");
  
  
}
