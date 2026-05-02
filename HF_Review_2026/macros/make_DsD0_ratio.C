#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

void make_DsD0_ratio(const Int_t mPlotTh = 1)
{
  gROOT->Reset();
  gSystem->mkdir("fig", kTRUE);

  const TString kHFInputDir = "/Users/starsdong/Work/work/HF";
  const TString kDsModelDir = "/Users/starsdong/Work/Paper/MINE/HFT_Review/macros/Signal/Ds/DsD0_ratio_Theory";
  const TString kDsTheoryFile = "/Users/starsdong/Work/Paper/MINE/2021_HFT_Ds/macros_20190322/root/ALICE_pp7TeV_PbPb5TeV_Pythia.root";
  const Int_t kPythia200Color = kGray+2;
  const Int_t kTsinghuaColor = kAzure+10;
  const Int_t kTamuColor = kRed-4;
  const Int_t kCaoKoColor = kRed+1;

  const Int_t NConfig = 5;
  const Int_t NMax = 8;
  const Int_t N[NConfig] = {4, 4, 5, 8, 4}; // ALICE_pp, ALICE_PbPb, STAR_AuAu
  const Char_t *ExpName[NConfig] = {"ALICE","ALICE","STAR","CMS","CMS"};
  const Char_t *CollName[NConfig] = {"pp_7TeV","PbPb_30_50_5p02TeV","AuAu_10_40_200GeV","pp_5p02TeV","PbPb_0_100_5p02TeV"};
  double x[NConfig][NMax], y[NConfig][NMax], ye[NConfig][NMax], yes[NConfig][NMax];
  TGraphErrors *gr[NConfig];

  ifstream inData;
  for(int i=0;i<NConfig;i++) {
    inData.open(Form("%s/dat/%s_DsD0_%s.txt", kHFInputDir.Data(), ExpName[i], CollName[i]));
    for(int j=0;j<N[i];j++) {
      inData >> x[i][j] >> y[i][j] >> ye[i][j] >> yes[i][j];
    }
    inData.close();
    gr[i] = new TGraphErrors(N[i], x[i], y[i], 0, ye[i]);
  }
  
  // Theory calculations
  // pythia
  double pt1[100], r1[100];
  int n = 0;
  TFile *fin = new TFile(kDsTheoryFile.Data());
  TGraph *gr_pythia_RHIC_tmp = new TGraph((TH1D *)fin->Get("Pythia_Monash_tune_200GeV"));
  double sum_pt_sub = 0;
  double sum_r_sub = 0;
  int n_sub = 0; // average over 10 GeV/c
  for(int i=0;i<gr_pythia_RHIC_tmp->GetN();i++) {
    double xx = gr_pythia_RHIC_tmp->GetX()[i];
    if(xx<10.0) {
      pt1[n] = xx;
      r1[n] = gr_pythia_RHIC_tmp->GetY()[i];
      n++;
    } else if(xx<15.0){
      sum_pt_sub += gr_pythia_RHIC_tmp->GetX()[i];
      sum_r_sub += gr_pythia_RHIC_tmp->GetY()[i];
      n_sub++;
    }
  }
  pt1[n] = sum_pt_sub/n_sub;
  r1[n] = sum_r_sub/n_sub;
  n++;
  TGraph *gr_pythia_RHIC = new TGraph(n, pt1, r1);
  gr_pythia_RHIC->SetLineWidth(2);
  gr_pythia_RHIC->SetLineColor(kPythia200Color);
  gr_pythia_RHIC->SetLineStyle(2);
  gr_pythia_RHIC->SetFillColor(kPythia200Color);
  

  // TFile *fin = new TFile("root/pythia_D_ratio.root");
  // TGraphErrors *gr_pythia = (TGraphErrors *)fin->Get("Ratio_Ds_D0");
  // gr_pythia->SetLineWidth(2);
  // gr_pythia->SetLineColor(8);

  // TFile *fin = new TFile("root/PYTHIA8_D_ratio_7TeV_20181029_0.root");
  // TDirectoryFile *fDs = (TDirectoryFile *)fin->Get("Ds");
  // TH1D *h_pythia8 = (TH1D *)fDs->Get("DstoD0_pt_whole_rebin");
  // TGraph *gr_pythia8 = new TGraph(h_pythia8);
  // gr_pythia8->SetLineWidth(1);
  // gr_pythia8->SetLineStyle(1);

  TGraphErrors *gr_Tsinghua = new TGraphErrors((kDsModelDir + "/Tsinghua/DsD0_1040_RHIC_sequential.dat").Data(), "%lg %lg %lg %lg");
  gr_Tsinghua->SetFillColor(kAzure-9);
  gr_Tsinghua->SetFillStyle(3004);
  gr_Tsinghua->SetLineColor(kTsinghuaColor);
  gr_Tsinghua->SetLineWidth(3);
  gr_Tsinghua->SetLineStyle(1);

  TGraphErrors *gr_TAMU = new TGraphErrors((kDsModelDir + "/TMMU_0-20%_DsD0_Ratio.dat").Data(), "%lg %lg %lg %lg");
  gr_TAMU->SetFillColor(kRed-9);
  gr_TAMU->SetFillStyle(3004);
  gr_TAMU->SetLineColor(kTamuColor);
  gr_TAMU->SetLineWidth(2);
  gr_TAMU->SetLineStyle(1);

  TGraphErrors *gr_CaoKo = new TGraphErrors((kDsModelDir + "/ratio_Ds_D0_1040_Shanshan.dat").Data(), "%lg %lg %lg %lg");
  gr_CaoKo->SetFillColor(kRed-9);
  gr_CaoKo->SetFillStyle(3005);
  gr_CaoKo->SetLineColor(kCaoKoColor);
  gr_CaoKo->SetLineWidth(2);
  gr_CaoKo->SetLineStyle(7);
  

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
  c1->SetTickx();
  c1->SetTicky();
  c1->Draw();
  c1->cd();


  double x1 = 0.01;
  double x2 = 8.5;
  double y1 = 0.;
  double y2 = 0.7;
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
  h0->GetYaxis()->SetTitle("D_{s}^{+}/D^{0}");
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

  gr_pythia_RHIC->Draw("c");

  if(mPlotTh) {
    gr_Tsinghua->Draw("e3");
    gr_TAMU->Draw("e3");
    gr_CaoKo->Draw("e3");
    gr_Tsinghua->Draw("c");
    gr_TAMU->Draw("c");
    gr_CaoKo->Draw("c");
  }

  // plotting the data points
  const Int_t kStyle[NConfig] = {25, 21, 20, 26, 22};
  const Int_t kColor[NConfig] = {1, 1, 4, 4, 4};
  const Double_t size[NConfig] = {2.0, 2.2, 1.8, 2.0, 2.5};
  const double xo = fabs(x2-x1)/100.;
  const double yo = fabs(y2-y1)/100.;
  for(int i=0;i<NConfig;i++) {

    if(i!=2) continue; // STAR Au+Au only
    for(int j=0;j<N[i];j++) {
      double x1 = x[i][j] - xo;
      double x2 = x[i][j] + xo;
      double y1 = y[i][j] - yes[i][j];
      double y2 = y[i][j] + yes[i][j];
      
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
    gr[i]->SetMarkerSize(size[i]);
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
   leg->AddEntry(gr[2], "Au+Au@200GeV, 10-40%", "pl");
   leg->Draw();
  
   TLegend *leg2 = new TLegend(0.69, 0.75, 0.99, 0.95);
   leg2->SetFillColor(10);
   leg2->SetFillStyle(10);
   leg2->SetLineStyle(4000);
   leg2->SetLineColor(10);
   leg2->SetLineWidth(0.);
   leg2->SetTextFont(42);
   leg2->SetTextSize(0.032);
   leg2->SetMargin(0.28);
   leg2->AddEntry(gr_Tsinghua, "Tsinghua 10-40%", "l");
   leg2->AddEntry(gr_TAMU, "He,Rapp 0-20%", "l");
   leg2->AddEntry(gr_CaoKo, "Cao,Ko 10-40%", "l");
   leg2->AddEntry(gr_pythia_RHIC, "PYTHIA8", "l");
   leg2->Draw();

   c1->Update();
  c1->cd();
  
  c1->SaveAs("fig/DsD0_ratio.pdf");
  c1->SaveAs("fig/DsD0_ratio.png");
  
  
}
