#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

void make_D0_v2(const int flag=0)
{
  gROOT->Reset();
  gSystem->mkdir("fig", kTRUE);

  const TString kV2InputDir = "/Users/starsdong/Work/Paper/MINE/2017_HFT_D0v2PRL/paper/figure4_v2CompareModel";

  const Int_t n_data = 6;
  const Int_t n_TAMU = 30;
  const Int_t n_Duke = 13;
  const Int_t n_SUB = 18; //50;
  const Int_t n_hydro = 15;
  const Int_t n_LBT = 15;
  const Int_t n_phi1 = 7;
  const Int_t n_phi = 15;
  const Float_t scale_phi = 1.028;
  const Int_t n_ks = 20;
  const Int_t n_omega = 11;
  const Double_t sc_omega = 1.177; // correction factor for omega centrality bias
  const Int_t n_scan = 7;
  const Int_t kDukeColor = kAzure+10;
  const Int_t kDukeStyle = 10;
  const Int_t kLBTColor = kBlue+1;
  const Int_t kLBTStyle = 4;
  const Int_t kTAMUColor = kGreen+2;
  const Int_t kTAMUStyle = 1;
  const Int_t kTAMUNoDiffStyle = 3;
  const Int_t kPHSDColor = kMagenta+1;
  const Int_t kPHSDStyle = 8;
  const Int_t kSubatechColor = kRed+1;
  const Int_t kSubatechStyle = 7;
  const Int_t kCataniaColor = kOrange+7;
  const Int_t kCataniaStyle = 5;
  const Int_t kHydroColor = kGray+2;
  const Int_t kHydroStyle = 9;


  Double_t x_data[n_data], y_data[n_data], ye_data[n_data], yes_data[n_data], yesL_data[n_data];
  Double_t yeL_data[n_data], yeU_data[n_data];
  //  Double_t x_data[n_data], y_data[n_data], ye_data[n_data];
  Double_t x_TAMU[n_TAMU], y_TAMU[n_TAMU], ye_TAMU[n_TAMU], y_0_TAMU[n_TAMU];
  Double_t x_Duke[n_Duke], y_Duke[n_Duke], ye_Duke[n_Duke];
  Double_t x_SUB[n_SUB], y_SUB[n_SUB];
  Double_t x_hydro[n_hydro], y_hydro[n_hydro];
						
  Double_t x_phi1[n_phi1], y_phi1[n_phi1], ye_phi1[n_phi1], yes_phi1[n_phi1];
  Double_t x_phi[n_phi], y_phi[n_phi], ye_phi[n_phi], yes_phi[n_phi];
  Double_t x_ks[n_ks], y_ks[n_ks], ye_ks[n_ks], yes_ks[n_ks];
  Double_t x_omega[n_omega], y_omega[n_omega], ye_omega[n_omega], yes_omega[n_omega];
  
  const Double_t scale_w_Ks = 1.00;  // when comparing MB data, centrality bias - no correction default
  ifstream inData;
  inData.open(Form("%s/Run14/D0_v2_04262016.txt", kV2InputDir.Data()));
  cout << " ======== Data ========= " << endl;
  for(int i=0;i<n_data;i++) {
    inData >> x_data[i] >> y_data[i] >> ye_data[i] >> yes_data[i] >> yesL_data[i];
    yeL_data[i] = sqrt(ye_data[i]*ye_data[i]+yes_data[i]*yes_data[i]+yesL_data[i]*yesL_data[i]);
    yeU_data[i] = sqrt(ye_data[i]*ye_data[i]+yes_data[i]*yes_data[i]);
    
    cout <<  x_data[i] << "\t" << y_data[i] << "\t" << ye_data[i] << "\t" << yes_data[i] << "\t" << yesL_data[i] << endl;
  }
  inData.close();

  TFile *fin = new TFile(Form("%s/Systematics_D0vn_SL16d_2016-10-14.ME.root", kV2InputDir.Data()));
  TGraphErrors *gr_data_0_80 = (TGraphErrors *)(fin->Get("vnStat_0_80"));
  gr_data_0_80->RemovePoint(0);
  TGraphErrors *gr_data_0_80_sys = (TGraphErrors *)(fin->Get("vnSyst_0_80"));
  gr_data_0_80_sys->RemovePoint(0);

  const Int_t n_data_new = 8;
  Double_t nonflow[n_data_new] = { 0.032629 , 0.0336555 , 0.0336555 , 0.033947 , 0.0346236 , 0.0353009 , 0.0361988 , 0.0382869};
                           //    { 0.032629 , 0.0336555 , 0.0336555 , 0.033947 , 0.0346236 , 0.0353009 , 0.035973 , 0.0366486 , 0.0382869};
  Double_t x_data_new[n_data_new], y_data_new[n_data_new], ye_data_new[n_data_new], yes_data_new[n_data_new], yesL_data_new[n_data_new];
  Double_t yeL_data_new[n_data_new], yeU_data_new[n_data_new];
  for(int i=0;i<gr_data_0_80->GetN();i++) {
    cout << i << endl;
    x_data_new[i] = gr_data_0_80->GetX()[i];
    y_data_new[i] = gr_data_0_80->GetY()[i];
    ye_data_new[i] = gr_data_0_80->GetEY()[i];
    yes_data_new[i] = gr_data_0_80_sys->GetEY()[i];
    yesL_data_new[i] = nonflow[i];
  
    yeL_data_new[i] = sqrt(ye_data_new[i]*ye_data_new[i]+yes_data_new[i]*yes_data_new[i]+yesL_data_new[i]*yesL_data_new[i]) * scale_w_Ks;
    yeU_data_new[i] = sqrt(ye_data_new[i]*ye_data_new[i]+yes_data_new[i]*yes_data_new[i]) * scale_w_Ks;
  }

  inData.open(Form("%s/Run14/TAMU_v2_0_80_09222015.txt", kV2InputDir.Data()));
  cout << " ======== TAMU ========= " << endl;
  for(int i=0;i<n_TAMU;i++) {
    double a, b, c;
    inData >> a >> b >> c;
    x_TAMU[i] = a;
    y_TAMU[i] = (b+c)/2.;
    ye_TAMU[i] = fabs(b-c)/2.;
    cout <<  x_TAMU[i] << "\t" << y_TAMU[i] << endl;
  }
  inData.close();

  inData.open(Form("%s/Run14/TAMU_v2_0_80_nodiff_09222015.txt", kV2InputDir.Data()));
  cout << " ======== TAMU ========= " << endl;
  for(int i=0;i<n_TAMU;i++) {
    inData >> x_TAMU[i] >> y_0_TAMU[i];
    cout <<  x_TAMU[i] << "\t" << y_0_TAMU[i] << endl;
  }
  inData.close();
  
  // inData.open("/Users/starsdong/work/datapoints/Run14/Duke_D0v2_20_40_1505.01413.txt");
  // cout << " ======== Duke ========= " << endl;
  // for(int i=0;i<n_Duke;i++) {
  //   double y1, y2;
  //   inData >> x_Duke[i] >> y1 >> y2;
  //   y_Duke[i] = (y1+y2)/2.;
  //   ye_Duke[i] = fabs(y2-y1)/2.;
  //   cout <<  x_Duke[i] << "\t" << y_Duke[i] << "\t" << ye_Duke[i] << endl;
  // }
  // inData.close();

  inData.open(Form("%s/Run14/SUBATECH_v2_0_80_EPOS3.txt", kV2InputDir.Data()));
  cout << " ======== SUBATECH ========= " << endl;
  for(int i=0;i<n_SUB;i++) {
    inData >> x_SUB[i] >> y_SUB[i];
    cout <<  x_SUB[i] << "\t" << y_SUB[i] << endl;
  }
  inData.close();

  inData.open(Form("%s/Run14/hydro_v2_b_7_8_new.txt", kV2InputDir.Data()));
  cout << " ======== hydro ========= " << endl;
  for(int i=0;i<n_hydro;i++) {
    inData >> x_hydro[i] >> y_hydro[i];
    cout <<  x_hydro[i] << "\t" << y_hydro[i] << endl;
  }
  inData.close();
  
  inData.open(Form("%s/Run14/phi_v2_0_80_PRC.txt", kV2InputDir.Data()));
  cout << " ======== Phi ========= " << endl;
  for(int i=0;i<n_phi1;i++) {
    inData >> x_phi1[i] >> y_phi1[i] >> ye_phi1[i] >> yes_phi1[i];
    cout <<  x_phi1[i] << "\t" << y_phi1[i] << "\t" << ye_phi1[i] << "\t" << yes_phi1[i] << endl;
  }
  inData.close();

  inData.open(Form("%s/Run14/phi_v2_0_80_Run1011.txt", kV2InputDir.Data()));
  cout << " ======== Phi ========= " << endl;
  for(int i=0;i<n_phi;i++) {
    inData >> x_phi[i] >> y_phi[i] >> ye_phi[i] >> yes_phi[i];
    y_phi[i] *= scale_phi;
    ye_phi[i] *= scale_phi;
    yes_phi[i] *= scale_phi;
    cout <<  x_phi[i] << "\t" << y_phi[i] << "\t" << ye_phi[i] << "\t" << yes_phi[i] << endl;    
  }
  inData.close();

  inData.open(Form("%s/Run14/ks_v2_0_80_PRC77.txt", kV2InputDir.Data()));
  cout << " ======== Ks ========= " << endl;
  for(int i=0;i<n_ks;i++) {
    inData >> x_ks[i] >> y_ks[i] >> ye_ks[i] >> yes_ks[i];
    cout <<  x_ks[i] << "\t" << y_ks[i] << "\t" << ye_ks[i] << "\t" << yes_ks[i] << endl;
  }
  inData.close();

  inData.open(Form("%s/Run14/omega_v2_0_80.txt", kV2InputDir.Data()));
  cout << " ======== Omega ========= " << endl;
  for(int i=0;i<n_omega;i++) {
    inData >> x_omega[i] >> y_omega[i] >> ye_omega[i] >> yes_omega[i];
    y_omega[i] *= sc_omega;
    ye_omega[i] *= sc_omega;
    yes_omega[i] *= sc_omega;
    //    x_omega[i] = sqrt(x_omega[i]*x_omega[i]+MassOmega*MassOmega)-MassOmega;
    cout <<  x_omega[i] << "\t" << y_omega[i] << "\t" << ye_omega[i] << "\t" << yes_omega[i] << endl;
  }
  inData.close();

  TGraph *gr_Duke_scan[n_scan]; 
  for(int i=0;i<n_scan;i++) {
    gr_Duke_scan[i] = new TGraph(Form("%s/Run14/Duke_/data-cen-00-80/v2_cen-00-80_D0%d0.dat", kV2InputDir.Data(), i+1),"%lg %lg");
    //    gr_Duke_scan[i] = new TGraph(Form("Run14/Duke/data-v2/data-cen-00-80/v2_cen-00-80_D0%d0.dat",i+1),"%lg %lg");
  }

  TGraph *gr_PHSD = new TGraph(Form("%s/Run14/PHSD_v2_0_80.txt", kV2InputDir.Data()),"%lg %lg");

  TGraph *gr_LBT = new TGraph(Form("%s/Run14/LBT_D0_v2.dat", kV2InputDir.Data()),"%lg %lg");

  TGraph *gr_Catania = new TGraph(Form("%s/Run14/D0_v2_Catania_PRC95.txt", kV2InputDir.Data()),"%lg %lg");
  
  TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //  gStyle->SetEndErrorSize(0);
  c1->SetFillColor(10);
  c1->SetFillStyle(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetFrameFillColor(10);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  //c1->SetLogy();
  c1->SetGridx(0);
  c1->SetGridy(0);
  c1->SetLeftMargin(0.16);
  c1->SetBottomMargin(0.16);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);

  double x1 = 0.0;
  double x2 = 7;
  double y1 = -0.05;
  double y2 = 0.292;
  TH1 *h0 = new TH1D("h0","",1,x1, x2);
  h0->SetMinimum(y1);
  h0->SetMaximum(y2);
  h0->GetXaxis()->SetNdivisions(208);
  h0->GetXaxis()->CenterTitle();
  h0->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h0->GetXaxis()->SetTitleOffset(1.1);
  h0->GetXaxis()->SetTitleSize(0.06);
  h0->GetXaxis()->SetLabelOffset(0.01);
  h0->GetXaxis()->SetLabelSize(0.05);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetTitleFont(42);
  h0->GetYaxis()->SetNdivisions(505);
  h0->GetYaxis()->CenterTitle();
  h0->GetYaxis()->SetTitle("Anisotropy Parameter, v_{2}");
  h0->GetYaxis()->SetTitleOffset(1.0);
  h0->GetYaxis()->SetTitleSize(0.06);
  h0->GetYaxis()->SetLabelOffset(0.015);
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

  TLine *l0 = new TLine(x1, 0, x2, 0);
  l0->SetLineWidth(1);
  l0->SetLineStyle(2);
  l0->Draw("same");

  TGraph *gr_hydro_ref = new TGraph(Form("%s/Run14/hydro.txt", kV2InputDir.Data()),"%lg %lg");
  gr_hydro_ref->SetLineStyle(2);
  gr_hydro_ref->SetLineWidth(2);
  //  gr_hydro_ref->Draw("c");
    
  //  TGraphErrors *gr_TAMU = new TGraphErrors(n_TAMU, x_TAMU, y_TAMU, 0, ye_TAMU);
  TGraph *gr_TAMU = new TGraph(n_TAMU, x_TAMU, y_TAMU);
  gr_TAMU->SetLineWidth(2);
  gr_TAMU->SetLineColor(kTAMUColor);
  gr_TAMU->SetLineStyle(kTAMUStyle);
  gr_TAMU->Draw("c");

  TGraph *gr_TAMU_0 = new TGraph(n_TAMU, x_TAMU, y_0_TAMU);
  gr_TAMU_0->SetLineWidth(2);
  gr_TAMU_0->SetLineColor(kTAMUColor);
  gr_TAMU_0->SetLineStyle(kTAMUNoDiffStyle);
  gr_TAMU_0->Draw("c");
  
  TGraph *gr_SUB = new TGraph(n_SUB, x_SUB, y_SUB);
  gr_SUB->SetLineWidth(2);
  gr_SUB->SetLineColor(kSubatechColor);
  gr_SUB->SetLineStyle(kSubatechStyle);
  gr_SUB->Draw("c");

  TGraph *gr_hydro = new TGraph(n_hydro, x_hydro, y_hydro);
  gr_hydro->SetLineWidth(2);
  gr_hydro->SetLineColor(kHydroColor);
  gr_hydro->SetLineStyle(kHydroStyle);
  gr_hydro->Draw("c");

  // D2piT = 7
  gr_Duke_scan[6]->SetLineWidth(2);
  gr_Duke_scan[6]->SetLineColor(kDukeColor);
  gr_Duke_scan[6]->SetLineStyle(kDukeStyle);
  gr_Duke_scan[6]->Draw("c");

  gr_PHSD->SetLineWidth(2);
  gr_PHSD->SetLineColor(kPHSDColor);
  gr_PHSD->SetLineStyle(kPHSDStyle);
  gr_PHSD->Draw("c");

  gr_LBT->SetLineWidth(2);
  gr_LBT->SetLineColor(kLBTColor);
  gr_LBT->SetLineStyle(kLBTStyle);
  gr_LBT->Draw("c");

  gr_Catania->SetLineWidth(2);
  gr_Catania->SetLineColor(kCataniaColor);
  gr_Catania->SetLineStyle(kCataniaStyle);
  gr_Catania->Draw("c");
  
  if(0) {
  for(int i=1;i<n_data;i++) {
    double x1 = x_data[i]-0.1;
    double x2 = x_data[i]+0.1;
    double y1 = y_data[i]-yes_data[i];
    double y2 = y_data[i]+yes_data[i];

    double y3 = y_data[i] - yesL_data[i];
    double y4 = y_data[i];
    TBox *box = new TBox(x1, y3, x2, y4);
    box->SetLineColor(16);
    box->SetFillColor(16);
    box->Draw("same");
    
    TLine *la = new TLine(x1, y1, x1, y1+0.003);
    la->Draw("same");
    TLine *lb = new TLine(x2, y1, x2, y1+0.003);
    lb->Draw("same");
    TLine *lc = new TLine(x1, y2, x1, y2-0.003);
    lc->Draw("same");
    TLine *ld = new TLine(x2, y2, x2, y2-0.003);
    ld->Draw("same");
    TLine *le = new TLine(x1, y1, x2, y1);
    le->SetLineWidth(2);
    le->Draw("same");
    TLine *lf = new TLine(x1, y2, x2, y2);
    lf->SetLineWidth(2);
    lf->Draw("same");
  }
  
  TGraphErrors *gr_data = new TGraphErrors(n_data-1, x_data+1, y_data+1, 0, ye_data+1);
  gr_data->SetMarkerStyle(24);
  gr_data->SetMarkerColor(1);
  gr_data->SetMarkerSize(1.8);
  gr_data->SetLineWidth(2);
  gr_data->Draw("p");
  }  

  for(int i=0;i<n_data_new;i++) {
    double x1 = x_data_new[i]-0.1;
    double x2 = x_data_new[i]+0.1;
    double y1 = y_data_new[i]-yes_data_new[i];
    double y2 = y_data_new[i]+yes_data_new[i];

    double y3 = y_data_new[i] - yesL_data_new[i];
    double y4 = y_data_new[i];
    TBox *box = new TBox(x1, y3, x2, y4);
    box->SetLineColor(16);
    box->SetFillColor(16);
    box->Draw("same");
    
    TLine *la = new TLine(x1, y1, x1, y1+0.003);
    la->SetLineColor(4);
    la->Draw("same");
    TLine *lb = new TLine(x2, y1, x2, y1+0.003);
    lb->SetLineColor(4);
    lb->Draw("same");
    TLine *lc = new TLine(x1, y2, x1, y2-0.003);
    lc->SetLineColor(4);
    lc->Draw("same");
    TLine *ld = new TLine(x2, y2, x2, y2-0.003);
    ld->SetLineColor(4);
    ld->Draw("same");
    TLine *le = new TLine(x1, y1, x2, y1);
    le->SetLineColor(4);
    le->SetLineWidth(2);
    le->Draw("same");
    TLine *lf = new TLine(x1, y2, x2, y2);
    lf->SetLineWidth(2);
    lf->SetLineColor(4);
    lf->Draw("same");
  }

  gr_data_0_80->SetMarkerStyle(20);
  gr_data_0_80->SetMarkerColor(4);
  gr_data_0_80->SetMarkerSize(1.8);
  gr_data_0_80->SetLineWidth(2);
  gr_data_0_80->SetLineColor(4);
  gr_data_0_80->Draw("p");

  
  TGraphErrors *gr_phi = new TGraphErrors(n_phi, x_phi, y_phi, 0, ye_phi);
  gr_phi->SetMarkerStyle(24);
  gr_phi->SetMarkerColor(2);
  gr_phi->SetMarkerSize(1.5);
  gr_phi->SetLineColor(2);
  gr_phi->SetLineWidth(2);
  //  gr_phi->Draw("p");

  TGraphErrors *gr_ks = new TGraphErrors(n_ks, x_ks, y_ks, 0, ye_ks);
  gr_ks->SetMarkerStyle(24);
  gr_ks->SetMarkerColor(1);
  gr_ks->SetMarkerSize(1.5);
  gr_ks->SetLineColor(1);
  gr_ks->SetLineWidth(2);
  //  gr_ks->Draw("p");

  TGraphErrors *gr_omega = new TGraphErrors(n_omega, x_omega, y_omega, 0, ye_omega);
  gr_omega->SetMarkerStyle(24);
  gr_omega->SetMarkerColor(4);
  gr_omega->SetMarkerSize(1.5);
  gr_omega->SetLineColor(1);
  gr_omega->SetLineWidth(2);
  //  gr_omega->Draw("p");
  
  TLegend *leg = new TLegend(0.18, 0.59, 0.57, 0.96);
  leg->SetFillStyle(0);
  leg->SetLineStyle(0);
  leg->SetLineColor(10);
  leg->SetLineWidth(0.);
  leg->SetTextSize(0.034);
  leg->SetTextFont(42);
  leg->SetMargin(0.27);
  leg->AddEntry(gr_Duke_scan[6], "Duke", "l");
  leg->AddEntry(gr_LBT, "LBT", "l");
  leg->AddEntry(gr_Catania, "Catania", "l");
  leg->AddEntry(gr_TAMU, "TAMU c-quark diff.", "l");
  leg->AddEntry(gr_TAMU_0, "TAMU no c-quark diff.", "l");
  leg->AddEntry(gr_PHSD, "PHSD", "l");
  leg->AddEntry(gr_SUB, "SUBATECH", "l");
  leg->AddEntry(gr_hydro, "3D viscous hydro", "l");
  leg->Draw();

  TLatex *tex = new TLatex(3.05, 0.264, "Au+Au #sqrt{s_{NN}} = 200 GeV");
  tex->SetTextFont(42);
  tex->SetTextSize(0.055);
  tex->Draw("same");

  tex = new TLatex(5.7, 0.234, "0-80%");
  tex->SetTextFont(42);
  tex->SetTextSize(0.055);
  tex->Draw("same");
  
  c1->SaveAs("fig/D0_v2.pdf");
  c1->SaveAs("fig/D0_v2.png");
}
