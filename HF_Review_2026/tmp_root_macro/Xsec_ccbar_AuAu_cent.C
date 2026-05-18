#include "style.C"

void Xsec_ccbar_AuAu_cent()
{
  gROOT->ProcessLine("set_style()");

  const Int_t nPhenAuAu = 5;
  double npartPhenAuAu[nPhenAuAu] = {352.2, 234.6, 140.4, 59.95, 14.5};
  double enpartPhenAuAu[nPhenAuAu] = {3.3, 4.7, 4.9, 3.6, 2.5};
  double sigPhenAuAu[nPhenAuAu] = {0.597, 0.596, 0.731, 0.841, 0.504};
  double errPhenAuAu[nPhenAuAu] = {0.1816, 0.1954, 0.2308, 0.3096, 0.4231};

  // PHENIX: direct total charm cross section per binary collision,
  // sigma_ccbar^bin = N_ccbar/T_AA from PRL 94, 082301 (2005), HEPData Table 17.
  // N_part values use the matching PHENIX Glauber table from PRC 90, 034903 (2014).

  const Int_t nPhenPP = 1;
  double npartPhenPP[nPhenPP] = {2.0};
  double enpartPhenPP[nPhenPP] = {0.0};
  double sigPhenPP[nPhenPP] = {0.567};
  double errPhenPP[nPhenPP] = {0.2012};

  const Int_t nStarPP = 1;
  double npartStarPP[nStarPP] = {1.};
  double enpartStarPP[nStarPP] = {0.0};
  double sigStarPP[nStarPP] = {0.797};
  double errStarPPd[nStarPP] = {0.3620};
  double errStarPPu[nStarPP] = {0.2956};

  const Int_t nStarDAu = 1;
  double npartStarDAu[nStarDAu] = {9.1};
  double enpartStarDAu[nStarDAu] = {0.4};
  double sigStarDAu[nStarDAu] = {1.400};
  double errStarDAu[nStarDAu] = {0.4472};

  const Int_t nFit = 8;
  double npartFit[nFit] = {2.0, 2.0, 9.1, 352.2, 234.6, 140.4, 59.95, 14.5};
  double enpartFit[nFit] = {0., 0., 0., 0., 0., 0., 0., 0.};
  double sigFit[nFit] = {0.567, 0.797, 1.400, 0.597, 0.596, 0.731, 0.841, 0.504};
  double errFit[nFit] = {0.2012, 0.3288, 0.4472, 0.1816, 0.1954, 0.2308, 0.3096, 0.4231};

  const Int_t phenixColor = kBlue+1;
  const Int_t starColor = kRed+1;
  const Int_t phenixMarker = 21;
  const Int_t starMarker = 29;

  TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.01);
  gStyle->SetGridColor(16);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.18);
  c1->SetTopMargin(0.03);
  c1->SetRightMargin(0.03);
  c1->SetTickx();
  c1->SetTicky();
  c1->Draw();
  c1->cd();

  double xx1 = -10.;
  double xx2 = 380.;
  double yy1 = 0.;
  double yy2 = 1.8;
  TH1D *d0 = new TH1D("d0", "", 1, xx1, xx2);
  d0->SetMinimum(yy1);
  d0->SetMaximum(yy2);
  d0->GetXaxis()->SetNdivisions(505);
  d0->GetXaxis()->CenterTitle();
  d0->GetXaxis()->SetTitle("N_{part}");
  d0->GetXaxis()->SetTitleOffset(1.08);
  d0->GetXaxis()->SetTitleSize(0.068);
  d0->GetXaxis()->SetLabelSize(0.050);
  d0->GetXaxis()->SetTitleFont(42);
  d0->GetXaxis()->SetLabelFont(42);
  d0->GetYaxis()->SetNdivisions(505);
  d0->GetYaxis()->SetTitle("#sigma_{c#bar{c}}^{NN} (mb)");
  d0->GetYaxis()->SetTitleOffset(0.88);
  d0->GetYaxis()->SetTitleSize(0.068);
  d0->GetYaxis()->SetLabelSize(0.050);
  d0->GetYaxis()->SetTitleFont(42);
  d0->GetYaxis()->SetLabelFont(42);
  d0->Draw();

  TLine *l1 = new TLine(xx1, yy1, xx2, yy1);
  l1->SetLineWidth(3);
  l1->Draw("same");
  TLine *l2 = new TLine(xx1, yy2, xx2, yy2);
  l2->SetLineWidth(3);
  l2->Draw("same");
  TLine *l3 = new TLine(xx1, yy1, xx1, yy2);
  l3->SetLineWidth(3);
  l3->Draw("same");
  TLine *l4 = new TLine(xx2, yy1, xx2, yy2);
  l4->SetLineWidth(3);
  l4->Draw("same");

  TGraphErrors *grFit = new TGraphErrors(nFit, npartFit, sigFit, enpartFit, errFit);
  TF1 *fitConst = new TF1("fitConst", "pol0", xx1, xx2);
  fitConst->SetParameter(0, 0.7);
  grFit->Fit(fitConst, "Q0");
  fitConst->SetLineColor(kBlack);
  fitConst->SetLineStyle(7);
  fitConst->SetLineWidth(3);
  fitConst->Draw("same");

  TGraphAsymmErrors *grPhenAuAu = new TGraphAsymmErrors(nPhenAuAu, npartPhenAuAu, sigPhenAuAu,
                                                        enpartPhenAuAu, enpartPhenAuAu,
                                                        errPhenAuAu, errPhenAuAu);
  grPhenAuAu->SetMarkerStyle(phenixMarker);
  grPhenAuAu->SetMarkerSize(1.8);
  grPhenAuAu->SetMarkerColor(phenixColor);
  grPhenAuAu->SetLineColor(phenixColor);
  grPhenAuAu->SetLineWidth(2);
  grPhenAuAu->Draw("p");

  TGraphAsymmErrors *grPhenPP = new TGraphAsymmErrors(nPhenPP, npartPhenPP, sigPhenPP,
                                                      enpartPhenPP, enpartPhenPP,
                                                      errPhenPP, errPhenPP);
  grPhenPP->SetMarkerStyle(phenixMarker);
  grPhenPP->SetMarkerSize(1.8);
  grPhenPP->SetMarkerColor(phenixColor);
  grPhenPP->SetLineColor(phenixColor);
  grPhenPP->SetLineWidth(2);
  grPhenPP->Draw("p");

  TGraphAsymmErrors *grStarPP = new TGraphAsymmErrors(nStarPP, npartStarPP, sigStarPP,
                                                      enpartStarPP, enpartStarPP,
                                                      errStarPPd, errStarPPu);
  grStarPP->SetMarkerStyle(starMarker);
  grStarPP->SetMarkerSize(2.5);
  grStarPP->SetMarkerColor(starColor);
  grStarPP->SetLineColor(starColor);
  grStarPP->SetLineWidth(2);
  grStarPP->Draw("p");

  TGraphAsymmErrors *grStarDAu = new TGraphAsymmErrors(nStarDAu, npartStarDAu, sigStarDAu,
                                                       enpartStarDAu, enpartStarDAu,
                                                       errStarDAu, errStarDAu);
  grStarDAu->SetMarkerStyle(starMarker);
  grStarDAu->SetMarkerSize(2.5);
  grStarDAu->SetMarkerColor(starColor);
  grStarDAu->SetLineColor(starColor);
  grStarDAu->SetLineWidth(2);
  grStarDAu->Draw("p");

  TLatex *tex = new TLatex(0.22, 0.88, "#sqrt{s_{NN}} = 200 GeV");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.050);
  tex->Draw("same");

  TLegend *leg = new TLegend(0.6, 0.68, 0.96, 0.94);
  leg->SetFillColor(10);
  leg->SetLineStyle(4000);
  leg->SetLineColor(10);
  leg->SetLineWidth(0.);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->AddEntry(grStarPP, "  STAR p+p", "p");
  leg->AddEntry(grStarDAu, "  STAR d+Au", "p");
  leg->AddEntry(grPhenPP, "  PHENIX p+p", "p");
  leg->AddEntry(grPhenAuAu, "  PHENIX Au+Au", "p");
  leg->AddEntry(fitConst, Form("  avg. = %.2f #pm %.2f mb",
                               fitConst->GetParameter(0),
                               fitConst->GetParError(0)), "l");
  leg->Draw();

  printf("Constant fit to all data points: %.6f +/- %.6f mb, chi2/ndf = %.3f/%d\n",
         fitConst->GetParameter(0), fitConst->GetParError(0),
         fitConst->GetChisquare(), fitConst->GetNDF());

  c1->Update();
  c1->SaveAs("Xsec_ccbar_Npart.pdf");
  c1->SaveAs("Xsec_ccbar_Npart.png");
}
