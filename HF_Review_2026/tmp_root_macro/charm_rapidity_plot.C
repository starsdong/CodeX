#include "style.C"

TGraph *makeCharmBand(const Int_t n, double *x, double *ylow, double *yhigh)
{
  TGraph *gr = new TGraph(2*n);
  for(int i=0; i<n; i++) {
    gr->SetPoint(i, x[i], yhigh[i]);
    gr->SetPoint(n+i, x[n-1-i], ylow[n-1-i]);
  }
  gr->SetFillColorAlpha(kGray+1, 0.25);
  gr->SetLineColor(kGray+1);
  gr->SetLineWidth(1);
  return gr;
}

void drawCharmFrame(double xx1, double xx2, double yy1, double yy2)
{
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
}

void charm_rapidity_plot()
{
  gROOT->ProcessLine("set_style()");
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.01);

  // FONLL charm rapidity distribution digitized from PHENIX PRC 86, 024909 Fig. 12.
  // The curves are based on the FONLL calculation from Cacciari, Nason, and Vogt,
  // Phys. Rev. Lett. 95, 122001 (2005), quoted as Ref. [39] in PRC 86, 024909.
  const Int_t nCurve = 13;
  double x[nCurve] = {-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
  double yCentral[nCurve] = {0.0100, 0.0186, 0.0303, 0.0431, 0.0543, 0.0619, 0.0646, 0.0619, 0.0543, 0.0431, 0.0303, 0.0186, 0.0100};
  double yLow[nCurve] = {0.0048, 0.0065, 0.0111, 0.0155, 0.0191, 0.0209, 0.0217, 0.0209, 0.0191, 0.0155, 0.0111, 0.0065, 0.0048};
  double yHigh[nCurve] = {0.0280, 0.0512, 0.0801, 0.1113, 0.1375, 0.1551, 0.1614, 0.1551, 0.1375, 0.1113, 0.0801, 0.0512, 0.0280};

  double ex0[2] = {0., 0.};
  double yStarX[1] = {0.0};
  // STAR PRL 94, 062301: d sigma_ccbar/dy(y=0) = 0.30 +/- 0.04 +/- 0.09 mb.
  double yStar[1] = {0.300};
  double yStarErr[1] = {0.0985};

  double yStarNewX[1] = {0.0};
  // STAR arXiv:1204.4244 / PRL 108, 202301:
  // d sigma_ccbar/dy(y=0) = 170 +/- 45(stat) +38/-59(syst) microbarn.
  double yStarNew[1] = {0.170};
  double yStarNewErrLow[1] = {0.0742};
  double yStarNewErrHigh[1] = {0.0597};

  double yPhenixSingleEX[1] = {0.0};
  // PHENIX PRL 97, 252002, HEPData Table 2:
  // d sigma_ccbar/dy = 123 +/- 12(stat) +/- 45(syst) microbarn at |y| < 0.35.
  double yPhenixSingleE[1] = {0.123};
  double yPhenixSingleEErr[1] = {0.0466};

  double yMuPrdX[2] = {-1.60, 1.60};
  // PHENIX PRD 76, 092002: d sigma_ccbar/dy(y=1.6) =
  // 0.243 +/- 0.013 +/- 0.105 -0.087/+0.049 mb.
  double yMuPrd[2] = {0.243, 0.243};
  double yMuPrdErrLow[2] = {0.137, 0.137};
  double yMuPrdErrHigh[2] = {0.117, 0.117};

  double yMuPrcX[2] = {-1.65, 1.65};
  // PHENIX PRC 86, 024909: d sigma_ccbar/dy(y=1.65) =
  // 0.139 +/- 0.029 -0.058/+0.051 mb.
  double yMuPrc[2] = {0.139, 0.139};
  double yMuPrcErrLow[2] = {0.0649, 0.0649};
  double yMuPrcErrHigh[2] = {0.0587, 0.0587};

  TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.14);
  c1->SetBottomMargin(0.17);
  c1->SetTopMargin(0.03);
  c1->SetRightMargin(0.03);
  c1->SetTickx();
  c1->SetTicky();
  c1->Draw();
  c1->cd();

  const double xx1 = -3.0;
  const double xx2 = 3.0;
  const double yy1 = 0.;
  const double yy2 = 0.45;
  TH1D *d0 = new TH1D("d0", "", 1, xx1, xx2);
  d0->SetMinimum(yy1);
  d0->SetMaximum(yy2);
  d0->GetXaxis()->SetNdivisions(505);
  d0->GetXaxis()->CenterTitle();
  d0->GetXaxis()->SetTitle("Rapidity y");
  d0->GetXaxis()->SetTitleOffset(1.10);
  d0->GetXaxis()->SetTitleSize(0.060);
  d0->GetXaxis()->SetLabelSize(0.045);
  d0->GetXaxis()->SetTitleFont(42);
  d0->GetXaxis()->SetLabelFont(42);
  d0->GetYaxis()->SetNdivisions(505);
  d0->GetYaxis()->SetTitle("d#sigma_{c#bar{c}}/dy (mb)");
  d0->GetYaxis()->SetTitleOffset(1.05);
  d0->GetYaxis()->SetTitleSize(0.060);
  d0->GetYaxis()->SetLabelSize(0.045);
  d0->GetYaxis()->SetTitleFont(42);
  d0->GetYaxis()->SetLabelFont(42);
  d0->Draw();

  drawCharmFrame(xx1, xx2, yy1, yy2);

  TGraph *grBand = makeCharmBand(nCurve, x, yLow, yHigh);
  grBand->Draw("f");

  TGraph *grHigh = new TGraph(nCurve, x, yHigh);
  grHigh->SetLineColor(kBlack);
  grHigh->SetLineStyle(7);
  grHigh->SetLineWidth(2);
  grHigh->Draw("c");

  TGraph *grLow = new TGraph(nCurve, x, yLow);
  grLow->SetLineColor(kBlack);
  grLow->SetLineStyle(7);
  grLow->SetLineWidth(2);
  grLow->Draw("c");

  TGraph *grCentral = new TGraph(nCurve, x, yCentral);
  grCentral->SetLineColor(kBlack);
  grCentral->SetLineStyle(1);
  grCentral->SetLineWidth(3);
  grCentral->Draw("c");

  TGraphErrors *grStar = new TGraphErrors(1, yStarX, yStar, ex0, yStarErr);
  grStar->SetMarkerStyle(21);
  grStar->SetMarkerSize(1.8);
  grStar->SetMarkerColor(kBlue+1);
  grStar->SetLineColor(kBlue+1);
  grStar->SetLineWidth(2);
  grStar->Draw("p");

  TGraphAsymmErrors *grStarNew = new TGraphAsymmErrors(1, yStarNewX, yStarNew, ex0, ex0, yStarNewErrLow, yStarNewErrHigh);
  grStarNew->SetMarkerStyle(25);
  grStarNew->SetMarkerSize(1.9);
  grStarNew->SetMarkerColor(kAzure+2);
  grStarNew->SetLineColor(kAzure+2);
  grStarNew->SetLineWidth(2);
  grStarNew->Draw("p");

  TGraphErrors *grPhenixSingleE = new TGraphErrors(1, yPhenixSingleEX, yPhenixSingleE, ex0, yPhenixSingleEErr);
  grPhenixSingleE->SetMarkerStyle(20);
  grPhenixSingleE->SetMarkerSize(1.9);
  grPhenixSingleE->SetMarkerColor(kRed);
  grPhenixSingleE->SetLineColor(kRed);
  grPhenixSingleE->SetLineWidth(2);
  grPhenixSingleE->Draw("p");

  TGraphAsymmErrors *grMuPrd = new TGraphAsymmErrors(2, yMuPrdX, yMuPrd, ex0, ex0, yMuPrdErrLow, yMuPrdErrHigh);
  grMuPrd->SetMarkerStyle(22);
  grMuPrd->SetMarkerSize(1.9);
  grMuPrd->SetMarkerColor(kGreen+2);
  grMuPrd->SetLineColor(kGreen+2);
  grMuPrd->SetLineWidth(2);
  grMuPrd->Draw("p");

  TGraphAsymmErrors *grMuPrc = new TGraphAsymmErrors(2, yMuPrcX, yMuPrc, ex0, ex0, yMuPrcErrLow, yMuPrcErrHigh);
  grMuPrc->SetMarkerStyle(23);
  grMuPrc->SetMarkerSize(1.9);
  grMuPrc->SetMarkerColor(kMagenta+2);
  grMuPrc->SetLineColor(kMagenta+2);
  grMuPrc->SetLineWidth(2);
  grMuPrc->Draw("p");

  TLatex *tex = new TLatex(0.18, 0.86, "p+p  #sqrt{s} = 200 GeV");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.046);
  tex->Draw("same");

  TLegend *leg = new TLegend(0.50, 0.62, 0.965, 0.91);
  leg->SetFillColor(10);
  leg->SetLineStyle(4000);
  leg->SetLineColor(10);
  leg->SetLineWidth(0.);
  leg->SetTextFont(42);
  leg->SetTextSize(0.030);
  leg->AddEntry(grBand, "  FONLL band", "f");
  leg->AddEntry(grCentral, "  FONLL central", "l");
  leg->AddEntry(grHigh, "  FONLL upper/lower", "l");
  leg->AddEntry(grStar, "  STAR D^{0}+e", "p");
  leg->AddEntry(grStarNew, "  STAR D^{0},D^{*}", "p");
  leg->AddEntry(grPhenixSingleE, "  PHENIX single e", "p");
  leg->AddEntry(grMuPrd, "  PHENIX #mu (PRD 76)", "p");
  leg->AddEntry(grMuPrc, "  PHENIX #mu (PRC 86)", "p");
  leg->Draw();

  c1->Update();
  c1->SaveAs("charm_rapidity_plot.pdf");
  c1->SaveAs("charm_rapidity_plot.png");
}
