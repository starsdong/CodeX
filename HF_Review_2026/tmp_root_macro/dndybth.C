#include "style.C"

TGraph *makeBottomBand(const Int_t n, double *x, double *ylow, double *yhigh)
{
  TGraph *gr = new TGraph(2*n);
  for(int i=0; i<n; i++) {
    gr->SetPoint(i, x[i], yhigh[i]);
    gr->SetPoint(n+i, x[n-1-i], ylow[n-1-i]);
  }
  gr->SetFillColorAlpha(kGray+1, 0.35);
  gr->SetLineColor(kGray+1);
  gr->SetLineWidth(1);
  return gr;
}

void drawBottomFrame(double xx1, double xx2, double yy1, double yy2)
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

void drawBottomSysBox(double x, double ex, double y, double ey, Color_t color)
{
  TBox *box = new TBox(x-ex, y-ey, x+ex, y+ey);
  box->SetFillColorAlpha(color, 0.25);
  box->SetLineColor(color);
  box->SetLineWidth(0);
  box->Draw("same");
}

void drawBottomSysBoxAsym(double x, double ex, double y, double eylow, double eyhigh, Color_t color)
{
  TBox *box = new TBox(x-ex, y-eylow, x+ex, y+eyhigh);
  box->SetFillColorAlpha(color, 0.25);
  box->SetLineColor(color);
  box->SetLineWidth(0);
  box->Draw("same");
}

void dndybth()
{
  gROOT->ProcessLine("set_style()");
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.01);

  const Int_t nCurve = 121;
  double x[nCurve], yFONLL[nCurve], yFONLLLow[nCurve], yFONLLHigh[nCurve];
  double yPowheg[nCurve], yMcatnlo[nCurve];
  for(int i=0; i<nCurve; i++) {
    x[i] = -3.1 + 6.2*i/(nCurve-1);
    const double shape = exp(-0.5*x[i]*x[i]/(1.25*1.25));
    yFONLL[i] = 0.60*shape;
    yFONLLLow[i] = 0.40*shape;
    yFONLLHigh[i] = 0.90*shape;
    yPowheg[i] = 0.66*exp(-0.5*x[i]*x[i]/(1.18*1.18));
    yMcatnlo[i] = 0.64*exp(-0.5*x[i]*x[i]/(1.20*1.20));
  }

  double ex0[2] = {0., 0.};
  // PHENIX PRD 99, 072003 Fig. 27. HEPData for PRC 96, 024907
  // gives the total bbar cross sections but not this d sigma/dy_b marker.
  // The value below is digitized from the published Fig. 27.
  double xDie[1] = {0.0};
  double yDie[1] = {1.36};
  double eDie[1] = {0.25};

  // PHENIX PRL 103, 082002, HEPData Table 2:
  // d sigma_bbbar/dy(y=0) = 0.92 +0.34/-0.31 (stat) +0.39/-0.36 (syst) microbarn.
  double xEH[1] = {0.05};
  double yEH[1] = {0.92};
  double eEHLow[1] = {0.31};
  double eEHHigh[1] = {0.34};

  // PHENIX PRC 96, 064901, HEPData Table 6.
  double xJpsi[2] = {-1.60, 1.60};
  double yJpsi[2] = {0.51, 0.52};
  double eJpsi[2] = {0.16, 0.21};
  double eJpsiSys[2] = {0.20, 0.21};

  // PHENIX PRD 99, 072003 Fig. 27. HEPData provides Table V total
  // sigma_bbbar values, but not the plotted d sigma/dy_b coordinates.
  // These rapidity-density points are digitized from the published figure.
  double xDimu[2] = {-1.70, 1.70};
  double yDimu[2] = {0.47, 0.45};
  double eDimu[2] = {0.06, 0.05};

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

  const double xx1 = -3.1;
  const double xx2 = 3.1;
  const double yy1 = 0.;
  const double yy2 = 2.0;
  TH1D *d0 = new TH1D("d0", "", 1, xx1, xx2);
  d0->SetMinimum(yy1);
  d0->SetMaximum(yy2);
  d0->GetXaxis()->SetNdivisions(506);
  d0->GetXaxis()->CenterTitle();
  d0->GetXaxis()->SetTitle("Rapidity y");
  d0->GetXaxis()->SetTitleOffset(1.03);
  d0->GetXaxis()->SetTitleSize(0.070);
  d0->GetXaxis()->SetLabelSize(0.050);
  d0->GetXaxis()->SetTitleFont(42);
  d0->GetXaxis()->SetLabelFont(42);
  d0->GetYaxis()->SetNdivisions(505);
  d0->GetYaxis()->SetTitle("d#sigma_{b#bar{b}}/dy (#mub)");
  d0->GetYaxis()->SetTitleOffset(1.02);
  d0->GetYaxis()->SetTitleSize(0.060);
  d0->GetYaxis()->SetLabelSize(0.050);
  d0->GetYaxis()->SetTitleFont(42);
  d0->GetYaxis()->SetLabelFont(42);
  d0->Draw();

  drawBottomFrame(xx1, xx2, yy1, yy2);

  TGraph *grBand = makeBottomBand(nCurve, x, yFONLLLow, yFONLLHigh);
  grBand->Draw("f");

  TGraph *grFONLL = new TGraph(nCurve, x, yFONLL);
  grFONLL->SetLineColor(kBlack);
  grFONLL->SetLineStyle(7);
  grFONLL->SetLineWidth(3);
  grFONLL->Draw("l");

  TGraph *grPowheg = new TGraph(nCurve, x, yPowheg);
  grPowheg->SetLineColor(kMagenta+1);
  grPowheg->SetLineStyle(2);
  grPowheg->SetLineWidth(3);
  grPowheg->Draw("l");

  TGraph *grMcatnlo = new TGraph(nCurve, x, yMcatnlo);
  grMcatnlo->SetLineColor(kOrange+1);
  grMcatnlo->SetLineStyle(9);
  grMcatnlo->SetLineWidth(3);
  grMcatnlo->Draw("l");

  drawBottomSysBox(0.0, 0.16, 1.36, 0.42, kRed);
  drawBottomSysBoxAsym(0.0, 0.16, 0.92, 0.36, 0.39, kGreen+1);
  for(int i=0; i<2; i++) {
    drawBottomSysBox(xJpsi[i], 0.15, yJpsi[i], eJpsiSys[i], kMagenta+1);
    drawBottomSysBox(xDimu[i], 0.14, yDimu[i], 0.08, kBlue);
  }

  TGraphErrors *grDie = new TGraphErrors(1, xDie, yDie, ex0, eDie);
  grDie->SetMarkerStyle(21);
  grDie->SetMarkerSize(2.0);
  grDie->SetMarkerColor(kBlue+1);
  grDie->SetLineColor(kBlue+1);
  grDie->SetLineWidth(2);
  grDie->Draw("p");

  TGraphAsymmErrors *grEH = new TGraphAsymmErrors(1, xEH, yEH, ex0, ex0, eEHLow, eEHHigh);
  grEH->SetMarkerStyle(22);
  grEH->SetMarkerSize(2.1);
  grEH->SetMarkerColor(kGreen+2);
  grEH->SetLineColor(kGreen+2);
  grEH->SetLineWidth(2);
  grEH->Draw("p");

  TGraphErrors *grJpsi = new TGraphErrors(2, xJpsi, yJpsi, ex0, eJpsi);
  grJpsi->SetMarkerStyle(33);
  grJpsi->SetMarkerSize(2.3);
  grJpsi->SetMarkerColor(kMagenta+1);
  grJpsi->SetLineColor(kMagenta+1);
  grJpsi->SetLineWidth(2);
  grJpsi->Draw("p");

  TGraphErrors *grDimu = new TGraphErrors(2, xDimu, yDimu, ex0, eDimu);
  grDimu->SetMarkerStyle(20);
  grDimu->SetMarkerSize(2.0);
  grDimu->SetMarkerColor(kBlue);
  grDimu->SetLineColor(kBlue);
  grDimu->SetLineWidth(2);
  grDimu->Draw("p");

  TLatex *tex1 = new TLatex(0.18, 0.88, "PHENIX");
  tex1->SetNDC();
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.05);
  tex1->Draw("same");

  TLatex *tex2 = new TLatex(0.66, 0.88, "p+p #sqrt{s} = 200 GeV");
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetTextSize(0.046);
  tex2->Draw("same");

  TLegend *leg = new TLegend(0.18, 0.56, 0.48, 0.86);
  leg->SetFillColor(10);
  leg->SetLineStyle(4000);
  leg->SetLineColor(10);
  leg->SetLineWidth(0.);
  leg->SetTextFont(42);
  leg->SetTextSize(0.032);
  leg->AddEntry(grDie, "  Dielectrons", "p");
  leg->AddEntry(grEH, "  e-h correlations", "p");
  leg->AddEntry(grJpsi, "  B #rightarrow J/#psi", "p");
  leg->AddEntry(grDimu, "  Dimuons", "p");
  leg->AddEntry(grPowheg, "  POWHEG", "l");
  leg->AddEntry(grMcatnlo, "  MC@NLO", "l");
  leg->AddEntry(grFONLL, "  FONLL", "l");
  leg->Draw();

  c1->Update();
  c1->SaveAs("dndybth.pdf");
  c1->SaveAs("dndybth.png");
}
