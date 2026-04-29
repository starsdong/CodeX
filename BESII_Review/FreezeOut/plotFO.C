////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A. Andronic et al, Nature 561 (2018) 321-330, arXiv:1710.09425
// TCF = TlimCF /{1 + exp[2.60 − ln(√sNN)/0.45]}
// muB =a/(1 + 0.288√sNN), with √sNN in GeV and TlimCF = 158.4MeV and a = 1307.5 MeV
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// J. Cleymans et al, PRC 73 (2006) 034905
// T(µB) = a − bµB^2 − cµB^4, where a = 0.166 ± 0.002 GeV, b = 0.139 ± 0.016 GeV−1, and c = 0.053 ± 0.021 GeV−3
//  µB(√s) = d/(1 + e√s), with d = 1.308 ± 0.028 GeV and e = 0.273 ± 0.008 GeV−1
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "draw.C+"
#include "style.C+"
#include "TBox.h"
#include "TSystem.h"

#include <cmath>

void plotFO()
{
  style();

  const Double_t EMIN = 2.0;
  const Double_t EMAX = 5e3;
  const Double_t MUMIN = 0.0;
  const Double_t MUMAX = 900.0;
  const Double_t TMIN = 30.;
  const Double_t TMAX = 340.;

  TF1 *f1_T = new TF1("f1_T", "158.4/(1.+exp(2.60-log(x)/0.45))", EMIN, EMAX);
  TF1 *f1_mu = new TF1("f1_mu", "1307.5/(1.+0.288*x)", EMIN, EMAX);

  TF1 *f2_T_mu = new TF1("f2_T_mu", "0.166 - 0.139*x*x - 0.053*x*x*x*x", 0., 1.3);
  TF1 *f2_mu = new TF1("f2_mu", "1.308/(1.+0.273*x)", EMIN, EMAX);

  const Int_t NM = 2;
  const Int_t NP = 5000;
  double T[NM][NP], mu[NM][NP];
  double logs[NP], sNN[NP];
  const Double_t MINLOGS = log(EMIN);
  const Double_t MAXLOGS = log(EMAX);
  const Char_t *Label[NM] = {"Andronic 2018", "Cleymans 2006"};

  for (int i = 0; i < NP; i++) {
    logs[i] = MINLOGS + (MAXLOGS - MINLOGS) / NP * i;
    sNN[i] = exp(logs[i]);

    T[0][i] = f1_T->Eval(sNN[i]);
    mu[0][i] = f1_mu->Eval(sNN[i]);

    mu[1][i] = f2_mu->Eval(sNN[i]);
    T[1][i] = f2_T_mu->Eval(mu[1][i]);
    mu[1][i] *= 1000.;
    T[1][i] *= 1000.;
  }

  TGraph *gr[NM];
  TGraph *grLog[NM];
  const Int_t kcolor[NM] = {kRed + 1, kBlue + 1};
  for (int im = 0; im < NM; im++) {
    gr[im] = new TGraph(NP, mu[im], T[im]);
    grLog[im] = new TGraph(NP, mu[im], T[im]);
    gr[im]->SetLineWidth(1);
    gr[im]->SetLineColor(kcolor[im]);
    grLog[im]->SetLineWidth(1);
    grLog[im]->SetLineColor(kcolor[im]);
  }

  const int NFOPI = 2;
  double mu_fopi[NFOPI] = {808, 810};
  double t_fopi[NFOPI] = {54, 52};
  double emu_fopi[NFOPI] = {5, 13};
  double et_fopi[NFOPI] = {2, 2};

  const int NHADES = 1;
  double mu_hades[NHADES] = {748};
  double t_hades[NHADES] = {70};
  double emu_hades[NHADES] = {8};
  double et_hades[NHADES] = {3};

  const int NAGS = 3;
  double mu_ags[NAGS] = {558, 578, 554};
  double t_ags[NAGS] = {123, 119, 118};
  double emu_ags[NAGS] = {15, 15, 13};
  double et_ags[NAGS] = {5, 5, 3};

  const int NSPS = 11;
  double mu_sps[NSPS] = {472, 467, 406, 414, 367, 381, 382, 284, 294, 298, 246};
  double t_sps[NSPS] = {136, 131, 144, 140, 148, 143, 146, 155, 150, 154, 155};
  double emu_sps[NSPS] = {14, 13, 19, 16, 14, 9, 9, 15, 11, 10, 10};
  double et_sps[NSPS] = {5, 5, 5, 3, 5, 3, 3, 5, 5, 4, 3};

  const int NSTAR = 7;
  double s_star[NSTAR] = {7.7, 11.5, 19.6, 27.0, 39.0, 62.4, 200.0};
  double mu_star[NSTAR] = {406, 297, 193, 151, 106, 72, 27};
  double t_star[NSTAR] = {146, 151, 158, 161, 162, 166, 166};
  double emu_star[NSTAR] = {14, 13, 10, 10, 12, 15, 12};
  double et_star[NSTAR] = {3, 3, 3, 3, 4, 6, 3};

  const int NALICE = 2;
  double mu_alice[NALICE] = {1, 1};
  double t_alice[NALICE] = {156, 157};
  double emu_alice[NALICE] = {4, 4};
  double et_alice[NALICE] = {2, 2};

  const int NLMR = 2;
  double mu_lmr[NLMR] = {83.235268, 151.0};
  double t_lmr[NLMR] = {177.6, 165.39};
  double emu_lmr[NLMR] = {0., 0.};
  double et_lmr[NLMR] = {15.268, 20.274};
  double esys_lmr[NLMR] = {12.964, 20.622};

  const int NIMR = 2;
  double mu_imr[NIMR] = {83.235268, 151.0};
  double t_imr[NIMR] = {287.3, 273.57};
  double emu_imr[NIMR] = {0., 0.};
  double et_imr[NIMR] = {70.298, 65.119};
  double esys_imr[NIMR] = {34.869, 9.9861};

  const int NNA60LMR = 1;
  double mu_na60_lmr[NNA60LMR] = {246.0};
  double t_na60_lmr[NNA60LMR] = {172.0};
  double emu_na60_lmr[NNA60LMR] = {0.0};
  double et_na60_lmr[NNA60LMR] = {6.0};

  const int NNA60IMR = 1;
  double mu_na60_imr[NNA60IMR] = {246.0};
  double t_na60_imr[NNA60IMR] = {245.0};
  double emu_na60_imr[NNA60IMR] = {0.0};
  double et_na60_imr[NNA60IMR] = {17.0};

  const int NTKIN = 9;
  const int NTKINSTAR = 8;
  const int NTKINALICE = 1;
  double s_tkin[NTKIN] = {7.7, 11.5, 19.6, 27.0, 39.0, 62.4, 130.0, 200.0, 2760.0};
  double mu_tkin[NTKIN];
  double t_tkin[NTKIN] = {116.0, 118.0, 113.0, 117.0, 117.0, 98.7, 96.5, 89.0, 90.0};
  double emu_tkin[NTKIN] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double et_tkin[NTKIN] = {10.0, 11.0, 11.0, 11.0, 11.0, 10.2, 8.0, 12.0, 0.0};
  for (int i = 0; i < NTKIN; ++i) {
    mu_tkin[i] = f1_mu->Eval(s_tkin[i]);
  }
  double mu_tkin_star[NTKINSTAR], t_tkin_star[NTKINSTAR], emu_tkin_star[NTKINSTAR], et_tkin_star[NTKINSTAR];
  double mu_tkin_alice[NTKINALICE], t_tkin_alice[NTKINALICE], emu_tkin_alice[NTKINALICE], et_tkin_alice[NTKINALICE];
  for (int i = 0; i < NTKINSTAR; ++i) {
    mu_tkin_star[i] = mu_tkin[i];
    emu_tkin_star[i] = emu_tkin[i];
    // Use the STAR chemical-freezeout fit muB when a matching energy is listed.
    // The 130 GeV kinetic point has no corresponding STAR chemical point here.
    for (int j = 0; j < NSTAR; ++j) {
      if (std::fabs(s_tkin[i] - s_star[j]) < 0.05) {
        mu_tkin_star[i] = mu_star[j];
        emu_tkin_star[i] = emu_star[j];
        break;
      }
    }
    t_tkin_star[i] = t_tkin[i];
    et_tkin_star[i] = et_tkin[i];
  }
  mu_tkin_alice[0] = mu_tkin[NTKIN - 1];
  t_tkin_alice[0] = t_tkin[NTKIN - 1];
  emu_tkin_alice[0] = emu_tkin[NTKIN - 1];
  et_tkin_alice[0] = et_tkin[NTKIN - 1];

  const int NCPREF = 1;
  double mu_cp_ref[NCPREF] = {602.1};
  double t_cp_ref[NCPREF] = {114.3};
  double emu_cp_ref[NCPREF] = {0.0};
  double et_cp_ref[NCPREF] = {0.0};

  const int NCPFUNC = 3;
  double mu_cp_func[NCPFUNC] = {635.0, 610.0, 600.0};
  double t_cp_func[NCPFUNC] = {107.0, 109.0, 117.0};

  const int NCPHOLO = 2;
  double mu_cp_holo[NCPHOLO] = {589.0, 626.0};
  double t_cp_holo[NCPHOLO] = {104.0, 119.0};

  TGraphErrors *g_fopi = new TGraphErrors(NFOPI, mu_fopi, t_fopi, emu_fopi, et_fopi);
  TGraphErrors *g_hades = new TGraphErrors(NHADES, mu_hades, t_hades, emu_hades, et_hades);
  TGraphErrors *g_ags = new TGraphErrors(NAGS, mu_ags, t_ags, emu_ags, et_ags);
  TGraphErrors *g_sps = new TGraphErrors(NSPS, mu_sps, t_sps, emu_sps, et_sps);
  TGraphErrors *g_star = new TGraphErrors(NSTAR, mu_star, t_star, emu_star, et_star);
  TGraphErrors *g_alice = new TGraphErrors(NALICE, mu_alice, t_alice, emu_alice, et_alice);
  TGraphErrors *g_lmr = new TGraphErrors(NLMR, mu_lmr, t_lmr, emu_lmr, et_lmr);
  TGraphErrors *g_imr = new TGraphErrors(NIMR, mu_imr, t_imr, emu_imr, et_imr);
  TGraphErrors *g_na60_lmr = new TGraphErrors(NNA60LMR, mu_na60_lmr, t_na60_lmr, emu_na60_lmr, et_na60_lmr);
  TGraphErrors *g_na60_imr = new TGraphErrors(NNA60IMR, mu_na60_imr, t_na60_imr, emu_na60_imr, et_na60_imr);
  TGraphErrors *g_tkin_star = new TGraphErrors(NTKINSTAR, mu_tkin_star, t_tkin_star, emu_tkin_star, et_tkin_star);
  TGraphErrors *g_tkin_alice = new TGraphErrors(NTKINALICE, mu_tkin_alice, t_tkin_alice, emu_tkin_alice, et_tkin_alice);
  TGraphErrors *g_cp_ref = new TGraphErrors(NCPREF, mu_cp_ref, t_cp_ref, emu_cp_ref, et_cp_ref);
  TGraph *g_cp_func = new TGraph(NCPFUNC, mu_cp_func, t_cp_func);
  TGraph *g_cp_holo = new TGraph(NCPHOLO, mu_cp_holo, t_cp_holo);

  TGraphErrors *glog_fopi = new TGraphErrors(NFOPI, mu_fopi, t_fopi, emu_fopi, et_fopi);
  TGraphErrors *glog_hades = new TGraphErrors(NHADES, mu_hades, t_hades, emu_hades, et_hades);
  TGraphErrors *glog_ags = new TGraphErrors(NAGS, mu_ags, t_ags, emu_ags, et_ags);
  TGraphErrors *glog_sps = new TGraphErrors(NSPS, mu_sps, t_sps, emu_sps, et_sps);
  TGraphErrors *glog_star = new TGraphErrors(NSTAR, mu_star, t_star, emu_star, et_star);
  TGraphErrors *glog_alice = new TGraphErrors(NALICE, mu_alice, t_alice, emu_alice, et_alice);
  TGraphErrors *glog_lmr = new TGraphErrors(NLMR, mu_lmr, t_lmr, emu_lmr, et_lmr);
  TGraphErrors *glog_imr = new TGraphErrors(NIMR, mu_imr, t_imr, emu_imr, et_imr);
  TGraphErrors *glog_na60_lmr = new TGraphErrors(NNA60LMR, mu_na60_lmr, t_na60_lmr, emu_na60_lmr, et_na60_lmr);
  TGraphErrors *glog_na60_imr = new TGraphErrors(NNA60IMR, mu_na60_imr, t_na60_imr, emu_na60_imr, et_na60_imr);
  TGraphErrors *glog_tkin_star = new TGraphErrors(NTKINSTAR, mu_tkin_star, t_tkin_star, emu_tkin_star, et_tkin_star);
  TGraphErrors *glog_tkin_alice = new TGraphErrors(NTKINALICE, mu_tkin_alice, t_tkin_alice, emu_tkin_alice, et_tkin_alice);
  TGraphErrors *glog_cp_ref = new TGraphErrors(NCPREF, mu_cp_ref, t_cp_ref, emu_cp_ref, et_cp_ref);
  TGraph *glog_cp_func = new TGraph(NCPFUNC, mu_cp_func, t_cp_func);
  TGraph *glog_cp_holo = new TGraph(NCPHOLO, mu_cp_holo, t_cp_holo);

  TBox *box_lmr[NLMR];
  TBox *box_imr[NIMR];
  TBox *box_lmr_log[NLMR];
  TBox *box_imr_log[NIMR];
  const double kBoxHalfWidth = 7.0;
  for (int i = 0; i < NLMR; ++i) {
    box_lmr[i] = new TBox(mu_lmr[i] - kBoxHalfWidth, t_lmr[i] - esys_lmr[i], mu_lmr[i] + kBoxHalfWidth, t_lmr[i] + esys_lmr[i]);
    box_lmr_log[i] = new TBox(mu_lmr[i] - kBoxHalfWidth, t_lmr[i] - esys_lmr[i], mu_lmr[i] + kBoxHalfWidth, t_lmr[i] + esys_lmr[i]);
    box_lmr[i]->SetFillColorAlpha(kAzure - 3, 0.28);
    box_lmr[i]->SetLineColor(kAzure - 3);
    box_lmr[i]->SetLineWidth(1);
    box_lmr_log[i]->SetFillColorAlpha(kAzure - 3, 0.28);
    box_lmr_log[i]->SetLineColor(kAzure - 3);
    box_lmr_log[i]->SetLineWidth(1);
  }
  for (int i = 0; i < NIMR; ++i) {
    box_imr[i] = new TBox(mu_imr[i] - kBoxHalfWidth, t_imr[i] - esys_imr[i], mu_imr[i] + kBoxHalfWidth, t_imr[i] + esys_imr[i]);
    box_imr_log[i] = new TBox(mu_imr[i] - kBoxHalfWidth, t_imr[i] - esys_imr[i], mu_imr[i] + kBoxHalfWidth, t_imr[i] + esys_imr[i]);
    box_imr[i]->SetFillColorAlpha(kOrange + 1, 0.28);
    box_imr[i]->SetLineColor(kOrange + 1);
    box_imr[i]->SetLineWidth(1);
    box_imr_log[i]->SetFillColorAlpha(kOrange + 1, 0.28);
    box_imr_log[i]->SetLineColor(kOrange + 1);
    box_imr_log[i]->SetLineWidth(1);
  }

  g_fopi->SetLineWidth(2);
  g_hades->SetLineWidth(2);
  g_ags->SetLineWidth(2);
  g_sps->SetLineWidth(2);
  g_star->SetLineWidth(2);
  g_alice->SetLineWidth(2);
  g_lmr->SetLineWidth(2);
  g_imr->SetLineWidth(2);
  g_na60_lmr->SetLineWidth(2);
  g_na60_imr->SetLineWidth(2);
  g_tkin_star->SetLineWidth(2);
  g_tkin_alice->SetLineWidth(2);
  g_cp_ref->SetLineWidth(2);
  glog_fopi->SetLineWidth(2);
  glog_hades->SetLineWidth(2);
  glog_ags->SetLineWidth(2);
  glog_sps->SetLineWidth(2);
  glog_star->SetLineWidth(2);
  glog_alice->SetLineWidth(2);
  glog_lmr->SetLineWidth(2);
  glog_imr->SetLineWidth(2);
  glog_na60_lmr->SetLineWidth(2);
  glog_na60_imr->SetLineWidth(2);
  glog_tkin_star->SetLineWidth(2);
  glog_tkin_alice->SetLineWidth(2);
  glog_cp_ref->SetLineWidth(2);

  g_fopi->SetMarkerStyle(33);  g_fopi->SetMarkerSize(2.6); g_fopi->SetMarkerColor(kRed + 1);      g_fopi->SetLineColor(kRed + 1);
  g_hades->SetMarkerStyle(33); g_hades->SetMarkerSize(2.6); g_hades->SetMarkerColor(kRed + 1);      g_hades->SetLineColor(kRed + 1);
  g_ags->SetMarkerStyle(33);   g_ags->SetMarkerSize(2.6); g_ags->SetMarkerColor(kRed + 1);          g_ags->SetLineColor(kRed + 1);
  g_sps->SetMarkerStyle(33);   g_sps->SetMarkerSize(2.6); g_sps->SetMarkerColor(kRed + 1);          g_sps->SetLineColor(kRed + 1);
  g_star->SetMarkerStyle(29);  g_star->SetMarkerSize(2.8); g_star->SetMarkerColor(kRed + 2);        g_star->SetLineColor(kRed + 2);
  g_alice->SetMarkerStyle(33); g_alice->SetMarkerSize(2.6); g_alice->SetMarkerColor(kRed + 1);      g_alice->SetLineColor(kRed + 1);
  g_lmr->SetMarkerStyle(29);   g_lmr->SetMarkerSize(1.9); g_lmr->SetMarkerColor(kAzure - 3);         g_lmr->SetLineColor(kAzure - 3);
  g_imr->SetMarkerStyle(29);   g_imr->SetMarkerSize(1.9); g_imr->SetMarkerColor(kOrange + 1);       g_imr->SetLineColor(kOrange + 1);
  g_na60_lmr->SetMarkerStyle(21); g_na60_lmr->SetMarkerSize(1.7); g_na60_lmr->SetMarkerColor(kAzure - 3); g_na60_lmr->SetLineColor(kAzure - 3);
  g_na60_imr->SetMarkerStyle(22); g_na60_imr->SetMarkerSize(1.7); g_na60_imr->SetMarkerColor(kOrange + 1); g_na60_imr->SetLineColor(kOrange + 1);
  g_tkin_star->SetMarkerStyle(30); g_tkin_star->SetMarkerSize(2.0); g_tkin_star->SetMarkerColor(kGray + 2);         g_tkin_star->SetLineColor(kGray + 2);
  g_tkin_alice->SetMarkerStyle(24); g_tkin_alice->SetMarkerSize(1.8); g_tkin_alice->SetMarkerColor(kGray + 2);         g_tkin_alice->SetLineColor(kGray + 2);
  g_cp_ref->SetMarkerStyle(20); g_cp_ref->SetMarkerSize(2.0); g_cp_ref->SetMarkerColor(kBlack);      g_cp_ref->SetLineColor(kBlack);
  g_cp_func->SetMarkerStyle(20); g_cp_func->SetMarkerSize(2.0); g_cp_func->SetMarkerColor(kBlack);   g_cp_func->SetLineColor(kBlack);
  g_cp_holo->SetMarkerStyle(20); g_cp_holo->SetMarkerSize(2.0); g_cp_holo->SetMarkerColor(kBlack);   g_cp_holo->SetLineColor(kBlack);
  glog_fopi->SetMarkerStyle(33);  glog_fopi->SetMarkerSize(2.6); glog_fopi->SetMarkerColor(kRed + 1);     glog_fopi->SetLineColor(kRed + 1);
  glog_hades->SetMarkerStyle(33); glog_hades->SetMarkerSize(2.6); glog_hades->SetMarkerColor(kRed + 1);    glog_hades->SetLineColor(kRed + 1);
  glog_ags->SetMarkerStyle(33);   glog_ags->SetMarkerSize(2.6); glog_ags->SetMarkerColor(kRed + 1);      glog_ags->SetLineColor(kRed + 1);
  glog_sps->SetMarkerStyle(33);   glog_sps->SetMarkerSize(2.6); glog_sps->SetMarkerColor(kRed + 1);      glog_sps->SetLineColor(kRed + 1);
  glog_star->SetMarkerStyle(29);  glog_star->SetMarkerSize(2.8); glog_star->SetMarkerColor(kRed + 2);     glog_star->SetLineColor(kRed + 2);
  glog_alice->SetMarkerStyle(33); glog_alice->SetMarkerSize(2.6); glog_alice->SetMarkerColor(kRed + 1);   glog_alice->SetLineColor(kRed + 1);
  glog_lmr->SetMarkerStyle(29);   glog_lmr->SetMarkerSize(1.9); glog_lmr->SetMarkerColor(kAzure - 3);      glog_lmr->SetLineColor(kAzure - 3);
  glog_imr->SetMarkerStyle(29);   glog_imr->SetMarkerSize(1.9); glog_imr->SetMarkerColor(kOrange + 1);    glog_imr->SetLineColor(kOrange + 1);
  glog_na60_lmr->SetMarkerStyle(21); glog_na60_lmr->SetMarkerSize(1.7); glog_na60_lmr->SetMarkerColor(kAzure - 3); glog_na60_lmr->SetLineColor(kAzure - 3);
  glog_na60_imr->SetMarkerStyle(22); glog_na60_imr->SetMarkerSize(1.7); glog_na60_imr->SetMarkerColor(kOrange + 1); glog_na60_imr->SetLineColor(kOrange + 1);
  glog_tkin_star->SetMarkerStyle(30); glog_tkin_star->SetMarkerSize(2.0); glog_tkin_star->SetMarkerColor(kGray + 2);     glog_tkin_star->SetLineColor(kGray + 2);
  glog_tkin_alice->SetMarkerStyle(24); glog_tkin_alice->SetMarkerSize(1.8); glog_tkin_alice->SetMarkerColor(kGray + 2);     glog_tkin_alice->SetLineColor(kGray + 2);
  glog_cp_ref->SetMarkerStyle(20); glog_cp_ref->SetMarkerSize(2.0); glog_cp_ref->SetMarkerColor(kBlack);  glog_cp_ref->SetLineColor(kBlack);
  glog_cp_func->SetMarkerStyle(20); glog_cp_func->SetMarkerSize(2.0); glog_cp_func->SetMarkerColor(kBlack); glog_cp_func->SetLineColor(kBlack);
  glog_cp_holo->SetMarkerStyle(20); glog_cp_holo->SetMarkerSize(2.0); glog_cp_holo->SetMarkerColor(kBlack); glog_cp_holo->SetLineColor(kBlack);


  TCanvas *c1 = new TCanvas("c1", "", 900, 700);
  c1->SetLeftMargin(0.13);
  c1->SetBottomMargin(0.13);
  c1->Draw();

  TH1D *d0 = new TH1D("d0", "", 1, MUMIN, MUMAX);
  d0->SetMinimum(TMIN);
  d0->SetMaximum(TMAX);
  d0->GetXaxis()->CenterTitle();
  d0->GetXaxis()->SetTitle("Baryon Chemical Potential #mu_{B} (MeV)");
  d0->GetXaxis()->SetLabelSize(0.045);
  d0->GetXaxis()->SetTickLength(0.03);
  d0->GetXaxis()->SetTitleOffset(1.1);
  d0->GetXaxis()->SetTitleSize(0.055);
  d0->GetYaxis()->SetTitle("Temperature T (MeV)");
  d0->GetYaxis()->SetTitleOffset(1.0);
  d0->GetYaxis()->SetTitleSize(0.06);
  d0->GetYaxis()->SetLabelSize(0.045);
  d0->Draw();

  for (int im = 0; im < NM; im++) {
    grLog[im]->Draw("c same");
  }

  glog_fopi->Draw("P SAME");
  glog_hades->Draw("P SAME");
  glog_ags->Draw("P SAME");
  glog_sps->Draw("P SAME");
  glog_star->Draw("P SAME");
  glog_alice->Draw("P SAME");
  for (int i = 0; i < NLMR; ++i) box_lmr_log[i]->Draw("same");
  for (int i = 0; i < NIMR; ++i) box_imr_log[i]->Draw("same");
  glog_lmr->Draw("P SAME");
  glog_imr->Draw("P SAME");
  glog_na60_lmr->Draw("P SAME");
  glog_na60_imr->Draw("P SAME");
  glog_tkin_star->Draw("P SAME");
  glog_tkin_alice->Draw("P SAME");
  glog_cp_ref->Draw("P SAME");
  glog_cp_func->Draw("P SAME");
  glog_cp_holo->Draw("P SAME");

  drawText(550, 315, "Chemical", 42, 0.028);
  drawText(550, 300, "Freezeout", 42, 0.028);
  
  TLegend *legCurve = new TLegend(0.74, 0.84, 0.90, 0.94);
  legCurve->SetFillStyle(4000);
  legCurve->SetBorderSize(0);
  legCurve->SetTextSize(0.03);
  legCurve->AddEntry(gr[0], Label[0], "l");
  legCurve->AddEntry(gr[1], Label[1], "l");
  legCurve->Draw();

  TLegend *legData1 = new TLegend(0.67, 0.56, 0.68, 0.82);
  legData1->SetFillStyle(4000);
  legData1->SetBorderSize(0);
  legData1->SetTextSize(0.03);
  legData1->AddEntry(g_alice, " ", "p");
  legData1->AddEntry(g_tkin_alice, " ", "p");
  legData1->AddEntry(g_na60_lmr, " ", "p");
  legData1->AddEntry(g_na60_imr, " ", "p");
  legData1->AddEntry(" ", " ", "");
  legData1->Draw();

  TLegend *legData0 = new TLegend(0.68, 0.56, 0.94, 0.82);
  legData0->SetFillStyle(4000);
  legData0->SetBorderSize(0);
  legData0->SetTextSize(0.03);
  legData0->AddEntry(g_star, "Chemical Freezeout", "p");
  legData0->AddEntry(g_tkin_star, "Kinetic Freezeout", "p");
  legData0->AddEntry(g_lmr, "LMR(e^{+}e^{-})", "p");
  legData0->AddEntry(g_imr, "IMR(e^{+}e^{-})", "p");
  legData0->AddEntry(g_cp_ref, "CP (theory)", "p");
  legData0->Draw();

  drawHistBox(MUMIN, MUMAX, TMIN, TMAX);

  gSystem->mkdir("fig", kTRUE);
  c1->Update();
  c1->SaveAs("fig/FreezeOut_Comp.pdf");
  c1->SaveAs("fig/FreezeOut_Comp.png");

  /*
  const Double_t MUMINLOG = 15.0;
  const Double_t MUMAXLOG = 1000.0;
  const Double_t TMINLOG = 45.;

  TCanvas *c2 = new TCanvas("c2", "", 900, 700);
  c2->SetLeftMargin(0.13);
  c2->SetBottomMargin(0.13);
  c2->SetLogx(1);
  c2->SetLogy(1);
  c2->Draw();

  TH1D *d1 = new TH1D("d1", "", 1, MUMINLOG, MUMAXLOG);
  d1->SetMinimum(TMINLOG);
  d1->SetMaximum(TMAX);
  d1->GetXaxis()->CenterTitle();
  d1->GetXaxis()->SetTitle("#mu_{B} (MeV)");
  d1->GetXaxis()->SetLabelSize(0.045);
  d1->GetXaxis()->SetTickLength(0.03);
  d1->GetXaxis()->SetTitleOffset(1.1);
  d1->GetXaxis()->SetTitleSize(0.055);
  d1->GetYaxis()->SetTitle("T (MeV)");
  d1->GetYaxis()->SetTitleOffset(1.0);
  d1->GetYaxis()->SetTitleSize(0.06);
  d1->GetYaxis()->SetLabelSize(0.045);
  d1->Draw();

  for (int im = 0; im < NM; im++) {
    gr[im]->Draw("c same");
  }

  g_fopi->Draw("P SAME");
  g_hades->Draw("P SAME");
  g_ags->Draw("P SAME");
  g_sps->Draw("P SAME");
  g_star->Draw("P SAME");
  g_alice->Draw("P SAME");
  for (int i = 0; i < NLMR; ++i) box_lmr[i]->Draw("same");
  for (int i = 0; i < NIMR; ++i) box_imr[i]->Draw("same");
  g_lmr->Draw("P SAME");
  g_imr->Draw("P SAME");
  g_na60_lmr->Draw("P SAME");
  g_na60_imr->Draw("P SAME");
  g_tkin_star->Draw("P SAME");
  g_tkin_alice->Draw("P SAME");
  g_cp_ref->Draw("P SAME");
  g_cp_func->Draw("P SAME");
  g_cp_holo->Draw("P SAME");

  drawText(550, 315, "Chemical", 42, 0.028);
  drawText(550, 300, "Freezeout", 42, 0.028);
  
  TLegend *legCurve2 = new TLegend(0.74, 0.84, 0.90, 0.94);
  legCurve2->SetFillStyle(4000);
  legCurve2->SetBorderSize(0);
  legCurve2->SetTextSize(0.03);
  legCurve2->AddEntry(gr[0], Label[0], "l");
  legCurve2->AddEntry(gr[1], Label[1], "l");
  legCurve2->Draw();

  TLegend *legData3 = new TLegend(0.67, 0.56, 0.68, 0.82);
  legData3->SetFillStyle(4000);
  legData3->SetBorderSize(0);
  legData3->SetTextSize(0.03);
  legData3->AddEntry(g_alice, " ", "p");
  legData3->AddEntry(g_tkin_alice, " ", "p");
  legData3->AddEntry(g_na60_lmr, "", "p");
  legData3->AddEntry(g_na60_imr, "", "p");
  legData3->AddEntry(" ", " ", "");
  legData3->Draw();

  TLegend *legData2 = new TLegend(0.68, 0.56, 0.94, 0.82);
  legData2->SetFillStyle(4000);
  legData2->SetBorderSize(0);
  legData2->SetTextSize(0.03);
  legData2->AddEntry(g_star, "Chemical Freezeout", "p");
  legData2->AddEntry(g_tkin_star, "Kinetic Freezeout", "p");
  legData2->AddEntry(g_lmr, "LMR(e^{+}e^{-})", "p");
  legData2->AddEntry(g_imr, "IMR(e^{+}e^{-})", "p");
  legData2->AddEntry(g_cp_ref, "CP (theory)", "p");
  legData2->Draw();

  drawHistBox(MUMINLOG, MUMAXLOG, TMINLOG, TMAX);

  c2->Update();
  c2->SaveAs("fig/FreezeOut_Comp_logxy.pdf");
  c2->SaveAs("fig/FreezeOut_Comp_logxy.png");
  */
}
