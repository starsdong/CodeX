#include "style.C"

void Xsec_pp() 
{
  gROOT->ProcessLine("set_style()");

  // data poionts of other experiments
  // float sdata[] = {21.87, 26.15, 27.65, 39.11, 39.11, 130.0, 1800.0};
  // float sedata[] = {0.,0.,0.,0.,0.,0.,0.};
  // float Xdata[] = {11.06, 18.57, 18.20, 28.22, 42.04, 424.7, 10731.};
  // float Xedata[] = {3.20, 10.0, 3.58, 15.5, 24.2, 262., 6.2e3};
  float sdata[] = {21.87, 26.15, 27.65, 39.11, 39.11, 1800.0};
  float sedata[] = {0.,0.,0.,0.,0.,0.};
  float Xdata[] = {11.06, 18.57, 18.20, 28.22, 42.04, 10731.};
  float Xedata[] = {3.20, 10.0, 3.58, 15.5, 24.2, 6.2e3};

  // data point of STAR
  // float sSTAR[] = {196.};
  // float XSTAR[] = {1400.};
  // float XeSTAR[] = {470.};
  float sSTAR[] = {195.};
  float XSTAR[] = {797.};
  float XeSTARu[] = {296};
  float XeSTARd[] = {362};

  // data point of PHENIX
  float sPHENIX[] = {130, 205.};
  float XPHENIX[] = {424.7, 709.};
  float XePHENIX[] = {262, 343.};

  // data points from Worhi's talk -- seems not correct
  float s1data[] = {21.87, 26.15, 27.65, 39.11, 39.11, 130.0, 200.0};
  float s1edata[] = {0.,0.,0.,0.,0.,0.,0.};
  float X1data[] = {11.06, 18.57, 18.20, 28.22, 42.04, 424.7, 1400.};
  float X1edata[] = {3.20, 10.0, 3.58, 15.5, 24.2, 262., 470.};

  // data points from FNAL PRL paper
  // pion beam
  float s2data[] = {200, 210, 230, 250, 360, 600};
  for(int i=0;i<6;i++) {
    s2data[i] = sqrt(2.*0.938*s2data[i]);
  }
  // xF>0 only D0 and Dch
  float X2data[] = {6.2, 8.0, 8.8, 11.7, 15.4, 29.5};
  float X2edata[] = {0.9, 1.0, 1.4, 1.0, 3.4, 5.4};
  //  const Double_t scale_pion = 1.6*1.24;
  const Double_t scale_pion = 1.6*1.5/2.0;
  for(int i=0;i<6;i++) {
    X2data[i] = X2data[i]*scale_pion;
     X2edata[i] = X2edata[i]*scale_pion;
  }

  // proton beam
  float s3data[] = {200, 250, 360, 400, 800, 800};
  for(int i=0;i<6;i++) {
    s3data[i] = sqrt(2.0*0.938*s3data[i]);
  }
  float X3data[] = {1.7, 8.7, 15.4, 14.8, 23.7, 34.5};
  float X3edata[] = {0.8, 1.6, 8.3, 2.0, 8.0, 9.7};
  //  const Double_t scale_pro = 2.0*1.24;
  const Double_t scale_pro = 2.0*1.5/2.0;
  for(int i=0;i<6;i++) {
    X3data[i] = X3data[i]*scale_pro;
    X3edata[i] = X3edata[i]*scale_pro;
  }


  // additional low energy data points
  float s3pdata[] = {200, 280, 350, 350, 360, 400};
  for(int i=0;i<6;i++) {
    s3pdata[i] = sqrt(2.0*0.938*s3pdata[i]);
  }
  float X3pdata[] = {3.9, 5.6, 22., 11.3, 24.6, 31.};
  float X3epdata[] = {2.5, 1.7, 9., 2.0, 12.0, 18.};
  for(int i=0;i<6;i++) {
    X3pdata[i] = X3pdata[i]*scale_pro;
    X3epdata[i] = X3epdata[i]*scale_pro;
  }

  //UA1, X-ray cosmic
  float s4data[] = {2.12e5, 5.5e4, 7.5e4};
  for(int i=0;i<3;i++) {
    s4data[i] = sqrt(2.0*0.938*s4data[i]);
  }
  float X4data[] = {696, 3300, 1700};
  float X4edata[] = {610, 1110, 530};

  // ISR ? from Romonar's paper
  float s5data[] = {57.5, 57.5, 57.5, 62, 62, 63, 63};
  float s5edata[] = {4.5, 4.5, 4.5, 0., 0., 0., 0};
  float X5data[] = {73., 70., 1390, 128, 650, 298, 840};
  float X5edata[] = {22, 38, 195, 75, 228, 153, 328};

  // LHC data
  // ALICE
  float sALICE[] = {2760, 7000*0.97};
  float XALICE[] = {4800, 8500};
  float XeALICEu[] = {2902, 5137};
  float XeALICEd[] = {1585, 2512};

  // ATLAS
  float sATLAS[] = {7000*1.03};
  float XATLAS[] = {8600};
  float XeATLASu[] = {3892};
  float XeATLASd[] = {3503};

  // Total bbar production cross sections in pp and pA.
  // pA values are normalized per target nucleon; energies are slightly offset
  // where measurements overlap so the points remain visible on the log axis.
  float sbbar_e0[] = {0., 0., 0., 0., 0., 0., 0.};
  float sbbar_PHENIX[] = {200.*0.97, 200.*1.03, 510.};
  float Xbbar_PHENIX[] = {3.20, 3.75, 13.1};
  float Xebbar_PHENIXd[] = {1.70, 0.62, 3.2};
  float Xebbar_PHENIXu[] = {1.84, 0.71, 3.2};

  float sbbar_ALICE[] = {2760., 5020.*0.97, 7000.*0.97, 13000.*0.97};
  float Xbbar_ALICE[] = {130., 218., 281., 502.};
  float Xebbar_ALICEd[] = {52.4, 49.1, 64.3, 53.7};
  float Xebbar_ALICEu[] = {45.2, 49.0, 63.4, 53.5};

  float sbbar_LHCb[] = {7000.*1.03, 13000.*1.03};
  float Xbbar_LHCb[] = {284., 495.};
  float Xebbar_LHCbd[] = {52.9, 52.0};
  float Xebbar_LHCbu[] = {52.9, 52.0};

  float sbbar_E789[] = {38.8*0.97};
  float Xbbar_E789[] = {0.0057};
  float Xebbar_E789d[] = {0.0020};
  float Xebbar_E789u[] = {0.0020};

  float sbbar_E771[] = {38.8*1.03};
  float Xbbar_E771[] = {0.0430};
  float Xebbar_E771d[] = {0.0184};
  float Xebbar_E771u[] = {0.0279};

  float sbbar_HERAB[] = {41.6};
  float Xbbar_HERAB[] = {0.0149};
  float Xebbar_HERABd[] = {0.0033};
  float Xebbar_HERABu[] = {0.0033};

  float sbbar_ALICEpA[] = {5020.*1.03};
  float Xbbar_ALICEpA[] = {162.5};
  float Xebbar_ALICEpAd[] = {19.1};
  float Xebbar_ALICEpAu[] = {19.1};
  
  
  
  ////////////
  // PYTHIA
  ///////////
//   //  low limit MSEL=0
//   float s1[14] = {18.5, 1.85e+01,2.59e+01,3.62e+01,5.06e+01,7.07e+01,9.89e+01,1.64e+02,3.20e+02,6.25e+02,1.22e+03,2.02e+03,3.96e+03, 3.96e3};
//   float x1[14] = {1.00, 2.51e+00,8.99e+00,2.86e+01,7.45e+01,1.69e+02,3.33e+02,8.01e+02,2.13e+03,5.02e+03,1.05e+04,1.73e+04,3.03e+04, 100.};

//   //  high limit MSEL=0
//   float s2[13] = {15.6,1.56e+01,2.19e+01,3.06e+01,4.28e+01,5.98e+01,8.36e+01,1.38e+02,2.70e+02,6.25e+02,1.45e+03,3.96e+03, 3.96e3};
//   float x2[13] = {1.00, 2.18e+00,7.66e+00,2.16e+01,5.75e+01,1.30e+02,2.73e+02,8.01e+02,2.55e+03,7.95e+03,1.91e+04,4.60e+04, 100.};
 
//   // low limit MSEL=4
//   float s3[12] = {15.6, 1.56e+01,2.19e+01,3.06e+01,5.06e+01,8.36e+01,1.64e+02,3.78e+02,8.75e+02,1.71e+03,3.87e+03, 3.87e3};
//   float x3[12] = {1.00, 2.41e+00,8.81e+00,2.30e+01,6.48e+01,1.50e+02,3.53e+02,7.70e+02,1.32e+03,1.85e+03,2.50e+03, 100.0};

//   // high limit MSEL=4
//   float s4[13] = {15.8, 1.58e+01,2.21e+01,3.09e+01,4.32e+01,7.15e+01,1.40e+02,2.74e+02,5.35e+02,1.05e+03,2.05e+03,3.87e+03, 3.87e3};
//   float x4[13] = {1.00, 4.13e+00,1.03e+01,2.49e+01,5.53e+01,1.47e+02,5.59e+02,1.37e+03,2.76e+03,4.92e+03,7.34e+03,1.03e+04, 100.};
  //
  // CTEQ6L 2002 MSEL=0
  //  float s1[13] = {1.79e+01,2.50e+01,3.50e+01,4.89e+01,6.84e+01,9.56e+01,1.34e+02,1.87e+02,3.66e+02,7.15e+02,1.40e+03,2.74e+03,3.83e+03};
  //  float x1[13] = {2.11e+00,7.88e+00,2.61e+01,7.20e+01,1.73e+02,3.77e+02,7.42e+02,1.38e+03,3.96e+03,9.15e+03,1.88e+04,3.21e+04,4.16e+04};

  // CTEQ6L 2002 MSEL=4
  //  float s2[11] = {1.58e+01,2.21e+01,3.09e+01,4.32e+01,7.15e+01,1.40e+02,2.74e+02,5.35e+02,1.05e+03,2.05e+03,3.78e+03};
  //  float x2[11] = {2.97e+00,9.24e+00,2.50e+01,5.24e+01,1.73e+02,5.61e+02,1.40e+03,2.77e+03,4.83e+03,7.50e+03,1.03e+04};
  // CTEQ5M1 MSEL=4
  float s1[10] = {10,20,40,80,130,200,500,1000, 1800, 5500};
  float x1[10] = {0.135, 3.05, 16.9, 52.9, 99.7, 160, 384, 688, 1083, 2307};

  // CTEQ5M1 MSEL=1
  float s2[10] = {10,20,40,80,130,200,500,1000, 1800, 5500};
  float x2[10] = {1.77, 12.4, 57.5, 171, 312, 498, 1082, 1712, 2349, 4320};

  // R. Vogt's prediction
  // MRST HO mu = m    m = 1.4
  float s5[13] = {10, 20, 40, 60, 80, 100, 200, 500, 1000, 2000, 5500, 9000, 140000};
  float x5[13] = {0.50, 7.81, 37.8, 72.6, 107, 142, 298, 669, 1120, 1766, 3175, 4110, 5114};

  // MRST HO mu = 2m    m=1.2
  float s6[13] = {10, 20, 40, 60, 80, 100, 200, 500, 1000, 2000, 5500, 9000, 140000};
  float x6[13] = {0.61, 8.60, 43.0, 85.3, 129.2, 173, 382, 923, 1667, 2873, 5833, 7921, 10235};

  // CTEQ mu = m   m = 1.4
  float s7[13] = {10, 20, 40, 60, 80, 100, 200, 500, 1000, 2000, 5500, 9000, 14000};
  float x7[13] = {0.65, 8.74, 40.3, 79.1, 120, 162, 366, 880, 1536, 2496, 4526, 5781, 7044};

  // CTEQ mu = 2m   m = 1.2
  float s8[13] = {10, 20, 40, 60, 80, 100, 200, 500, 1000, 2000, 5500, 9000, 14000};
  float x8[13] = {0.72, 9.20, 44.5, 89.8, 139, 190, 445, 1147, 2136, 3740, 7385, 9679, 12014};

  // GRV mu = m    m = 1.3
  float s9[13] = {10, 20, 40, 60, 80, 100, 200, 500, 1000, 2000, 5500, 9000, 14000};
  float x9[13] = {1.17, 9.49, 34.9, 64.8, 96.4, 129, 289, 721, 1317, 2266, 4593, 6286, 8210};

  // Romona's additional predictions
  // MRST HO mu_R = 0.5m, mu_F = 2m, m=1.2
  float ss6[10] = {10, 30., 63, 200, 500, 1000, 2000, 5500, 9000, 14000};
  float xx6[10] = {22.3, 824., 3075., 12640., 30970, 5.67e4, 9.90e4, 2.04e5, 2.80e5, 3.5e5};

  // MRST HO mu_R = 1m, mu_F = 2m, m=1.2
  float ss7[10] = {10, 30., 63, 200, 500, 1000, 2000, 5500, 9000, 14000};
  float xx7[10] = {1.98, 75.6, 287., 1188., 2886., 5243., 9085., 1.86e4, 2.53e4, 3.28e4};

  // Mangano 1993  mu_F = 2mu, mu_R = 0.5mu, mu = m_c = 1.2
  float s10[] = {0.101, 0.177, 0.349, 0.689, 1.357, 2.674, 6.605, 18.27, 20.61};
  float x10[] = {51.14, 135.0, 277., 479., 762., 1210., 1850., 2820., 2820.};
  for(int i=0;i<9;i++) {
    s10[i] = sqrt(2.0*0.938*s10[i]*1000.);
  }

  // FONLL calculation dsigma/dy
  // mu_R = mu_F = m_T  m = 1.2
  float sa[] = {20, 40, 80, 200, 500, 1000, 1960};
  float xa[] = {10.1, 26.1, 60.1, 118., 189., 262., 359.};
  float ra[] = {2.47, 2.83, 2.86, 3.75, 4.46, 4.77, 5.01};
  for(int i=0;i<7;i++) {
    xa[i] = xa[i]*ra[i];
  }

  // mu_R = 2m_T mu_F = m_T  m = 1.2
  float sb[] = {20, 40, 80, 200, 500, 1000, 1960};
  float xb[] = {11.0, 28.5, 75.6, 191., 354., 523., 750.};
  for(int i=0;i<7;i++) {
    xb[i] = xb[i]*ra[i];
  }

  // mu_R = 2m_T mu_F = m_T  m = 1.2
  float sc[] = {20, 40, 80, 200, 500, 1000, 1960};
  float xc[] = {8.3, 15.3, 26.3, 35.7, 46.8, 59.0, 76.4};
  for(int i=0;i<7;i++) {
    xc[i] = xc[i]*ra[i];
  }

   const Int_t charmColor = kBlue+2;
   const Int_t bottomColor = kMagenta+2;
   const Int_t lowEnergyMarker = 24;   // below RHIC
   const Int_t rhicMarker = 20;        // RHIC
   const Int_t highEnergyMarker = 26;  // above RHIC
   const Int_t starColor = kRed+1;
   const Int_t phenixColor = kBlue+1;
   const Int_t starMarker = 29;
   const Int_t phenixMarker = 21;

   TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetEndErrorSize(0.01);
   gStyle->SetGridColor(16);
   //   c1->SetLogy();
   c1->SetFillColor(10);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameFillColor(0);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   c1->SetLeftMargin(0.15);
   c1->SetBottomMargin(0.16);
   c1->SetTopMargin(0.03);
   c1->SetRightMargin(0.03);
   c1->SetLogx();
   c1->SetLogy();
   c1->SetTickx();
   c1->SetTicky();
   c1->Draw();
   c1->cd();
  
   double xx1 = 10;
   double xx2 = 20000;
   double yy1 = 1e-3;
   double yy2 = 5e4;
   TH1D *d0 = new TH1D("d0","",1,xx1, xx2);
   d0->SetMinimum(yy1);
   d0->SetMaximum(yy2);
   d0->GetXaxis()->SetNdivisions(201);
   d0->GetXaxis()->CenterTitle();
   d0->GetXaxis()->SetTitle("Collision Energy #sqrt{s} (GeV)");
   d0->GetXaxis()->SetTitleOffset(1.12);
   d0->GetXaxis()->SetTitleSize(0.052);
   d0->GetXaxis()->SetLabelOffset(999.);
   d0->GetXaxis()->SetLabelSize(0.0001);
   d0->GetXaxis()->SetTitleFont(42);
   d0->GetXaxis()->SetLabelFont(42);
   d0->GetYaxis()->SetNdivisions(402);
   d0->GetYaxis()->SetTitle("#sigma_{Q#bar{Q}}^{NN} (#mub)");
   d0->GetYaxis()->SetTitleOffset(1.15);
   d0->GetYaxis()->SetTitleSize(0.052);
   d0->GetYaxis()->SetLabelOffset(999.);
   d0->GetYaxis()->SetLabelSize(0.0001);
   d0->GetYaxis()->SetTitleFont(42);
   d0->GetYaxis()->SetLabelFont(42);
   d0->Draw();


   TLine *l1 = new TLine(xx1,yy1,xx2,yy1);
   l1->SetLineWidth(3);
   l1->Draw("same");
   TLine *l2 = new TLine(xx1,yy2,xx2,yy2);
   l2->SetLineWidth(3);
   l2->Draw("same");
   TLine *l3 = new TLine(xx1,yy1,xx1,yy2);
   l3->SetLineWidth(3);
   l3->Draw("same");
   TLine *l4 = new TLine(xx2,yy1,xx2,yy2);
   l4->SetLineWidth(3);
   l4->Draw("same");

   for(int i=-3;i<=4;i++) {
     TLatex *tex = new TLatex(8.4, pow(10,i), Form("10^{%d}",i));
     tex->SetTextFont(42);
     tex->SetTextAlign(32);
     tex->SetTextSize(0.040);
     tex->Draw("same");
   }

   for(int i=1;i<=4;i++) {
     double xpos = pow(10,i);
     if(i==1) xpos = 14.;
     if(i==4) xpos = 8500.;
     TLatex *tex = new TLatex(xpos, yy1*2.5,  Form("10^{%d}",i));
     tex->SetTextFont(42);
     tex->SetTextAlign(22);
     tex->SetTextSize(0.040);
     tex->Draw("same");
   }

   /*
   double f=1;
   TLine *la1 = new TLine(370*f, 0.9, 440*f, 0.9);
   la1->SetLineWidth(2);
   la1->Draw("same");
   TLine *la2 = new TLine(330*f, 0.55, 370*f, 0.9);
   la2->SetLineWidth(2);
   la2->Draw("same");
   TLine *la3 = new TLine(310*f, 0.70, 330*f, 0.55);
   la3->SetLineWidth(2);
   la3->Draw("same");
   TLine *la4 = new TLine(295*f, 0.60, 310*f, 0.70);
   la4->SetLineWidth(2);
   la4->Draw("same");
   */
  TGraph *gr2 = new TGraph(10, s2, x2);
  gr2->SetLineColor(4);
  //  gr2->SetFillColor(7);
  gr2->SetLineStyle(2);
  gr2->SetLineWidth(2);
  gr2->Draw("c");
  TGraph *gr1 = new TGraph(10, s1, x1);
  //  gr1->SetFillColor(10);
  gr1->SetLineColor(2);
  gr1->SetLineStyle(2);
  gr1->SetLineWidth(2);
  //  gr1->Draw("c");
  /*
  TGraph *gr4 = new TGraph(11, s4+1, x4+1);
  gr4->SetFillColor(8);
  gr4->SetLineColor(1);
  gr4->SetLineStyle(3);
  gr4->Draw("c");
  TGraph *gr3 = new TGraph(10, s3+1, x3+1);
  gr3->SetFillColor(10);
  gr3->SetLineColor(1);
  gr3->SetLineStyle(3);
  gr3->Draw("c");
  */
  /*
  TGraph *gr5 = new TGraph(13, s5, x5);
  gr5->SetLineColor(1);
  gr5->SetLineStyle(1);
  gr5->SetLineWidth(2);
  gr5->Draw("c");
  */
  TGraph *gr6 = new TGraph(13, s6, x6);
  gr6->SetLineColor(2);
  gr6->SetLineStyle(4);
  gr6->SetLineWidth(2);
  gr6->Draw("c");
  /*
  TGraph *gr7 = new TGraph(13, s7, x7);
  gr7->SetLineColor(2);
  gr7->SetLineStyle(4);
  gr7->SetLineWidth(2);
  gr7->Draw("c");

  TGraph *gr8 = new TGraph(13, s8, x8);
  gr8->SetLineColor(1);
  gr8->SetLineStyle(1);
  gr8->SetLineWidth(2);
  gr8->Draw("c");

  TGraph *gr9 = new TGraph(13, s9, x9);
  gr9->SetLineColor(5);
  gr9->SetLineStyle(5);
  gr9->Draw("c");
  */
  
  TGraph *grr6 = new TGraph(10, ss6, xx6);
  grr6->SetLineColor(2);
  grr6->SetLineStyle(4);
  grr6->SetLineWidth(2);
  //  grr6->Draw("c");

  TGraph *grr7 = new TGraph(10, ss7, xx7);
  grr7->SetLineColor(1);
  grr7->SetLineStyle(1);
  grr7->SetLineWidth(2);
  grr7->Draw("c");

  // FONLL
  TGraph *gra = new TGraph(7, sa, xa);
  gra->SetLineColor(2);
  gra->SetLineStyle(4);
  gra->SetLineWidth(2);
  //  gra->Draw("c");
  TGraph *grb = new TGraph(7, sb, xb);
  grb->SetLineColor(2);
  grb->SetLineStyle(4);
  grb->SetLineWidth(2);
  //  grb->Draw("c");
  TGraph *grc = new TGraph(7, sc, xc);
  grc->SetLineColor(2);
  grc->SetLineStyle(4);
  grc->SetLineWidth(2);
  //  grc->Draw("c");

  

  //  gStyle->SetLineWidth(2);
  TF1 *fitfun = new TF1("fitfun","[1]*pow(x,[0])");
  //  fitfun->SetParameters(1.78, 0.0985);
  //  fitfun->SetParameters(1.13, 0.845);
  fitfun->SetParameters(1.93, 0.0348);
  fitfun->SetLineStyle(3);
  fitfun->SetLineWidth(3);
  fitfun->SetLineColor(6);
  fitfun->SetRange(15,800);
  //  fitfun->Draw("same");


//   TGraphErrors *grdata1 = new TGraphErrors(5, sdata, Xdata, sedata, Xedata);
//   grdata1->SetMarkerSize(1.5);
//   grdata1->SetMarkerStyle(24);
//   grdata1->SetMarkerColor(1);
//   grdata1->SetLineColor(1);
//   grdata1->Draw("p");
  TGraphErrors *grpion = new TGraphErrors(6, s2data, X2data, sedata, X2edata);
  grpion->SetMarkerSize(1.5);
  grpion->SetMarkerStyle(lowEnergyMarker);
  grpion->SetMarkerColor(charmColor);
  grpion->SetLineColor(charmColor);
  grpion->SetLineWidth(2);
  //  grpion->Draw("p");

  TGraphErrors *grpro = new TGraphErrors(6, s3data, X3data, sedata, X3edata);
  grpro->SetMarkerSize(1.7);
  grpro->SetMarkerStyle(lowEnergyMarker);
  grpro->SetMarkerColor(charmColor);
  grpro->SetLineColor(charmColor);
  grpro->SetLineWidth(2);
  grpro->Draw("p");

  TGraphErrors *grpro_p = new TGraphErrors(6, s3pdata, X3pdata, sedata, X3epdata);
  grpro_p->SetMarkerSize(1.7);
  grpro_p->SetMarkerStyle(lowEnergyMarker);
  grpro_p->SetMarkerColor(charmColor);
  grpro_p->SetLineColor(charmColor);
  grpro_p->SetLineWidth(2);
  grpro_p->Draw("p");

  TGraphErrors *grISR = new TGraphErrors(7, s5data, X5data, s5edata, X5edata);
  grISR->SetMarkerSize(1.5);
  grISR->SetMarkerStyle(highEnergyMarker);
  grISR->SetMarkerColor(charmColor);
  grISR->SetLineColor(charmColor);
  grISR->SetLineWidth(2);
  //  grISR->Draw("p");

  TGraphErrors *grUA1 = new TGraphErrors(1, s4data, X4data, sedata, X4edata);
  grUA1->SetMarkerSize(1.5);
  grUA1->SetMarkerStyle(highEnergyMarker);
  grUA1->SetMarkerColor(charmColor);
  grUA1->SetLineColor(charmColor);
  grUA1->SetLineWidth(2);
  grUA1->Draw("p");

  // TLine *lb1 = new TLine((580),(1900),(680),(1900));
  // lb1->SetLineWidth(2);
  // lb1->Draw("same");
  // TLine *lb2 = new TLine((580),(1600),(580),(1900));
  // lb2->SetLineWidth(1);
  // lb2->Draw("same");
  // TLine *lb3 = new TLine((680),(1600),(680),(1900));
  // lb3->SetLineWidth(1);
  // lb3->Draw("same");


  TGraphErrors *grHi = new TGraphErrors(2, s4data+1, X4data+1, sedata, X4edata+1);
  grHi->SetMarkerSize(1.9);
  grHi->SetMarkerStyle(highEnergyMarker);
  grHi->SetMarkerColor(charmColor);
  grHi->SetLineColor(charmColor);
  grHi->SetLineWidth(2);
  grHi->Draw("p");
  

  TGraphErrors *grdata4 = new TGraphErrors(1, sdata+5, Xdata+5, sedata+5, Xedata+5);
  grdata4->SetMarkerSize(1.5);
  grdata4->SetMarkerStyle(highEnergyMarker);
  grdata4->SetMarkerColor(charmColor);
  grdata4->SetLineColor(charmColor);
  grdata4->SetLineWidth(2);
  //  grdata4->Draw("p");

//   TGraphErrors *grdata5 = new TGraphErrors(1, sdata+6, Xdata+6, sedata+6, Xedata+6);
//   grdata5->SetMarkerSize(1.5);
//   grdata5->SetMarkerStyle(24);
//   grdata5->SetMarkerColor(6);
//   grdata5->SetLineColor(6);
//   grdata5->Draw("p");
  
  
  TGraphErrors *grphenix1 = new TGraphErrors(2, sPHENIX, XPHENIX, sedata, XePHENIX);
  grphenix1->SetMarkerSize(1.8);
  grphenix1->SetMarkerStyle(phenixMarker);
  grphenix1->SetMarkerColor(phenixColor);
  grphenix1->SetLineColor(phenixColor);
  grphenix1->SetLineWidth(2);
  grphenix1->Draw("p");

//  TGraphErrors *grstar1 = new TGraphErrors(1, sSTAR, XSTAR, sedata, XeSTAR);
  TGraphAsymmErrors *grstar1 = new TGraphAsymmErrors(1, sSTAR, XSTAR, sedata, sedata, XeSTARd, XeSTARu);
  grstar1->SetMarkerSize(2.5);
  grstar1->SetMarkerStyle(starMarker);
  grstar1->SetMarkerColor(starColor);
  grstar1->SetLineColor(starColor);
  grstar1->SetLineWidth(2);
  grstar1->Draw("p");
  
  //  LHC
  TGraphAsymmErrors *gralice = new TGraphAsymmErrors(2, sALICE, XALICE, 0, 0, XeALICEd, XeALICEu);
  gralice->SetMarkerSize(1.8);
  gralice->SetMarkerStyle(highEnergyMarker);
  gralice->SetMarkerColor(charmColor);
  gralice->SetLineColor(charmColor);
  gralice->SetLineWidth(2);
  gralice->Draw("p");
  
  TGraphAsymmErrors *gratlas = new TGraphAsymmErrors(1, sATLAS, XATLAS, 0, 0, XeATLASd, XeATLASu);
  gratlas->SetMarkerSize(1.8);
  gratlas->SetMarkerStyle(highEnergyMarker);
  gratlas->SetMarkerColor(charmColor);
  gratlas->SetLineColor(charmColor);
  gratlas->SetLineWidth(2);
  gratlas->Draw("p");

  TGraphAsymmErrors *grbbar_PHENIX = new TGraphAsymmErrors(3, sbbar_PHENIX, Xbbar_PHENIX, sbbar_e0, sbbar_e0, Xebbar_PHENIXd, Xebbar_PHENIXu);
  grbbar_PHENIX->SetMarkerSize(1.8);
  grbbar_PHENIX->SetMarkerStyle(phenixMarker);
  grbbar_PHENIX->SetMarkerColor(phenixColor);
  grbbar_PHENIX->SetLineColor(phenixColor);
  grbbar_PHENIX->SetLineWidth(2);
  grbbar_PHENIX->Draw("p");

  TGraphAsymmErrors *grbbar_ALICE = new TGraphAsymmErrors(4, sbbar_ALICE, Xbbar_ALICE, sbbar_e0, sbbar_e0, Xebbar_ALICEd, Xebbar_ALICEu);
  grbbar_ALICE->SetMarkerSize(2.0);
  grbbar_ALICE->SetMarkerStyle(highEnergyMarker);
  grbbar_ALICE->SetMarkerColor(bottomColor);
  grbbar_ALICE->SetLineColor(bottomColor);
  grbbar_ALICE->SetLineWidth(2);
  grbbar_ALICE->Draw("p");

  TGraphAsymmErrors *grbbar_LHCb = new TGraphAsymmErrors(2, sbbar_LHCb, Xbbar_LHCb, sbbar_e0, sbbar_e0, Xebbar_LHCbd, Xebbar_LHCbu);
  grbbar_LHCb->SetMarkerSize(2.0);
  grbbar_LHCb->SetMarkerStyle(highEnergyMarker);
  grbbar_LHCb->SetMarkerColor(bottomColor);
  grbbar_LHCb->SetLineColor(bottomColor);
  grbbar_LHCb->SetLineWidth(2);
  grbbar_LHCb->Draw("p");

  TGraphAsymmErrors *grbbar_E789 = new TGraphAsymmErrors(1, sbbar_E789, Xbbar_E789, sbbar_e0, sbbar_e0, Xebbar_E789d, Xebbar_E789u);
  grbbar_E789->SetMarkerSize(1.8);
  grbbar_E789->SetMarkerStyle(lowEnergyMarker);
  grbbar_E789->SetMarkerColor(bottomColor);
  grbbar_E789->SetLineColor(bottomColor);
  grbbar_E789->SetLineWidth(2);
  grbbar_E789->Draw("p");

  TGraphAsymmErrors *grbbar_E771 = new TGraphAsymmErrors(1, sbbar_E771, Xbbar_E771, sbbar_e0, sbbar_e0, Xebbar_E771d, Xebbar_E771u);
  grbbar_E771->SetMarkerSize(1.8);
  grbbar_E771->SetMarkerStyle(lowEnergyMarker);
  grbbar_E771->SetMarkerColor(bottomColor);
  grbbar_E771->SetLineColor(bottomColor);
  grbbar_E771->SetLineWidth(2);
  grbbar_E771->Draw("p");

  TGraphAsymmErrors *grbbar_HERAB = new TGraphAsymmErrors(1, sbbar_HERAB, Xbbar_HERAB, sbbar_e0, sbbar_e0, Xebbar_HERABd, Xebbar_HERABu);
  grbbar_HERAB->SetMarkerSize(1.8);
  grbbar_HERAB->SetMarkerStyle(lowEnergyMarker);
  grbbar_HERAB->SetMarkerColor(bottomColor);
  grbbar_HERAB->SetLineColor(bottomColor);
  grbbar_HERAB->SetLineWidth(2);
  grbbar_HERAB->Draw("p");

  TGraphAsymmErrors *grbbar_ALICEpA = new TGraphAsymmErrors(1, sbbar_ALICEpA, Xbbar_ALICEpA, sbbar_e0, sbbar_e0, Xebbar_ALICEpAd, Xebbar_ALICEpAu);
  grbbar_ALICEpA->SetMarkerSize(1.8);
  grbbar_ALICEpA->SetMarkerStyle(highEnergyMarker);
  grbbar_ALICEpA->SetMarkerColor(bottomColor);
  grbbar_ALICEpA->SetLineColor(bottomColor);
  grbbar_ALICEpA->SetLineWidth(2);
  grbbar_ALICEpA->Draw("p");
  
  //   TGraphErrors *grall = new TGraphErrors(7, s1data, X1data, s1edata, X1edata);
  //   TF1 *fun1 = new TF1("fun1","[0]*pow(x,[1])");
  //   fun1->SetParameters(0.02, 2.);
  //   fun1->SetRange(1,300);
  //   grall->Fit("fun1","R");
  /*
    TLatex *tex = new TLatex(log10(700), 2.5, "UA2");
    tex->SetTextSize(0.035);
    tex->Draw("same");
  */
  // TLatex *tex1 = new TLatex(300, pow(10.,3.7), "Pamir");
  // tex1->SetTextSize(0.045);
  // tex1->Draw("same");
  
  // TLatex *tex2 = new TLatex(400, pow(10, 3.2), "Muon");
  // tex2->SetTextSize(0.045);
  // tex2->Draw("same");

  TGraph *legCharmLow = new TGraph();
  legCharmLow->SetMarkerStyle(lowEnergyMarker);
  legCharmLow->SetMarkerColor(charmColor);
  legCharmLow->SetLineColor(charmColor);
  legCharmLow->SetMarkerSize(1.7);
  TGraph *legBottomLow = new TGraph();
  legBottomLow->SetMarkerStyle(lowEnergyMarker);
  legBottomLow->SetMarkerColor(bottomColor);
  legBottomLow->SetLineColor(bottomColor);
  legBottomLow->SetMarkerSize(1.7);
  TGraph *legCharmRHIC = new TGraph();
  legCharmRHIC->SetMarkerStyle(starMarker);
  legCharmRHIC->SetMarkerColor(starColor);
  legCharmRHIC->SetLineColor(starColor);
  legCharmRHIC->SetMarkerSize(2.0);
  TGraph *legBottomRHIC = new TGraph();
  legBottomRHIC->SetMarkerStyle(phenixMarker);
  legBottomRHIC->SetMarkerColor(phenixColor);
  legBottomRHIC->SetLineColor(phenixColor);
  legBottomRHIC->SetMarkerSize(1.7);
  TGraph *legCharmHigh = new TGraph();
  legCharmHigh->SetMarkerStyle(highEnergyMarker);
  legCharmHigh->SetMarkerColor(charmColor);
  legCharmHigh->SetLineColor(charmColor);
  legCharmHigh->SetMarkerSize(1.7);
  TGraph *legBottomHigh = new TGraph();
  legBottomHigh->SetMarkerStyle(highEnergyMarker);
  legBottomHigh->SetMarkerColor(bottomColor);
  legBottomHigh->SetLineColor(bottomColor);
  legBottomHigh->SetMarkerSize(1.7);

  TLegend *leg1 = new TLegend(0.50, 0.25, 0.96, 0.47);
  leg1->SetFillColor(10);
  leg1->SetLineStyle(4000);
  leg1->SetLineColor(10);
  leg1->SetLineWidth(0.);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.030);
  leg1->SetNColumns(2);
  leg1->AddEntry(legCharmLow, "  c#bar{c}: below RHIC","p");
  leg1->AddEntry(legBottomLow, "  b#bar{b}: below RHIC","p");
  leg1->AddEntry(legCharmRHIC, "  STAR","p");
  leg1->AddEntry(legBottomRHIC, "  PHENIX","p");
  leg1->AddEntry(legCharmHigh, "  c#bar{c}: above RHIC","p");
  leg1->AddEntry(legBottomHigh, "  b#bar{b}: above RHIC","p");
  leg1->Draw();
  


  TGraph *gr10 = new TGraph(9, s10, x10);
  gr10->SetLineColor(2);
  gr10->SetLineStyle(4);
  gr10->SetLineWidth(2);

  TLegend *leg2 = new TLegend(0.50, 0.78, 0.94, 0.94);
  leg2->SetFillColor(10);
  leg2->SetLineStyle(4000);
  leg2->SetLineColor(10);
  leg2->SetLineWidth(0.);
  leg2->SetTextFont(12);
  leg2->SetTextSize(0.038);
  //  leg->AddEntry(fitfun, " (   s )^{n} (n=1.9#pm0.2) ", "l");
//   leg->AddEntry(gr1, " MSEL=4 CTEQ5M1", "l");
//   leg->AddEntry(gr2, " MSEL=1 CTEQ5M1", "l");
//   leg->AddEntry(gr8, " NLO CTEQ #mu=2m m=1.2", "l");
  //  leg->AddEntry(gr2, " PYTHIA II", "l");
  leg2->AddEntry(grr7, " NLO pQCD (#mu_{R}=m_{c})", "l");
  leg2->AddEntry(gr6, " NLO pQCD (#mu_{R}=2m_{c})", "l");
  leg2->AddEntry(gr2, " PYTHIA", "l");
  leg2->Draw();

/*   TLine *l1 = new TLine(620, 45, 750, 45); */
/*   l1->SetLineWidth(1); */
/*   l1->Draw("same"); */
/*   TLine *l1 = new TLine(570, 32, 620, 45); */
/*   l1->SetLineWidth(1); */
/*   l1->Draw("same"); */
/*   TLine *l1 = new TLine(555, 38, 570, 32); */
/*   l1->SetLineWidth(1); */
/*   l1->Draw("same"); */
/*   TLine *l1 = new TLine(540, 34, 555, 38); */
/*   l1->SetLineWidth(1); */
/*   l1->Draw("same"); */

  //  p1->Modified();

 
 c1->Update();
 c1->SaveAs("Xsec_pp.pdf");
 c1->SaveAs("Xsec_pp.png");
}
