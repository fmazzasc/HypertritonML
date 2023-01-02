#include "Riostream.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>

constexpr double kPreliminaryTau{253};
constexpr double kPreliminaryTauStat[2]{11, 11}; // first + then -
constexpr double kPreliminaryTauSyst[2]{6, 6};   // first + then -
constexpr bool kWOPrel{true};

void CollectionPlotNew()
{

  gStyle->SetOptStat(0);

  constexpr float kOffset = -0.5;
  constexpr int nMeasures{13};
  constexpr int nPublished{9};

  const Int_t N = nPublished + 1; // number of lifetime values
  Float_t point[N] = {9, 8, 7, 6, 5, 4, 3, 2, 1};
  Float_t lifetau[N] = { 232, 264, 246,  183, 181, 142, 242, 221, 253};
  Float_t err_y[N] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  Float_t err_x_low[N] = { 34, 52, 41,  32, 39, 21, 38, 15, 11};
  Float_t err_x_high[N] = { 45, 84, 62,  42, 54, 24, 34, 15, 11};
  Float_t errsyst_x[N] = {0, 0, 0, 37, 33, 29, 17, 19, 6};
  Float_t errx_sumq_h[N] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  Float_t errx_sumq_l[N] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  Float_t errsyst_y[N] = {0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

  double w[N] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  double v[N] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  double s[N] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  double sum_w = 0;
  double sum_v = 0;
  double chi2 = 0;

  const Int_t kMarkTyp = 20;              // marker type
  const Int_t kMarkCol = 1;               //...and color
  const Float_t kTitSize = 0.055;         // axis title size
  const Float_t kAsize = 0.85 * kTitSize; //...and label size
  const Float_t kToffset = 0.8;
  const Int_t kFont = 42;

  for (int i = 0; i < N; i++)
  {
    std::cout << "Lifetime: " << lifetau[i] << std::endl;
    s[i] = std::hypot(err_x_high[i], err_x_low[i]);
    errx_sumq_h[i] = std::hypot(err_x_high[i], errsyst_x[i]);
    errx_sumq_l[i] = std::hypot(err_x_low[i], errsyst_x[i]);
    std::cout << "Sum of errors: " << s[i] << std::endl;
    std::cout << "Sum of errors high: " << errx_sumq_h[i] << std::endl;
    std::cout << "Sum of errors low: " << errx_sumq_l[i] << std::endl;
  }
  TGraphAsymmErrors *gSpect3;
  gSpect3 = new TGraphAsymmErrors(N, point, lifetau, 0x0, 0x0, errx_sumq_l, errx_sumq_h);
  auto fitPtr = gSpect3->Fit("pol0", "SE");
  std::cout << "Average lifetime: " << fitPtr->Parameter(0) << fitPtr->LowerError(0) << " + " << fitPtr->UpperError(0) << std::endl;

  TCanvas *cfit = new TCanvas("cfit", "cfit");
  cfit->cd();
  gSpect3->Draw("ape");

  std::cout << "Chi2: " << fitPtr->Chi2() << std::endl;
  std::cout << "NDF: " << fitPtr->Ndf() << std::endl;

  cout << "Sum of weights after all it : " << sum_w << endl;
  cout << "Sum of weighted measurements after all it : " << sum_v << endl;
  cout << "Final value : " << sum_v/sum_w << endl;
  cout << "Uncertainty : " << 1/sqrt(sum_w) << endl;
  cout << "Chi2 : " << chi2 << endl;
  cout << "Chi2/(N-1) : " << chi2/(N-1) << endl;

  TCanvas *cv = new TCanvas("cv", "lifetime collection", 800, 800);
  // cv->SetMargin(0.340961, 0.0514874, 0.17, 0.070162);
  cv->SetMargin(0.390961, 0.0514874, 0.121294, 0.140162);
  TH2D *frame = new TH2D("frame", ";{}^{3}_{#Lambda}H lifetime (ps);", 1000, 0, 510, nMeasures, kOffset, kOffset + nMeasures);
  std::string names[nMeasures]{"PR 136 (1964) B1803", "PRL 20 (1968) 819", "PR 180 (1969) 1307", "NPB 16 (1970) 46", "PRD 1 (1970) 66", "NPB 67 (1973) 269", "STAR, Science 328 (2010) 58", "HypHI, NPA 913 (2013) 170", "ALICE, PLB 754 (2016) 360", "STAR, PRC 97 (2018) 054909", "ALICE, PLB 797 (2019) 134905", "STAR, PRL 128 (2021) 202301"};
  if (kWOPrel == false)
    names[nMeasures - 1] = "ALICE, Pb#minusPb 5.02 TeV";

  std::reverse(std::begin(names), std::end(names));
  for (int i{0}; i < nMeasures; ++i)
  {
    if (i == 0)
      frame->GetYaxis()->SetBinLabel(i + 1, Form("#color[%d]{%s}", 2, names[i].data()));
    else
      frame->GetYaxis()->SetBinLabel(i + 1, names[i].data());
  }

  frame->Draw("col");

  TGraphAsymmErrors *gSpect;
  gSpect = new TGraphAsymmErrors(N, lifetau, point, err_x_low, err_x_high, err_y, err_y);

  TGraphAsymmErrors *gSpect2; // syst. err.
  gSpect2 = new TGraphAsymmErrors(N, lifetau, point, errsyst_x, errsyst_x, errsyst_y, errsyst_y);

  Float_t point_a[1] = {0};
  Float_t lifetau_a[1] = {kPreliminaryTau};
  Float_t err_y_a[1] = {0};
  Float_t err_x_low_a[1] = {kPreliminaryTauStat[0]};
  Float_t err_x_high_a[1] = {kPreliminaryTauStat[1]};
  Float_t errsyst_x_a[1] = {kPreliminaryTauSyst[0]};
  Float_t errsyst_y_a[1] = {0.1};
  TGraphAsymmErrors *gSpect_alice = new TGraphAsymmErrors(1, lifetau_a, point_a, err_x_low_a, err_x_high_a, err_y_a, err_y_a);

  TGraphAsymmErrors *gSpect2_alice = new TGraphAsymmErrors(1, lifetau_a, point_a, errsyst_x_a, errsyst_x_a, errsyst_y_a, errsyst_y_a);

  TLine *Fcn = new TLine(263.2, kOffset, 263.2, nMeasures + kOffset);
  Fcn->SetLineWidth(2);
  Fcn->SetLineColor(kBlack);
  TBox *boxLambda = new TBox(263.2 - 2, kOffset, 263.2 + 2, nMeasures + kOffset);
  boxLambda->SetFillColorAlpha(kBlack, 0.35);
  boxLambda->SetFillStyle(3001);
  boxLambda->SetLineWidth(1);
  boxLambda->SetLineStyle(2);
  boxLambda->SetLineColor(kBlack);
  boxLambda->Draw("same");

  TLine *Fcn1 = new TLine(206.335, kOffset, 206.335, nMeasures + kOffset);
  Fcn1->SetLineWidth(1);
  Fcn1->SetLineStyle(2);
  Fcn1->SetLineColor(kOrange - 3);

  TLine *Fcn2 = new TLine(239.3, kOffset, 239.3, nMeasures + kOffset);
  Fcn2->SetLineWidth(2);
  Fcn2->SetLineStyle(7);
  Fcn2->SetLineColor(kMagenta + 2);

  // 1.97085e+02   1.31749e+01
  // 2.08360e+02   1.42092e+01

  TLine *Fcn3 = new TLine(255.5, kOffset, 255.5, nMeasures + kOffset);
  Fcn3->SetLineWidth(2);
  Fcn3->SetLineStyle(5);
  Fcn3->SetLineColor(kBlue);

  TLine *Fcn4 = new TLine(268, kOffset, 268, nMeasures + kOffset);
  Fcn4->SetLineWidth(2);
  Fcn4->SetLineStyle(10);
  Fcn4->SetLineColor(kCyan - 7);
  // TLine *Fcn5 = new TLine(213, kOffset, 213, nMeasures + kOffset);
  // Fcn5->SetLineWidth(2);
  // Fcn5->SetLineStyle(7);
  // Fcn5->SetLineColor(kMagenta + 2);
  TLine *Fcn6 = new TLine(232, kOffset, 232, nMeasures + kOffset);
  Fcn6->SetLineWidth(2);
  Fcn6->SetLineStyle(9);
  Fcn6->SetLineColor(kGreen + 2);

  TLine *Fcn7 = new TLine(173, kOffset, 173, nMeasures + kOffset);
  Fcn7->SetLineWidth(2);
  Fcn7->SetLineStyle(6);
  Fcn7->SetLineColor(kPink - 5);
  const int kBCT = TColor::GetColor("#b8d4ff");
  TBox *boxFcn7 = new TBox(234 - 27, kOffset, 234 + 27, nMeasures + kOffset);
  boxFcn7->SetFillColorAlpha(kOrange - 3, 0.6);
  boxFcn7->SetFillStyle(3004);
  boxFcn7->SetLineWidth(1);
  boxFcn7->SetLineStyle(2);
  boxFcn7->SetLineColor(kOrange - 3);
  boxFcn7->Draw("same");

  double y[] = {kOffset, nMeasures + kOffset};
  double x[] = {206.335, 206.335};
  double ey[] = {0, 0};
  double ex[] = {13.498, 12.806};
  TGraphErrors *ge = new TGraphErrors(2, x, y, ex, ey);
  ge->SetFillColorAlpha(kOrange - 3, 0.35);
  ge->SetFillStyle(3004);
  ge->Draw("3SAME");

  TLine *Fcn1_up = new TLine(219.83, kOffset, 219.83, nMeasures + kOffset);
  Fcn1_up->SetLineWidth(1);
  Fcn1_up->SetLineColor(kOrange - 3);
  TLine *Fcn1_low = new TLine(193.529, kOffset, 193.529, nMeasures + kOffset);
  Fcn1_low->SetLineWidth(1);
  Fcn1_low->SetLineColor(kOrange - 3);

  TBox *boxFcn1 = new TBox(193.529, kOffset, 219.83, nMeasures + kOffset);
  boxFcn1->SetFillColorAlpha(kOrange - 3, 0.35);
  boxFcn1->SetFillStyle(3004);
  boxFcn1->SetLineWidth(1);
  boxFcn1->SetLineStyle(2);
  boxFcn1->SetLineColor(kOrange - 3);
  // boxFcn1->Draw("same");

  gSpect->SetMarkerStyle(kMarkTyp);
  gSpect->SetMarkerSize(1.6);
  gSpect->SetMarkerColor(kMarkCol);
  gSpect->SetLineColor(kMarkCol);

  gSpect2->SetMarkerStyle(0);
  gSpect2->SetMarkerColor(kMarkCol);
  gSpect2->SetMarkerSize(0.1);
  gSpect2->SetLineStyle(1);
  gSpect2->SetLineColor(kMarkCol);
  gSpect2->SetFillColor(0);
  gSpect2->SetFillStyle(0);

  gSpect_alice->SetMarkerStyle(kFullSquare);
  gSpect_alice->SetMarkerSize(1.);
  gSpect_alice->SetMarkerColor(kRed);
  gSpect_alice->SetLineColor(kRed);

  gSpect2_alice->SetMarkerStyle(0);
  gSpect2_alice->SetMarkerColor(kRed);
  gSpect2_alice->SetMarkerSize(0.1);
  gSpect2_alice->SetLineStyle(1);
  gSpect2_alice->SetLineColor(kRed);
  gSpect2_alice->SetFillColor(0);
  gSpect2_alice->SetFillStyle(0);

  Fcn->Draw("SAME");
  Fcn2->Draw("SAME");
  // Fcn1_up->Draw("SAME");
  // Fcn1_low->Draw("SAME");
  Fcn3->Draw("SAME");
  Fcn4->Draw("SAME");
  // Fcn5->Draw("SAME");
  Fcn6->Draw("SAME");
  Fcn7->Draw("SAME");
  gSpect->Draw("pzsame");
  gSpect2->Draw("spe2");
  if (kWOPrel == false)
  {
    gSpect_alice->Draw("pzsame");
    gSpect2_alice->Draw("spe2");
  }

  // TLegend *leg1 = new TLegend(0.341533, 0.95, 0.99085, 0.99);
  TLegend *leg1 = new TLegend(.68, 0.26, .78, 0.43);
  leg1->SetFillStyle(0);
  leg1->SetMargin(0.2); // separation symbol-text
  leg1->SetBorderSize(0);
  // leg1->SetNColumns(2);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.025);
  leg1->AddEntry(Fcn, "#Lambda lifetime - PDG value", "l");
  // leg1->AddEntry(boxFcn1, "Published world average", "fl");

  TLegend *leg2 = new TLegend(0.340961, 0.87062, 0.949085, 0.995175);
  leg2->SetFillStyle(0);
  leg2->SetMargin(0.16); // separation symbol-text
  leg2->SetBorderSize(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.022);
  leg2->SetNColumns(2);
  leg2->SetHeader("Theoretical predictions");
  leg2->AddEntry(Fcn2, "Nuo. Cim. 46 (1966) 786", "l");
  leg2->AddEntry(Fcn7, "Nuo. Cim. 51 (1979) 180-186", "l");
  leg2->AddEntry(Fcn6, "J.Phys. G18 (1992) 339-357", "l");
  leg2->AddEntry(Fcn3, "PRC 57 (1998) 1595", "l");
  leg2->AddEntry(Fcn4, "PRC 102 (2020) 064002", "l");
  leg2->AddEntry(boxFcn7, "PLB 811 (2020) 135916", "f");

  leg1->Draw();
  leg2->Draw();

  cv->SaveAs("../Results/2Body/Collection.eps");
  cv->SaveAs("../Results/2Body/Collection.pdf");
  cv->SaveAs("../Results/2Body/Collection.png");
}
