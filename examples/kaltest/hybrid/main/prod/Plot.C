#include <iomanip>
#include <sstream>

void Plot(int layer=259, const char *fname="p-100.0t000.10k.root")
{
  using namespace std;

  gStyle->SetOptStat(0);

  TFile *filep = TFile::Open(fname);
  //--
  // Residual including the hit in fit
  //--
  TH2D *hin = new TH2D("hin","",15,0.,250.,10,-0.05,+0.05);
  stringstream str1;
  str1 << "dxin"      << setw(3) << setfill('0') << layer << ":"
       << "250-abs(z" << setw(3) << setfill('0') << layer << ")"
       << ">>hin" << ends;
  track->Draw(str1.str().data(),"");
  hin->FitSlicesY();
  //--
  // Residual excluding the hit in fit
  //--
  TH2D *hot = new TH2D("hot","",15,0.,250.,10,-0.05,+0.05);
  stringstream str2;
  str2 << "dxot"      << setw(3) << setfill('0') << layer << ":"
       << "250-abs(z" << setw(3) << setfill('0') << layer << ")"
       << ">>hot" << ends;
  track->Draw(str2.str().data(),"");
  hot->FitSlicesY();
  //--
  // Plot them together
  //--
  hin_2->SetMinimum(0.);
  hin_2->SetMaximum(0.03);
  hin_2->SetMarkerStyle(4);
  hin_2->SetMarkerSize(1);
  hin_2->Draw();

  hot_2->SetMarkerStyle(5);
  hot_2->SetMarkerSize(1);
  hot_2->Draw("same");
  //--
  // Calculate and plot geometric mean
  //--
  TH1D *hgm = new TH1D((*hin_2) * (*hot_2));
  int nbins = hgm->GetNbinsX();
  for (int i=0; i<=nbins; i++) {
    hgm->SetBinContent(i,TMath::Sqrt(hgm->GetBinContent(i)));
    hgm->SetBinError  (i,(hin_2->GetBinError(i)+hot_2->GetBinError(i))/2);
  }

  stringstream titlestr;
  titlestr << "GM Resolutin (Row" << layer << ")" << ends;
  hin_2->SetTitle(titlestr.str().data());
  hin_2->GetXaxis()->SetTitle("Drift Length [cm]");
  hin_2->GetYaxis()->SetTitle("#sigma_{x} [cm]");
  hin_2->GetYaxis()->SetTitleOffset(1.24);

  hgm->SetMarkerStyle(20);
  hgm->SetMarkerSize(1);
  hgm->SetMarkerColor(2);
  hgm->Draw("same");
}
