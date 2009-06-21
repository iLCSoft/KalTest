#include <iomanip>
#include <sstream>

void Plot(int layer=259, const char *fname="p-100.0t000.10k.root")
{
  using namespace std;

  TFile *filep = TFile::Open(fname);

  TH2D *h1 = new TH2D("h1","",15,0.,250.,10,-0.05,+0.05);
  stringstream str1;
  str1 << "dxin"      << setw(3) << setfill('0') << layer << ":"
       << "250-abs(z" << setw(3) << setfill('0') << layer << ")"
       << ">>h1" << ends;
  track->Draw(str1.str().data(),"");
  h1->FitSlicesY();

  TH2D *h2 = new TH2D("h2","",15,0.,250.,10,-0.05,+0.05);
  stringstream str2;
  str2 << "dxot"      << setw(3) << setfill('0') << layer << ":"
       << "250-abs(z" << setw(3) << setfill('0') << layer << ")"
       << ">>h2" << ends;
  track->Draw(str2.str().data(),"");
  h2->FitSlicesY();

  h1_2->SetMinimum(0.);
  h1_2->SetMaximum(0.03);
  h1_2->SetMarkerStyle(4);
  h1_2->SetMarkerSize(1);
  h1_2->Draw();

  h2_2->SetMarkerStyle(5);
  h2_2->SetMarkerSize(1);
  h2_2->Draw("same");
}
