#include <TProfile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TPaveLabel.h>

#include <sstream>
#include <iostream>

TFile *f;
TH2I *h_tw[102];
TH2I *h_tw_ind[5][4];

TString rootFileFolder = "root/";
TString rootFilePrefix = "hd_root-r";

// datasets
// TString rootFileFolder = "/d/grid17/sdobbs/2017-01-mon/mon_ver34/rootfiles/"; //03XXX
// TString rootFileFolder = "/d/grid17/sdobbs/2018-01-mon/mon_ver21/rootfiles/"; //04XXX
// TString rootFileFolder = "/d/grid17/sdobbs/2018-08-mon/mon_ver15/rootfiles/"; //05XXX
// TString rootFileFolder = "/d/grid17/sdobbs/2019-01-mon/mon_ver15/rootfiles/"; //06XXX
// TString rootFileFolder = "/d/grid17/sdobbs/2019-11-mon/mon_ver17/"; //07XXX

// TString rootFileFolder = "";
// TString rootFilePrefix = "hd_root_";

TString resultFolder = "print_corrected/";

void print_corrected(TString runNumber = "72369-29June2021") {

   TString inputFile=rootFileFolder+rootFilePrefix+runNumber+".root";
   
   TCanvas *c1 = new TCanvas();
   c1->Divide(3,3);

   f = new TFile(inputFile);
   std::cout << "file: " << inputFile << std::endl;

   TString makeDir = "mkdir -p ";
   makeDir += resultFolder; 
   system(makeDir);
   TString resultSubFolder = resultFolder+runNumber+"/";
   makeDir = "mkdir -p ";
   makeDir += resultSubFolder;
   system(makeDir);

   Int_t idx = 0;
   for (Int_t i = 0; i < 102; ++i)
   {
      c1->cd(idx+1);
      gPad->SetLogz();
      h_tw[i] = (TH2I*)f->Get(Form("TAGM_TW/t-rf/h_dt_vs_pp_%i",i+1));
      h_tw[i]->SetTitle(Form("Timewalk Col %i",i+1));
      h_tw[i]->Draw("colz");
      idx++;
      // cout << i+1 << " " << idx << endl;
      if (idx == 9){
         TString pdfName;
         pdfName += resultSubFolder;
         pdfName += "h_dt_vs_pp_";
         pdfName += i-7;
         pdfName += "_";
         pdfName += i+1;
         pdfName += ".pdf";
         c1->Print(pdfName);
         idx = 0;
      }
   }
}