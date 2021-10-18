
using namespace RooFit;

void WriteTWFitResults(ofstream &fouttw, TGraphErrors *gr, Int_t col);

const Int_t fcol=1;//first col
TString rootFileOutputPrefix = "TW_";
TString resultFolder = "resultTAGMTWExtractor/";

void TAGMTWFit(TString runNumber="72369"){
    TString rootFileOutputName = resultFolder+rootFileOutputPrefix+runNumber+".root";
    TString resultSubFolder = resultFolder+runNumber+"/";
    TString makeDir = "mkdir -p ";
    makeDir += resultSubFolder; 
    system(makeDir);
    TFile *fOutput = new TFile(rootFileOutputName,"read");
    ofstream fouttw; fouttw.open(resultFolder+"TWFit_"+runNumber+".csv");
    const Int_t ncol = fOutput->GetNkeys();
    TGraphErrors *g[ncol];
    TCanvas *c0 = new TCanvas();
    for (int j=fcol;j<fcol+ncol;j++){
        TString grName = "gr_pp_vs_dt_fitted_";
        grName += j;
        g[j-fcol] = (TGraphErrors*)fOutput->Get(grName);
        WriteTWFitResults(fouttw,g[j-fcol],j);
        g[j-fcol]->Draw("ap");
        TString pdfName;
        pdfName = resultSubFolder;
        pdfName += grName;
        pdfName += "_twfit.pdf";
        c0->Print(pdfName);
    }
    fouttw.close();
    fOutput->Close();
}

void WriteTWFitResults(ofstream &fouttw, TGraphErrors *gr, Int_t col){
    TString separator = " ";
    Double_t xLowLims,y1;
    gr->GetPoint(0,xLowLims,y1);
    Double_t xUpLims,y2;
    gr->GetPoint(gr->GetN()-1,xUpLims,y2);
    // cout << xLowLims << " " << xUpLims << endl;
    TF1 *twFitFunction = new TF1("twFitFunction","[0] + [1]*((1/(x+[3]))**[2])",xLowLims,xUpLims);
    twFitFunction->SetParameter(0,0);
    twFitFunction->SetParameter(1,100);
    twFitFunction->SetParameter(2,2);
    twFitFunction->SetParameter(3,-1*gr->GetXaxis()->GetXmin());
    gr->Fit("twFitFunction","EQR");
    twFitFunction->FixParameter(2,twFitFunction->GetParameter(2));
    TFitResultPtr TWFitResultPtr = gr->Fit("twFitFunction","EQSR");
    Double_t p0 = TWFitResultPtr->Value(0);
    Double_t p1 = TWFitResultPtr->Value(1);
    Double_t p2 = TWFitResultPtr->Value(2);
    Double_t p3 = TWFitResultPtr->Value(3);
    Double_t chi2tw = TWFitResultPtr->Chi2();
    Int_t fitStatus = TWFitResultPtr;
    if((col == 9) || (col == 27) || (col == 81) || (col == 99)){
        for (Int_t i=0;i<6;i++) fouttw << i << separator << col << separator << p0 << separator << p1 << separator << p2 << separator << p3 << separator << "0." << endl;
    }
    else fouttw << "0" << separator << col << separator << p0 << separator << p1 << separator << p2 << separator << p3 << separator << "0." << endl;
    delete twFitFunction;
}