
using namespace RooFit;

void WriteTWFitResults(ofstream &fouttw, TGraphErrors *gr, Int_t col);

const Int_t col=2;//first col
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
    for (int j=col;j<col+ncol;j++){
        TString grName = "gr_pp_vs_dt_fitted_";
        grName += j;
        g[j-col] = (TGraphErrors*)fOutput->Get(grName);
        WriteTWFitResults(fouttw,g[j-col],j);
        g[j-col]->Draw("ap");
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
    TString separator = ",";

    Float_t xLowLims = gr->GetXaxis()->GetXmin();
    Float_t xUpLims = gr->GetXaxis()->GetXmax();

    TF1 *twFitFunction = new TF1("twFitFunction","[0] + [1]*((1/([4]*x+[3]))^[2])",xLowLims,xUpLims);
    twFitFunction->SetParameter(0,-0.0661505);
    twFitFunction->SetParameter(1,5.55105e+08);
    twFitFunction->SetParameter(2,42.7598);
    twFitFunction->SetParameter(3,1.43117);
    twFitFunction->SetParameter(4,0.000330912);

    TFitResultPtr TWFitResultPtr = gr->Fit("twFitFunction","MS");
    Double_t p0 = TWFitResultPtr->Value(0);
    Double_t p1 = TWFitResultPtr->Value(1);
    Double_t p2 = TWFitResultPtr->Value(2);
    Double_t p3 = TWFitResultPtr->Value(3);
    Double_t p4 = TWFitResultPtr->Value(4);
    Double_t chi2tw = TWFitResultPtr->Chi2();
    Int_t fitStatus = TWFitResultPtr;

    fouttw << col << separator << p0 << separator << p1 << separator << p2 << separator << p3 << separator << p4 << separator << chi2tw << separator << fitStatus << endl;
    delete twFitFunction;
}