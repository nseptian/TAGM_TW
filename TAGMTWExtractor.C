
using namespace RooFit;

struct gaussianFitResults
{
    int ph_bin;
    float mean;
};

enum model{
    singleGaussian,doubleGaussian,tripleGaussian
};

TH1* GetHistogram(TFile *f, int col, int ph_bin);
gaussianFitResults WriteGaussianFitResults(ofstream &fout, TH1 *h, int col, int ph_bin, double bkg2MeanMax, TString resultSubFolder);
double GetMode(TH1 *h);
double GetDipBinCenter(TH1 *h);
vector<int> GetNumberOfPP(TFile *f);
void WriteTWFitResults(ofstream &fouttw, TGraph *gr, Int_t col);

//some parameters to configure
const Int_t col=2;//first column to fit
const Int_t ncol=99;//number of column to fit
Int_t ppLowLims[ncol];
                        //  /*20,*/23//,24,20,23,//1-5
                        // 15,19,21,25,26,//6-10
                        // 24,25,24,22,25,//11-15
                        // 23,24,23,25,23,//16-20
                        // 27,26,25,26,21,//21-25
                        // 23,25,23,27,23,//26-30
                        // 24,31,23,22,22,//31-35
                        // 21,27,27,25,28,//36-40
                        // 25,23,13,11,15,//41-45
                        // 14,13,11,14,13,//46-50
                        // 11,16,15,12,14,//51-55
                        // 14,15,16,14,14,//56-60
                        // 13,17,19,17,15,//61-65
                        // 13,17,19,16,17,//66-70
                        // 15,17,11,17,12,//71-76
                        // 15,15,15,13,18,//76-80
                        // 15,15,15,17,16,//81-85
                        // 14,10,18,17,15,//86-90
                        // 17,17,15,18,14,//91-95
                        // 15,12,17,14,16
                        // };//96-100

const int numbOfEntriesLims = 3500;
Int_t ppUpLims[ncol];
const Float_t dPPbin=16;
const Int_t dtLowLims=-3;
const Int_t dtUpLims=7;

const double sgnSigMin = 0.01;
const double sgnSigMax = 1.0;
const double bkg1SigMin = 0.8;
const double bkg1SigMax = 10.0;
const double bkg2SigMin = 0.3;
const double bkg2SigMax = 0.6;

const double deltaMean = 2;
const double sgnSig = 0.6;
const double bkg1Sig = 2.0;
const double bkg2Sig = 0.2;
const double c1 = 0.1;
const double c2 = 0.4;
// TString resultFolder = "resultTAGMTWExtractor/";
// TString rootFileFolder = "root/";
// TString rootFilePrefix = "hd_root-r";

TString resultFolder = "resultTAGMTWExtractor/";
TString rootFileFolder = "/d/grid17/sdobbs/2019-11-mon/mon_ver17/";
TString rootFilePrefix = "hd_root_";

// model fitModel = tripleGaussian;
// const bool useModelBasedOnChi2 = kTRUE;//if TRUE set ch2Thres and model options;
const bool useMinos = kTRUE; 
const int nModel = 2;
const model fitModelOptions[nModel] = {doubleGaussian,tripleGaussian}; //fitModelOptions[0] is the primary model
const double chi2Thres[2] = {1.0,3.0};
const double chi2PPLims = 0.8; 

int TAGMTWExtractor(TString runNumber="72369") {

    TString rootFile=rootFileFolder+rootFilePrefix+runNumber+".root";

    gROOT->SetBatch(kTRUE);

    TFile *f = new TFile(rootFile,"read");
    vector<int> numberOfPP = GetNumberOfPP(f);
    TString makeDir = "mkdir -p ";
    makeDir += resultFolder; 
    system(makeDir);
    TCanvas *c0 = new TCanvas();
    TString resultSubFolder = resultFolder+runNumber+"/";
    TGraph *g[ncol];
    ofstream fouttw; fouttw.open(resultFolder+"TWFit_"+runNumber+".csv");
    for(int j=col;j<col+ncol;j++){
        stringstream ss; ss << j;
        ofstream fout; fout.open(resultSubFolder+"chi2Fit_col_"+TString(ss.str())+".csv");
        const int n = numberOfPP[j-col];
        Float_t meanGraph[n];
        Float_t ppGraph[n];
        double bkg2MeanMax;
        for (int i=ppLowLims[j-col];i<ppUpLims[j-col];i++) {
            TH1 *h = GetHistogram(f,j,i);
            h->Draw();
            if (i==ppLowLims[j-col]) bkg2MeanMax = GetDipBinCenter(h);
            gaussianFitResults fR = WriteGaussianFitResults(fout,h,j,i,bkg2MeanMax,resultSubFolder);
            meanGraph[i-ppLowLims[j-col]]=fR.mean;
            ppGraph[i-ppLowLims[j-col]]=(i-0.5)*(dPPbin);
            delete h;
        }
        g[j-col] = new TGraph(n,ppGraph,meanGraph);
        g[j-col]->SetMarkerSize(1);
        g[j-col]->SetMarkerStyle(kStar);
        g[j-col]->GetXaxis()->SetTitle("Pulse peak (adc counts)");
        g[j-col]->GetYaxis()->SetTitle("T - RF (ns)");
        TString namegrTitle = "Timewalk col ";
        namegrTitle+=j;
        g[j-col]->SetTitle(namegrTitle);
        WriteTWFitResults(fouttw,g[j-col],j);
        g[j-col]->Draw("ap");
        TString pdfName;
        pdfName = resultSubFolder;
        pdfName += "gr_pp_vs_dt_fitted_";
        pdfName += j;
        pdfName += ".pdf";
        c0->Print(pdfName);
        // if (j-col>0) delete g[j-col-1];
        // if (j==col+ncol-1) delete g[j-col];
        // cout << "bkg2MeanMax = " << bkg2MeanMax << endl;
        fout.close();
    }
    fouttw.close();
    return 0;
}

TH1* GetHistogram(TFile *f, int col, int ph_bin) {
    stringstream ss_c; ss_c << col;
    TH2I *h2 = (TH2I*)f->Get("TAGM_TW/t-rf/h_dt_vs_pp_"+TString(ss_c.str()));
    // h2->Draw("colz");
    stringstream ss_ph; ss_ph << h2->GetXaxis()->GetBinCenter(ph_bin);
    TH1 *h;
    h = h2->ProjectionY(TString(h2->GetName())+"_"+TString(ss_ph.str()),ph_bin,ph_bin);
    h->GetXaxis()->SetRangeUser(dtLowLims,dtUpLims);
    // h->Draw();
    h->SetTitle("TAGM column "+TString(ss_c.str())+", Pulse-peak "+TString(ss_ph.str()));
    h->GetXaxis()->SetTitle("time(TDC) - time(RF) [ns]");
    return h;
}

gaussianFitResults WriteGaussianFitResults(ofstream &fout, TH1 *h, int col, int ph_bin, double bkg2MeanMax,TString resultSubFolder) {
    TString sep = ",";
    const double sgnMean = GetMode(h); //get maximum bin center

    double chi2Th = chi2Thres[0];
    if (TMath::Abs(sgnMean) < chi2PPLims) chi2Th = chi2Thres[1];

    RooRealVar x("TDC time difference","TDC time difference [ns]",dtLowLims,dtUpLims);
    RooDataHist data("data","data",RooArgList(x),h);
    
    stringstream ss; ss << col;
    TString makeDirFit = "mkdir -p ";
    makeDirFit += resultSubFolder;
    makeDirFit += "column_";
    makeDirFit += TString(ss.str());
    system(makeDirFit);

    gaussianFitResults fResults;
    for (int iModel=0;iModel<nModel;iModel++){
        TCanvas *canvas = new TCanvas("c","c",800,500);
        RooPlot *plot = x.frame();

        plot->SetTitle(h->GetTitle());
        plot->SetXTitle(h->GetXaxis()->GetTitle());
        // plot->SetYTitle("TAGM hits");
        plot->SetTitleOffset(1.1,"Y");

        RooRealVar meanSgn("meanSgn","meanSgn",sgnMean,sgnMean-0.5,sgnMean+0.5);
        RooRealVar sigmaSgn("sigmaSgn","sigmaSgn",sgnSig,sgnSigMin,sgnSigMax);

        double chi2Fit=99999;

        if (fitModelOptions[iModel]==singleGaussian){
                RooGaussian fitFunction("gauss1","gauss1",x,meanSgn,sigmaSgn);
                fitFunction.fitTo(data,RooFit::Minos(useMinos));//
                data.plotOn(plot);
                fitFunction.plotOn(plot,LineColor(kRed));
                chi2Fit = plot->chiSquare();
                fitFunction.paramOn(plot,Layout(0.9,0.6,0.9),Format("NEU",AutoPrecision(1)));
                plot->Draw();
                fResults.mean = meanSgn.getVal();
            }
            else{
                if (fitModelOptions[iModel]==doubleGaussian){
                    RooRealVar meanBkg1("meanBkg1","meanBkg1",sgnMean,dtLowLims,dtUpLims);
                    RooGaussian gauss1("gauss1","gauss1",x,meanSgn,sigmaSgn);
                    RooRealVar sigmaBkg1("sigmaBkg1","sigmaBkg1",bkg1Sig,bkg1SigMin,bkg1SigMax);
                    RooGaussian gauss2("gauss2","gauss2",x,meanBkg1,sigmaBkg1);
                    RooRealVar cBkg1("cBkg1","cBkg1",c1,0.00,0.8);
                    RooAddPdf fitFunction("doubleGauss","doubleGauss",RooArgList(gauss2,gauss1),RooArgList(cBkg1));
                    fitFunction.fitTo(data,RooFit::PrintLevel(-1),RooFit::Minos(useMinos));//
                    data.plotOn(plot);
                    fitFunction.plotOn(plot,LineColor(kRed));
                    chi2Fit = plot->chiSquare();
                    fitFunction.plotOn(plot,Components("gauss1"),LineColor(kBlue),LineStyle(kDashed));
                    fitFunction.plotOn(plot,Components("gauss2"),LineColor(kGreen),LineStyle(kDashed));
                    fitFunction.paramOn(plot,Layout(0.9,0.6,0.9),Format("NEU",AutoPrecision(1)));
                    plot->Draw();
                    fResults.mean = meanSgn.getVal();
                    fout << col << sep << (ph_bin-0.5)*dPPbin << sep << "doubleGaussian" << sep << chi2Fit << endl;
                }
                else
                {
                    if (fitModelOptions[iModel]==tripleGaussian){
                        RooRealVar meanBkg1("meanBkg1","meanBkg1",sgnMean,sgnMean-0.5,sgnMean+0.5);
                        RooRealVar meanBkg2("meanBkg2","meanBkg2",sgnMean-deltaMean,dtLowLims,bkg2MeanMax);
                        RooGaussian gauss1("gauss1","gauss1",x,meanSgn,sigmaSgn);
                        RooRealVar sigmaBkg2("sigmaBkg2","sigmaBkg2",bkg2Sig,bkg2SigMin,bkg2SigMax);
                        RooGaussian gauss2("gauss2","gauss2",x,meanBkg2,sigmaBkg2);
                        RooRealVar sigmaBkg1("sigmaBkg1","sigmaBkg1",sgnSig,bkg1SigMin,bkg1SigMax);
                        RooGaussian gauss3("gauss3","gauss3",x,meanBkg1,sigmaBkg1);
                        RooRealVar cBkg1("cBkg1","cBkg1",c1,0.01,0.4);
                        RooRealVar cBkg2("cBkg2","cBkg2",c2,0.01,0.4);
                        RooAddPdf fitFunction("tripleGaussian","tripleGaussian",RooArgList(gauss2,gauss3,gauss1),RooArgList(cBkg2,cBkg1));
                        fitFunction.fitTo(data,RooFit::PrintLevel(-1),RooFit::Minos(useMinos));//
                        data.plotOn(plot);
                        fitFunction.plotOn(plot,LineColor(kRed));
                        chi2Fit = plot->chiSquare();
                        fitFunction.plotOn(plot,Components("gauss1"),LineColor(kBlue),LineStyle(kDashed));
                        fitFunction.plotOn(plot,Components("gauss2"),LineColor(kGreen),LineStyle(kDashed));
                        fitFunction.plotOn(plot,Components("gauss3"),LineColor(kYellow),LineStyle(kDashed));
                        fitFunction.paramOn(plot,Layout(0.9,0.6,0.9),Format("NEU",AutoPrecision(1)));
                        plot->Draw();
                        fResults.mean = meanSgn.getVal();
                        fout << col << sep << (ph_bin-0.5)*dPPbin << sep << "tripleGaussian" << sep << chi2Fit << endl;
                    }
                    else
                    {  
                        cout << "Error! Incorrect chosen model.";
                        EXIT_FAILURE;
                    }
                }
            }
        
        // cout << endl << "chi2 = " << chi2Fit << endl;
        canvas->Print(resultSubFolder+"column_"+TString(ss.str())+"/"+TString(h->GetName())+".pdf");
        delete plot;    
        delete canvas;

        if (chi2Fit < chi2Th) break;
    
    }
    
    fResults.ph_bin = ph_bin;
    return fResults;
}

double GetMode(TH1 *h) {
    if (h->GetEntries() == 0.0) return 999.0;
    int max_bin = h->GetMaximumBin();
    double max = h->GetBinContent(max_bin);
    if (max < 7.0) return 999.0;
    return h->GetBinCenter(max_bin);
}

vector<int> GetNumberOfPP(TFile *f){ //get number of PP
    vector<int> vecIterator;
    for(int j=col;j<col+ncol;j++){
        int iterator=0;
        int ppLLims=-1;
        vector<double> vecMode;
        for (int i=0;i<99999;i++) {
            TH1 *h = GetHistogram(f,j,i);
            int entries = h->GetEntries();
            // cout << entries << " " << GetMode(h) <<  endl;
            if ((entries > numbOfEntriesLims) && (ppLLims==-1)){
                ppLLims = i;
            }
            if (ppLLims != -1){
                if (entries < numbOfEntriesLims-(0.1*numbOfEntriesLims)) break;
                else {
                    iterator++;
                    vecMode.push_back(GetMode(h));
                }
            }
        }
        int maxModeIdx = max_element(vecMode.begin(),vecMode.end()) - vecMode.begin();
        ppLowLims[j-col] = ppLLims+maxModeIdx;
        ppUpLims[j-col] = ppLLims+iterator;
        iterator -= maxModeIdx;
        vecIterator.push_back(iterator);
    }
    return vecIterator;
}

double GetDipBinCenter(TH1 *h){
    int upperBin = h->GetMaximumBin();
    const int nBinUsed = 15;
    int binUsed[nBinUsed];
    double prevAvgBinUsed;
    double temp = 9999;
    double avgBinUsed;
    do{
        avgBinUsed=0;
        prevAvgBinUsed = temp;
        for (int i=0;i<nBinUsed;i++){
            avgBinUsed+=h->GetBinContent(upperBin-i);
        }
        avgBinUsed/=nBinUsed;
        temp = avgBinUsed;
        upperBin-=1;
    } while(avgBinUsed < prevAvgBinUsed);
    return h->GetBinCenter(upperBin-9);
}

void WriteTWFitResults(ofstream &fouttw, TGraph *gr, Int_t col){
    TString separator = ",";

    Float_t xLowLims = gr->GetXaxis()->GetXmin();
    Float_t xUpLims = gr->GetXaxis()->GetXmax();

    TF1 *twFitFunction = new TF1("twFitFunction","[0] + [1]*((1/([4]*x+[3]))^[2])",xLowLims,xUpLims);
    twFitFunction->SetParameter(0,-0.0661505);
    twFitFunction->SetParameter(1,5.55105e+08);
    twFitFunction->SetParameter(2,42.7598);
    twFitFunction->SetParameter(3,1.43117);
    twFitFunction->SetParameter(4,0.000330912);

    TFitResultPtr TWFitResultPtr = gr->Fit("twFitFunction","EMS");
    Double_t p0 = TWFitResultPtr->Value(0);
    Double_t p1 = TWFitResultPtr->Value(1);
    Double_t p2 = TWFitResultPtr->Value(2);
    Double_t p3 = TWFitResultPtr->Value(3);
    Double_t p4 = TWFitResultPtr->Value(4);
    Double_t chi2tw = TWFitResultPtr->Chi2();
    Int_t fitStatus = TWFitResultPtr;

    fouttw << col << separator << p0 << separator << p1 << separator << p2 << separator << p3 << separator << p4 << separator << chi2tw << separator << fitStatus << endl;
}