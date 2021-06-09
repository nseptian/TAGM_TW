

using namespace RooFit;

struct gaussianFitResults
{
    int ph_bin;
    float mean1;
    float mean2;
    float error1;
    float error2;
    float cMean2;
};

struct parsForTripleGaussian
{
    double mean1;
    double mean2;
    double dip;
};

// enum model{
//     singleGaussian,doubleGaussian,tripleGaussian
// };

TH1* GetHistogram(TFile *f, int col, int ph_bin);
gaussianFitResults WriteGaussianFitResults(ofstream &fout, TH1 *h, int col, int ph_bin, TString resultSubFolder);
double GetMode(TH1 *h);
parsForTripleGaussian GetParsForTripleGaussian(TH1 *h);
double GetDipBinCenter(TH1 *h, double mode);
vector<vector<int>> GetPP(TFile *f);
// void WriteTWFitResults(ofstream &fouttw, TGraphErrors *gr, Int_t col);

//some parameters to configure
const Int_t col=1;//first column to fit
const Int_t ncol=102;//number of column to fit
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

// const int numbOfEntriesLims = 3500; 07XXX
// const int minModeBinContentSelPP = 100;

int numbOfEntriesLims;
const int minModeBinContentSelPP = 0;

// const int numbOfEntriesLims = 50;
// const int minModeBinContentSelPP = 10;

Int_t ppUpLims[ncol];
const Float_t dPPbin=16;
const Int_t dtLowLims=-3;
const Int_t dtUpLims=7;

const double sgnSigMin = 0.01;
const double sgnSigMax = 2.0;
const double bkg1SigMin = 2.0;
const double bkg1SigMax = 3.0;
const double bkg2SigMin = 0.2;
const double bkg2SigMax = 0.5;

double deltaMean = 2.5;
const double sgnSig = 0.3;
const double bkg1Sig = 2.0;
const double bkg2Sig = 0.2;
const double c1 = 0.1;
const double c2 = 0.4;

TString resultFolder = "resultTAGMTWExtractor/";

// local sample
TString rootFileFolder = "root/";
TString rootFilePrefix = "hd_root-r";

// datasets
// TString rootFileFolder = "/d/grid17/sdobbs/2017-01-mon/mon_ver34/rootfiles/"; //03XXX
// TString rootFileFolder = "/d/grid17/sdobbs/2018-01-mon/mon_ver21/rootfiles/"; //04XXX
// TString rootFileFolder = "/d/grid17/sdobbs/2018-08-mon/mon_ver15/rootfiles/"; //05XXX
// TString rootFileFolder = "/d/grid17/sdobbs/2019-01-mon/mon_ver15/rootfiles/"; //06XXX
// TString rootFileFolder = "/d/grid17/sdobbs/2019-11-mon/mon_ver17/"; //07XXX
// TString rootFilePrefix = "hd_root_";

TString rootFileOutputPrefix = "TW_";

// model fitModel = tripleGaussian;
// const bool useModelBasedOnChi2 = kTRUE;//if TRUE set ch2Thres and model options;
const bool useMinos = kTRUE; 
// const int nModel = 2;
// const model fitModelOptions[nModel] = {doubleGaussian,tripleGaussian}; //fitModelOptions[0] is the primary model
const double chi2DiffThres = 0.6;
const double maxModeBinContentTrpGauss = 1000;

int TAGMTWExtractor(TString runNumber="72369") {

    TString rootFile=rootFileFolder+rootFilePrefix+runNumber+".root";

    gROOT->SetBatch(kTRUE);

    TFile *f = new TFile(rootFile,"read");
    // vector<int> numberOfPP = GetNumberOfPP(f);
    
    vector<vector<int>> SelectedPP = GetPP(f);
    // return 1;
    TString makeDir = "mkdir -p ";
    makeDir += resultFolder; 
    system(makeDir);
    TCanvas *c0 = new TCanvas();
    TString resultSubFolder = resultFolder+runNumber+"/";
    TGraphErrors *g[ncol];
    // ofstream fouttw; fouttw.open(resultFolder+"TWFit_"+runNumber+".txt");
    TString rootFileOutputName = resultFolder+rootFileOutputPrefix+runNumber+".root";
    TFile *fOutput = new TFile(rootFileOutputName,"recreate");

    for(int j=col;j<col+ncol;j++){
        const int n = SelectedPP[j-col].size();
        // cout << "n = " << n << endl;
        stringstream ss; ss << j;
        ofstream fout; fout.open(resultSubFolder+"chi2Fit_col_"+TString(ss.str())+".csv");
        if (n<=0) {
            //dummy graph for empty pulse peak
            Double_t xempty[6] = {200,300,400,500,600,700};
            Double_t yempty[6] = {0.,0.,0.,0.,0.,0.};
            g[j-col] = new TGraphErrors(6,xempty,yempty);
        }
        else {
            Float_t meanGraph[n];
            Float_t meanErrorGraph[n];
            Float_t ppGraph[n];
            double prevMean;
            for (int i=0;i<n;i++) {
                int idxPP = SelectedPP[j-col][i];
                TH1 *h = GetHistogram(f,j,idxPP);
                h->Draw();
                // if (i==0) bkg2MeanMax = GetDipBinCenter(h);
                gaussianFitResults fR = WriteGaussianFitResults(fout,h,j,idxPP,resultSubFolder);
                if (i==0){
                    if (fR.cMean2 > 0.5) {
                        meanGraph[i]=max(fR.mean1,fR.mean2);
                        if (fR.mean1 > fR.mean2) meanErrorGraph[0]=fR.error1;
                        else meanErrorGraph[i]=fR.error2;
                    }
                    else {
                        meanGraph[i]=fR.mean1;
                        meanErrorGraph[i]=fR.error1;
                    }
                } 
                else{
                    double mean1Diff = abs(fR.mean1-prevMean);
                    double mean2Diff = abs(fR.mean2-prevMean);
                    if ((mean1Diff > mean2Diff) && (fR.cMean2 > 0.5)) {
                        meanGraph[i]=fR.mean2;
                        meanErrorGraph[i]=fR.error2;
                    }
                    else {
                        meanGraph[i]=fR.mean1;
                        meanErrorGraph[i]=fR.error1;
                    }
                }       
                ppGraph[i]=(idxPP-0.5)*(dPPbin);
                prevMean = meanGraph[i];
                delete h;
            }
            g[j-col] = new TGraphErrors(n,ppGraph,meanGraph,0,meanErrorGraph);
        }
        g[j-col]->SetMarkerSize(1);
        g[j-col]->SetMarkerStyle(kStar);
        g[j-col]->GetXaxis()->SetTitle("Pulse peak (adc counts)");
        g[j-col]->GetYaxis()->SetTitle("T - RF (ns)");
        TString namegrTitle = "Timewalk col ";
        namegrTitle+=j;
        g[j-col]->SetTitle(namegrTitle);
        // WriteTWFitResults(fouttw,g[j-col],j);
        g[j-col]->Draw("ap");
        TString pdfName;
        pdfName = resultSubFolder;
        pdfName += "gr_pp_vs_dt_fitted_";
        TString grName = "gr_pp_vs_dt_fitted_";
        pdfName += j;
        grName += j;
        g[j-col]->Write(grName);
        pdfName += ".pdf";
        c0->Print(pdfName);
        // if (j-col>0) delete g[j-col-1];
        // if (j==col+ncol-1) delete g[j-col];
        // cout << "bkg2MeanMax = " << bkg2MeanMax << endl;
        fout.close();
    }
    // fouttw.close();
    return 0;
}

TH1* GetHistogram(TFile *f, int col, int ph_bin) {
    stringstream ss_c; ss_c << col;
    TH2I *h2 = (TH2I*)f->Get("TAGM_TW/t-rf/h_dt_vs_pp_"+TString(ss_c.str()));
    // h2->Draw("colz");
    numbOfEntriesLims = 0.005*h2->GetEntries();
    stringstream ss_ph; ss_ph << h2->GetXaxis()->GetBinCenter(ph_bin);
    TH1 *h;
    h = h2->ProjectionY(TString(h2->GetName())+"_"+TString(ss_ph.str()),ph_bin,ph_bin);
    h->GetXaxis()->SetRangeUser(dtLowLims,dtUpLims);
    
    // h->Draw();
    h->SetTitle("TAGM column "+TString(ss_c.str())+", Pulse-peak "+TString(ss_ph.str()));
    h->GetXaxis()->SetTitle("time(TDC) - time(RF) [ns]");
    return h;
}

gaussianFitResults WriteGaussianFitResults(ofstream &fout, TH1 *h, int col, int ph_bin,TString resultSubFolder) {
    TString sep = ",";

    RooRealVar x("TDC time difference","TDC time difference [ns]",dtLowLims,dtUpLims);
    RooDataHist data("data","data",RooArgList(x),h);
    
    stringstream ss; ss << col;
    TString makeDirFit = "mkdir -p ";
    makeDirFit += resultSubFolder;
    makeDirFit += "column_";
    makeDirFit += TString(ss.str());
    system(makeDirFit);

    gaussianFitResults fResults;

    //double gaussian fit
    double sgnMean = GetMode(h); //get maximum bin center
    // if (bkg2MeanMax >= sgnMean) bkg2MeanMax = sgnMean-0.3;
    TCanvas *canvas = new TCanvas("c","c",800,500);
    RooPlot *plotDouble = x.frame();
    plotDouble->SetTitle(h->GetTitle());
    plotDouble->SetXTitle(h->GetXaxis()->GetTitle());
    plotDouble->SetTitleOffset(1.1,"Y");
    RooRealVar meanSgnDouble("meanSgnDouble","meanSgnDouble",sgnMean,sgnMean-3.0,sgnMean+3.0);
    RooRealVar sigmaSgnDouble("sigmaSgnDouble","sigmaSgnDouble",sgnSig,sgnSigMin,sgnSigMax);

    RooRealVar meanBkg1Double("meanBkg1Double","meanBkg1Double",sgnMean,sgnMean-3.0,sgnMean+3.0);
    RooGaussian gauss1Double("gauss1Double","gauss1Double",x,meanSgnDouble,sigmaSgnDouble);
    RooRealVar sigmaBkg1Double("sigmaBkg1Double","sigmaBkg1Double",bkg1Sig,bkg1SigMin,bkg1SigMax);
    RooGaussian gauss2Double("gauss2Double","gauss2Double",x,meanBkg1Double,sigmaBkg1Double);
    RooRealVar cBkg1Double("cBkg1Double","cBkg1Double",c1,0.01,0.6);
    RooAddPdf fitFunctionDouble("doubleGauss","doubleGauss",RooArgList(gauss2Double,gauss1Double),RooArgList(cBkg1Double));
    fitFunctionDouble.fitTo(data,RooFit::PrintLevel(-1),RooFit::Minos(useMinos));//
    data.plotOn(plotDouble);
    fitFunctionDouble.plotOn(plotDouble,LineColor(kRed));
    double chi2FitDouble = plotDouble->chiSquare();
    fitFunctionDouble.plotOn(plotDouble,Components("gauss1Double"),LineColor(kBlue),LineStyle(kDashed));
    fitFunctionDouble.plotOn(plotDouble,Components("gauss2Double"),LineColor(kGreen),LineStyle(kDashed));
    fitFunctionDouble.paramOn(plotDouble,Layout(0.9,0.6,0.9),Format("NEU",AutoPrecision(1)));
    plotDouble->Draw();
    fResults.mean1 = meanSgnDouble.getVal();
    fResults.error1 = meanSgnDouble.getError();
    fResults.mean2 = meanBkg1Double.getVal();
    fResults.error2 = meanBkg1Double.getError();
    fResults.cMean2 = cBkg1Double.getVal();
    fout << col << sep << (ph_bin-0.5)*dPPbin << sep << "doubleGaussian" << sep << chi2FitDouble << endl;

    // triple gaussian fit
    parsForTripleGaussian params = GetParsForTripleGaussian(h);
    sgnMean = params.mean1;
    double bkg2Mean = params.mean2;
    double bkg2MeanMax = params.dip;
    RooPlot *plotTriple = x.frame();
    plotTriple->SetTitle(h->GetTitle());
    plotTriple->SetXTitle(h->GetXaxis()->GetTitle());
    plotTriple->SetTitleOffset(1.1,"Y");
    RooRealVar meanSgnTriple("meanSgnTriple","meanSgnTriple",sgnMean,sgnMean-2.0,sgnMean+2.0);
    RooRealVar sigmaSgnTriple("sigmaSgnTriple","sigmaSgnTriple",sgnSig,sgnSigMin,sgnSigMax);

    RooRealVar meanBkg1Triple("meanBkg1Triple","meanBkg1Triple",sgnMean,sgnMean-3.0,sgnMean+3.0);
    RooRealVar meanBkg2Triple("meanBkg2Triple","meanBkg2Triple",bkg2Mean,bkg2Mean,bkg2MeanMax);
    RooGaussian gauss1Triple("gauss1Triple","gauss1Triple",x,meanSgnTriple,sigmaSgnTriple);
    RooRealVar sigmaBkg2Triple("sigmaBkg2Triple","sigmaBkg2Triple",bkg2Sig,bkg2SigMin,bkg2SigMax);
    RooGaussian gauss2Triple("gauss2Triple","gauss2Triple",x,meanBkg2Triple,sigmaBkg2Triple);
    RooRealVar sigmaBkg1Triple("sigmaBkg1Triple","sigmaBkg1Triple",sgnSig,bkg1SigMin,bkg1SigMax);
    RooGaussian gauss3Triple("gauss3Triple","gauss3Triple",x,meanBkg1Triple,sigmaBkg1Triple);
    RooRealVar cBkg1Triple("cBkg1Triple","cBkg1Triple",c1,0.02,0.4);
    RooRealVar cBkg2Triple("cBkg2Triple","cBkg2Triple",c2,0.1,0.6);
    RooAddPdf fitFunctionTriple("tripleGaussian","tripleGaussian",RooArgList(gauss2Triple,gauss3Triple,gauss1Triple),RooArgList(cBkg2Triple,cBkg1Triple));
    fitFunctionTriple.fitTo(data,RooFit::PrintLevel(-1),RooFit::Minos(useMinos));//
    data.plotOn(plotTriple);
    fitFunctionTriple.plotOn(plotTriple,LineColor(kRed));
    double chi2FitTriple = plotTriple->chiSquare();
    fout << col << sep << (ph_bin-0.5)*dPPbin << sep << "tripleGaussian" << sep << chi2FitTriple << endl;
    double diffChi2Fit = abs(chi2FitDouble-chi2FitTriple);//(diffChi2Fit > chi2DiffThres) //|| ((cBkg1Double.getVal() <= 0.1) && (diffChi2Fit > 2.0)
    double diffMeanTriple = abs(meanSgnTriple.getVal()-meanBkg2Triple.getVal());
    double modeBinContent = h->GetBinContent(h->GetBin(meanSgnTriple.getVal()));
    bool isTripleFunctionUsed = (diffMeanTriple > 1.0) && (cBkg2Triple.getVal() < 0.4) && (diffChi2Fit > chi2DiffThres) && (((chi2FitTriple < chi2FitDouble) && (modeBinContent < maxModeBinContentTrpGauss) && (cBkg1Double.getVal() >= 0.04)) || ((cBkg1Double.getVal() < 0.04) && (chi2FitTriple < chi2FitDouble)) || ((chi2FitDouble > 3.0)));
    if (isTripleFunctionUsed) {
        fitFunctionTriple.plotOn(plotTriple,Components("gauss1Triple"),LineColor(kBlue),LineStyle(kDashed));
        fitFunctionTriple.plotOn(plotTriple,Components("gauss2Triple"),LineColor(kGreen),LineStyle(kDashed));
        fitFunctionTriple.plotOn(plotTriple,Components("gauss3Triple"),LineColor(kYellow),LineStyle(kDashed));
        fitFunctionTriple.paramOn(plotTriple,Layout(0.9,0.6,0.9),Format("NEU",AutoPrecision(1)));
        plotTriple->Draw();
        fResults.mean1 = meanSgnTriple.getVal();
        fResults.error1 = meanSgnTriple.getError();
        fResults.mean2 = meanBkg2Triple.getVal();
        fResults.error2 = meanBkg2Triple.getError();
        fResults.cMean2 = cBkg2Triple.getVal();
    }
    // cout << endl << "chi2 = " << chi2Fit << endl;
    canvas->Print(resultSubFolder+"column_"+TString(ss.str())+"/"+TString(h->GetName())+".pdf");
    
    delete plotDouble;
    delete plotTriple;    
    delete canvas;

    fResults.ph_bin = ph_bin;
    return fResults;
}

// gaussianFitResults WriteGaussianFitResults(ofstream &fout, TH1 *h, int col, int ph_bin, double bkg2MeanMax,TString resultSubFolder) {
//     TString sep = ",";
//     const double sgnMean = GetMode(h); //get maximum bin center

//     if (bkg2MeanMax >= sgnMean) bkg2MeanMax = sgnMean-0.3;

//     double chi2Th = chi2Thres[0];
//     if (TMath::Abs(sgnMean) < chi2PPLims) chi2Th = chi2Thres[1];

//     RooRealVar x("TDC time difference","TDC time difference [ns]",dtLowLims,dtUpLims);
//     RooDataHist data("data","data",RooArgList(x),h);
    
//     stringstream ss; ss << col;
//     TString makeDirFit = "mkdir -p ";
//     makeDirFit += resultSubFolder;
//     makeDirFit += "column_";
//     makeDirFit += TString(ss.str());
//     system(makeDirFit);

//     gaussianFitResults fResults;
//     for (int iModel=0;iModel<nModel;iModel++){
//         TCanvas *canvas = new TCanvas("c","c",800,500);
//         RooPlot *plot = x.frame();

//         plot->SetTitle(h->GetTitle());
//         plot->SetXTitle(h->GetXaxis()->GetTitle());
//         // plot->SetYTitle("TAGM hits");
//         plot->SetTitleOffset(1.1,"Y");

//         RooRealVar meanSgn("meanSgn","meanSgn",sgnMean,sgnMean-0.5,sgnMean+0.5);
//         RooRealVar sigmaSgn("sigmaSgn","sigmaSgn",sgnSig,sgnSigMin,sgnSigMax);

//         double chi2Fit=99999;

//         if (fitModelOptions[iModel]==singleGaussian){
//                 RooGaussian fitFunction("gauss1","gauss1",x,meanSgn,sigmaSgn);
//                 fitFunction.fitTo(data,RooFit::Minos(useMinos));//
//                 data.plotOn(plot);
//                 fitFunction.plotOn(plot,LineColor(kRed));
//                 chi2Fit = plot->chiSquare();
//                 fitFunction.paramOn(plot,Layout(0.9,0.6,0.9),Format("NEU",AutoPrecision(1)));
//                 plot->Draw();
//                 fResults.mean1 = meanSgn.getVal();
//                 fResults.mean2 = 0;
//                 fResults.cMean2 = 0;
//             }
//             else{
//                 if (fitModelOptions[iModel]==doubleGaussian){
//                     RooRealVar meanBkg1("meanBkg1","meanBkg1",sgnMean,dtLowLims,dtUpLims);
//                     RooGaussian gauss1("gauss1","gauss1",x,meanSgn,sigmaSgn);
//                     RooRealVar sigmaBkg1("sigmaBkg1","sigmaBkg1",bkg1Sig,bkg1SigMin,bkg1SigMax);
//                     RooGaussian gauss2("gauss2","gauss2",x,meanBkg1,sigmaBkg1);
//                     RooRealVar cBkg1("cBkg1","cBkg1",c1,0.01,0.6);
//                     RooAddPdf fitFunction("doubleGauss","doubleGauss",RooArgList(gauss2,gauss1),RooArgList(cBkg1));
//                     fitFunction.fitTo(data,RooFit::PrintLevel(-1),RooFit::Minos(useMinos));//
//                     data.plotOn(plot);
//                     fitFunction.plotOn(plot,LineColor(kRed));
//                     chi2Fit = plot->chiSquare();
//                     fitFunction.plotOn(plot,Components("gauss1"),LineColor(kBlue),LineStyle(kDashed));
//                     fitFunction.plotOn(plot,Components("gauss2"),LineColor(kGreen),LineStyle(kDashed));
//                     fitFunction.paramOn(plot,Layout(0.9,0.6,0.9),Format("NEU",AutoPrecision(1)));
//                     plot->Draw();
//                     fResults.mean1 = meanSgn.getVal();
//                     fResults.error1 = meanSgn.getError();
//                     fResults.mean2 = meanBkg1.getVal();
//                     fResults.error2 = meanBkg1.getError();
//                     fResults.cMean2 = cBkg1.getVal();
//                     fout << col << sep << (ph_bin-0.5)*dPPbin << sep << "doubleGaussian" << sep << chi2Fit << endl;
//                 }
//                 else
//                 {
//                     if (fitModelOptions[iModel]==tripleGaussian){
//                         RooRealVar meanBkg1("meanBkg1","meanBkg1",sgnMean,sgnMean-0.5,sgnMean+0.5);
//                         RooRealVar meanBkg2("meanBkg2","meanBkg2",sgnMean-deltaMean,dtLowLims,bkg2MeanMax);
//                         RooGaussian gauss1("gauss1","gauss1",x,meanSgn,sigmaSgn);
//                         RooRealVar sigmaBkg2("sigmaBkg2","sigmaBkg2",bkg2Sig,bkg2SigMin,bkg2SigMax);
//                         RooGaussian gauss2("gauss2","gauss2",x,meanBkg2,sigmaBkg2);
//                         RooRealVar sigmaBkg1("sigmaBkg1","sigmaBkg1",sgnSig,bkg1SigMin,bkg1SigMax);
//                         RooGaussian gauss3("gauss3","gauss3",x,meanBkg1,sigmaBkg1);
//                         RooRealVar cBkg1("cBkg1","cBkg1",c1,0.01,0.4);
//                         RooRealVar cBkg2("cBkg2","cBkg2",c2,0.01,0.6);
//                         RooAddPdf fitFunction("tripleGaussian","tripleGaussian",RooArgList(gauss2,gauss3,gauss1),RooArgList(cBkg2,cBkg1));
//                         fitFunction.fitTo(data,RooFit::PrintLevel(-1),RooFit::Minos(useMinos));//
//                         data.plotOn(plot);
//                         fitFunction.plotOn(plot,LineColor(kRed));
//                         chi2Fit = plot->chiSquare();
//                         fitFunction.plotOn(plot,Components("gauss1"),LineColor(kBlue),LineStyle(kDashed));
//                         fitFunction.plotOn(plot,Components("gauss2"),LineColor(kGreen),LineStyle(kDashed));
//                         fitFunction.plotOn(plot,Components("gauss3"),LineColor(kYellow),LineStyle(kDashed));
//                         fitFunction.paramOn(plot,Layout(0.9,0.6,0.9),Format("NEU",AutoPrecision(1)));
//                         plot->Draw();
//                         fResults.mean1 = meanSgn.getVal();
//                         fResults.error1 = meanSgn.getError();
//                         fResults.mean2 = meanBkg2.getVal();
//                         fResults.error2 = meanBkg2.getError();
//                         fResults.cMean2 = cBkg2.getVal();
//                         fout << col << sep << (ph_bin-0.5)*dPPbin << sep << "tripleGaussian" << sep << chi2Fit << endl;
//                     }
//                     else
//                     {  
//                         cout << "Error! Incorrect chosen model.";
//                         EXIT_FAILURE;
//                     }
//                 }
//             }
        
//         // cout << endl << "chi2 = " << chi2Fit << endl;
//         canvas->Print(resultSubFolder+"column_"+TString(ss.str())+"/"+TString(h->GetName())+".pdf");
//         delete plot;    
//         delete canvas; 

//         if (chi2Fit < chi2Th) break;
    
//     }
    
//     fResults.ph_bin = ph_bin;
//     return fResults;
// }

double GetMode(TH1 *h) {
    if (h->GetEntries() == 0.0) return 9999.0;
    int max_bin = h->GetMaximumBin();
    double max = h->GetBinContent(max_bin);
    // if (max < minModeBinContentSelPP) return 9999.0;
    double maxBinCenter = h->GetBinCenter(max_bin);
    return maxBinCenter;
}

parsForTripleGaussian GetParsForTripleGaussian(TH1 *h) {
    TH1 *h0 = (TH1*)h->Clone("h0");
    int entries = h0->GetEntries();
    parsForTripleGaussian pars;
    if (entries == 0.0) {
        pars.mean1 = 9999.0;
        pars.mean2 = 0;
        pars.dip = 0;
        delete h0;
        return pars;
    };
    double xmin = dtLowLims;//h->GetXaxis()->GetXmin();
    double xmax = dtUpLims;//h->GetXaxis()->GetXmax();
    double xmid = (xmax + xmin)/2;
    // double binCenterSeparator = (xmax+xmin)/2.0;
    // double binWidth = h->GetXaxis()->GetBinWidth(0);
    // int upBin = h0->GetXaxis;
    // int lowBin = 1;
    // int midBin = upBin/2;
    int max_bin1,max_bin2,max1,max2;
    int iter = 0;
    do {
        h0->GetXaxis()->SetRangeUser(xmid,xmax);
        max_bin1 = h0->GetMaximumBin();
        max1 = h0->GetBinContent(max_bin1);
    
        h0->GetXaxis()->SetRangeUser(xmin,xmid);
        max_bin2 = h0->GetMaximumBin();
        max2 = h0->GetBinContent(max_bin2);

        xmax = xmid;
        xmid = (xmax + xmin)/2;

        cout << entries << " " << (max1-max2) << " " << max1 << " " << max2 << " " << h0->GetBinCenter(max_bin1) << " " << h0->GetBinCenter(max_bin2) << endl;
        iter++;
    } while (iter < 2 && ((((max1-max2) < -0.5*max2)) || (max_bin1-max_bin2 < 30)));
    
    double mode1 = h0->GetBinCenter(max_bin1);
    double mode2 = h0->GetBinCenter(max_bin2);
    double dipBinCenter = GetDipBinCenter(h,mode1);

    pars.mean1 = mode1;
    if ((mode1-mode2)>3.0) pars.mean2 = mode1 - 3.0;
    else pars.mean2 = mode2;
    if ((dipBinCenter <= pars.mean2) || (dipBinCenter >= (pars.mean1+pars.mean2)*0.8)) pars.dip = (pars.mean1+pars.mean2)/2;    
    else pars.dip = dipBinCenter;

    delete h0;
    return pars;
}

vector<vector<int>> GetPP(TFile *f){ //get PP
    vector<vector<int>> vecIdxPPAllCol;
    for(int j=col;j<col+ncol;j++){
        vector<int> vecIdxPPcol;
        for (int i=1;i<125;i++) {
            TH1 *h = GetHistogram(f,j,i);
            int entries = h->GetEntries();
            double mode = GetMode(h);
            // cout << i << " ";
            if ((entries > numbOfEntriesLims) && mode<9000.0){    
                cout << i << " " << entries << " " << GetMode(h) <<  endl;
                vecIdxPPcol.push_back(i);
            }
        }
        vecIdxPPAllCol.push_back(vecIdxPPcol);
    }
    return vecIdxPPAllCol;
}

// vector<int> GetNumberOfPP(TFile *f){ //get number of PP
//     vector<int> vecIterator;
//     for(int j=col;j<col+ncol;j++){
//         int iterator=0;
//         int ppLLims=-1;
//         vector<double> vecMode;
//         for (int i=10;i<125;i++) {
//             // if (iterator == 8) break;
//             TH1 *h = GetHistogram(f,j,i);
//             int entries = h->GetEntries();
//             double mode = GetMode(h);
//             cout << entries << " " << GetMode(h) <<  endl;
//             if ((entries > numbOfEntriesLims) && (ppLLims==-1) && mode<9000.0){
//                 ppLLims = i+1;
//             }
//             if (ppLLims != -1){
//                 if (entries < numbOfEntriesLims-(0.1*numbOfEntriesLims)) break;
//                 else {
//                     if (mode<9000.0) {
//                         iterator++;
//                         if (iterator<11) vecMode.push_back(mode);
//                     }
//                     else {
//                         if (iterator<11){
//                             iterator=0;
//                             ppLLims=i+1;
//                             vecMode.clear();
//                         }
//                         else break;
//                     }
//                 }
//             }
//         }
//         iterator--;
//         int maxModeIdx = max_element(vecMode.begin(),vecMode.end()) - vecMode.begin();
//         ppLowLims[j-col] = ppLLims+maxModeIdx;
//         ppUpLims[j-col] = ppLLims+iterator;
//         iterator -= maxModeIdx;
//         // iterator--;
//         vecIterator.push_back(iterator);
//     }
//     return vecIterator;
// }

double GetDipBinCenter(TH1 *h, double mode){
    int upperBin = 1+(mode - h->GetXaxis()->GetXmin())/h->GetXaxis()->GetBinWidth(0);
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

// void WriteTWFitResults(ofstream &fouttw, TGraph *gr, Int_t col){
//     TString separator = ",";

//     Float_t xLowLims = gr->GetXaxis()->GetXmin();
//     Float_t xUpLims = gr->GetXaxis()->GetXmax();

//     TF1 *twFitFunction = new TF1("twFitFunction","[0] + [1]*((1/([4]*x+[3]))^[2])",xLowLims,xUpLims);
//     twFitFunction->SetParameter(0,-0.0661505);
//     twFitFunction->SetParameter(1,5.55105e+08);
//     twFitFunction->SetParameter(2,42.7598);
//     twFitFunction->SetParameter(3,1.43117);
//     twFitFunction->SetParameter(4,0.000330912);

//     TFitResultPtr TWFitResultPtr = gr->Fit("twFitFunction","MS");
//     Double_t p0 = TWFitResultPtr->Value(0);
//     Double_t p1 = TWFitResultPtr->Value(1);
//     Double_t p2 = TWFitResultPtr->Value(2);
//     Double_t p3 = TWFitResultPtr->Value(3);
//     Double_t p4 = TWFitResultPtr->Value(4);
//     Double_t chi2tw = TWFitResultPtr->Chi2();
//     Int_t fitStatus = TWFitResultPtr;

//     fouttw << col << separator << p0 << separator << p1 << separator << p2 << separator << p3 << separator << chi2tw << separator << fitStatus << endl;
//     delete twFitFunction;
// }
//  separator << p4