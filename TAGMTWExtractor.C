

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

TH1* GetHistogram(TFile *f, int col, int ph_bin);
gaussianFitResults DoGaussianFit(ofstream &fout, TH1 *h, int col, int ph_bin, TString resultSubFolder);
double GetMode(TH1 *h);
parsForTripleGaussian GetParsForTripleGaussian(TH1 *h);
double GetDipBinCenter(TH1 *h, double mode);
vector<vector<int>> GetPP(TFile *f);
vector<int> IdentifyBadFit(Int_t n, Float_t* means,Float_t* errors);
int GetGraphLowerLimsIdx(const int n, Float_t* means, Float_t* pps);

//some parameters to configure
const Int_t col=1;//first column to fit
const Int_t ncol=102;//number of column to fit
Int_t ppLowLims[ncol];

int numbOfEntriesLims;
const int minModeBinContentSelPP = 0;

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

const bool useMinos = kTRUE;
const double chi2DiffThres = 0.6;
const double maxModeBinContentTrpGauss = 1000;

int TAGMTWExtractor(TString runNumber="72369") {

    TString rootFile=rootFileFolder+rootFilePrefix+runNumber+".root";

    gROOT->SetBatch(kTRUE);

    TFile *f = new TFile(rootFile,"read");
    
    //select pulse peak to fit for all columns
    vector<vector<int>> SelectedPP = GetPP(f);

    TString makeDir = "mkdir -p ";
    makeDir += resultFolder; 
    system(makeDir);
    TCanvas *c0 = new TCanvas();
    TString resultSubFolder = resultFolder+runNumber+"/";
    TGraphErrors *g[ncol];
    TString rootFileOutputName = resultFolder+rootFileOutputPrefix+runNumber+".root";
    TFile *fOutput = new TFile(rootFileOutputName,"recreate");

    //loop over columns
    for(int j=col;j<col+ncol;j++){
        const int n = SelectedPP[j-col].size(); //number of pulse peak to fit on each column
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
            //loop over selected pulse peak
            for (int i=0;i<n;i++) {
                int idxPP = SelectedPP[j-col][i];
                //Get histogram of selected pulse peak
                TH1 *h = GetHistogram(f,j,idxPP);
                h->Draw();

                //Do gaussian fit
                gaussianFitResults fR = DoGaussianFit(fout,h,j,idxPP,resultSubFolder);

                //Choose gaussian mean and standard deviation
                if (i==0){
                    if (fR.cMean2 > 0.6) {
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

            //Discard bad fits
            vector<int> badFit = IdentifyBadFit(n,meanGraph,meanErrorGraph);
            const int cleanedN = n - badFit.size();
            Float_t cleanedMeanGraph[cleanedN];
            Float_t cleanedMeanErrorGraph[cleanedN];
            Float_t cleanedPPGraph[cleanedN];

            int nDifference = 0;
            for (int z=0;z<n;z++){
                bool isDiscarded = false;
                for (int y=0;y<badFit.size();y++){
                    if (z==badFit[y]){
                        isDiscarded = true;
                        nDifference++;
                    }
                }
                if (!isDiscarded) {
                    cleanedMeanGraph[z-nDifference] = meanGraph[z];
                    cleanedMeanErrorGraph[z-nDifference] = meanErrorGraph[z];
                    cleanedPPGraph[z-nDifference] = ppGraph[z];
                }
            }
            // cout << cleanedN << " " << badFit.size() << endl;

            if (cleanedN<=0) {
                //dummy graph for empty pulse peak after cleaning
                Double_t xempty[6] = {200,300,400,500,600,700};
                Double_t yempty[6] = {0.,0.,0.,0.,0.,0.};
                g[j-col] = new TGraphErrors(6,xempty,yempty);
            }
            else{
                //Re-set lower limit after discarding bad fits
                int newLowerLims = GetGraphLowerLimsIdx(cleanedN,cleanedMeanGraph,cleanedPPGraph);
                if (newLowerLims != 0) {
                    const int finalN = cleanedN-newLowerLims;
                    Float_t finalMeanGraph[finalN];
                    Float_t finalMeanErrorGraph[finalN];
                    Float_t finalPPGraph[finalN];
                    for (int z=newLowerLims;z<cleanedN;z++){
                        finalMeanGraph[z-newLowerLims] = cleanedMeanGraph[z];
                        finalMeanErrorGraph[z-newLowerLims] = cleanedMeanErrorGraph[z];
                        finalPPGraph[z-newLowerLims] = cleanedPPGraph[z];
                    }
                    g[j-col] = new TGraphErrors(finalN,finalPPGraph,finalMeanGraph,0,finalMeanErrorGraph);
                }
                else g[j-col] = new TGraphErrors(cleanedN,cleanedPPGraph,cleanedMeanGraph,0,cleanedMeanErrorGraph);
            }
        }
        g[j-col]->SetMarkerSize(1);
        g[j-col]->SetMarkerStyle(kStar);
        g[j-col]->GetXaxis()->SetTitle("Pulse peak (adc counts)");
        g[j-col]->GetYaxis()->SetTitle("T - RF (ns)");
        TString namegrTitle = "Timewalk col ";
        namegrTitle+=j;
        g[j-col]->SetTitle(namegrTitle);
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
        fout.close();
    }
    // fouttw.close();
    return 0;
}

//Get histogram
TH1* GetHistogram(TFile *f, int col, int ph_bin) {
    stringstream ss_c; ss_c << col;
    TH2I *h2 = (TH2I*)f->Get("TAGM_TW/tdc-rf/h_dt_vs_pp_tdc_"+TString(ss_c.str()));
    // h2->Draw("colz");
    numbOfEntriesLims = 0.003*h2->GetEntries();
    stringstream ss_ph; ss_ph << h2->GetXaxis()->GetBinCenter(ph_bin);
    TH1 *h;
    h = h2->ProjectionY(TString(h2->GetName())+"_"+TString(ss_ph.str()),ph_bin,ph_bin);
    h->GetXaxis()->SetRangeUser(dtLowLims,dtUpLims);
    
    // h->Draw();
    h->SetTitle("TAGM column "+TString(ss_c.str())+", Pulse-peak "+TString(ss_ph.str()));
    h->GetXaxis()->SetTitle("time(TDC) - time(RF) [ns]");
    return h;
}

//Gaussian fitter
gaussianFitResults DoGaussianFit(ofstream &fout, TH1 *h, int col, int ph_bin,TString resultSubFolder) {
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

    double sgnMean = GetMode(h); //get maximum bin center
    TCanvas *canvas = new TCanvas("c","c",800,500);
    RooPlot *plotDouble = x.frame();
    plotDouble->SetTitle(h->GetTitle());
    plotDouble->SetXTitle(h->GetXaxis()->GetTitle());
    plotDouble->SetTitleOffset(1.1,"Y");

    //double gaussian fit
    RooRealVar meanSgnDouble("meanSgnDouble","meanSgnDouble",sgnMean,sgnMean-0.5,sgnMean+0.5);
    RooRealVar sigmaSgnDouble("sigmaSgnDouble","sigmaSgnDouble",sgnSig,sgnSigMin,sgnSigMax);
    RooRealVar meanBkg1Double("meanBkg1Double","meanBkg1Double",sgnMean,sgnMean-2.0,sgnMean+2.0);
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

    //triple gaussian fit
    //determine initial parameters of tripe gaussian fitting
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
    //logic for choosing choosing appropriate result (double/triple), need more evaluation
    bool isTripleFunctionUsed = (diffMeanTriple > 1.0) && (cBkg2Triple.getVal() < 0.4) && (diffChi2Fit > chi2DiffThres) && (((chi2FitTriple < chi2FitDouble) && (modeBinContent < maxModeBinContentTrpGauss) && (cBkg1Double.getVal() >= 0.04)) || ((cBkg1Double.getVal() < 0.04) && (chi2FitTriple < chi2FitDouble))); //7XXXX  || ((chi2FitDouble > 3.0))
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

double GetMode(TH1 *h) {
    if (h->GetEntries() == 0.0) return 9999.0;
    int max_bin = h->GetMaximumBin();
    double max = h->GetBinContent(max_bin);
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
    double xmin = dtLowLims;
    double xmax = dtUpLims;
    double xmid = (xmax + xmin)/2;
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

        // cout << entries << " " << (max1-max2) << " " << max1 << " " << max2 << " " << h0->GetBinCenter(max_bin1) << " " << h0->GetBinCenter(max_bin2) << endl;
        iter++;
    } while (iter < 2 && ((((max1-max2) < -0.5*max2)) || (max_bin1-max_bin2 < 30))); //this logic need more evaluation, test on different kind of dataset
    
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
            if (((entries > numbOfEntriesLims) && (mode < 9000.0) && (entries > 300)) || entries > 4000.0){    
                cout << i << " " << entries << " " << GetMode(h) <<  endl;
                vecIdxPPcol.push_back(i);
            }
        }
        vecIdxPPAllCol.push_back(vecIdxPPcol);
    }
    return vecIdxPPAllCol;
}

//determine bin center between two peaks
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

//identify bad fits
vector<int> IdentifyBadFit(Int_t n, Float_t* means,Float_t* errors){
    vector<int> badFitIdx;
    for (Int_t i=0;i<n;i++) {
        if (errors[i] > 0.1){
            badFitIdx.push_back(i);
        }
    }
    return badFitIdx;
}

//Set new lower limit at discontinuity
int GetGraphLowerLimsIdx(const int n, Float_t* means, Float_t* pps){
    Float_t grads[n-1];
    vector<int> badGrads;
    
    if (n<10) return 0;
    for (int i=0;i<n-1;i++) {
        grads[i] = (means[i+1]-means[i])/(pps[i+1]-pps[i]);
        // cout << i << " " << grads[i] << endl;
    }
    for (int i=0;i<10;i++) if ((grads[i]-grads[i+1])>0.003 && (grads[i+1] < -0.0001) && (grads[i]>-0.006)) badGrads.push_back(i+1);
    if (badGrads.size() != 0) return badGrads[badGrads.size()-1];
    else return 0;
}