
using namespace RooFit;

struct fitResults
{
    int ph_bin;
    float mean;
};

TH1* GetHistogram(TFile *f, int col, int ph_bin);
fitResults WriteGaussianFitResults(ofstream &fout, TH1 *h, int col, int ph_bin,
                                    double deltaMean,
                                    double mean1val,
                                    double sig1,
                                    double sig2,
                                    double c2);
double GetMode(TH1 *h);

//some parameters to configure
const Int_t col=40;//first column to fit
const Int_t ncol=1;//number of column to fit
const Int_t ph_binlow[ncol]= {28
                        // 20,23,24,20,23,//1-5
                        // 15,19,20,25,26,//6-10
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
                        };//96-100
const Int_t n = 30;
Int_t ph_binup[ncol];
const Float_t dph_bin=16;
const Int_t dt_low=-3;
const Int_t dt_up=10;

bool doVariation = false;
const int nVar = 10;
double deltaMeanVar[4] = {0,1,2,3};
double mean1Var[nVar];
double sig1Var[nVar];
double sig2Var[nVar];
double c2Var[nVar];


double mean1Min = 0.;
double mean1Max = 2.;
double sig1min = 0.;
double sig1max = 2.;
double sig2min = 0.;
double sig2max = 5.;
double c2min = 0.;
double c2max = 0.5;

int doVariationGaussianFit(TString rootFile="root/hd_root-r72369.root") {
    
    for(int icol=0;icol<ncol;icol++) ph_binup[icol] = ph_binlow[icol]+n;
    
    double deltaMean1 = (mean1Max-mean1Min)/nVar;
    double deltaSig1 = (sig1max-sig1min)/nVar;
    double deltaSig2 = (sig2max-sig2min)/nVar;
    double deltaC2 = (c2max-c2min)/nVar;

    for(int i=0;i<nVar;i++){
        sig1Var[i] = sig1min + i*deltaSig1;
        sig2Var[i] = sig2min + i*deltaSig2;
        mean1Var[i] = mean1Min + i*deltaMean1;
        c2Var[i] = c2min + i*deltaC2;
    }

    TFile *f = new TFile(rootFile,"read");

    for(int j=col;j<col+ncol;j++){
        Float_t mean[n];
        Float_t ph_bina[n];
        stringstream ss; ss << j;
        for (int i=ph_binlow[j-col];i<=ph_binup[j-col];i++)
        {
            stringstream sph; sph << (i-0.5)*dph_bin;
            system("mkdir -p variation-gaussian-fits-csv/col_"+TString(ss.str()));
            ofstream fout; fout.open("variation-gaussian-fits-csv/col_"+TString(ss.str())+"/"+ TString(sph.str()) +".csv");
            for (size_t ideltaMeanVar = 0; ideltaMeanVar < 4; ideltaMeanVar++)
            {
                for (size_t isig1Var = 0; isig1Var < nVar; isig1Var++)
                {
                    for (size_t isig2Var = 0; isig2Var < nVar; isig2Var++)
                    {
                        for (size_t imean1Var = 0; imean1Var < nVar; imean1Var++)
                        {
                            for (size_t ic2Var = 0; ic2Var < nVar; ic2Var++)
                            {
                                TH1 *h = GetHistogram(f,j,i);
                                fitResults fR = WriteGaussianFitResults(fout,h,j,i,
                                deltaMeanVar[ideltaMeanVar],
                                mean1Var[imean1Var],
                                sig1Var[isig1Var],
                                sig2Var[isig2Var],
                                c2Var[ic2Var]
                                );
                                delete h;
                            }
                        }
                    }
                }
            }
            fout.close();
        }
        
    }
    return 0;
}

TH1* GetHistogram(TFile *f, int col, int ph_bin){
    stringstream ss_c; ss_c << col;
    TH2I *h2 = (TH2I*)f->Get("TAGM_TW/t-rf/h_dt_vs_pp_"+TString(ss_c.str()));
    stringstream ss_ph; ss_ph << h2->GetXaxis()->GetBinCenter(ph_bin);
    TH1 *h;
    h = h2->ProjectionY(TString(h2->GetName())+"_"+TString(ss_ph.str()),ph_bin,ph_bin);
    h->GetXaxis()->SetRangeUser(dt_low,dt_up);
    h->SetTitle("TAGM column "+TString(ss_c.str())+", Pulse-peak "+TString(ss_ph.str()));
    h->GetXaxis()->SetTitle("time(TDC) - time(RF) [ns]");
    delete h2;
    return h;
}

fitResults WriteGaussianFitResults(ofstream &fout, TH1 *h, int col, int ph_bin,
                                    double deltaMean,
                                    double mean1val,
                                    double sig1,
                                    double sig2,
                                    double c2)
{
    TString sep = ",";
    
    RooRealVar x("TDC time difference","TDC time difference [ns]",dt_low,dt_up);
    RooDataHist data("data","data",RooArgList(x),h);
    // model: double gaussian
    
    RooRealVar mean1("mean1","mean1",mean1val,dt_low,dt_up);
    RooRealVar mean2("mean2","mean2",mean1val-deltaMean,dt_low,dt_up);
    RooRealVar sigma1("sigma1","sigma1",sig1,sig1min,sig1max);
    RooGaussian gauss1("gauss1","gauss1",x,mean1,sigma1);
    RooRealVar sigma2("sigma2","sigma2",sig2,sig2min,sig2max);
    RooGaussian gauss2("gauss2","gauss2",x,mean2,sigma2);
    // f1: fraction of entries in first gaussian
    RooRealVar f1("f1","f1",0.9,0.1,1.0);
    RooRealVar f2("f2","f2",0.1,0.01,c2);
    RooAddPdf doubleGauss("doubleGauss","doubleGauss",RooArgList(gauss1,gauss2),RooArgList(f1,f2));
    
    doubleGauss.fitTo(data,PrintLevel(-1));//
    RooPlot *plot = x.frame();
    data.plotOn(plot);
    doubleGauss.plotOn(plot);
    fout << col << sep << ph_bin << sep << deltaMean << sep << mean1val << sep << sig1 << sep << sig2 << sep << c2 << sep << plot->chiSquare() << endl;
    
    // delete canvas;
    delete plot;

    fitResults fResults;
    
    fResults.mean = mean1.getVal();
    fResults.ph_bin = ph_bin;
    return fResults;
}