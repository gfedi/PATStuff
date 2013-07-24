{
gROOT->Reset();
gSystem->Load("libRooFit");
using namespace RooFit ;
#include </usr/users/ulmerk/stuff/tdrstyle.C>
//#include "RooKeysPdf.h"
gROOT->Reset();
setTDRStyle();
tdrStyle->SetPadRightMargin(0.03);
tdrStyle->SetMarkerSize(0.5);

bool doKsPsiPDFs = false;
bool doSXF = false;
bool usePerEventErrors = false;
bool fitData = true;
bool fitMC = false;
bool doToys = false;
bool doPtBins = true;
bool doYBins = true;
bool floatKsMass = false;
bool useDavidBinning = false;
int nToys = 1;

TFile seleBB("seleJPsiKs_BJPsiXCombo_Fall10_7TeV_386.root");
TTree *treeBB=tofit;

TFile seleSig("seleJPsiKs_B0_sigMC-all.root");
TTree *treeSigPrivate=tofit;


RooRealVar* BMass = new RooRealVar("bMass","m_{B}",4.9,5.7,"GeV");
RooRealVar* BMass_sub = new RooRealVar("bMass","m_{B}",5.4,5.7,"GeV");
//RooRealVar* Bct = new RooRealVar("bctauMPVRf","bctauMPVRf",-0.05,0.4,"cm");
//RooRealVar* BctE = new RooRealVar("bctauMPVRfE","bctauMPVRfE",0.001,0.10,"cm");  
//RooRealVar* Bct = new RooRealVar("bctauRf","bctauRf",-0.05,0.4,"cm");
//RooRealVar* BctE = new RooRealVar("bctauRfE","bctauRfE",0.0,0.02,"cm"); 
RooRealVar* Bct = new RooRealVar("bctau2D","ct",-0.05,0.4,"cm");
RooRealVar* BctE = new RooRealVar("bctau2DE","bctau2DE",0.0,0.02,"cm"); 
RooRealVar* PsiMass = new RooRealVar("JMass","JPsi mass",2.946916,3.246916,"GeV");
RooRealVar* VMass = new RooRealVar("VMass","Ks mass",0.478,0.518,"GeV");
RooRealVar* genKsPsi = new RooRealVar("genKsPsi","genKsPsi",-2,2);
RooRealVar* genKsPsi2 = new RooRealVar("genKsPsi2","genKsPsi2",-2,2);
RooRealVar* truthMatch = new RooRealVar("truthMatch","truthMatch",-2,2);
RooRealVar* feedup = new RooRealVar("feedup","feedup",-2,20);
RooRealVar* feeddown = new RooRealVar("feeddown","feeddown",-2,20);
RooRealVar* hlt_2mu3 = new RooRealVar("hlt_2mu3","hlt_2mu3",-2,2);
RooRealVar* hlt_2mu0 = new RooRealVar("hlt_2mu0","hlt_2mu0",-2,2);
RooRealVar* hlt_mu9 = new RooRealVar("hlt_mu9","hlt_mu9",-2,2);
RooRealVar* hlt_2muOpen = new RooRealVar("hlt_L12muOpen","hlt_L12muOpen",-2,2);
RooRealVar* hlt_2muOpenTight = new RooRealVar("hlt_L12muOpenTight","hlt_L12muOpenTight",-2,2);
RooRealVar* hlt_mu0trk0 = new RooRealVar("hlt_mu0trk0","hlt_mu0trk0",-2,2);
RooRealVar* hlt_mu0trkmu0 = new RooRealVar("hlt_mu0trkmu0","hlt_mu0trkmu0",-2,2);
RooRealVar* muTrigmu0trk0 = new RooRealVar("muTrigmu0trk0","muTrigmu0trk0",-2,2);
RooRealVar* muTrigmu0trkmu0 = new RooRealVar("muTrigmu0trkmu0","muTrigmu0trkmu0",-2,2);
RooRealVar* mumTrig2mu0 = new RooRealVar("mumTrig2mu0","mumTrig2mu0",-2,2);
RooRealVar* mupTrig2mu0 = new RooRealVar("mupTrig2mu0","mupTrig2mu0",-2,2);
RooRealVar* mumTrigL1Open = new RooRealVar("mumTrigL1Open","mumTrigL1Open",-2,2);
RooRealVar* mupTrigL1Open = new RooRealVar("mupTrigL1Open","mupTrigL1Open",-2,2);  
RooRealVar* bPx = new RooRealVar("bPx","bPx",-1000.,1000.);
RooRealVar* bPy = new RooRealVar("bPy","bPy",-1000.,1000.);
RooRealVar* bY = new RooRealVar("bY","bY",-1000.,1000.);

RooArgSet Vars;
Vars.add(*BMass);
Vars.add(*Bct);
Vars.add(*BctE);
Vars.add(*PsiMass);
Vars.add(*VMass);
Vars.add(*genKsPsi);
Vars.add(*genKsPsi2);
Vars.add(*truthMatch);
Vars.add(*feedup);
Vars.add(*feeddown);
Vars.add(*hlt_2mu3);
Vars.add(*hlt_2mu0);
Vars.add(*hlt_mu9);
Vars.add(*hlt_2muOpen);
Vars.add(*hlt_2muOpenTight);
Vars.add(*hlt_mu0trk0);
Vars.add(*hlt_mu0trkmu0);
Vars.add(*muTrigmu0trk0);
Vars.add(*muTrigmu0trkmu0);
Vars.add(*mumTrig2mu0);
Vars.add(*mupTrig2mu0);
Vars.add(*mumTrigL1Open);
Vars.add(*mupTrigL1Open);
Vars.add(*bPx);
Vars.add(*bPy);
Vars.add(*bY);


TString trigger = TString("hlt_2mu0_nominal_plots");
gSystem->mkdir("JPsiKsFit_dataRes/"+trigger);

TString triggerCut = TString("mumTrig2mu0==1&&mupTrig2mu0==1");

ofstream outTextFile;
outTextFile.open ("JPsiKsFit_dataRes/"+trigger+"/JPsiKs_Results.txt");

cout << "Using trigger = " << trigger << endl;
outTextFile << "Using trigger = " << trigger << endl;


TString bptCut = TString("&&sqrt(bPx*bPx+bPy*bPy)>5&&abs(bY)<2.2");

RooDataSet MC_BB_noHLT("MC_BB_noHLT", "MC_BB_noHLT", treeBB, Vars);
RooDataSet* MC_BB = MC_BB_noHLT->reduce(triggerCut+bptCut);

RooDataSet* MC_Sig = MC_BB->reduce("genKsPsi==1 && truthMatch==1");
RooDataSet* MC_Sxf = MC_BB->reduce("genKsPsi==1 && truthMatch!=1");
RooDataSet* MC_BkgBB = MC_BB->reduce("((genKsPsi!=1 && feedup!=2 && feeddown!=2 && feeddown!=3 && feeddown!=6)||((genKsPsi==1||feedup==2||feeddown==2||feeddown==3||feeddown==6)&&truthMatch!=1))");

RooDataSet* MC_Feed = MC_BB->reduce("(feedup==2||feeddown==2||feeddown==3||feeddown==6)&&truthMatch==1");

RooDataSet MC_SigPrivate_noHLT("MC_SigPrivate_noHLT", "MC_SigPrivate_noHLT", treeSigPrivate, Vars);
RooDataSet* MC_SigPrivate = MC_SigPrivate_noHLT->reduce(triggerCut+bptCut+"&&genKsPsi2==1&&truthMatch==1");


RooRealVar* fracBCore = new RooRealVar("fracBCore","fracBCore",0.8,0.5,1.0);
RooRealVar* fracBOut = new RooRealVar("fracBOut","fracBOut",0.05,0.00,0.10);

RooRealVar* meanBCore = new RooRealVar("meanBCore","meanBCore",5.279,5.2,5.5);
RooRealVar* widthBCore = new RooRealVar("widthBCore","widthBCore",0.012,0.001,0.4);
RooRealVar* meanBOut = new RooRealVar("meanBOut","meanBOut",5.279,5.2,5.5);
RooRealVar* widthBOut = new RooRealVar("widthBOut","widthBOut",0.10,0.001,0.4);
RooRealVar* meanBTail = new RooRealVar("meanBTail","meanBTail",5.279,5.2,5.5);
RooRealVar* widthBTail = new RooRealVar("widthBTail","widthBTail",0.2,0.001,0.4);

RooGaussian bSigCore("bSigCore","bSigCore",*BMass,*meanBCore,*widthBCore);
RooGaussian bSigOut("bSigOut","bSigOut",*BMass,*meanBOut,*widthBOut);
RooGaussian bSigTail("bSigTail","bSigTail",*BMass,*meanBTail,*widthBTail);
RooAddPdf bSig("bSig","bSig",RooArgList(bSigCore,bSigTail),RooArgList(*fracBCore));

TCanvas *PDFs = new TCanvas("PDFs","PDFs",800,800);
PDFs->Divide(2,4);

TCanvas *MCPDFs = new TCanvas("MCPDFs","MCPDFs",800,800);
MCPDFs->Divide(2,5);

TCanvas *MCPDFs_all = new TCanvas("MCPDFs_all","MCPDFs_all",800,300);
MCPDFs_all->Divide(2,1);

RooFitResult* fitres= bSig.fitTo(*MC_SigPrivate);

TCanvas *bmass = new TCanvas("bmass","bmass",800,800);
bmass->Divide(2,3);
bmass->cd(1);
RooPlot* frame1=BMass->frame();
MC_SigPrivate->plotOn(frame1,(200));
bSig.plotOn(frame1,LineColor(4));
frame1->Draw();
bSig.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
frame1->GetXaxis()->SetRange(35,65);
frame1->Draw();
TPaveText* pbox = (TPaveText*) frame1->findObject("bSig_paramBox");
if (pbox) {
  pbox->AddText(Form("%s%.3f", "#chi^{2}/DOF = ", frame1->chiSquare()));
}
PDFs->cd(1);
frame1->Draw();

MCPDFs_all->cd(1);
frame1->Draw();

fracBCore.setConstant(kTRUE);
fracBOut.setConstant(kTRUE);
meanBCore.setConstant(kTRUE);
meanBTail.setConstant(kTRUE);
meanBOut.setConstant(kTRUE);
widthBCore.setConstant(kTRUE);
widthBTail.setConstant(kTRUE);
widthBOut.setConstant(kTRUE);

//  make signal lifetime pdf
// Decay (x) res model 

cout << "Starting ctau signal PDF" << endl;  

if (usePerEventErrors) {
  //////////////////////////////////////////////////////////////////////////
  // Use per event errors for signal lifetime PDF

  RooRealVar* meanB = new RooRealVar("meanB","meanB",0.0,-1.0,1.0);
  RooRealVar* scaleB = new RooRealVar("scaleB","scaleB",1.0,0.05,2.2);
  RooGaussModel resB("resB","resB",*Bct,*meanB,*scaleB,RooConst(1.0),*BctE) ; 

  RooRealVar* ctauB= new RooRealVar("ctauB","ctauB",0.046,0,0.1) ; 
  RooDecay decaySig("decaySig","decaySig",*Bct,*ctauB,resB,RooDecay::SingleSided) ;  

  RooFitResult* fitres= decaySig.fitTo(*MC_SigPrivate,ConditionalObservables(*BctE));

  /////////////////////////////////////////////////////////////////////////
} else {
  RooRealVar* meanResCore = new RooRealVar("meanResCore","meanResCore",0.0,-0.02,0.02);
  RooRealVar* sigmaResCore = new RooRealVar("sigmaResCore","sigmaResCore",0.002,0.0001,0.1);
  RooGaussModel resCore("resCore","resCore",*Bct,*meanResCore,*sigmaResCore); 

  RooRealVar* meanResTail = new RooRealVar("meanResTail","meanResTail",0.0,-0.01,0.01);
  RooRealVar* sigmaResTail = new RooRealVar("sigmaResTail","sigmaResTail",0.002,0.0000001,0.1);
  RooGaussModel resTail("resTail","resTail",*Bct,*meanResTail,*sigmaResTail); 

  RooRealVar* fracResCore = new RooRealVar("fracResCore","fracResCore",0.80,0.1,1.0);

  RooAddModel resDG("resDG","resDG",RooArgList(resCore, resTail), RooArgList(*fracResCore));

  RooRealVar* ctauB= new RooRealVar("ctauB","ctauB",0.046,0,0.1) ; 

  RooDecay decaySig("decaySig","decaySig",*Bct,*ctauB,resDG,RooDecay::SingleSided) ;  

  RooFitResult* fitres= decaySig.fitTo(*MC_SigPrivate);
}

TCanvas *ctau = new TCanvas("ctau","ctau",800,800);
ctau->Divide(2,3);
ctau->cd(1).SetLogy();  
//ctau->cd(1);
RooPlot* frame1=Bct->frame();
MC_SigPrivate->plotOn(frame1,Binning(45));
if (usePerEventErrors) {
  decaySig.plotOn(frame1,ProjWData(*MC_SigPrivate),LineColor(4));
} else {
  decaySig.plotOn(frame1,LineColor(4));
}
frame1->SetMinimum(0.1);
frame1->Draw();
decaySig.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
frame1->Draw();
TPaveText* pbox = (TPaveText*) frame1->findObject("decaySig_paramBox");
if (pbox) {
  pbox->AddText(Form("%s%.3f", "#chi^{2}/DOF = ", frame1->chiSquare()));
}

PDFs->cd(2).SetLogy();
frame1->Draw();

//////////////////////////
// make background pdfs for feed
//////////////////////////
//

//get b mass background pdf


RooRealVar* fracBCoreFeed = new RooRealVar("fracBCoreFeed","fracBCoreFeed",0.75,0.3,1.);
RooRealVar* fracBOutFeed = new RooRealVar("fracBOutFeed","fracBOutFeed",0.007,0.,0.1);
  
RooRealVar* meanBCoreFeed = new RooRealVar("meanBCoreFeed","meanBCoreFeed",5.00,4.9,6.0);
RooRealVar* widthBCoreFeed = new RooRealVar("widthBCoreFeed","widthBCoreFeed",0.06,0.001,2.0);
RooGaussian gausBCoreFeed("gausBCoreFeed","gausBCoreFeed",*BMass,*meanBCoreFeed,*widthBCoreFeed);

RooRealVar* meanBTailFeed = new RooRealVar("meanBTailFeed","meanBTailFeed",5.11,4.9,6.0);
RooRealVar* widthBTailFeed = new RooRealVar("widthBTailFeed","widthBTailFeed",0.03,0.001,2.0);
RooGaussian gausBTailFeed("gausBTailFeed","gausBTailFeed",*BMass,*meanBTailFeed,*widthBTailFeed);

RooRealVar* meanBOutFeed = new RooRealVar("meanBOutFeed","meanBOutFeed",5.28,4.9,6.0);
RooRealVar* widthBOutFeed = new RooRealVar("widthBOutFeed","widthBOutFeed",0.025,0.001,2.0);
RooGaussian gausBOutFeed("gausBOutFeed","gausBOutFeed",*BMass,*meanBOutFeed,*widthBOutFeed);

RooAddPdf bFeed("bFeed","bFeed",RooArgList(gausBCoreFeed,gausBOutFeed,gausBTailFeed),RooArgList(*fracBCoreFeed, *fracBOutFeed));
//  RooFitResult* fitres= bFeed.fitTo(*MC_Feed,Hesse(true),Minos(true));    
RooFitResult* fitres= bFeed.fitTo(*MC_Feed,Hesse(true));  
  
//try using RooKeysPDFs?
//RooKeysPDF bFeed("bFeed","bFeed",*BMass, *MC_Feed);


RooFitResult* fitres= bFeed.fitTo(*MC_Feed,Hesse(true),Minos(true));  

bmass->cd(2);
RooPlot* frame1=BMass->frame();
MC_Feed->plotOn(frame1,Binning(40));
bFeed.plotOn(frame1,LineColor(4));
frame1->Draw();
bFeed.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
frame1->Draw();
TPaveText* pbox = (TPaveText*) frame1->findObject("bFeed_paramBox");
if (pbox) {
  pbox->AddText(Form("%s%.3f", "#chi^{2}/DOF = ", frame1->chiSquare()));
}

PDFs->cd(3);
frame1->Draw();

MCPDFs_all->cd(2);
frame1->Draw();

meanBOutFeed.setConstant(kTRUE);
fracBOutFeed.setConstant(kTRUE);  
fracBCoreFeed.setConstant(kTRUE);
meanBCoreFeed.setConstant(kTRUE);
widthBCoreFeed.setConstant(kTRUE);
meanBTailFeed.setConstant(kTRUE);
widthBTailFeed.setConstant(kTRUE);
widthBOutFeed.setConstant(kTRUE);

//  make Feed background lifetime pdf
  
//Gaussian convoluted decay component

RooRealVar* ctauFeed = new RooRealVar("ctauFeed","ctauFeed",0.05,0,1.0) ; 
RooDecay decayFeed("decayFeed","decayFeed",*Bct,*ctauFeed,resDG,RooDecay::SingleSided) ; 

RooFitResult* fitres= decayFeed.fitTo(*MC_Feed,Hesse(true),Minos(true));

ctau->cd(2).SetLogy();
RooPlot* frame1=Bct->frame();
MC_Feed->plotOn(frame1,Binning(90));
decayFeed.plotOn(frame1,LineColor(4));
frame1->Draw();
decayFeed.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
frame1->SetMinimum(0.1);
frame1->Draw();
TPaveText* pbox = (TPaveText*) frame1->findObject("decayFeed_paramBox");
if (pbox) {
  pbox->AddText(Form("%s%.3f", "#chi^{2}/DOF = ", frame1->chiSquare()));
}

PDFs->cd(4).SetLogy();
frame1->Draw();

ctauFeed.setConstant(kTRUE);







///////////////////////////
// get b mass background pdf for prompt JPsi
///////////////////////////
  
TFile bkgPro("seleJPsiKs_prompt_Fall10_7TeV_386.root");
TTree *treeBkgPro=tofit;
  
RooDataSet MC_BkgPro_noHLT("MC_BkgPro_noHLT","MC_BkgPro_noHLT",treeBkgPro,Vars);
//RooDataSet* MC_BkgPro = MC_BkgPro_noHLT->reduce("hlt_2mu3==1 || hlt_mu9==1");
RooDataSet* MC_BkgPro = MC_BkgPro_noHLT->reduce(triggerCut+bptCut);

//get b mass background pdf

RooRealVar* shapePro = new RooRealVar("shapePro","shapePro",-1.,-10000,10000);

RooExponential bPro("bPro","bPro",*BMass,*shapePro);
  
RooFitResult* fitres= bPro.fitTo(*MC_BkgPro,Hesse(true),Minos(true));

bmass->cd(4);
RooPlot* frame1=BMass->frame();
MC_BkgPro->plotOn(frame1,Binning(16));
bPro.plotOn(frame1,LineColor(4));
frame1->Draw();
bPro.paramOn(frame1, Layout(0.65,0.95,0.45), Format("NELU",AutoPrecision(3)) ) ;
frame1->Draw();
TPaveText* pbox = (TPaveText*) frame1->findObject("bPro_paramBox");
if (pbox) {
  pbox->AddText(Form("%s%.3f", "#chi^{2}/DOF = ", frame1->chiSquare()));
}

PDFs->cd(7);
frame1->Draw();

shapePro.setConstant(kTRUE);


//  make prompt JPsi background lifetime pdf  
RooAddPdf decayPro("decayPro","decayPro",RooArgList(resCore, resTail), RooArgList(*fracResCore));

fracResCore.setVal(0.87);
meanResCore.setVal(0.0);
meanResTail.setVal(0.0);
sigmaResCore.setVal(0.0045);
sigmaResTail.setVal(0.013);

RooFitResult* fitres= resDG.fitTo(*MC_BkgPro,Hesse(true),Minos(true));

ctau->cd(4).SetLogy();
RooPlot* frame1=Bct->frame();
MC_BkgPro->plotOn(frame1,Binning(150));
decayPro.plotOn(frame1,LineColor(4));
frame1->Draw();
decayPro.paramOn(frame1, Layout(0.6,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
frame1->GetXaxis()->SetRange(0,40);
frame1->SetMinimum(0.1);
frame1->Draw();
TPaveText* pbox = (TPaveText*) frame1->findObject("decayPro_paramBox");
if (pbox) {
  pbox->AddText(Form("%s%.3f", "#chi^{2}/DOF = ", frame1->chiSquare()));
}

PDFs->cd(8).SetLogy();
frame1->Draw();  

fracResCore.setConstant(kTRUE);
meanResCore.setConstant(kTRUE);
meanResTail.setConstant(kTRUE);
sigmaResCore.setConstant(kTRUE);
sigmaResTail.setConstant(kTRUE);






//////////////////////////
// make background pdfs for non-peaking BBbar
//////////////////////////
//

//get b mass background pdf

RooRealVar* shapeBB = new RooRealVar("shapeBB","shapeBB",0.,-5,5);

RooExponential bBB("bBB","bBB",*BMass,*shapeBB);

RooFitResult* fitres= bBB.fitTo(*MC_BkgBB,Hesse(true),Minos(true));

bmass->cd(3);
RooPlot* frame1=BMass->frame();
MC_BkgBB->plotOn(frame1,Binning(20));
bBB.plotOn(frame1,LineColor(4));
frame1->Draw();
bBB.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
frame1->Draw();
TPaveText* pbox = (TPaveText*) frame1->findObject("bBB_paramBox");
if (pbox) {
  pbox->AddText(Form("%s%.3f", "#chi^{2}/DOF = ", frame1->chiSquare()));
}

PDFs->cd(5);
frame1->Draw();

shapeBB.setConstant(kTRUE);

//  make BB background lifetime pdf



RooRealVar* ctauBB = new RooRealVar("ctauBB","ctauBB",0.010,0,1.0) ; 
RooDecay decayBBConv("decayBBConv","decayBBConv",*Bct,*ctauBB,resDG,RooDecay::SingleSided) ; 

RooRealVar* ctauBB2 = new RooRealVar("ctauBB2","ctauBB2",0.044,0,1.0) ; 
RooDecay decayBBConv2("decayBBConv2","decayBBConv2",*Bct,*ctauBB2,resDG,RooDecay::SingleSided) ; 

RooRealVar* fracBBDecay = new RooRealVar("fracBBDecay","frac",0.61,0.0,1.);
RooAddPdf decayBB("decayBB","decayBB",RooArgSet(decayBBConv,decayBBConv2),*fracBBDecay);

//take resolution from prompt JPsi component to simplify the fit
meanResCore.setConstant(kTRUE); 
sigmaResCore.setConstant(kTRUE); 
meanResTail.setConstant(kTRUE); 
sigmaResTail.setConstant(kTRUE); 
fracResCore.setConstant(kTRUE); 

RooFitResult* fitres= decayBB.fitTo(*MC_BkgBB,Hesse(true),Minos(true));

ctau->cd(3).SetLogy();
RooPlot* frame1=Bct->frame();
MC_BkgBB->plotOn(frame1,Binning(90));
decayBB.plotOn(frame1,LineColor(4));
frame1->Draw();
decayBB.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
frame1->SetMinimum(0.1);
frame1->Draw();
TPaveText* pbox = (TPaveText*) frame1->findObject("decayBB_paramBox");
if (pbox) {
  pbox->AddText(Form("%s%.3f", "#chi^{2}/DOF = ", frame1->chiSquare()));
}

PDFs->cd(6).SetLogy();
frame1->Draw();






bmass->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKs_bMass_PDF.eps");
ctau->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKs_ctau_PDF.eps");
bmass->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKs_bMass_PDF.pdf");
ctau->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKs_ctau_PDF.pdf");


PDFs->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKs_PDFs.pdf");
MCPDFs_all->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKs_MCPDFs_all.pdf");


// Now construct 2D model of the PDFs  with Bmass and lifetime

RooProdPdf modelSig("modelSig","modelSig",RooArgList(bSig, decaySig));
RooProdPdf modelFeed("modelFeed","modelFeed",RooArgList(bFeed, decayFeed));
RooProdPdf modelBkgBB("modelBkgBB","modelBkgBB",RooArgList(bBB, decayBB));
RooProdPdf modelBkgPrompt("modelBkgPrompt","modelBkgPrompt",RooArgList(bPro, decayPro));

RooRealVar* nsig = new RooRealVar("nsig","nsig",900,0.0,50000);
RooRealVar* nfeed = new RooRealVar("nfeed","nfeed",100,0.0,50000);
RooRealVar* nbkgBB = new RooRealVar("nbkgBB","nbkgBB",300,0.0,50000);
RooRealVar* nbkgPro = new RooRealVar("nbkgPro","nbkgPro",50000,0.0,200000);

RooExtendPdf sig2De("sig2De", "extended gauss", modelSig, *nsig);
RooExtendPdf feed2De("feed2De", "extended gauss", modelFeed, *nfeed);
RooExtendPdf bkgBB2De("bkgBB2De", "extended gauss", modelBkgBB, *nbkgBB);
RooExtendPdf bkgPro2De("bkgPro2De", "extended gauss", modelBkgPrompt, *nbkgPro);

RooAddPdf model2D("model2D","model2D",RooArgList(sig2De,feed2De,bkgBB2De,bkgPro2De),RooArgList(*nsig,*nfeed,*nbkgBB,*nbkgPro));
RooAddPdf model2DSB("model2DSB","model2DSB",RooArgList(bkgBB2De,bkgPro2De),RooArgList(*nbkgBB,*nbkgPro));

BMass->setRange("HighMassSB",5.4,5.70);

TFile *outFile = new TFile("JPsiKsFit_dataRes/"+trigger+"/JPsiKsRooFit_PDFs.root","RECREATE");

model2D->Write();
outFile->Close();

//set yields to expected values to obtain a good fit
nsig->setVal( MC_Sig->numEntries() );
nfeed->setVal( MC_Feed->numEntries() );
nbkgBB->setVal( MC_BkgBB->numEntries() );
nbkgPro->setVal( 7*MC_BkgPro->numEntries() );


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//	 PDFs are finished. Now moving to fits
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////












if(fitMC) {

  cout << "doing MC fit" << endl;

  // construct a subsample dataset of the MC with the event yields from data
  int nSigData = 753;
  int nFeedData = 364;
  int nBBData = 4006;
  int nPromptData = 16696;

  RooDataSet* subSample(0);
  RooArgSet subSampleSet;
  subSampleSet = &Vars;
  subSample = new RooDataSet("subSample", "subSample", subSampleSet);

  int eventsSigSubsample = 0; int eventsFeedSubsample = 0; int eventsBBSubsample = 0; int eventsPromptSubsample = 0;

  //first select nSigData signal MC events
  for (Int_t j=0; j<(nSigData); j++) {
    Int_t I=RooRandom::randomGenerator()->Integer(MC_Sig->numEntries());
    const RooArgSet* row=MC_Sig->get(I);
    if(!row) {
      cout<<"skipping null row entry I="<<I<<endl;
      continue;
    }
    subSample->add(*row);
    eventsSigSubsample++;
  }
  
  outTextFile << "size of subsample with just signal = " << subSample->numEntries() << endl;

  //next select nFeedData feed MC events
  for (Int_t j=0; j<(nFeedData); j++) {
    Int_t I=RooRandom::randomGenerator()->Integer(MC_Feed->numEntries());
    const RooArgSet* row=MC_Feed->get(I);
    if(!row) {
      cout<<"skipping null row entry I="<<I<<endl;
      continue;
    }
    subSample->add(*row);
    eventsFeedSubsample++;
  }
  
  outTextFile << "size of subSample after feed = " << subSample->numEntries() << "  feed size = " << eventsFeedSubsample << endl;


  //next get fraction of non-peaking B events
  for (Int_t j=0; j<(nBBData); j++) {
    Int_t I=RooRandom::randomGenerator()->Integer(MC_BkgBB->numEntries());
    const RooArgSet* row=MC_BkgBB->get(I);
    if(!row) {
      cout<<"skipping null row entry I="<<I<<endl;
      continue;
    }
    subSample->add(*row);
    eventsBBSubsample++;
  }
  
  outTextFile << "size of subsample after BkgBB = " << subSample->numEntries() << "  non-peaking BB size = " << eventsBBSubsample << endl;

  //next get fraction of prompt events
  for (Int_t j=0; j<(nPromptData); j++) {
    Int_t I=RooRandom::randomGenerator()->Integer(MC_BkgPro->numEntries());
    const RooArgSet* row=MC_BkgPro->get(I);
    if(!row) {
      cout<<"skipping null row entry I="<<I<<endl;
      continue;
    }
    subSample->add(*row);
    eventsPromptSubsample++;
  }
  
  outTextFile << "size of subsample after BkgPro = " << subSample->numEntries() << "  size of prompt = " << eventsPromptSubsample << endl;

  cout << "After making subsample MC dataset representing yields from data" << endl;




  //next fit the subsample of MC events representing 40 pb-1 
  
  //set yields close to expected values to obtain a good fit
  nsig->setVal( MC_Sig->numEntries()/1.5 );
  nfeed->setVal( MC_Feed->numEntries()/1.5 );
  nbkgBB->setVal( MC_BkgBB->numEntries()/1.5 );
  nbkgPro->setVal( MC_BkgPro->numEntries()*3 );

  outTextFile << " size of subsample = " << subSample->numEntries() << endl;

  RooDataSet* Data_HighMassSB = subSample->reduce("bMass>5.4");

  // first fit high B mass sideband to get ctau values for the BB background

  nbkgBB->setConstant(kFALSE);
  nbkgPro->setConstant(kFALSE);
  fracBBDecay->setConstant(kFALSE);
  ctauBB->setConstant(kTRUE);       // would like to see if I can float these here
  ctauBB2->setConstant(kTRUE);

  meanResCore->setConstant(kFALSE);
  meanResTail->setConstant(kFALSE);
  sigmaResCore->setConstant(kFALSE);
  sigmaResTail->setConstant(kFALSE);
  fracResCore->setConstant(kFALSE);

  ctauB->setConstant(kTRUE);

  model2DSB.fitTo(*Data_HighMassSB,Hesse(true),Minos(true),Timer(true),Save(true));
  
  TCanvas c3a("c3a","c3a",10,10,700,400);
  c3a.Divide(2,1);
  c3a.cd(1);
  RooPlot* frame1=BMass->frame();
  Data_HighMassSB->plotOn(frame1, Binning(32));
  model2DSB.plotOn(frame1,LineColor(3),Range("HighMassSB"), Components(RooArgSet(bkgPro2De)),LineStyle(7),LineWidth(2.0));
  model2DSB.plotOn(frame1,LineColor(4),Range("HighMassSB"), Components(RooArgSet(bkgPro2De,bkgBB2De)),LineWidth(2.0));
  frame1->Draw();
  model2DSB.paramOn(frame1, Layout(0.15,0.95,0.60), Format("NELU",AutoPrecision(3)) ) ;
  frame1->getAttText()->SetTextSize(0.04); 
  frame1->Draw();

  c3a.cd(2).SetLogy();
  RooPlot* frame1=Bct->frame();
  Data_HighMassSB->plotOn(frame1, CutRange("HighMassSB"));
  model2DSB.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7),LineWidth(2.0));
  model2DSB.plotOn(frame1,LineColor(4));
  frame1->Draw();

  frame1->SetMinimum(0.01);
  frame1->Draw();

  c3a->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsMCMassSB.eps");
  c3a->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsMCMassSB.pdf");
  c3a->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsMCMassSB.png");


  // next float only signal yeilds in full sample

  nbkgBB->setConstant(kFALSE);
  nbkgPro->setConstant(kFALSE);
  nfeed->setConstant(kFALSE);
  nsig->setConstant(kFALSE);

  meanResCore->setConstant(kTRUE);
  meanResTail->setConstant(kTRUE);
  sigmaResCore->setConstant(kTRUE);
  sigmaResTail->setConstant(kTRUE);
  fracResCore->setConstant(kTRUE);

  ctauB->setConstant(kTRUE);
  ctauBB.setConstant(kTRUE);
  ctauBB2.setConstant(kTRUE);
  fracBBDecay->setConstant(kTRUE);
  

  if (usePerEventErrors) {
    model2D.fitTo(*subSample,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true),RooFit::ConditionalObservables(RooArgSet(*BctE)));
  } else {
    model2D.fitTo(*subSample,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));
  }
  

  // finally float yields, ctau resolution and signal lifetime together

  nbkgBB->setConstant(kFALSE);
  nbkgPro->setConstant(kFALSE);
  nfeed->setConstant(kFALSE);
  nsig->setConstant(kFALSE);

  meanResCore->setConstant(kFALSE);
  meanResTail->setConstant(kFALSE);
  sigmaResCore->setConstant(kFALSE);
  sigmaResTail->setConstant(kFALSE);
  fracResCore->setConstant(kFALSE);

  ctauB->setConstant(kFALSE);

  //would also like to float B mass shapes for prompt and combinatorial if possible here
  shapeBB.setConstant(kFALSE);
  shapePro.setConstant(kFALSE);


  ctauBB->setConstant(kTRUE);
  ctauBB2->setConstant(kTRUE);
  fracBBDecay->setConstant(kTRUE);

  if (usePerEventErrors) {
    model2D.fitTo(*subSample,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true),RooFit::ConditionalObservables(RooArgSet(*BctE)));
  } else {
    model2D.fitTo(*subSample,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));
  }



  
  TCanvas c4("c4","c4",10,10,800,400);
  c4.Divide(2,1);
  c4.cd(1);
  RooPlot* frame1=BMass->frame();
  subSample->plotOn(frame1,Binning(32));
  model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
  model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
  model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
  model2D.plotOn(frame1,LineColor(4),Components(RooArgSet(bkgBB2De,bkgPro2De,sig2De,feed2De)));
  frame1->Draw();
  model2D.paramOn(frame1, Layout(0.25,0.95,0.45), Format("NELU",AutoPrecision(3)) ) ;
  frame1->getAttText()->SetTextSize(0.04); 
  frame1->Draw();
  
  c4.cd(2).SetLogy();
  RooPlot* frame1=Bct->frame();
  subSample->plotOn(frame1);
  model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
  model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
  model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
  if (usePerEventErrors) {
    model2D.plotOn(frame1,LineColor(4),Components(RooArgSet(bkgBB2De,bkgPro2De,sig2De,feed2De)),ProjWData(*subSample));
  } else {
    model2D.plotOn(frame1,LineColor(4),Components(RooArgSet(bkgBB2De,bkgPro2De,sig2De,feed2De)));
  }
  
  frame1->Draw();
  frame1->SetMinimum(0.1);
  frame1->Draw();

  c4->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitMCResults_subsample.eps");
  c4->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitMCResults_subsample.pdf");
  c4->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitMCResults_subsample.png");



  //save total events from full fit for toys later
  float Ntot=nsig->getVal()+nfeed->getVal()+nbkgBB->getVal()+nbkgPro->getVal();
 
  outTextFile << "====================================================================================================================================================" << endl;
  outTextFile << "Sum of yields from MC subsample with yields from 40 pb-1 in data = " << Ntot;
  outTextFile << " Results summary" << endl;
  outTextFile << " nSig = " << nsig->getVal() << "+-" << nsig->getError() << " embedded " << eventsSigSubsample << endl;
  outTextFile << " nFeed = " << nfeed->getVal() << "+-" << nfeed->getError() << " embedded " << eventsFeedSubsample << endl;
  outTextFile << " nBkgBB = " << nbkgBB->getVal() << "+-" << nbkgBB->getError() << " embedded " << eventsBBSubsample << endl;
  outTextFile << " nBkgPro = " << nbkgPro->getVal() << "+-" << nbkgPro->getError() << " embedded " << eventsPromptSubsample << endl;
  outTextFile << " ctau = " << ctauB->getVal() << "+-" << ctauB->getError() << " embedded " << 0.0460 << endl; 
  outTextFile << "====================================================================================================================================================" << endl;
   
   
  cout  << "====================================================================================================================================================" << endl;
  cout  << "Sum of yields from MC subsample with yields from 40 pb-1 in data = " << Ntot;
  cout << " Results summary" << endl;
  cout  << " nSig = " << nsig->getVal() << "+-" << nsig->getError() << " embedded " << eventsSigSubsample << endl;
  cout << " nFeed = " << nfeed->getVal() << "+-" << nfeed->getError() << " embedded " << eventsFeedSubsample << endl;
  cout  << " nBkgBB = " << nbkgBB->getVal() << "+-" << nbkgBB->getError() << " embedded " << eventsBBSubsample << endl;
  cout  << " nBkgPro = " << nbkgPro->getVal() << "+-" << nbkgPro->getError() << " embedded " << eventsPromptSubsample << endl;
  cout  << " ctau = " << ctauB->getVal() << "+-" << ctauB->getError() << " embedded " << 0.0460 << endl; 
  cout << "====================================================================================================================================================" << endl;
   




  ////////////////////////////////////////////////////////////////////////////////////
  // Now split MC into different pT bins

  int nBins = 5;

  float SigYield[nBins];
  float ErrYield[nBins];
  float SigGen[nBins];

  TCanvas cBinPDFs[nBins];
  TCanvas *cBinProj = new TCanvas("cBinProj","cBinProj",800,800);
  TCanvas *cBinProj2 = new TCanvas("cBinProj2","cBinProj2",800,600);
  cBinProj->Divide(3,3);
  cBinProj2->Divide(3,2);

  for(int i=0;i<nBins;i++){

    float cut1(0);
    float cut2(0);

    if(i==0) {cut1=5;cut2=10;}
    if(i==1) {cut1=10;cut2=13;}
    if(i==2) {cut1=13;cut2=17;}
    if(i==3) {cut1=17;cut2=24;}
    if(i==4) {cut1=24;cut2=40;} 

    //create datasets for each pT cut range for each component

    char ptCut[150];
    sprintf(ptCut,"sqrt((bPx*bPx)+(bPy*bPy))>=%g && sqrt((bPx*bPx)+(bPy*bPy))<%g",cut1,cut2);
  
    outTextFile << "pT cut = " << ptCut << endl;
    RooAbsData* MC_Sig_sub = MC_Sig->reduce(ptCut);
    RooAbsData* MC_BB_sub = MC_BkgBB->reduce(ptCut);
    RooAbsData* MC_Feed_sub = MC_Feed->reduce(ptCut);
    RooAbsData* MC_Pro_sub = MC_BkgPro->reduce(ptCut);
    RooAbsData* MC_HighMassSB_subsample_sub = Data_HighMassSB->reduce(ptCut);

    outTextFile << "size of each component from full sample after pT cut: nSig = " << MC_Sig_sub->numEntries() <<
            " nFeed = " << MC_Feed_sub->numEntries() << 
	    " nBB = " << MC_BB_sub->numEntries() << 
            " nPro = " << MC_Pro_sub->numEntries() << endl;


    outTextFile << "size of reduced pt cut full sample = " << subSample->reduce(ptCut)->numEntries() << endl;

    RooDataSet* MC_subsample_sub(0);
    MC_subsample_sub = new RooDataSet("MC_subsample_sub", "MC_subsample_sub", subSampleSet);

    eventsSigSubsample = 0; eventsFeedSubsample = 0; eventsBBSubsample = 0; eventsPromptSubsample = 0;

    //first get fraction of signal events
    // For the pt bins, I'm going to just generate events in proportion to the expectations from MC
    // instead of the actual yields from data
    float sub_lumi = 40.0;
    float BB_lumi = 64.3;
    float prompt_lumi = 16.5;

    sub_frac = sub_lumi/BB_lumi;
    for (Int_t j=0; j<(MC_Sig_sub->numEntries()*sub_frac); j++) {
      Int_t I=RooRandom::randomGenerator()->Integer(MC_Sig_sub->numEntries());
      const RooArgSet* row=MC_Sig_sub->get(I);
      if(!row) continue;
      MC_subsample_sub->add(*row);
      eventsSigSubsample++;
    }

    //next get fraction of peaking B events
    for (Int_t j=0; j<(MC_Feed_sub->numEntries()*sub_frac); j++) {
      Int_t I=RooRandom::randomGenerator()->Integer(MC_Feed_sub->numEntries());
      const RooArgSet* row=MC_Feed_sub->get(I);
      if(!row) continue;
      MC_subsample_sub->add(*row);
      eventsFeedSubsample++;
    }

    //next get fraction of non-peaking B events
    for (Int_t j=0; j<(MC_BB_sub->numEntries()*sub_frac); j++) {
      Int_t I=RooRandom::randomGenerator()->Integer(MC_BB_sub->numEntries());
      const RooArgSet* row=MC_BB_sub->get(I);
      if(!row) continue;
      MC_subsample_sub->add(*row);
      eventsBBSubsample++;
    }
  
    //next get fraction of prompt events
    sub_frac = sub_lumi/prompt_lumi;
    for (Int_t j=0; j<(MC_Pro_sub->numEntries()*sub_frac); j++) {
      Int_t I=RooRandom::randomGenerator()->Integer(MC_Pro_sub->numEntries());
      const RooArgSet* row=MC_Pro_sub->get(I);
      if(!row) continue;
      MC_subsample_sub->add(*row);
      eventsPromptSubsample++;
    }

    outTextFile << "size of reduce pt cut from each component summed = " << MC_subsample_sub->numEntries() << endl;

    //  RooDataSet* Form("%s%i","MC_subsample_sub_",i) = MC_subsample_sub->clone();
    MC_subsample_sub->Write(Form("%s%i","MC_subsample_sub_",i));
    //  sprintf(file,Form("%s%i%s","/JPsiKsFitPDFs_bin", i, ".pdf") );

    // set some parameters constand and get others from fits to the PDFs
         
    //signal and peaking BB B mass params will be taken directly from these MC fits
    //signal B mass params
    fracBCore.setConstant(kFALSE);
    fracBOut.setConstant(kFALSE);
    meanBCore.setConstant(kFALSE);
    meanBTail.setConstant(kFALSE);
    meanBOut.setConstant(kFALSE);
    widthBCore.setConstant(kFALSE);
    widthBTail.setConstant(kFALSE);
    widthBOut.setConstant(kFALSE);
  
    // peaking BB B mass params
    // There's not enough MC bin by bin, so use the value from the single fit
    fracBOutFeed.setConstant(kFALSE);  
    fracBCoreFeed.setConstant(kFALSE);
    meanBCoreFeed.setConstant(kFALSE);
    widthBCoreFeed.setConstant(kFALSE);
    meanBTailFeed.setConstant(kFALSE);
    widthBTailFeed.setConstant(kFALSE);
    meanBOutFeed.setConstant(kFALSE);
    widthBOutFeed.setConstant(kFALSE);
    
    // non-peaking BB and prompt B mass params will float in final fit, but get
    // MC values here as a good starting point
    // non-peaking BB B mass params
    shapeBB.setConstant(kFALSE);
    // prompt B mass params
    shapePro.setConstant(kFALSE);
    
    // lifetime resolution will come from data, but fit MC signal here to
    // get good starting values
    //signal lifetime parameters
    meanResCore.setConstant(kFALSE);
    sigmaResCore.setConstant(kFALSE);
    meanResTail.setConstant(kFALSE);
    sigmaResTail.setConstant(kFALSE);
    fracResCore.setConstant(kFALSE);
    
    //keep signal lifetime value from fit to full range
    ctauB.setConstant(kTRUE);


  
    //peaking BB lifetime params
    //take this from MC. Should check if different in different bins or not.
    ctauFeed.setConstant(kTRUE);
  
    //non-peaking BB lifetime params
    //will get from data. Should check if it's important to do it bin by bin or not.
    // for now assume it's ok and use the values from the single fit above in data.
    ctauBB.setConstant(kTRUE);
    ctauBB2.setConstant(kTRUE);
    fracBBDecay.setConstant(kTRUE);


    cBinPDFs[i].Divide(2,4);

    //first fit for prompt background, so the lifetime resolution can be fixed for the other fits
  
  
    RooFitResult* fitres= bPro.fitTo(*MC_Pro_sub,Hesse(true),Minos(true));
    cBinPDFs[i].cd(7);
    RooPlot* frame1=BMass->frame();
    MC_Pro_sub->plotOn(frame1,Binning(32));
    bPro.plotOn(frame1,LineColor(4));
    frame1->Draw();
    bPro.paramOn(frame1,Layout(0.7,0.95,0.95),Format("NELU",AutoPrecision(3)));
    frame1->Draw();
  
    RooFitResult* fitres= decayPro.fitTo(*MC_Pro_sub,Hesse(true),Minos(true));
    cBinPDFs[i].cd(8).SetLogy();
    RooPlot* frame1=Bct->frame();
    MC_Pro_sub->plotOn(frame1,Binning(225));
    decayPro.plotOn(frame1,LineColor(4));
    frame1->Draw();
    decayPro.paramOn(frame1,Layout(0.5,0.95,0.95),Format("NELU",AutoPrecision(3)));
    frame1->GetXaxis()->SetRange(0,40);
    frame1->SetMinimum(0.1);
    frame1->Draw();
  
    meanResCore->setConstant(kTRUE);
    meanResTail->setConstant(kTRUE);
    sigmaResCore->setConstant(kTRUE);
    sigmaResTail->setConstant(kTRUE);
    fracResCore->setConstant(kTRUE);
  
    cBinPDFs[i].cd(1);
    RooFitResult* fitres= bSig.fitTo(*MC_Sig_sub,Hesse(true),Minos(true));
    RooPlot* frame1=BMass->frame();
    MC_Sig_sub->plotOn(frame1,Binning(80));
    bSig.plotOn(frame1,LineColor(4));
    frame1->Draw();
    bSig.paramOn(frame1,Layout(0.65,0.95,0.95),Format("NELU",AutoPrecision(3)));
    frame1->GetXaxis()->SetRange(35,65);
    frame1->Draw();
  
    cBinPDFs[i].cd(2).SetLogy();
    if(usePerEventErrors) {
      RooFitResult* fitres= decaySig.fitTo(*MC_Sig_sub,Hesse(true),Minos(true),ConditionalObservables(*BctE));
    } else {
      RooFitResult* fitres= decaySig.fitTo(*MC_Sig_sub,Hesse(true),Minos(true));
    }
    
    RooPlot* frame1=Bct->frame();
    MC_Sig_sub->plotOn(frame1,Binning(45));
    if (usePerEventErrors) {
      decaySig.plotOn(frame1,ProjWData(*MC_Sig),LineColor(4));
    } else {
      decaySig.plotOn(frame1,LineColor(4));
    }
    frame1->SetMinimum(0.1);
    frame1->Draw();
    decaySig.paramOn(frame1,Layout(0.65,0.95,0.95),Format("NELU",AutoPrecision(3)));
    frame1->Draw();
    
    cBinPDFs[i].cd(3);
    // Don't fit for the peaking BB PDF shape in each bin because too few 
    // to get a good fit.
    //RooFitResult* fitres= bFeed.fitTo(*MC_Feed_sub,Hesse(true),Minos(true));
    RooPlot* frame1=BMass->frame();
    MC_Feed_sub->plotOn(frame1,Binning(32));
    bFeed.plotOn(frame1,LineColor(4));
    frame1->Draw();
    bFeed.paramOn(frame1,Layout(0.7,0.90,0.90),Format("NELU",AutoPrecision(3)));
    TPaveText* pbox = (TPaveText*) frame1->findObject("bFeed_paramBox");
    if (pbox) {
      pbox->AddText("Fixed parameters");
    }
    frame1->Draw();
  
    cBinPDFs[i].cd(4).SetLogy();
    RooPlot* frame1=Bct->frame();
    MC_Feed_sub->plotOn(frame1,Binning(90));
    decayFeed.plotOn(frame1,LineColor(4));
    frame1->Draw();
    decayFeed.paramOn(frame1,Layout(0.65,0.90,0.90),Format("NELU",AutoPrecision(3)));
    TPaveText* pbox = (TPaveText*) frame1->findObject("decayFeed_paramBox");
    if (pbox) {
      pbox->AddText("Fixed parameters");
    }
    frame1->SetMinimum(0.1);
    frame1->Draw();
  
    RooFitResult* fitres= bBB.fitTo(*MC_BB_sub,Hesse(true),Minos(true));
    cBinPDFs[i].cd(5);
    RooPlot* frame1=BMass->frame();
    MC_BB_sub->plotOn(frame1,Binning(32));
    bBB.plotOn(frame1,LineColor(4));
    frame1->Draw();
    bBB.paramOn(frame1,Layout(0.7,0.95,0.95),Format("NELU",AutoPrecision(3)));
    frame1->Draw();
  
    RooFitResult* fitres= decayBB.fitTo(*MC_BB_sub,Hesse(true),Minos(true));
    cBinPDFs[i].cd(6).SetLogy();
    RooPlot* frame1=Bct->frame();
    MC_BB_sub->plotOn(frame1,Binning(90));
    decayBB.plotOn(frame1,LineColor(4));
    frame1->Draw();
    decayBB.paramOn(frame1,Layout(0.5,0.95,0.95),Format("NELU",AutoPrecision(3)));
    frame1->SetMinimum(0.1);
    frame1->Draw();
  
  
    // after PDF fits now fix signal B mass shapes for signal and peaking BB bkg
    fracBCore.setConstant(kTRUE);
    fracBOut.setConstant(kTRUE);
    meanBCore.setConstant(kTRUE);
    meanBTail.setConstant(kTRUE);
    meanBOut.setConstant(kTRUE);
    widthBCore.setConstant(kTRUE);
    widthBTail.setConstant(kTRUE);
    widthBOut.setConstant(kTRUE);
  
    fracBOutFeed.setConstant(kTRUE);  
    fracBCoreFeed.setConstant(kTRUE);
    meanBCoreFeed.setConstant(kTRUE);
    widthBCoreFeed.setConstant(kTRUE);
    meanBTailFeed.setConstant(kTRUE);
    widthBTailFeed.setConstant(kTRUE);
    meanBOutFeed.setConstant(kTRUE);
    widthBOutFeed.setConstant(kTRUE);
  
    // fix signal B lifetime
    ctauB->setConstant(kTRUE);
  
    // fix non-peaking and prompt B mass shapes, for now
    shapeBB.setConstant(kTRUE);
    shapePro.setConstant(kTRUE);
    
  
    // Do fit as above for full sample, except without floating the b0 lifetime
    // first fit high B mass sideband to get ctau values for the BB background
  
    nbkgBB->setConstant(kFALSE);
    nbkgPro->setConstant(kFALSE);
    ctauBB.setConstant(kTRUE);
    ctauBB2.setConstant(kTRUE);
    fracBBDecay->setConstant(kFALSE);
  
    meanResCore->setConstant(kFALSE);
    meanResTail->setConstant(kFALSE);
    sigmaResCore->setConstant(kFALSE);
    sigmaResTail->setConstant(kFALSE);
    fracResCore->setConstant(kFALSE);

    cout << "starting fit to MC for bin " << i << endl;
    
    // first fit to sideband
    model2DSB.fitTo(*MC_HighMassSB_subsample_sub,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));

    // next float only signal yeilds
  
    nbkgBB->setConstant(kFALSE);
    nbkgPro->setConstant(kFALSE);
    nfeed->setConstant(kFALSE);
    nsig->setConstant(kFALSE);
  
    meanResCore->setConstant(kTRUE);
    meanResTail->setConstant(kTRUE);
    sigmaResCore->setConstant(kTRUE);
    sigmaResTail->setConstant(kTRUE);
    fracResCore->setConstant(kTRUE);
  
    ctauB->setConstant(kTRUE);
    ctauBB->setConstant(kTRUE);
    ctauBB2->setConstant(kTRUE);
    fracBBDecay->setConstant(kTRUE);
  
    if (usePerEventErrors) {
      RooFitResult* fitres= model2D.fitTo(*MC_subsample_sub,"mhe",RooFit::ConditionalObservables(RooArgSet(*BctE)));
    } else {
      RooFitResult* fitres= model2D.fitTo(*MC_subsample_sub,"mhe");
    }


   // finally try floating yields, resolution and B mass shapes for prompt and non-peaking B
  
    nbkgBB->setConstant(kFALSE);
    nbkgPro->setConstant(kFALSE);
    nfeed->setConstant(kFALSE);
    nsig->setConstant(kFALSE);
  
    meanResCore->setConstant(kFALSE);
    meanResTail->setConstant(kFALSE);
    sigmaResCore->setConstant(kFALSE);
    sigmaResTail->setConstant(kFALSE);
    fracResCore->setConstant(kFALSE);
  
    shapeBB.setConstant(kFALSE);
    shapePro.setConstant(kFALSE);
   
    ctauB->setConstant(kTRUE);
    ctauBB->setConstant(kTRUE);
    ctauBB2->setConstant(kTRUE);
    fracBBDecay->setConstant(kTRUE);
  

    cout << " before fit to MC for bin " << i << endl;
    if (usePerEventErrors) {
      RooFitResult* fitres= model2D.fitTo(*MC_subsample_sub,"mhe",RooFit::ConditionalObservables(RooArgSet(*BctE)));
    } else {
      RooFitResult* fitres= model2D.fitTo(*MC_subsample_sub,"mhe");
    }
    cout << " after fit to MC for bin " << i << endl;

    if(i<3) cBinProj->cd(i*3+1);
    else cBinProj2->cd((i-3)*3+1);
    RooPlot* frame1=BMass->frame();
    MC_subsample_sub->plotOn(frame1,Binning(32));
    model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(4));
    frame1->Draw();
    model2D.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
    frame1->Draw();
  
    if(i<3) cBinProj->cd(i*3+2);
    else cBinProj2->cd((i-3)*3+2);
    RooPlot* frame1=BMass->frame();
    MC_subsample_sub->plotOn(frame1,Binning(32),CutRange("ct0p1"));
    model2D.plotOn(frame1,LineColor(3),ProjectionRange("ct0p1"),Components(RooArgSet(bkgPro2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(2),ProjectionRange("ct0p1"),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(6),ProjectionRange("ct0p1"),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(4),ProjectionRange("ct0p1"));
    frame1->Draw();
    model2D.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
    frame1->Draw();
  
    if(i<3) cBinProj->cd(i*3+3).SetLogy();
    else cBinProj2->cd((i-3)*3+3).SetLogy();
    RooPlot* frame1=Bct->frame();
    MC_subsample_sub->plotOn(frame1);
    if (usePerEventErrors) {
      model2D.plotOn(frame1,LineColor(4),ProjWData(*MC_subsample_sub));  
    } else {
      model2D.plotOn(frame1,LineColor(4)); 
    }
    model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
    frame1->Draw();
    frame1->SetMinimum(0.1);
    model2D.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
    frame1->Draw();
  
    SigYield[i]=nsig->getVal();
    ErrYield[i]=nsig->getError();
    SigGen[i]=eventsSigSubsample;

    outTextFile << "====================================================" << endl;
    outTextFile << "Summary of fit results from bin" << i << endl;
    outTextFile << " nSig = " << nsig->getVal() << "+-" << nsig->getError() << " embedded " << eventsSigSubsample << endl;
    outTextFile << " nFeed = " << nfeed->getVal() << "+-" << nfeed->getError() << " embedded " << eventsFeedSubsample << endl;
    outTextFile << " nBkgBB = " << nbkgBB->getVal() << "+-" << nbkgBB->getError() << " embedded " << eventsBBSubsample << endl;
    outTextFile << " nBkgPro = " << nbkgPro->getVal() << "+-" << nbkgPro->getError() << " embedded " << eventsPromptSubsample << endl;
    outTextFile << cut1 << "-" << cut2 << " & $" << Form("%.1f",nsig->getVal()) << "\\pm" << Form ("%.1f",nsig->getError()) << "$ & $" << eventsSigSubsample << 
      "$ & $" << Form("%.1f",nfeed->getVal()) << "\\pm" << Form ("%.1f",nfeed->getError()) << "$ & $" << eventsFeedSubsample <<
      "$ & $" << Form("%.1f",nbkgBB->getVal()) << "\\pm" << Form ("%.1f",nbkgBB->getError()) << "$ & $" << eventsBBSubsample <<
      "$ & $" << Form("%.1f",nbkgPro->getVal()) << "\\pm" << Form ("%.1f",nbkgPro->getError()) << "$ & $" << eventsPromptSubsample << "$ \\\\" << endl;
    outTextFile <<  "====================================================" << endl;

    char file[150];
    sprintf(file,Form("%s%i%s","/JPsiKsFitPDFs_MC_ptbin", i, ".pdf") );
    cout << "Saving pdfs for MC bin with file name " << file << endl;
    cBinPDFs[i]->SaveAs("JPsiKsFit_dataRes/"+trigger+file);

  } //end loop over bins

  cBinProj->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResults_MC_ptbins1.pdf");
  cBinProj2->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResults_MC_ptbins2.pdf");

  for(int i=0;i<nBins;i++){
    outTextFile << "Pt bin " << i << " yield: " << SigYield[i]  
    << " +/- " << ErrYield[i] << " Gen: " << SigGen[i] << "\n";
  }

  cout << "past the binned fit for MC" << endl;

}



if (fitData) {

  //   Now do the fit to data

  cout << "starting the data fit" << endl;
  TFile seleBB("seleJPsiKs_data_Nov4ReReco.root");
  TTree *treeData=tofit;

  RooDataSet Data_All_noHLT("Data_All_noHLT", "Data_All_noHLT", treeData, Vars);  

  RooDataSet* Data_All = Data_All_noHLT->reduce(triggerCut+bptCut);

  RooDataSet* Data_HighMassSB = Data_All->reduce("bMass>5.4");

  // first fit high B mass sideband to get ctau values for the BB background

  nbkgBB->setVal(1200);
  nbkgPro->setVal(6700);
  nbkgBB->setConstant(kFALSE);
  nbkgPro->setConstant(kFALSE);
  fracBBDecay->setConstant(kFALSE);
  ctauBB->setConstant(kFALSE);       // would like to see if I can float these here and then fix them for each bin.
  ctauBB2->setConstant(kFALSE);

  meanResCore->setConstant(kTRUE);
  meanResTail->setConstant(kTRUE);
  sigmaResCore->setConstant(kTRUE);
  sigmaResTail->setConstant(kTRUE);
  fracResCore->setConstant(kTRUE);

  ctauB->setConstant(kTRUE);

  model2DSB.fitTo(*Data_HighMassSB,Hesse(true),Minos(true),Timer(true),Save(true));

  TCanvas *SBCanv = new TCanvas("SBCanv","SBCanv",800,800);
  SBCanv->Divide(2,6);


  TCanvas c3a("c3a","c3a",10,10,700,400);
  c3a.Divide(2,1);
  c3a.cd(1);
  RooPlot* frame1=BMass->frame();
  Data_HighMassSB->plotOn(frame1, Binning(32));
  model2DSB.plotOn(frame1,LineColor(3),Range("HighMassSB"), Components(RooArgSet(bkgPro2De)),LineStyle(7));
  model2DSB.plotOn(frame1,LineColor(4),Range("HighMassSB"), Components(RooArgSet(bkgPro2De,bkgBB2De)));
  frame1->Draw();

  c3a.cd(2).SetLogy();
  RooPlot* frame1=Bct->frame();
  Data_HighMassSB->plotOn(frame1, CutRange("HighMassSB"));
  model2DSB.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
  model2DSB.plotOn(frame1,LineColor(4));
  frame1->Draw();

  frame1->SetMinimum(0.01);
  model2DSB.paramOn(frame1, Layout(0.45,0.95,0.90), Format("NELU",AutoPrecision(3)) ) ;
  frame1->getAttText()->SetTextSize(0.04); 

  frame1->Draw();
  
  c3a->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsDataMassSB.eps");
  c3a->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsDataMassSB.pdf");
  c3a->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsDataMassSB.png");

  //print out results from this fit:
outTextFile << endl;
outTextFile << "Results from SB fit to full pt and y range" << endl;
outTextFile << "==========================================" << endl;
outTextFile << "  Combinatorial BB background lifetime 1 & " << ctauBB->getVal() << "$\\pm$" << ctauBB->getError() << endl;
outTextFile << "  Combinatorial BB background lifetime 2 & " << ctauBB2->getVal() << "$\\pm$" << ctauBB2->getError() << endl;
outTextFile << "  Combinatorial BB background lifetime 1 fraction & " << fracBBDecay->getVal() << "$\\pm$" << fracBBDecay->getError() << endl;
outTextFile << "  Combinatorial BB background yield & " << nbkgBB->getVal() << "$\\pm$" << nbkgBB->getError() << endl;
outTextFile << "  Prompt yield & " << nbkgPro->getVal() << "$\\pm$" << nbkgPro->getError() << endl;


  // next float only signal yeilds in full sample

  nbkgBB->setConstant(kFALSE);
  nbkgPro->setConstant(kFALSE);
  nfeed->setConstant(kFALSE);
  nsig->setConstant(kFALSE);

  meanResCore->setConstant(kTRUE);
  meanResTail->setConstant(kTRUE);
  sigmaResCore->setConstant(kTRUE);
  sigmaResTail->setConstant(kTRUE);
  fracResCore->setConstant(kTRUE);

  ctauB->setConstant(kTRUE);
  ctauBB.setConstant(kTRUE);
  ctauBB2.setConstant(kTRUE);
  fracBBDecay->setConstant(kTRUE);




  if (usePerEventErrors) {
    model2D.fitTo(*Data_All,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true),ConditionalObservables(*BctE));
  } else {
    model2D.fitTo(*Data_All,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));
  }



  // finally float yields, ctau resolution and signal lifetime together

  nbkgBB->setConstant(kFALSE);
  nbkgPro->setConstant(kFALSE);
  nfeed->setConstant(kFALSE);
  nsig->setConstant(kFALSE);

  meanResCore->setConstant(kFALSE);
  meanResTail->setConstant(kFALSE);
  sigmaResCore->setConstant(kFALSE);
  sigmaResTail->setConstant(kFALSE);
  fracResCore->setConstant(kFALSE);

  ctauB->setConstant(kFALSE);

  //would also like to float B mass shapes for prompt and combinatorial if possible here
  shapeBB.setConstant(kFALSE);
  shapePro.setConstant(kFALSE);

  if(floatKsMass) {
    meanBCore->setConstant(kFALSE);
    widthBCore->setConstant(kFALSE);
  }

  ctauBB->setConstant(kTRUE);
  ctauBB2->setConstant(kTRUE);
  fracBBDecay->setConstant(kTRUE);

  //try getting NLL
cout << "  Try getting NLL value before fit = : " <<  model2D.createNLL(*Data_All).getVal() << endl;

  if (usePerEventErrors) {
    model2D.fitTo(*Data_All,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true),ConditionalObservables(*BctE));
  } else {
    model2D.fitTo(*Data_All,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));
  }
  
  //try getting NLL
cout << "  Try getting NLL value after fit = : " <<  model2D.createNLL(*Data_All).getVal() << endl;

  
  //make projection plots, including with cut on ct>0.1 for Jim
  Bct->setRange("ct0p1",0.01,0.4);

  
  TLatex* cms = new TLatex();
  cms->SetNDC();
  cms->SetTextSize(0.055);
  cms->SetTextAlign(22);
  cms->SetText(0.750,0.90,"CMS  #sqrt{s} = 7 TeV");
  TLatex* lumi = new TLatex();
  lumi->SetNDC();
  lumi->SetTextSize(0.055);
  lumi->SetTextAlign(22);
  lumi->SetText(0.750,0.83,"L = 40 pb^{-1}");

  TLatex* alab = new TLatex();
  alab->SetNDC();
  alab->SetTextSize(0.055);
  alab->SetTextAlign(22);
  alab->SetText(0.23,0.20,"(a)");
  
  TLatex* blab = new TLatex();
  blab->SetNDC();
  blab->SetTextSize(0.055);
  blab->SetTextAlign(22);
  blab->SetText(0.23,0.20,"(b)");
 
  TCanvas cPub("cPub","cPub",10,10,900,400);
  cPub.Divide(2,1);
  TCanvas cPubA("cPubA","cPubA",10,10,700,500);
  TCanvas cPubB("cPubB","cPubB",10,10,700,500);
  TCanvas c3("c3","c3",10,10,900,400);
  c3.Divide(3,1);
  c3.cd(1);
  RooPlot* frame1=BMass->frame();
  Data_All->plotOn(frame1, Binning(32));
  model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(3),LineWidth(2.5));
  model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,feed2De)),LineStyle(7),LineWidth(2.0));
  model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(5),LineWidth(2.0));
  model2D.plotOn(frame1,LineColor(4),Components(RooArgSet(bkgBB2De,bkgPro2De,sig2De,feed2De)),LineWidth(2.0));
  frame1->Draw();
  cms->Draw();
  lumi->Draw();
  
  cPub.cd(1);
  frame1->GetYaxis()->SetTitleOffset(1.3);
  frame1->Draw();
  cms->Draw();
  lumi->Draw();

  cPubA.cd();
  frame1->GetYaxis()->SetTitleOffset(1.3);
  frame1->GetYaxis()->SetTitle("Candidates / ( 0.025 GeV )");
  frame1->Draw();

  TH1F a;
  a.SetLineColor(3);
  a.SetLineStyle(3);
  TH1F b;
  b.SetLineColor(2);
  b.SetLineStyle(7);
  TH1F c;
  c.SetLineColor(6);
  c.SetLineStyle(5);
  TH1F d;
  d.SetLineColor(4);
  d.SetLineStyle(1);

  TLegend* legA = new TLegend(0.62,0.20,0.96,0.45);
  legA->AddEntry(Data_All, "CMS data", "lp");
  legA->AddEntry(&a,"Prompt J/#Psi","l");
  legA->AddEntry(&b,"  + peaking B","l");
  legA->AddEntry(&c,"  + non-peaking B","l");
  legA->AddEntry(&d,"  + signal","l");
  legA->SetFillColor(0);
  legA->SetBorderSize(0);
  legA->Draw();
  alab->Draw();
  cms->Draw();
  lumi->Draw();

  c3.cd(2);
  RooPlot* frame1=BMass->frame();
  Data_All->plotOn(frame1, Binning(32),CutRange("ct0p1"));
  model2D.plotOn(frame1,LineColor(3),ProjectionRange("ct0p1"),Components(RooArgSet(bkgPro2De)),LineStyle(3),LineWidth(2.5));
  model2D.plotOn(frame1,LineColor(2),ProjectionRange("ct0p1"),Components(RooArgSet(bkgPro2De,feed2De)),LineStyle(7),LineWidth(2.0));
  model2D.plotOn(frame1,LineColor(6),ProjectionRange("ct0p1"),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(5),LineWidth(2.0));
  model2D.plotOn(frame1,LineColor(4),ProjectionRange("ct0p1"),Components(RooArgSet(bkgBB2De,bkgPro2De,sig2De,feed2De)),LineWidth(2.0));
  frame1->Draw();
  cms->Draw();
  lumi->Draw();

  c3.cd(3).SetLogy();
  RooPlot* frame1=Bct->frame();
  Data_All->plotOn(frame1);
  model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(3),LineWidth(2.5));
  model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,feed2De)),LineStyle(7),LineWidth(2.0));
  model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(5),LineWidth(2.0));
  if (usePerEventErrors) {
    model2D.plotOn(frame1,LineColor(4),ProjWData(*Data_All),LineWidth(2.0));
  } else {
    model2D.plotOn(frame1,LineColor(4),LineWidth(2.0));
  }

  frame1->Draw();
  frame1->SetMinimum(0.1);
  frame1->SetMaximum(9000);
  frame1->Draw();
  cms->Draw();
  lumi->Draw();

  cPub.cd(2).SetLogy(); 
  frame1->Draw();
  cms->Draw();
  lumi->Draw();

  cPubB.cd().SetLogy(); 
  frame1->GetYaxis()->SetTitle("Candidates / ( 0.0045 cm )");
  frame1->Draw();
  cms->Draw();
  lumi->Draw();
  
  TLegend* legB = new TLegend(0.62,0.47,0.95,0.72);
  legB->AddEntry(Data_All, "CMS data", "lp");
  legB->AddEntry(&a,"Prompt J/#Psi","l");
  legB->AddEntry(&b,"  + peaking B","l");
  legB->AddEntry(&c,"  + non-peaking B","l");
  legB->AddEntry(&d,"  + signal","l");
  legB->SetFillColor(0);
  legB->SetBorderSize(0);
  legB->Draw();
  blab->Draw();
  
  c3->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsData2.eps");
  c3->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsData2.pdf");
  c3->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsData2.png");
  
  cPub->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsDataPub.eps");
  cPub->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsDataPub.pdf");
  
  cPubA->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsDataPubA.eps");
  cPubA->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsDataPubA.pdf");
  
  cPubB->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsDataPubB.eps");
  cPubB->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsDataPubB.pdf");
  
  
// now draw ct with the stat box for the AN
  model2D.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
  frame1->getAttText()->SetTextSize(0.02); 
  frame1->Draw();
  cms->Draw();
  lumi->Draw();
  
  c3->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsData.eps");
  c3->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsData.pdf");
  c3->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResultsData.png");


  outTextFile << "====================================================================================================================================================" << endl;
  outTextFile << "Sum of signal yields from fit to full data sample = " <<
 	   nsig->getVal()+nfeed->getVal()+nbkgBB->getVal()+nbkgPro->getVal() << endl;
  outTextFile << " Results summary" << endl;
  outTextFile << " nSig = " << nsig->getVal() << "+-" << nsig->getError() << endl;
  outTextFile << " nFeed = " << nfeed->getVal() << "+-" << nfeed->getError() << endl;
  outTextFile << " nBkgBB = " << nbkgBB->getVal() << "+-" << nbkgBB->getError()  << endl;
  outTextFile << " nBkgPro = " << nbkgPro->getVal() << "+-" << nbkgPro->getError() << endl;
  outTextFile << " ctau = " << ctauB->getVal() << "+-" << ctauB->getError() << " PDG = " << 0.04572 << endl; 
  outTextFile << "====================================================================================================================================================" << endl;
   
  outTextFile << endl;
  outTextFile << "Results from fit to full pt and y range" << endl;
  outTextFile << "==========================================" << endl;
  outTextFile << " Signal yield & " << nsig->getVal() << "$\\pm$" << nsig->getError() << endl;
  outTextFile << " Peaking BB yield & " << nfeed->getVal() << "$\\pm$" << nfeed->getError() << endl;
  outTextFile << " Combinatorial BB yield & " << nbkgBB->getVal() << "$\\pm$" << nbkgBB->getError()  << endl;
  outTextFile << " Prompt yield & " << nbkgPro->getVal() << "$\\pm$" << nbkgPro->getError() << endl;
  outTextFile << " Lifetime & " << ctauB->getVal() << "$\\pm$" << ctauB->getError() << endl; 
  outTextFile << " ct core mean & " << meanResCore->getVal() << "$\\pm$" << meanResCore->getError() << endl;
  outTextFile << " ct tail mean & " << meanResTail->getVal() << "$\\pm$" << meanResTail->getError() << endl;
  outTextFile << " ct core sigma & " << sigmaResCore->getVal() << "$\\pm$" << sigmaResCore->getError() << endl;
  outTextFile << " ct tail sigma & " <<  sigmaResTail->getVal() << "$\\pm$" << sigmaResTail->getError() << endl;
  outTextFile << " ct core fraction & " <<  fracResCore->getVal() << "$\\pm$" << fracResCore->getError() << endl;
  outTextFile << " Combinatorial BB shape & " <<  shapeBB->getVal() << "$\\pm$" << shapeBB->getError() << endl;
  outTextFile << " Prompt shape & " << shapePro->getVal() << "$\\pm$" << shapePro->getError() << endl;
  outTextFile << "====================================================================================================================================================" << endl;



  //next let's do some toy experiments
  if (doToys) {

    float Ntot=nsig->getVal()+nfeed->getVal()+nbkgBB->getVal()+nbkgPro->getVal();

    outTextFile << "================================================" << endl;
    outTextFile << "Starting pure toys" << endl;
    outTextFile << "===================================================" << endl;

    outTextFile << "Now Ntot = " << Ntot << "  nSig = " << nsig->getVal() << " nbkgBB peaking= " << nfeed->getVal() << "  nbkgBB nonpeaking = " <<
    nbkgBB->getVal() << "  nPrompt = " << nbkgPro->getVal() << endl;

// float all the same params
//    ctauB.setConstant(kFALSE);
//    fracResCore->setConstant(kTRUE);
//    meanResCore->setConstant(kTRUE);
//    meanResTail->setConstant(kTRUE);
//    sigmaResCore->setConstant(kTRUE);
//    sigmaResTail->setConstant(kTRUE);
//    shapeBB->setConstant(kTRUE);
//    shapePro->setConstant(kTRUE);
    

    if (usePerEventErrors) {
      RooMCStudy mgr(model2D,RooArgList(*BMass,*Bct,*BctE),FitModel(model2D),Extended(1),ConditionalObservables(RooArgSet(*BctE)),Minos(1),Hesse(1),SplitRange(false));
      mgr.generateAndFit(nToys ,Ntot,kTRUE);
      mgr.Write();
      RooDataSet* toySample = mgr.genData(0);

      model2D.fitTo(*toySample,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true),ConditionalObservables(RooArgSet(*BctE)));
    } else {
      RooMCStudy mgr(model2D,RooArgList(*BMass,*Bct),FitModel(model2D),Extended(1),Minos(1),Hesse(1),SplitRange(false));
      mgr.generateAndFit(nToys ,Ntot,kTRUE);
      mgr.Write();
      RooDataSet* toySample = mgr.genData(0);

      model2D.fitTo(*toySample,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));
      
      cout << " Single toy nll = " << model2D.createNLL(*toySample).getVal() << endl;
      
    }
    
    TCanvas *cToy1 = new TCanvas("cToy1","cToy1",800,500);
    cToy1->Divide(2,1);
    cToy1->cd(1);
    RooPlot* frame1=BMass->frame();
    toySample->plotOn(frame1,Binning(32));
    model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(4),Components(RooArgSet(bkgBB2De,bkgPro2De,sig2De,feed2De)));
    frame1->Draw();
 
    cToy1->cd(2).SetLogy();
    RooPlot* frame1=Bct->frame();
    toySample->plotOn(frame1);
    model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(4),Components(RooArgSet(bkgBB2De,bkgPro2De,sig2De,feed2De)));
    frame1->Draw();
    model2D.paramOn(frame1,toySample,"Fit parameters",1,"NELU");
    frame1->Draw();
    TPaveText* pbox = (TPaveText*) frame1->findObject("model2D_paramBox");

    cToy1->SaveAs("JPsiKsFit_dataRes/"+trigger+"/pToys_sample.eps");

    TCanvas *cToy2 = new TCanvas("cToy2","cToy2",800,800);
    cToy2->Divide(3,3);
    TCanvas *cToy3 = new TCanvas("cToy3","cToy3",800,800);
    cToy3->Divide(3,2);

    cToy2->cd(1);
    RooPlot* ctauMeanFrame = ctauB->frame(0.0400,0.0550,75);
    mgr.plotParamOn(ctauMeanFrame);
    ctauMeanFrame->Draw();

    cToy2->cd(2);
    RooPlot* ctauErrorFrame = mgr.plotError(*ctauB,0.0012,0.0018,60);
    ctauErrorFrame->Draw();

    cToy2->cd(3);
    RooPlot* ctauPullFrame = mgr.plotPull(*ctauB,-3.0,3.0,60,kTRUE);
    ctauPullFrame->Draw();
    cToy3->cd(1);
    ctauPullFrame->Draw();

    cToy2->cd(4);
    RooPlot* bmassMeanFrame = nsig->frame(nsig->getVal()-150,nsig->getVal()+150 ,30);
    mgr.plotParamOn(bmassMeanFrame);
    bmassMeanFrame->Draw();

    cToy2->cd(5);
    RooPlot* bmassErrorFrame = mgr.plotError(*nsig,30,45,60);

    bmassErrorFrame->Draw();
  
    cToy2->cd(6);
    RooPlot* bmassPullFrame = mgr.plotPull(*nsig,-3.0,3.0,60,kTRUE);
    bmassPullFrame->Draw();
    cToy3->cd(2);
    bmassPullFrame->Draw();
  
    cToy2->cd(7);
    RooPlot* frame = mgr.plotPull(*nfeed,-3.0,3.0,60,kTRUE);
    frame->Draw();
    cToy3->cd(3);
    frame->Draw();
  
    cToy2->cd(8);
    RooPlot* frame = mgr.plotPull(*nbkgBB,-3.0,3.0,60,kTRUE);
    frame->Draw();
    cToy3->cd(4);
    frame->Draw();
  
    cToy2->cd(9);
    RooPlot* frame = mgr.plotPull(*nbkgPro,-3.0,3.0,60,kTRUE);
    frame->Draw();
    cToy3->cd(5);
    frame->Draw();
  
    cToy2->SaveAs("JPsiKsFit_dataRes/"+trigger+"/pToys.eps");
    cToy3->SaveAs("JPsiKsFit_dataRes/"+trigger+"/pToys_pulls.eps");
    cToy3->SaveAs("JPsiKsFit_dataRes/"+trigger+"/pToys_pulls.pdf");
    
//    mgr->Write();
//    outFile2->Close();
  
  }


  if (doPtBins) {

  ////////////////////////////////////////////////////////////////////////////////////
  // Now split data into different pT bins. Still need MC to fit for PDFs in each bin. Would be better to get PDFs from data here, too.

  int nBins = 6;

  float SigYield[nBins];
  float ErrYield[nBins];

  TCanvas cBinPDFs[nBins];
  TCanvas *cBinProj = new TCanvas("cBinProj","cBinProj",800,800);
  TCanvas *cBinProj2 = new TCanvas("cBinProj2","cBinProj2",800,800);
  TCanvas *cBinProj3 = new TCanvas("cBinProj3","cBinProj3",800,800);
  cBinProj->Divide(3,3);
  cBinProj2->Divide(3,3);
  cBinProj3->Divide(2,3);

  //want to fix the ct resolution to the value from the 3rd bin for the 4th and 5th bins (and 6th)
  float meanResCoreValue;
  float sigmaResCoreValue;
  float meanResTailValue;
  float sigmaResTailValue;
  float fracResCoreValue;

  for(int i=0;i<nBins;i++){

    float cut1(0);
    float cut2(0);

    if(i==0) {cut1=5;cut2=10;}
    if(i==1) {cut1=10;cut2=13;}
    if(i==2) {cut1=13;cut2=17;}
    if(i==3) {cut1=17;cut2=24;}
    if(i==4) {cut1=24;cut2=40;}
    if(i==5) {cut1=24;cut2=999;}

    //create datasets for each pT cut range for each component

    char ptCut[150];
    sprintf(ptCut,"sqrt((bPx*bPx)+(bPy*bPy))>=%g && sqrt((bPx*bPx)+(bPy*bPy))<%g",cut1,cut2);
  
    outTextFile << "pT cut = " << ptCut << endl;
    RooAbsData* MC_Sig_sub = MC_SigPrivate->reduce(ptCut);
    RooAbsData* MC_BB_sub = MC_BkgBB->reduce(ptCut);
    RooAbsData* MC_Feed_sub = MC_Feed->reduce(ptCut);
    RooAbsData* MC_Pro_sub = MC_BkgPro->reduce(ptCut);
    RooAbsData* Data_subsample_sub = Data_All->reduce(ptCut);

    RooAbsData* Data_HighMassSB_subsample_sub = Data_HighMassSB->reduce(ptCut);

    // reset all params that need to float in the PDF fits. Done for each bin separately.

    //signal and peaking BB B mass params will be taken directly from these MC fits
    //signal B mass params
    fracBCore.setConstant(kFALSE);
    fracBOut.setConstant(kFALSE);
    meanBCore.setConstant(kFALSE);
    meanBTail.setConstant(kFALSE);
    meanBOut.setConstant(kFALSE);
    widthBCore.setConstant(kFALSE);
    widthBTail.setConstant(kFALSE);
    widthBOut.setConstant(kFALSE);

    // peaking BB B mass params
    // There's not enough MC bin by bin, so use the value from the single fit
    fracBOutFeed.setConstant(kFALSE);  
    fracBCoreFeed.setConstant(kFALSE);
    meanBCoreFeed.setConstant(kFALSE);
    widthBCoreFeed.setConstant(kFALSE);
    meanBTailFeed.setConstant(kFALSE);
    widthBTailFeed.setConstant(kFALSE);
    meanBOutFeed.setConstant(kFALSE);
    widthBOutFeed.setConstant(kFALSE);
    
    // non-peaking BB and prompt B mass params will float in final fit, but get
    // MC values here as a good starting point
    // non-peaking BB B mass params
    shapeBB.setConstant(kFALSE);
    // prompt B mass params
    shapePro.setConstant(kFALSE);
    
    // lifetime resolution will come from data, but fit MC signal here to
    // get good starting values
    //signal lifetime parameters
    meanResCore.setConstant(kFALSE);
    sigmaResCore.setConstant(kFALSE);
    meanResTail.setConstant(kFALSE);
    sigmaResTail.setConstant(kFALSE);
    fracResCore.setConstant(kFALSE);
    
    //keep signal lifetime value from fit to full range
    ctauB.setConstant(kTRUE);

    //peaking BB lifetime params
    //take this from MC. Should check if different in different bins or not.
    ctauFeed.setConstant(kTRUE);

    //non-peaking BB lifetime params
    //will get from data. Should check if it's important to do it bin by bin or not.
    // for now assume it's ok and use the values from the single fit above in data.
    ctauBB.setConstant(kTRUE);
    ctauBB2.setConstant(kTRUE);
    fracBBDecay.setConstant(kTRUE);


    cBinPDFs[i].Divide(2,4);

    //first fit for prompt background, so the lifetime resolution can be fixed for the other fits


    RooFitResult* fitres= bPro.fitTo(*MC_Pro_sub,Hesse(true),Minos(true));
    cBinPDFs[i].cd(7);
    RooPlot* frame1=BMass->frame();
    MC_Pro_sub->plotOn(frame1,Binning(32));
    bPro.plotOn(frame1,LineColor(4));
    frame1->Draw();
    bPro.paramOn(frame1,Layout(0.7,0.95,0.95),Format("NELU",AutoPrecision(3)));
    frame1->Draw();

    RooFitResult* fitres= decayPro.fitTo(*MC_Pro_sub,Hesse(true),Minos(true));
    cBinPDFs[i].cd(8).SetLogy();
    RooPlot* frame1=Bct->frame();
    MC_Pro_sub->plotOn(frame1,Binning(225));
    decayPro.plotOn(frame1,LineColor(4));
    frame1->Draw();
    decayPro.paramOn(frame1,Layout(0.5,0.95,0.95),Format("NELU",AutoPrecision(3)));
    frame1->GetXaxis()->SetRange(0,40);
    frame1->SetMinimum(0.1);
    frame1->Draw();

    meanResCore->setConstant(kTRUE);
    meanResTail->setConstant(kTRUE);
    sigmaResCore->setConstant(kTRUE);
    sigmaResTail->setConstant(kTRUE);
    fracResCore->setConstant(kTRUE);

    cBinPDFs[i].cd(1);
    RooFitResult* fitres= bSig.fitTo(*MC_Sig_sub,Hesse(true),Minos(true));
    RooPlot* frame1=BMass->frame();
    MC_Sig_sub->plotOn(frame1,Binning(80));
    bSig.plotOn(frame1,LineColor(4));
    frame1->Draw();
    bSig.paramOn(frame1,Layout(0.65,0.95,0.95),Format("NELU",AutoPrecision(3)));
    frame1->GetXaxis()->SetRange(35,65);
    frame1->Draw();
    MCPDFs->cd(i*2+1);
    frame1->Draw();

    cBinPDFs[i].cd(2).SetLogy();
    if(usePerEventErrors) {
      RooFitResult* fitres= decaySig.fitTo(*MC_Sig_sub,Hesse(true),Minos(true),ConditionalObservables(*BctE));
    } else {
      RooFitResult* fitres= decaySig.fitTo(*MC_Sig_sub,Hesse(true),Minos(true));
    }
    
    RooPlot* frame1=Bct->frame();
    MC_Sig_sub->plotOn(frame1,Binning(45));
    if (usePerEventErrors) {
      decaySig.plotOn(frame1,ProjWData(*MC_Sig_sub),LineColor(4));
    } else {
      decaySig.plotOn(frame1,LineColor(4));
    }
    frame1->SetMinimum(0.1);
    frame1->Draw();
    decaySig.paramOn(frame1,Layout(0.65,0.95,0.95),Format("NELU",AutoPrecision(3)));
    frame1->Draw();
  
    cBinPDFs[i].cd(3);
    // Don't fit for the peaking BB PDF shape in each bin because too few 
    // to get a good fit.
    //RooFitResult* fitres= bFeed.fitTo(*MC_Feed_sub,Hesse(true),Minos(true));
    RooPlot* frame1=BMass->frame();
    MC_Feed_sub->plotOn(frame1,Binning(32));
    bFeed.plotOn(frame1,LineColor(4));
    frame1->Draw();
    bFeed.paramOn(frame1,Layout(0.7,0.90,0.90),Format("NELU",AutoPrecision(3)));
    TPaveText* pbox = (TPaveText*) frame1->findObject("bFeed_paramBox");
    if (pbox) {
      pbox->AddText("Fixed parameters");
    }
    frame1->Draw();
    MCPDFs->cd(i*2+2);
    frame1->Draw();
    
    cBinPDFs[i].cd(4).SetLogy();
    RooPlot* frame1=Bct->frame();
    MC_Feed_sub->plotOn(frame1,Binning(90));
    decayFeed.plotOn(frame1,LineColor(4));
    frame1->Draw();
    decayFeed.paramOn(frame1,Layout(0.65,0.90,0.90),Format("NELU",AutoPrecision(3)));
    TPaveText* pbox = (TPaveText*) frame1->findObject("decayFeed_paramBox");
    if (pbox) {
      pbox->AddText("Fixed parameters");
    }
    frame1->SetMinimum(0.1);
    frame1->Draw();

    RooFitResult* fitres= bBB.fitTo(*MC_BB_sub,Hesse(true),Minos(true));
    cBinPDFs[i].cd(5);
    RooPlot* frame1=BMass->frame();
    MC_BB_sub->plotOn(frame1,Binning(32));
    bBB.plotOn(frame1,LineColor(4));
    frame1->Draw();
    bBB.paramOn(frame1,Layout(0.7,0.95,0.95),Format("NELU",AutoPrecision(3)));
    frame1->Draw();

    RooFitResult* fitres= decayBB.fitTo(*MC_BB_sub,Hesse(true),Minos(true));
    cBinPDFs[i].cd(6).SetLogy();
    RooPlot* frame1=Bct->frame();
    MC_BB_sub->plotOn(frame1,Binning(90));
    decayBB.plotOn(frame1,LineColor(4));
    frame1->Draw();
    decayBB.paramOn(frame1,Layout(0.5,0.95,0.95),Format("NELU",AutoPrecision(3)));
    frame1->SetMinimum(0.1);
    frame1->Draw();


    // after PDF fits now fix signal B mass shapes for signal and peaking BB bkg
    fracBCore.setConstant(kTRUE);
    fracBOut.setConstant(kTRUE);
    meanBCore.setConstant(kTRUE);
    meanBTail.setConstant(kTRUE);
    meanBOut.setConstant(kTRUE);
    widthBCore.setConstant(kTRUE);
    widthBTail.setConstant(kTRUE);
    widthBOut.setConstant(kTRUE);

    fracBOutFeed.setConstant(kTRUE);  
    fracBCoreFeed.setConstant(kTRUE);
    meanBCoreFeed.setConstant(kTRUE);
    widthBCoreFeed.setConstant(kTRUE);
    meanBTailFeed.setConstant(kTRUE);
    widthBTailFeed.setConstant(kTRUE);
    meanBOutFeed.setConstant(kTRUE);
    widthBOutFeed.setConstant(kTRUE);

    // fix signal B lifetime
    ctauB->setConstant(kTRUE);

    // fix non-peaking and prompt B mass shapes, for now
    shapeBB.setConstant(kTRUE);
    shapePro.setConstant(kTRUE);
   



    // next float only signal yeilds

    nbkgBB->setConstant(kFALSE);
    nbkgPro->setConstant(kFALSE);
    nfeed->setConstant(kFALSE);
    nsig->setConstant(kFALSE);

    meanResCore->setConstant(kTRUE);
    meanResTail->setConstant(kTRUE);
    sigmaResCore->setConstant(kTRUE);
    sigmaResTail->setConstant(kTRUE);
    fracResCore->setConstant(kTRUE);

    ctauB->setConstant(kTRUE);
    ctauBB->setConstant(kTRUE);
    ctauBB2->setConstant(kTRUE);
    fracBBDecay->setConstant(kTRUE);




    if (usePerEventErrors) {
      model2D.fitTo(*Data_subsample_sub,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true),ConditionalObservables(*BctE));
    } else {
      model2D.fitTo(*Data_subsample_sub,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));
    }

    // finally try floating yields, resolution and B mass shapes for prompt and non-peaking B

    nbkgBB->setConstant(kFALSE);
    nbkgPro->setConstant(kFALSE);
    nfeed->setConstant(kFALSE);
    nsig->setConstant(kFALSE);

    meanResCore->setConstant(kFALSE);
    meanResTail->setConstant(kFALSE);
    sigmaResCore->setConstant(kFALSE);
    sigmaResTail->setConstant(kFALSE);
    fracResCore->setConstant(kFALSE);

    shapeBB.setConstant(kFALSE);
    shapePro.setConstant(kFALSE);
 
    ctauB->setConstant(kTRUE);
    ctauBB->setConstant(kTRUE);
    ctauBB2->setConstant(kTRUE);
    fracBBDecay->setConstant(kTRUE);

    if(floatKsMass){
      meanBCore->setConstant(kFALSE);
      widthBCore->setConstant(kFALSE);
    }

    if (i==1) {
      meanResCore->setVal(0.00011);
      sigmaResCore->setVal(.00366);
      meanResTail->setVal(-0.00033);
      sigmaResTail->setVal(.01102);
      fracResCore->setVal(0.874);
    }

    if(i>2) {
      meanResCore->setVal(meanResCoreValue);
      sigmaResCore->setVal(sigmaResCoreValue);
      meanResTail->setVal(meanResTailValue);
      sigmaResTail->setVal(sigmaResTailValue);
      fracResCore->setVal(fracResCoreValue);
      meanResCore->setConstant(kTRUE);
      meanResTail->setConstant(kTRUE);
      sigmaResCore->setConstant(kTRUE);
      sigmaResTail->setConstant(kTRUE);
      fracResCore->setConstant(kTRUE);
    }

    cout << " before fit to data for bin " << i << endl;
    if (usePerEventErrors) {
      model2D.fitTo(*Data_subsample_sub,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true),ConditionalObservables(*BctE));
    } else {
      model2D.fitTo(*Data_subsample_sub,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));
    }
    cout << " after fit to data for bin " << i << endl;

    if (i==2) {
      meanResCoreValue = meanResCore->getVal();
      sigmaResCoreValue = sigmaResCore->getVal();
      meanResTailValue = meanResTail->getVal();
      sigmaResTailValue = sigmaResTail->getVal();
      fracResCoreValue = fracResCore->getVal();
    }

    if (i<3) {cBinProj->cd(i*3+1); cout << "cd to cbinproj = " << i*3+1 << endl;}
    else { cBinProj2->cd((i-3)*3+1); cout << "cd to cbinproj2 = " << (i-3)*3+1 << endl;}
    RooPlot* frame1=BMass->frame();
    Data_subsample_sub->plotOn(frame1,Binning(32));
    model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(4));
    frame1->Draw();
    if (i<3) cBinProj->cd(i*3+2);
    else cBinProj2->cd((i-3)*3+2);
    RooPlot* frame1=BMass->frame();
    Data_subsample_sub->plotOn(frame1,Binning(32),CutRange("ct0p1"));
    model2D.plotOn(frame1,LineColor(3),ProjectionRange("ct0p1"),Components(RooArgSet(bkgPro2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(2),ProjectionRange("ct0p1"),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(6),ProjectionRange("ct0p1"),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(4),ProjectionRange("ct0p1"));
    frame1->Draw();
    if (i<3) cBinProj->cd(i*3+3).SetLogy();
    else cBinProj2->cd((i-3)*3+3).SetLogy();
    RooPlot* frame1=Bct->frame();
    Data_subsample_sub->plotOn(frame1);
    model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
    if (usePerEventErrors) {
      model2D.plotOn(frame1,LineColor(4),ProjWData(*Data_subsample_sub));  
    } else {
      model2D.plotOn(frame1,LineColor(4));
    }
    frame1->Draw();
    frame1->SetMinimum(0.1);
    model2D.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
    frame1->Draw();

cout << "moving to cbinproj3 cd = " << i << endl;
    cBinProj3->cd(i+1);
    RooPlot* frame1=Bct->frame();
    Data_subsample_sub->plotOn(frame1);
    model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
    model2D.plotOn(frame1,LineColor(4));
    
    frame1->Draw();
    model2D.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
    frame1->Draw();

    SigYield[i]=nsig->getVal();
    ErrYield[i]=nsig->getError();

    outTextFile << "====================================================" << endl;
    outTextFile << "Summary of fit results from bin" << i << endl;
    outTextFile << " nSig = " << nsig->getVal() << "+-" << nsig->getError()  << endl;
    outTextFile << " nFeed = " << nfeed->getVal() << "+-" << nfeed->getError()  << endl;
    outTextFile << " nBkgBB = " << nbkgBB->getVal() << "+-" << nbkgBB->getError() << endl;
    outTextFile << " nBkgPro = " << nbkgPro->getVal() << "+-" << nbkgPro->getError()  << endl;
    outTextFile << cut1 << "-" << cut2 << " & $" << Form("%.1f",nsig->getVal()) << "\\pm" << Form ("%.1f",nsig->getError()) << 
      "$ & $" << Form("%.1f",nfeed->getVal()) << "\\pm" << Form ("%.1f",nfeed->getError()) << 
      "$ & $" << Form("%.1f",nbkgBB->getVal()) << "\\pm" << Form ("%.1f",nbkgBB->getError()) << 
      "$ & $" << Form("%.1f",nbkgPro->getVal()) << "\\pm" << Form ("%.1f",nbkgPro->getError()) << "$ \\\\" << endl;
    outTextFile <<  "====================================================" << endl;



  outTextFile << endl;
  outTextFile << "Results from fit to pt bin " << i << endl;
  outTextFile << "==========================================" << endl;
  outTextFile << " Signal yield & " << nsig->getVal() << "$\\pm$" << nsig->getError() << endl;
  outTextFile << " Peaking BB yield & " << nfeed->getVal() << "$\\pm$" << nfeed->getError() << endl;
  outTextFile << " Combinatorial BB yield & " << nbkgBB->getVal() << "$\\pm$" << nbkgBB->getError()  << endl;
  outTextFile << " Prompt yield & " << nbkgPro->getVal() << "$\\pm$" << nbkgPro->getError() << endl;
  outTextFile << " ct core mean & " << meanResCore->getVal() << "$\\pm$" << meanResCore->getError() << endl;
  outTextFile << " ct tail mean & " << meanResTail->getVal() << "$\\pm$" << meanResTail->getError() << endl;
  outTextFile << " ct core sigma & " << sigmaResCore->getVal() << "$\\pm$" << sigmaResCore->getError() << endl;
  outTextFile << " ct tail sigma & " <<  sigmaResTail->getVal() << "$\\pm$" << sigmaResTail->getError() << endl;
  outTextFile << " ct core fraction & " <<  fracResCore->getVal() << "$\\pm$" << fracResCore->getError() << endl;
  outTextFile << " Combinatorial BB shape & " <<  shapeBB->getVal() << "$\\pm$" << shapeBB->getError() << endl;
  outTextFile << " Prompt shape & " << shapePro->getVal() << "$\\pm$" << shapePro->getError() << endl;
  outTextFile << "=============================================================================================" << endl;

    // last fit high B mass sideband for validation plots. Only float yields.
    // BB lifetimes fixed to values from full fit.

    nbkgBB->setConstant(kFALSE);
    nbkgPro->setConstant(kFALSE);
    nsig->setConstant(kTRUE);
    nfeed->setConstant(kTRUE);
    
    meanResCore->setConstant(kTRUE);
    meanResTail->setConstant(kTRUE);
    sigmaResCore->setConstant(kTRUE);
    sigmaResTail->setConstant(kTRUE);
    fracResCore->setConstant(kTRUE);
    shapeBB.setConstant(kTRUE);
    shapePro.setConstant(kTRUE);
    
    cout<< "fitting sideband region for pt bin = " << i << endl;
    
    model2DSB.fitTo(*Data_HighMassSB_subsample_sub,Hesse(true),Minos(true),Timer(true),Save(true));
  

    SBCanv->cd(i*2+1);
    RooPlot* frame3=BMass->frame();
    Data_HighMassSB_subsample_sub->plotOn(frame3, Binning(32));
    model2DSB.plotOn(frame3,LineColor(3),Range("HighMassSB"), Components(RooArgSet(bkgPro2De)),LineStyle(7));
    model2DSB.plotOn(frame3,LineColor(4),Range("HighMassSB"), Components(RooArgSet(bkgPro2De,bkgBB2De)));
    frame3->Draw();

    SBCanv->cd(i*2+2).SetLogy();
    RooPlot* frame2=Bct->frame();
    Data_HighMassSB_subsample_sub->plotOn(frame2, CutRange("HighMassSB"));
    model2DSB.plotOn(frame2,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
    model2DSB.plotOn(frame2,LineColor(4));
    frame2->Draw();

    frame2->SetMinimum(0.01);
    model2DSB.paramOn(frame2, Layout(0.65,0.95,0.90), Format("NELU",AutoPrecision(3)) ) ;
    frame2->getAttText()->SetTextSize(0.04); 
    frame2->Draw();

    char file[150];
    sprintf(file,Form("%s%i%s","/JPsiKsFitPDFs_data_ptbin", i, ".pdf") );
    cout << "Saving pdfs for data bin with file name " << file << endl;
    cBinPDFs[i]->SaveAs("JPsiKsFit_dataRes/"+trigger+file);

    /////////////////////////////////////////////////////////////////
    //next let's do some toy experiments for this bin
    if (doToys) {


      outTextFile << "================================================" << endl;
      outTextFile << "Starting pure toys" << endl;
      outTextFile << "===================================================" << endl;

      outTextFile << "Now Ntot = " << Ntot << "  nSig = " << nsig->getVal() << " nbkgBB peaking= " << nfeed->getVal() << "  nbkgBB nonpeaking = " <<
      nbkgBB->getVal() << "  nPrompt = " << nbkgPro->getVal() << endl;

      Ntot=nsig->getVal()+nfeed->getVal()+nbkgBB->getVal()+nbkgPro->getVal();

      ctauB.setConstant(kTRUE);
      fracResCore->setConstant(kTRUE);
      meanResCore->setConstant(kTRUE);
      meanResTail->setConstant(kTRUE);
      sigmaResCore->setConstant(kTRUE);
      sigmaResTail->setConstant(kTRUE);
      shapeBB->setConstant(kTRUE);
      shapePro->setConstant(kTRUE);

      if (usePerEventErrors) {
        RooMCStudy mgrBin(model2D,RooArgList(*BMass,*Bct,*BctE),FitModel(model2D),Extended(1),ConditionalObservables(RooArgSet(*BctE)),Minos(1),Hesse(1),SplitRange(false));
        mgrBin.generateAndFit(nToys ,Ntot,kTRUE);
        mgrBin.Write();
        RooDataSet* toySample = mgrBin.genData(0);

        model2D.fitTo(*toySample,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true),ConditionalObservables(RooArgSet(*BctE)));
      } else {
        RooMCStudy mgrBin(model2D,RooArgList(*BMass,*Bct),FitModel(model2D),Extended(1),Minos(1),Hesse(1),SplitRange(false));
        mgrBin.generateAndFit(nToys ,Ntot,kTRUE);
        mgrBin.Write();
        RooDataSet* toySample = mgrBin.genData(0);

        model2D.fitTo(*toySample,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));
      }


      TCanvas *cToy1 = new TCanvas("cToy1","cToy1",800,800);
      cToy1->Divide(2,1);
      cToy1->cd(1);
      RooPlot* frame1=BMass->frame();
      toySample->plotOn(frame1,Binning(32));
      model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(4),Components(RooArgSet(bkgBB2De,bkgPro2De,sig2De,feed2De)));
      frame1->Draw();
      model2D.paramOn(frame1,toySample,"Fit parameters",1,"NELU");
      frame1->Draw();
      TPaveText* pbox = (TPaveText*) frame1->findObject("model2D_paramBox");
      if (pbox) {
        pbox->AddText(Form("%s%.3f", "#chi^{2}/DOF = ", frame1->chiSquare()));
      }
 
      cToy1->cd(2).SetLogy();
      RooPlot* frame1=Bct->frame();
      toySample->plotOn(frame1);
      model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(4),Components(RooArgSet(bkgBB2De,bkgPro2De,sig2De,feed2De)));
      frame1->Draw();
      model2D.paramOn(frame1,toySample,"Fit parameters",1,"NELU");
      frame1->Draw();
      TPaveText* pbox = (TPaveText*) frame1->findObject("model2D_paramBox");
      if (pbox) {
        pbox->AddText(Form("%s%.3f", "#chi^{2}/DOF = ", frame1->chiSquare()));
      }
    
    
      sprintf(file,Form("%s%i%s","/pToys_sample_ptbin", i, ".pdf") );
      cout << "Saving toy plot for bin with file name " << file << endl;
      cToy1->SaveAs("JPsiKsFit_dataRes/"+trigger+file);

      TCanvas *cToy3 = new TCanvas("cToy3","cToy3",800,800);
      cToy3->Divide(2,2);
  
      cToy3->cd(1);
      RooPlot* bmassPullFrame = mgrBin.plotPull(*nsig,-3.0,3.0,60,kTRUE);
      bmassPullFrame->Draw();
  
      cToy3->cd(2);
      RooPlot* frame = mgrBin.plotPull(*nfeed,-3.0,3.0,60,kTRUE);
      frame->Draw();
  
      cToy3->cd(3);
      RooPlot* frame = mgrBin.plotPull(*nbkgBB,-3.0,3.0,60,kTRUE);
      frame->Draw();
  
      cToy3->cd(4);
      RooPlot* frame = mgrBin.plotPull(*nbkgPro,-3.0,3.0,60,kTRUE);
      frame->Draw();
    
      char file2[150];
      sprintf(file,Form("%s%i%s","/pToys_pulls_ptbin", i, ".pdf") );
      cout << "Saving toy plot for bin with file name " << file << endl;
      cToy1->SaveAs("JPsiKsFit_dataRes/"+trigger+file);
      cToy3->SaveAs("JPsiKsFit_dataRes/"+trigger+file);

    }

outTextFile << " Done with bin " << i << endl;
outTextFile << "===============================================" << endl;

  } //end loop over bins



  cBinProj->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResults_data_ptbins1.pdf");
  cBinProj2->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResults_data_ptbins2.pdf");
  cBinProj3->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResults_data_ptbins3.pdf");
  MCPDFs->SaveAs("JPsiKsFit_dataRes/"+trigger+"/BMassPDFs_sigPeak_ptbins.pdf");
  SBCanv->SaveAs("JPsiKsFit_dataRes/"+trigger+"/SB_ptbins.pdf");

  for(int i=0;i<nBins;i++){
    outTextFile << "Pt bin " << i << " yield: " << SigYield[i]  
    << " +/- " << ErrYield[i] << "\n";
  }
  
  } // closing if (doPtBins) loop


  if (doYBins) {

    ////////////////////////////////////////////////////////////////////////////////////
    // Now split data into different Y bins. Still need MC to fit for PDFs in each bin. Would be better to get PDFs from data here, too.

    int nBins = 5;
    float SigYield[nBins];
    float ErrYield[nBins];

    TCanvas cBinPDFs[nBins];
    TCanvas *cBinProj = new TCanvas("cBinProj","cBinProj",800,800);
    TCanvas *cBinProj2 = new TCanvas("cBinProj2","cBinProj2",800,600);
    TCanvas *cBinProj3 = new TCanvas("cBinProj3","cBinProj3",800,800);
    cBinProj->Divide(3,3);
    cBinProj2->Divide(3,2);  
    cBinProj3->Divide(2,3);

    for(int i=0;i<nBins;i++){

      float cut1(0);
      float cut2(0);


    if (useDavidBinning) {
      if(i==0) {cut1=0.0;cut2=0.6;}
      if(i==1) {cut1=0.6;cut2=1.1;}
      if(i==2) {cut1=1.1;cut2=1.45;}
      if(i==3) {cut1=1.45;cut2=1.8;}
      if(i==4) {cut1=1.8;cut2=2.2;}
    } else {
      if(i==0) {cut1=0.0;cut2=0.5;}
      if(i==1) {cut1=0.5;cut2=1.0;}
      if(i==2) {cut1=1.0;cut2=1.4;}
      if(i==3) {cut1=1.4;cut2=1.8;}
      if(i==4) {cut1=1.8;cut2=2.2;}  
    }
  
      //create datasets for each Y cut range for each component

      char ptCut[150];
      sprintf(ptCut,"abs(bY)>=%g && abs(bY)<%g",cut1,cut2);
    
      outTextFile << "Y cut = " << ptCut << endl;
      RooAbsData* MC_Sig_sub = MC_SigPrivate->reduce(ptCut);
      RooAbsData* MC_BB_sub = MC_BkgBB->reduce(ptCut);
      RooAbsData* MC_Feed_sub = MC_Feed->reduce(ptCut);
      RooAbsData* MC_Pro_sub = MC_BkgPro->reduce(ptCut);
      RooAbsData* Data_subsample_sub = Data_All->reduce(ptCut);

      RooAbsData* Data_HighMassSB_subsample_sub = Data_HighMassSB->reduce(ptCut);

      // reset all params that need to float in the PDF fits. Done for each bin separately.
  
      //signal and peaking BB B mass params will be taken directly from these MC fits
      //signal B mass params
      fracBCore.setConstant(kFALSE);
      fracBOut.setConstant(kFALSE);
      meanBCore.setConstant(kFALSE);
      meanBTail.setConstant(kFALSE);
      meanBOut.setConstant(kFALSE);
      widthBCore.setConstant(kFALSE);
      widthBTail.setConstant(kFALSE);
      widthBOut.setConstant(kFALSE);
  
      // peaking BB B mass params
      // There's not enough MC bin by bin, so use the value from the single fit
      fracBOutFeed.setConstant(kFALSE);  
      fracBCoreFeed.setConstant(kFALSE);
      meanBCoreFeed.setConstant(kFALSE);
      widthBCoreFeed.setConstant(kFALSE);
      meanBTailFeed.setConstant(kFALSE);
      widthBTailFeed.setConstant(kFALSE);
      meanBOutFeed.setConstant(kFALSE);
      widthBOutFeed.setConstant(kFALSE);
      
      // non-peaking BB and prompt B mass params will float in final fit, but get
      // MC values here as a good starting point
      // non-peaking BB B mass params
      shapeBB.setConstant(kFALSE);
      // prompt B mass params
      shapePro.setConstant(kFALSE);
      
      // lifetime resolution will come from data, but fit MC signal here to
      // get good starting values
      //signal lifetime parameters
      meanResCore.setConstant(kFALSE);
      sigmaResCore.setConstant(kFALSE);
      meanResTail.setConstant(kFALSE);
      sigmaResTail.setConstant(kFALSE);
      fracResCore.setConstant(kFALSE);
      
      //keep signal lifetime value from fit to full range
      ctauB.setConstant(kTRUE);
  
      //peaking BB lifetime params
      //take this from MC. Should check if different in different bins or not.
      ctauFeed.setConstant(kTRUE);
  
      //non-peaking BB lifetime params
      //will get from data. Should check if it's important to do it bin by bin or not.
      // for now assume it's ok and use the values from the single fit above in data.
      ctauBB.setConstant(kTRUE);
      ctauBB2.setConstant(kTRUE);
      fracBBDecay.setConstant(kTRUE);
  
  
      cBinPDFs[i].Divide(2,4);
  
      //first fit for prompt background, so the lifetime resolution can be fixed for the other fits
  
  
      RooFitResult* fitres= bPro.fitTo(*MC_Pro_sub,Hesse(true),Minos(true));
      cBinPDFs[i].cd(7);
      RooPlot* frame1=BMass->frame();
      MC_Pro_sub->plotOn(frame1,Binning(32));
      bPro.plotOn(frame1,LineColor(4));
      frame1->Draw();
      bPro.paramOn(frame1,Layout(0.7,0.95,0.95),Format("NELU",AutoPrecision(3)));
      frame1->Draw();
  
      RooFitResult* fitres= decayPro.fitTo(*MC_Pro_sub,Hesse(true),Minos(true));
      cBinPDFs[i].cd(8).SetLogy();
      RooPlot* frame1=Bct->frame();
      MC_Pro_sub->plotOn(frame1,Binning(225));
      decayPro.plotOn(frame1,LineColor(4));
      frame1->Draw();
      decayPro.paramOn(frame1,Layout(0.5,0.95,0.95),Format("NELU",AutoPrecision(3)));
      frame1->GetXaxis()->SetRange(0,40);
      frame1->SetMinimum(0.1);
      frame1->Draw();
  
      meanResCore->setConstant(kTRUE);
      meanResTail->setConstant(kTRUE);
      sigmaResCore->setConstant(kTRUE);
      sigmaResTail->setConstant(kTRUE);
      fracResCore->setConstant(kTRUE);
  
      fracBCore->setVal(0.8);
      widthBCore->setVal(0.01);
      widthBTail->setVal(0.03);
      
      cBinPDFs[i].cd(1);
      RooFitResult* fitres= bSig.fitTo(*MC_Sig_sub,Hesse(true),Minos(true));
      RooPlot* frame1=BMass->frame();
      MC_Sig_sub->plotOn(frame1,Binning(80));
      bSig.plotOn(frame1,LineColor(4));
      frame1->Draw();
      bSig.paramOn(frame1,Layout(0.65,0.95,0.95),Format("NELU",AutoPrecision(3)));
      frame1->GetXaxis()->SetRange(35,65);
      frame1->Draw();
      MCPDFs->cd(i*2+1);
      frame1->Draw();

      cBinPDFs[i].cd(2).SetLogy();
      if(usePerEventErrors) {
  	RooFitResult* fitres= decaySig.fitTo(*MC_Sig_sub,Hesse(true),Minos(true),ConditionalObservables(*BctE));
      } else {
  	RooFitResult* fitres= decaySig.fitTo(*MC_Sig_sub,Hesse(true),Minos(true));
      }
      
      RooPlot* frame1=Bct->frame();
      MC_Sig_sub->plotOn(frame1,Binning(45));
      if (usePerEventErrors) {
  	decaySig.plotOn(frame1,ProjWData(*MC_Sig_sub),LineColor(4));

      } else {
  	decaySig.plotOn(frame1,LineColor(4));

      }
      frame1->SetMinimum(0.1);
      frame1->Draw();
      decaySig.paramOn(frame1,Layout(0.65,0.95,0.95),Format("NELU",AutoPrecision(3)));
      frame1->Draw();
    
      cBinPDFs[i].cd(3);
      // Don't fit for the peaking BB PDF shape in each bin because too few 
      // to get a good fit.
      //RooFitResult* fitres= bFeed.fitTo(*MC_Feed_sub,Hesse(true),Minos(true));
      RooPlot* frame1=BMass->frame();
      MC_Feed_sub->plotOn(frame1,Binning(32));
      bFeed.plotOn(frame1,LineColor(4));
      MCPDFs->cd(i*2+2);
      frame1->Draw();


      frame1->Draw();
      bFeed.paramOn(frame1,Layout(0.7,0.90,0.90),Format("NELU",AutoPrecision(3)));
      TPaveText* pbox = (TPaveText*) frame1->findObject("bFeed_paramBox");
      if (pbox) {
  	pbox->AddText("Fixed parameters");
      }
      frame1->Draw();
  
      cBinPDFs[i].cd(4).SetLogy();
      RooPlot* frame1=Bct->frame();
      MC_Feed_sub->plotOn(frame1,Binning(90));
      decayFeed.plotOn(frame1,LineColor(4));
      frame1->Draw();
      decayFeed.paramOn(frame1,Layout(0.65,0.90,0.90),Format("NELU",AutoPrecision(3)));
      TPaveText* pbox = (TPaveText*) frame1->findObject("decayFeed_paramBox");
      if (pbox) {
  	pbox->AddText("Fixed parameters");
      }
      frame1->SetMinimum(0.1);
      frame1->Draw();
  
      RooFitResult* fitres= bBB.fitTo(*MC_BB_sub,Hesse(true),Minos(true));
      cBinPDFs[i].cd(5);
      RooPlot* frame1=BMass->frame();
      MC_BB_sub->plotOn(frame1,Binning(32));
      bBB.plotOn(frame1,LineColor(4));
      frame1->Draw();
      bBB.paramOn(frame1,Layout(0.7,0.95,0.95),Format("NELU",AutoPrecision(3)));
      frame1->Draw();

      RooFitResult* fitres= decayBB.fitTo(*MC_BB_sub,Hesse(true),Minos(true));
      cBinPDFs[i].cd(6).SetLogy();
      RooPlot* frame1=Bct->frame();
      MC_BB_sub->plotOn(frame1,Binning(90));
      decayBB.plotOn(frame1,LineColor(4));
      frame1->Draw();
      decayBB.paramOn(frame1,Layout(0.5,0.95,0.95),Format("NELU",AutoPrecision(3)));
      frame1->SetMinimum(0.1);
      frame1->Draw();
  
  
      // after PDF fits now fix signal B mass shapes for signal and peaking BB bkg
      fracBCore.setConstant(kTRUE);
      fracBOut.setConstant(kTRUE);
      meanBCore.setConstant(kTRUE);
      meanBTail.setConstant(kTRUE);
      meanBOut.setConstant(kTRUE);
      widthBCore.setConstant(kTRUE);
      widthBTail.setConstant(kTRUE);
      widthBOut.setConstant(kTRUE);
  
      fracBOutFeed.setConstant(kTRUE);  
      fracBCoreFeed.setConstant(kTRUE);
      meanBCoreFeed.setConstant(kTRUE);
      widthBCoreFeed.setConstant(kTRUE);
      meanBTailFeed.setConstant(kTRUE);
      widthBTailFeed.setConstant(kTRUE);
      meanBOutFeed.setConstant(kTRUE);
      widthBOutFeed.setConstant(kTRUE);
  
      // fix signal B lifetime
      ctauB->setConstant(kTRUE);
  
      // fix non-peaking and prompt B mass shapes, for now
      shapeBB.setConstant(kTRUE);
      shapePro.setConstant(kTRUE);
     
  
      // first fit high B mass sideband only float yields. BB background lifetime params from
      // fit to full pT range. 
  
      nbkgBB->setConstant(kFALSE);
      nbkgPro->setConstant(kFALSE);
      ctauBB.setConstant(kTRUE);
      ctauBB2.setConstant(kTRUE);
      fracBBDecay->setConstant(kTRUE);
  
      meanResCore->setConstant(kTRUE);
      meanResTail->setConstant(kTRUE);
      sigmaResCore->setConstant(kTRUE);
      sigmaResTail->setConstant(kTRUE);
      fracResCore->setConstant(kTRUE);
  
      cout << "starting fit to data for bin " << i << endl;
  
      model2DSB.fitTo(*Data_HighMassSB_subsample_sub,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));
  
  
      // next float only signal yeilds
  
      nbkgBB->setConstant(kFALSE);
      nbkgPro->setConstant(kFALSE);
      nfeed->setConstant(kFALSE);
      nsig->setConstant(kFALSE);
  
      meanResCore->setConstant(kTRUE);
      meanResTail->setConstant(kTRUE);
      sigmaResCore->setConstant(kTRUE);
      sigmaResTail->setConstant(kTRUE);
      fracResCore->setConstant(kTRUE);
  
      ctauB->setConstant(kTRUE);
      ctauBB->setConstant(kTRUE);
      ctauBB2->setConstant(kTRUE);
      fracBBDecay->setConstant(kTRUE);
  
  
  
  
      if (usePerEventErrors) {
  	model2D.fitTo(*Data_subsample_sub,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true),ConditionalObservables(*BctE));
      } else {
  	model2D.fitTo(*Data_subsample_sub,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));
      }
  
      // finally try floating yields, resolution and B mass shapes for prompt and non-peaking B
  
      nbkgBB->setConstant(kFALSE);
      nbkgPro->setConstant(kFALSE);
      nfeed->setConstant(kFALSE);
      nsig->setConstant(kFALSE);
  
      meanResCore->setConstant(kFALSE);
      meanResTail->setConstant(kFALSE);
      sigmaResCore->setConstant(kFALSE);
      sigmaResTail->setConstant(kFALSE);
      fracResCore->setConstant(kFALSE);
  
      shapeBB.setConstant(kFALSE);
      shapePro.setConstant(kFALSE);

      if(floatKsMass) {
        meanBCore->setConstant(kFALSE);
        widthBCore->setConstant(kFALSE);
      }
   
      ctauB->setConstant(kTRUE);
      ctauBB->setConstant(kTRUE);
      ctauBB2->setConstant(kTRUE);
      fracBBDecay->setConstant(kTRUE);
  
  
      cout << " before fit to data for bin " << i << endl;
      if (usePerEventErrors) {
  	model2D.fitTo(*Data_subsample_sub,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true),ConditionalObservables(*BctE));
      } else {
  	model2D.fitTo(*Data_subsample_sub,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));
      }
      cout << " after fit to data for bin " << i << endl;
  
  
  
      if (i<3) cBinProj->cd(i*3+1);
      else cBinProj2->cd((i-3)*3+1);
      RooPlot* frame1=BMass->frame();
      Data_subsample_sub->plotOn(frame1,Binning(32));
      model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(4));
      frame1->Draw();
      if (i<3) cBinProj->cd(i*3+2);
      else cBinProj2->cd((i-3)*3+2);
      RooPlot* frame1=BMass->frame();
      Data_subsample_sub->plotOn(frame1,Binning(32),CutRange("ct0p1"));
      model2D.plotOn(frame1,LineColor(3),ProjectionRange("ct0p1"),Components(RooArgSet(bkgPro2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(2),ProjectionRange("ct0p1"),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(6),ProjectionRange("ct0p1"),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(4),ProjectionRange("ct0p1"));
      frame1->Draw();      
      if (i<3) cBinProj->cd(i*3+3).SetLogy();
      else cBinProj2->cd((i-3)*3+3).SetLogy();
      RooPlot* frame1=Bct->frame();
      Data_subsample_sub->plotOn(frame1);
      model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
      if (usePerEventErrors) {
  	model2D.plotOn(frame1,LineColor(4),ProjWData(*Data_subsample_sub));  
      } else {
  	model2D.plotOn(frame1,LineColor(4));
      }
      frame1->Draw();
      frame1->SetMinimum(0.1);
      model2D.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
      frame1->Draw();
      
      cBinProj3->cd(i+1);
      
      RooPlot* frame1=Bct->frame();
      Data_subsample_sub->plotOn(frame1);
      model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
      model2D.plotOn(frame1,LineColor(4));
      
      frame1->Draw();
      model2D.paramOn(frame1, Layout(0.65,0.95,0.95), Format("NELU",AutoPrecision(3)) ) ;
      frame1->Draw();      
    
      SigYield[i]=nsig->getVal();
      ErrYield[i]=nsig->getError();
  
      outTextFile << "====================================================" << endl;
      outTextFile << "Summary of fit results from y bin" << i << endl;
      outTextFile << " nSig = " << nsig->getVal() << "+-" << nsig->getError()  << endl;
      outTextFile << " nFeed = " << nfeed->getVal() << "+-" << nfeed->getError()  << endl;
      outTextFile << " nBkgBB = " << nbkgBB->getVal() << "+-" << nbkgBB->getError() << endl;
      outTextFile << " nBkgPro = " << nbkgPro->getVal() << "+-" << nbkgPro->getError()  << endl;
      outTextFile << cut1 << "-" << cut2 << " & $" << Form("%.1f",nsig->getVal()) << "\\pm" << Form ("%.1f",nsig->getError()) << 
  	"$ & $" << Form("%.1f",nfeed->getVal()) << "\\pm" << Form ("%.1f",nfeed->getError()) << 
  	"$ & $" << Form("%.1f",nbkgBB->getVal()) << "\\pm" << Form ("%.1f",nbkgBB->getError()) << 
  	"$ & $" << Form("%.1f",nbkgPro->getVal()) << "\\pm" << Form ("%.1f",nbkgPro->getError()) << "$ \\\\" << endl;
      outTextFile <<  "====================================================" << endl;
  
      outTextFile << endl;
      outTextFile << "Results from fit to y bin " << i << endl;
      outTextFile << "==========================================" << endl;
      outTextFile << " Signal yield & " << nsig->getVal() << "$\\pm$" << nsig->getError() << endl;
      outTextFile << " Peaking BB yield & " << nfeed->getVal() << "$\\pm$" << nfeed->getError() << endl;
      outTextFile << " Combinatorial BB yield & " << nbkgBB->getVal() << "$\\pm$" << nbkgBB->getError()  << endl;
      outTextFile << " Prompt yield & " << nbkgPro->getVal() << "$\\pm$" << nbkgPro->getError() << endl;
      outTextFile << " ct core mean & " << meanResCore->getVal() << "$\\pm$" << meanResCore->getError() << endl;
      outTextFile << " ct tail mean & " << meanResTail->getVal() << "$\\pm$" << meanResTail->getError() << endl;
      outTextFile << " ct core sigma & " << sigmaResCore->getVal() << "$\\pm$" << sigmaResCore->getError() << endl;
      outTextFile << " ct tail sigma & " <<  sigmaResTail->getVal() << "$\\pm$" << sigmaResTail->getError() << endl;
      outTextFile << " ct core fraction & " <<  fracResCore->getVal() << "$\\pm$" << fracResCore->getError() << endl;
      outTextFile << " Combinatorial BB shape & " <<  shapeBB->getVal() << "$\\pm$" << shapeBB->getError() << endl;
      outTextFile << " Prompt shape & " << shapePro->getVal() << "$\\pm$" << shapePro->getError() << endl;
      outTextFile << "=============================================================================================" << endl;



      char file[150];
      sprintf(file,Form("%s%i%s","/JPsiKsFitPDFs_data_ybin", i, ".pdf") );
      cout << "Saving pdfs for data bin with file name " << file << endl;
      cBinPDFs[i]->SaveAs("JPsiKsFit_dataRes/"+trigger+file);
  
      /////////////////////////////////////////////////////////////////
      //next let's do some toy experiments for this bin
      if (doToys) {
  
  
  	outTextFile << "================================================" << endl;
  	outTextFile << "Starting pure toys" << endl;
  	outTextFile << "===================================================" << endl;
  
  	outTextFile << "Now Ntot = " << Ntot << "  nSig = " << nsig->getVal() << " nbkgBB peaking= " << nfeed->getVal() << "  nbkgBB nonpeaking = " <<
  	nbkgBB->getVal() << "  nPrompt = " << nbkgPro->getVal() << endl;
  
  	Ntot=nsig->getVal()+nfeed->getVal()+nbkgBB->getVal()+nbkgPro->getVal();
  
  	ctauB.setConstant(kTRUE);
  	fracResCore->setConstant(kTRUE);
  	meanResCore->setConstant(kTRUE);
  	meanResTail->setConstant(kTRUE);
  	sigmaResCore->setConstant(kTRUE);
  	sigmaResTail->setConstant(kTRUE);
  	shapeBB->setConstant(kTRUE);
  	shapePro->setConstant(kTRUE);
  
  	if (usePerEventErrors) {
  	  RooMCStudy mgrBin(model2D,RooArgList(*BMass,*Bct,*BctE),FitModel(model2D),Extended(1),ConditionalObservables(RooArgSet(*BctE)),Minos(1),Hesse(1),SplitRange(false));
  	  mgrBin.generateAndFit(nToys ,Ntot,kTRUE);
  	  mgrBin.Write();
  	  RooDataSet* toySample = mgrBin.genData(0);
  
  	  model2D.fitTo(*toySample,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true),ConditionalObservables(RooArgSet(*BctE)));
  	} else {
  	  RooMCStudy mgrBin(model2D,RooArgList(*BMass,*Bct),FitModel(model2D),Extended(1),Minos(1),Hesse(1),SplitRange(false));
  	  mgrBin.generateAndFit(nToys ,Ntot,kTRUE);
  	  mgrBin.Write();
  	  RooDataSet* toySample = mgrBin.genData(0);
  
  	  model2D.fitTo(*toySample,Extended(true),SplitRange(false),Hesse(true),Minos(true),Timer(true),Save(true));
  	}
  
  
  	TCanvas *cToy1 = new TCanvas("cToy1","cToy1",800,800);
  	cToy1->Divide(2,1);
  	cToy1->cd(1);
  	RooPlot* frame1=BMass->frame();
  	toySample->plotOn(frame1,Binning(32));
  	model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
  	model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
  	model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
  	model2D.plotOn(frame1,LineColor(4),Components(RooArgSet(bkgBB2De,bkgPro2De,sig2De,feed2De)));
  	frame1->Draw();
  	model2D.paramOn(frame1,toySample,"Fit parameters",1,"NELU");
  	frame1->Draw();
  	TPaveText* pbox = (TPaveText*) frame1->findObject("model2D_paramBox");
  	if (pbox) {
  	  pbox->AddText(Form("%s%.3f", "#chi^{2}/DOF = ", frame1->chiSquare()));
  	}
   
  	cToy1->cd(2).SetLogy();
  	RooPlot* frame1=Bct->frame();
  	toySample->plotOn(frame1);
  	model2D.plotOn(frame1,LineColor(3),Components(RooArgSet(bkgPro2De)),LineStyle(7));
  	model2D.plotOn(frame1,LineColor(2),Components(RooArgSet(bkgPro2De,bkgBB2De)),LineStyle(7));
  	model2D.plotOn(frame1,LineColor(6),Components(RooArgSet(bkgPro2De,bkgBB2De,feed2De)),LineStyle(7));
  	model2D.plotOn(frame1,LineColor(4),Components(RooArgSet(bkgBB2De,bkgPro2De,sig2De,feed2De)));
  	frame1->Draw();
  	model2D.paramOn(frame1,toySample,"Fit parameters",1,"NELU");
  	frame1->Draw();
  	TPaveText* pbox = (TPaveText*) frame1->findObject("model2D_paramBox");
  	if (pbox) {
  	  pbox->AddText(Form("%s%.3f", "#chi^{2}/DOF = ", frame1->chiSquare()));
  	}
      
      
  	sprintf(file,Form("%s%i%s","/pToys_sample_ybin", i, ".pdf") );
  	cout << "Saving toy plot for bin with file name " << file << endl;
  	cToy1->SaveAs("JPsiKsFit_dataRes/"+trigger+file);

  	TCanvas *cToy3 = new TCanvas("cToy3","cToy3",800,800);
  	cToy3->Divide(2,2);
  
  	cToy3->cd(1);
  	RooPlot* bmassPullFrame = mgrBin.plotPull(*nsig,-3.0,3.0,60,kTRUE);
  	bmassPullFrame->Draw();
    
  	cToy3->cd(2);
  	RooPlot* frame = mgrBin.plotPull(*nfeed,-3.0,3.0,60,kTRUE);
  	frame->Draw();
    
  	cToy3->cd(3);
  	RooPlot* frame = mgrBin.plotPull(*nbkgBB,-3.0,3.0,60,kTRUE);
  	frame->Draw();
	
        cToy3->cd(4);    
        RooPlot* frame = mgrBin.plotPull(*nbkgPro,-3.0,3.0,60,kTRUE);
        frame->Draw();
    
        char file2[150];
        sprintf(file,Form("%s%i%s","/pToys_pulls_ybin", i, ".pdf") );
        cout << "Saving toy plot for bin with file name " << file << endl;
        cToy1->SaveAs("JPsiKsFit_dataRes/"+trigger+file);
        cToy3->SaveAs("JPsiKsFit_dataRes/"+trigger+file);

      }

  outTextFile << " Done with bin " << i << endl;
  outTextFile << "===============================================" << endl;

    } //end loop over bins

outTextFile << " After y bins loop" << endl;

    cBinProj->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResults_data_ybins1.pdf");
    cBinProj2->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResults_data_ybins2.pdf");
    cBinProj3->SaveAs("JPsiKsFit_dataRes/"+trigger+"/JPsiKsFitResults_data_ybins3.pdf");
    MCPDFs->SaveAs("JPsiKsFit_dataRes/"+trigger+"/BMassPDFs_sigPeak_ybins.pdf");

outTextFile << " nBins = " << nBins << endl;

    for(int i=0;i<nBins;i++){
      outTextFile << "Y bin " << i << " yield: " << SigYield[i]  
      << " +/- " << ErrYield[i] << "\n";
    }
  
  } // closing if (doYBins) loop


} // closing if (fitData) loop

cout << "got to the end" << endl;

}
