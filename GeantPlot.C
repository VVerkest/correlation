void GeantPlot(){
  const float pi = 3.141592;

  TString exportName = "geant";
  TString fileSaveName = "FIT/";
  fileSaveName += exportName;
  fileSaveName += ".root";
  TString canvasTitle = "Geant : ";
  
  // create a new Root file
  TFile *top = new TFile(fileSaveName,"recreate");
 
  // create a subdirectory "tof" in this file
  TDirectory *cdtof = top->mkdir("tof");
  cdtof->cd();    // make the "tof" directory the current directory
  
// create a new subdirectory for each plane
 const Int_t nfiles = 4;
 const Int_t nPtBins = 5;
 double ptBinLo[nPtBins] = { 3, 5, 9, 13, 17 };
 double ptBinHi[nPtBins] = { 4, 8, 12, 16, 24 };
 TString ptBinString[nPtBins] = { "0.5-1.0", "1.0-2.0", "2.0-3.0", "3.0-4.0", "4.0-6.0" };
 TString minSubPt[nfiles] = {"7",  "8",  "9", "10"};
 TString minLeadPt[nfiles] = {"14", "16", "18", "20"};
 TString CName[nPtBins] = {"assoc05to10", "assoc10to20", "assoc20to30", "assoc30to40", "assoc40to60"};
 TString CTitle[nPtBins] = {"0.5 < Pt_{assoc} < 1.0","1.0 < Pt_{assoc} < 2.0","2.0 < Pt_{assoc} < 3.0","3.0 < Pt_{assoc} < 4.0","4.0 < Pt_{assoc} < 6.0"};
 int color[nfiles] = { 2, 51, 4, 8 };
 Int_t g,h,i,j;
 TDirectory *cdplane[nfiles];
 TH1D *assocPtCorr[nfiles][nPtBins];
 // TCanvas *PtCanvas;
 TH1D* PtPlot[nfiles][nPtBins];
 TF1 *ppfit[nPtBins];
 
 for (i=0;i<nfiles;i++) {
   TString l = minLeadPt[i];
   TString s = minSubPt[i];
   TString directory = "lead";
   directory += l;
   directory += "_sub";
   directory += s;
   cdplane[i] = cdtof->mkdir(directory);
   cdplane[i]->cd();

   for (j=0;j<nPtBins;j++) {
     TString f = ptBinString[j];
     TString PtRange = " min ";
     PtRange += f;
     PtRange += " Pt_assoc";
     
     TString hname = "L_";
     hname += directory;
     hname += "_bin_";
     hname += j;
     TString htitle = "Lead Corr: ";
     htitle += directory;
     htitle += PtRange;
     // hn[i][j] = new TH1F(hname,htitle,100,0,100);

     TString hname = "S_";
     hname += directory;
     hname += "_bin_";
     hname += j;
     TString htitle = "Sub Corr: ";
     htitle += directory;
     htitle += PtRange;
     // hs[i][j] = new TH1F(hname,htitle,100,0,100);
   }
   cdtof->cd();    // change current directory to top
 }

 for (i=0;i<nfiles;i++) {    // import histograms here!
   
   for (j=0;j<nPtBins;j++) {
     l = minLeadPt[i];
     s = minSubPt[i];
     TString importName = "out/geant";
     importName += "_lead";
     importName += l;
     importName += "_sub";
     importName += s;
     importName += ".root";
     TFile* ppdijetFILE = new TFile( importName, "READ" );
     TH3D* leadJetCorr = (TH3D**) ppdijetFILE->Get("ppleadjetcorr");
     TH3D* subJetCorr = (TH3D**) ppdijetFILE->Get("ppsubjetcorr");
     TH3D* leadJetPt = (TH3D*) ppdijetFILE->Get("ppleadjetpt");
     TH2D* ppdijetEvents = (TH2D*) ppdijetFILE->Get("binvzdist");
     TH2D* pplead = (TH2D*) leadJetCorr->Project3D("ZY");
     TH2D* ppsub = (TH2D*) subJetCorr->Project3D("ZY");
     TString nameSet = "leadCorr_";        // NAME HISTOGRAM
     nameSet += l;
     nameSet += "_";
     nameSet += s;
     pplead->SetName( nameSet );
     nameSet = "subCorr_";        // NAME HISTOGRAM
     nameSet += l;
     nameSet += "_";
     nameSet += s;
     ppsub->SetName(nameSet);
     nameSet = "leadPt_";        // NAME HISTOGRAM
     nameSet += l;
     nameSet += "_";
     nameSet += s;
     leadJetPt->SetName(nameSet);
     
     cdtof->cd();    // change current directory to top
     cdplane[i]->cd();   // ENTER DIRECTORY TO SAVE CORRESPONDING HISTOGRAMS
     
     nameSet = "leadCorr_";
     nameSet += l;
     nameSet += "_";
     nameSet += s;
     nameSet += "__bin_";
     nameSet += j;

     // PROJECT
     assocPtCorr[i][j] = (TH1D*) pplead->ProjectionX( nameSet , ptBinLo[j], ptBinHi[j] );
     assocPtCorr[i][j]->Scale(1/assocPtCorr[i][j]->GetBinWidth(1));
     assocPtCorr[i][j]->Scale( 1/double(ppdijetEvents->Integral()) );
     
     // FIT
     TString fiteq = "[0]+[1]*exp(-0.5*((x-[2])/[3])**2)+[4]*exp(-0.5*((x-[5])/[6])**2)";
     double phiMin = -pi + pi/2.0;
     double phiMax = pi + pi/2.0;
     nameSet = "corrFit_";
     nameSet += l;
     nameSet += "_";
     nameSet += s;
     nameSet += "__bin_";
     nameSet += j;
     
     ppfit[j] = new TF1( nameSet, fiteq, phiMin, phiMax);
     ppfit[j]->FixParameter(2, 0);
     ppfit[j]->FixParameter(5, pi);
     ppfit[j]->SetParameter(3, 0.2);
     ppfit[j]->SetParameter(6, 0.2);
     ppfit[j]->SetLineColor(color[i]);
     ppfit[j]->SetLineWidth(2);
     assocPtCorr[i][j]->Fit(ppfit[j]);


     gROOT->SetEditHistograms();
     assocPtCorr[i][j]->SetLineColor(color[i]);
     assocPtCorr[i][j]->SetLineWidth(2);
     //assocPtCorr[i][j]->SetMaximum(2.5);
     assocPtCorr[i][j]->SetMinimum(0);
     gStyle->SetOptStat(0);
     
     // WRITE

     assocPtCorr[i][j]->Write();
     pplead->Write();
     ppsub->Write();
   }
     ppdijetEvents->Write();
     leadJetPt->Write();
 }
 
 cdtof->cd();
 
 for (j=0; j<nPtBins; j++){
   i=0;
   canvasTitle = "Geant : ";
   canvasTitle += CTitle[j];
   TString canvasName = exportName;
   canvasName += "_";
   canvasName += CName[j];   
   TCanvas *PtCanvas = new TCanvas( canvasName , CTitle[j] ,0 ,23 ,1280 ,709 );
   assocPtCorr[i][j]->Draw();
   assocPtCorr[i][j]->SetTitle( canvasTitle );
   i+=1;
   for (i=1; i<nfiles; i++){
     TString SaveName = CName[j];
     SaveName += ".png";
     assocPtCorr[i][j]->Draw("SAME");
   }
   cdtof->cd();
   PtCanvas->Write();
 }
     
     


 cdtof->cd(); 
 delete top;
}
