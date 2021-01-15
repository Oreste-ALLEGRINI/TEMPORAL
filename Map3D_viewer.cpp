#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <math.h>
#include <vector>
#include <Riostream.h>

// Root
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TLatex.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TF2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TGraph.h"
#include "TMath.h"

namespace {
Double_t my_transfer_function(const Double_t *x, const Double_t * /*param*/)
{
   // Bin values in our example range from -2 to 1.
   // Let's make values from -2. to -1.5 more transparent:
    if (*x <= 1)
       return 0.00;

    if (*x > 1 && *x <= 5)
       return 0.05;

    if (*x < 5 && *x <= 10)
       return 0.1;

    if (*x < 10 && *x <= 20)
       return 0.5;

    if (*x > 20 && *x < 50)
       return 0.75;

    if (*x > 50)
       return 0.9;

    return 0.2;
}

}

void Map3D_viewer(std::string path)
{
   TFile *rootfile = new TFile(path.c_str(),"READ");
   TH3F *CeBr1 = (TH3F*)rootfile->Get("3Dmap_T2");
   TH3F *CeBr2 = (TH3F*)rootfile->Get("3Dmap_T4");
   gStyle->SetCanvasPreferGL(kTRUE);
   TCanvas *c        = new TCanvas("glbox","TH3 Drawing", 100, 10, 850, 400);
   TPaveLabel *title = new TPaveLabel(0.04, 0.86, 0.96, 0.98,
                           "3D Maps of CeBr1 and CeBr2.");
   title->SetFillColor(32);
   title->Draw();

   TPad *boxPad  = new TPad("box", "box", 0.02, 0.02, 0.48, 0.82);
   TPad *box1Pad = new TPad("box1", "box1", 0.52, 0.02, 0.98, 0.82);
   boxPad->Draw();
   box1Pad->Draw();

   TList * lf = CeBr1->GetListOfFunctions();
      if (lf) {
         TF1 * tf = new TF1("TransferFunction", my_transfer_function); //IMPORTANT, the name of TF1 object MUST be "TransferFunction" to be taken into account by "glcol"to set the transparency
         lf->Add(tf);
      }

   boxPad->cd();
   CeBr1->Draw("glcol");

   TList * lf2 = CeBr2->GetListOfFunctions();
   if (lf2) {
       TF1 * tf2 = new TF1("TransferFunction", my_transfer_function);
       lf2->Add(tf2);
   }

   box1Pad->cd();
   CeBr2->Draw("glcol");
}
