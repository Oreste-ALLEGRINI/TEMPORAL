/**
 * @author      : o.allegrini@ipnl.in2p3.fr
 * @file        : main
 * @created     : jeudi oct. 15, 2020 16:46:04 CEST
 * @comment     : This program analyses the data recorded in the CeBr3 detectors of the classic temporal camera and the data from the detector trigger (diamond hodoscope)
 */
// TEMPORAL CAMERA LIBRARY
#include "quadratorevt.h"
#include <filereader.h>

//FUNCTIONS
#include <functions.h>

//C++ LIBRARIES
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <cstring>
#include <cstdlib>
#include <cstddef>
#include <cmath>
#include <dirent.h> //class used to manipulate and traverse directories
#include <stdio.h>
#include <mutex>
#include <vector>
#include <array>
#include <algorithm> //class used for the max_element function
#include <functional>

//ROOT LIBRARIES
#include "TLatex.h"
#include "TROOT.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TFile.h"
#include "TApplication.h"
#include "TRint.h"
#include "TAxis.h"
#include "TAttLine.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TGWindow.h"
#include "TGClient.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TH3.h"

int main(int argc,char ** argv)
{
    /////////////////// Starting ROOT in the c++ program //////////////////
    TApplication *theApp = new TApplication ("App", &argc, argv);
    if (gROOT->IsBatch()) {
        fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
        return 1;
    }
    /////////////////////// Variables declaration ////////////////////////
    /** Debug/Test variables (probably not currently used)**/
    double test1 = 0, test2 = 0, test3 = 0;
    /** Reading file variables */
    FileReader reader;
    TuileDam readTile;
    int it_Tile1, it_Tile2, it_Tile4; //iterator used to browse the vectors for the coincidence research and the CRT assessment.
    bool real_value; //boolean variable to check if a CRT value has been recorded in CRT1_4 and CRT1_4 variables.
    std::string filename;
    //Loading the data file
    std::cout<<"Enter <path/filename.type>"<<std::endl;
    std::cin>>filename;
    QString Qfilename = QString::fromStdString(filename); //Filereader variables uses QString input ==> This line converts string to QString

    /** data recovering variables */
    //Tile 1
    //// @brief: TuileEvt vectors recover the whole data in a Tile. TuileImg vectors recovers the pixel position of photons of a TuileEvt
    std::vector<TuileEvt> Tile1;
    std::vector<TuileImg> Tile1_Img;

    //Tile 2
    std::vector<TuileEvt> Tile2;
    std::vector<TuileImg> Tile2_Img;

    //Tile 4
    std::vector<TuileEvt> Tile4;
    std::vector<TuileImg> Tile4_Img;

    //Computed CRT vectors
    std::array<std::vector<double>,11> record_CRT1_4;
    std::array<std::vector<double>,11> record_CRT1_2;
    std::array<std::vector<double>,11> record_CRT2_4;
    std::array<std::vector<double>,11> record_CRTA_B;

    //Energy vectors
    std::vector<int> record_energy_Tile1;
    std::vector<int> record_energy_Tile2;
    std::vector<int> record_energy_Tile4;

    /** computation of CRT peak*/
    double CRT_peak_T1_T4, CRT_peak_T1_T2, CRT_peak_T2_T4;

    /** computation of reference energy peak*/
    char energy_filter_answer;
    bool energy_filtering;
    int energy_peak_Tile1, energy_peak_Tile2, energy_peak_Tile4;
    int peak_found;
    double* peaks_position;

    /** computation of temperature range filter*/
    char temper_filter_answer;
    bool temper_filtering;
    int temper_peak_Tile1, temper_peak_Tile2, temper_peak_Tile4;

    /** Output */
    std::string outfolder = "/Analysis_output/";
    std::string fileout = SplitFilename((const std::string&)filename, outfolder);
    std::stringstream path;
    std::string extension = ".root";
    path.str("");
    path<<fileout<<extension;
    std::cout<< "The output file will be : "<<path.str().c_str()<<std::endl;
    TFile outroot(path.str().c_str(), "RECREATE");
    //TFile outroot(cwd, "RECREATE");

    /////////////////////// Histo/Graphs definition ////////////////////////

    /** Temporary graphs for raw data filtering **/
    TH1F* TA_TB = new TH1F("","",83886080,-41943040,41943040);

    /** Coincidence Time Resolution */

    gStyle->SetLabelFont(22,"X");
    gStyle->SetTitleFont(22,"X");
    gStyle->SetLabelFont(22,"Y");
    gStyle->SetTitleFont(22,"Y");
    gStyle->SetHistLineWidth(2);
    int mean_CRT_T1_T2;
    int mean_CRT_T1_T4;
    int mean_CRT_T2_T4;
    TH1F* CRT_Tile1_Tile4;
    TH1F* CRT_Tile1_Tile2;
    TH1F* CRT_Tile2_Tile4;

    /** Photons spectrum */
    TH1F* Photons_Spectrum_Tile1 = new TH1F("Photons_Spectrum_Tile1","Photons_Spectrum_Tile1", 1000, 0, 20000);
    Photons_Spectrum_Tile1->SetTitle(";Energy (a.u);Counts");
    TH1F* Photons_Spectrum_Tile2 = new TH1F("Photons_Spectrum_Tile2","Photons_Spectrum_Tile2", 1000, 0, 20000);
    Photons_Spectrum_Tile2->SetTitle(";Energy (a.u);Counts");
    TH1F* Photons_Spectrum_Tile4 = new TH1F("Photons_Spectrum_Tile4","Photons_Spectrum_Tile4", 1000, 0, 20000);
    Photons_Spectrum_Tile4->SetTitle(";Energy (a.u);Counts");
    Photons_Spectrum_Tile4->SetLineColor(kRed);
    TSpectrum *spectrum = new TSpectrum();

    TH1F* Coinc_Photons_Spectrum_Tile12 = new TH1F("Coinc_Photons_Spectrum_Tile1","Coinc_Photons_Spectrum_Tile1", 1000, 0, 20000);
    Coinc_Photons_Spectrum_Tile12->SetTitle(";Energy (a.u);Counts");
    TH1F* Coinc_Photons_Spectrum_Tile14 = new TH1F("Coinc_Photons_Spectrum_Tile1","Coinc_Photons_Spectrum_Tile1", 1000, 0, 20000);
    Coinc_Photons_Spectrum_Tile14->SetTitle(";Energy (a.u);Counts");
    TH1F* Coinc_Photons_Spectrum_Tile21 = new TH1F("Coinc_Photons_Spectrum_Tile2","Coinc_Photons_Spectrum_Tile2", 1000, 0, 20000);
    Coinc_Photons_Spectrum_Tile21->SetTitle(";Energy (a.u);Counts");
    TH1F* Coinc_Photons_Spectrum_Tile24 = new TH1F("Coinc_Photons_Spectrum_Tile2","Coinc_Photons_Spectrum_Tile2", 1000, 0, 20000);
    Coinc_Photons_Spectrum_Tile24->SetTitle(";Energy (a.u);Counts");
    TH1F* Coinc_Photons_Spectrum_Tile41 = new TH1F("Coinc_Photons_Spectrum_Tile4","Coinc_Photons_Spectrum_Tile4", 1000, 0, 20000);
    Coinc_Photons_Spectrum_Tile41->SetTitle(";Energy (a.u);Counts");
    TH1F* Coinc_Photons_Spectrum_Tile42 = new TH1F("Coinc_Photons_Spectrum_Tile4","Coinc_Photons_Spectrum_Tile4", 1000, 0, 20000);
    Coinc_Photons_Spectrum_Tile42->SetTitle(";Energy (a.u);Counts");

    TH1F* X = new TH1F("X","X", 1000, 0, 10000);
    X->SetTitle(";Energy (a.u);Counts");

    TH1F* Y = new TH1F("X","X", 1000, 0, 10000);
    Y->SetTitle(";Energy (a.u);Counts");

    TH1F* Z = new TH1F("X","X", 1000, 0, 10000);
    Z->SetTitle(";Energy (a.u);Counts");

    TH2D* PhotonSpectrumT2vsT4 = new TH2D("", "", 1200, 0, 12000, 1200, 0, 12000);

    /** Temperature spectrum */
    TH1F* Temper_Tile1 = new TH1F("Temper_Tile1","Temper_Tile1", 125, -2000, 500);
    TH1F* Temper_Tile2 = new TH1F("Temper_Tile2","Temper_Tile2", 125, -2000, 500);
    TH1F* Temper_Tile4 = new TH1F("Temper_Tile4","Temper_Tile4", 125, -2000, 500);

    /** 1D depth profiles in CeBr detectors*/
    char axis;
    TH1F* DepthProfileT2coincT1_T2;
    TH1F* DepthProfileT4coincT1_T4;
    TH1F* DepthProfileT2coincT2_T4;
    TH1F* DepthProfileT4coincT2_T4;

    /** 2D mapping of the position of the photons*/
    int DimXY=32;
    TH2D* mapT2coincT1_T2 = new TH2D("", "", DimXY, 0.5, 8.5, DimXY, 0.5, 8.5);
    TH2D* mapT4coincT1_T4 = new TH2D("", "", DimXY, 0.5, 8.5, DimXY, 0.5, 8.5);
    TH2D* mapT2coincT2_T4 = new TH2D("", "", DimXY, 0.5, 8.5, DimXY, 0.5, 8.5);
    TH2D* mapT4coincT2_T4 = new TH2D("", "", DimXY, 0.5, 8.5, DimXY, 0.5, 8.5);
    // Raw floodmap study without any filtering
    TH2D* mapTAcoincTA_TB = new TH2D("", "", DimXY, 0.5, 8.5, DimXY, 0.5, 8.5);
    TH2D* mapTBcoincTA_TB = new TH2D("", "", DimXY, 0.5, 8.5, DimXY, 0.5, 8.5);

    /** 3D mapping of the position of the photons*/
    TH3F* map3DT2coincT1_T2 = new TH3F("", "", DimXY, 0.5, 8.5, DimXY, 0.5, 8.5, 18, 0.5, 5);
    map3DT2coincT1_T2->GetXaxis()->SetTitle("X");
    map3DT2coincT1_T2->GetYaxis()->SetTitle("Y");
    map3DT2coincT1_T2->GetZaxis()->SetTitle("Z");
    TH3F* map3DT4coincT1_T4 = new TH3F("", "", DimXY, 0.5, 8.5, DimXY, 0.5, 8.5, 18, 0.5, 5);
    map3DT4coincT1_T4->GetXaxis()->SetTitle("X");
    map3DT4coincT1_T4->GetYaxis()->SetTitle("Y");
    map3DT4coincT1_T4->GetZaxis()->SetTitle("Z");
    TH3F* map3DT2coincT2_T4 = new TH3F("", "", DimXY, 0.5, 8.5, DimXY, 0.5, 8.5, 18, 0.5, 5);
    map3DT2coincT2_T4->GetXaxis()->SetTitle("X");
    map3DT2coincT2_T4->GetYaxis()->SetTitle("Y");
    map3DT2coincT2_T4->GetZaxis()->SetTitle("Z");
    TH3F* map3DT4coincT2_T4 = new TH3F("", "", DimXY, 0.5, 8.5, DimXY, 0.5, 8.5, 18, 0.5, 5);
    map3DT4coincT2_T4->GetXaxis()->SetTitle("X");
    map3DT4coincT2_T4->GetYaxis()->SetTitle("Y");
    map3DT4coincT2_T4->GetZaxis()->SetTitle("Z");

    ///////////////////////   Canvas and PAD definition    ////////////////////////
    TCanvas* ESpectrum_Canvas = new TCanvas("ESpectrum","ESpectrum", 1200, 1200);
    TCanvas* CRT_Canvas = new TCanvas("CRT_Canvas","CRT_Canvas", 1200, 1200);
    CRT_Canvas->Close();

    if(reader.openReadFile(Qfilename) == FileReader::OPENNING_ERROR){
        std::cout << "Recherche du fichier: "<< filename <<" pas de fichier" << std::endl;
        return -1;
    }

    std::cout << "**** START ANALYSIS ****" << std::endl;

    std::cout << "Recovering data of each tiles ..." <<std::endl;
    while(reader.getNextEvent(&readTile) == FileReader::OK){

        if(readTile.tuile.tuile==1 || readTile.tuile.tuile==3){
            Tile1.push_back(readTile.tuile);
            Tile1_Img.push_back(readTile.image);
        }

        else if(readTile.tuile.tuile==2){
            Tile2.push_back(readTile.tuile);
            Tile2_Img.push_back(readTile.image);
        }

        else if(readTile.tuile.tuile==4){
            Tile4.push_back(readTile.tuile);
            Tile4_Img.push_back(readTile.image);
        }
        else {std::cout<<"ERREUR D'IDENTIFICATION DE LA TUILE: "<<readTile.tuile.tuile<<std::endl; return 0;}
    }

    std::cout<<"TuileEvt in tile 1: "<<Tile1.size()<<std::endl;
    std::cout<<"TuileEvt in tile 2: "<<Tile2.size()<<std::endl;
    std::cout<<"TuileEvt in tile 4: "<<Tile4.size()<<std::endl;

    std::cout<< "Do you want to apply an energy filtering (window 20% centered) ? (Y/N)"<<std::endl;
    std::cin>> energy_filter_answer;
    energy_filter_answer = towlower (energy_filter_answer);
    switch (energy_filter_answer)
    {
    case 'y': case 'yes':
        energy_filtering = true;
        break;

    case 'n': case 'no':
        energy_filtering = false;
        break;
    default: std::cout<<"Energy filter has not been correctly initialized. Please check that your choice is Y or N."<<std::endl;
        break;
    }

    if (energy_filtering == true){
        std::cout<<"Research of the energy peak positions in the spectrums"<<std::endl;

        /** Filling, drawing and printing the energy spectrums **/
        for(int i=0; i<Tile1.size(); i++){
            Photons_Spectrum_Tile1->Fill(Tile1.at(i).photons);
        }
        for(int i=0; i<Tile2.size(); i++){
            Photons_Spectrum_Tile2->Fill(Tile2.at(i).photons);
        }
        for(int i=0; i<Tile4.size(); i++){
            Photons_Spectrum_Tile4->Fill(Tile4.at(i).photons);
        }

        ESpectrum_Canvas->Divide(2,2);
        ESpectrum_Canvas->cd(1);
        Photons_Spectrum_Tile1->Draw();
        ESpectrum_Canvas->cd(2);
        Photons_Spectrum_Tile2->Draw();
        ESpectrum_Canvas->cd(3);
        Photons_Spectrum_Tile4->Draw();

        /** User's choice of energy filter positions per detector **/
        if(Photons_Spectrum_Tile1->GetBinContent(1) == Photons_Spectrum_Tile1->GetEntries()){
            std::cout<<"\n ******** Energy data are not provided for the detector connected on the Tile 1 ******** \n"<<std::endl;
        }
        else{
            peak_found = spectrum->Search(Photons_Spectrum_Tile1,10,"",0.05); //Search the peaks in the spectrum and displays them in the histogram
            peaks_position = spectrum->GetPositionX();
            std::cout<<"List of peak positions for Tile 1:"<<std::endl;
            std::cout<<PrintEnergyPeaks(peaks_position,peak_found)<<std::endl;
        }

        if(Photons_Spectrum_Tile2->GetBinContent(1) == Photons_Spectrum_Tile2->GetEntries()){
            std::cout<<"\n ******** Energy data are not provided for the detector connected on the Tile 2 ******** \n"<<std::endl;
        }
        else{
            peak_found = spectrum->Search(Photons_Spectrum_Tile2,10,"",0.05); //Search the peaks in the spectrum and displays them in the histogram
            peaks_position = spectrum->GetPositionX();
            std::cout<<"List of peak positions for Tile 2:"<<std::endl;
            std::cout<<PrintEnergyPeaks(peaks_position,peak_found)<<std::endl;
        }

        if(Photons_Spectrum_Tile4->GetBinContent(1) == Photons_Spectrum_Tile4->GetEntries()){
            std::cout<<"\n ******** Energy data are not provided for the detector connected on the Tile 4 ******** \n"<<std::endl;
        }
        else{
            peak_found = spectrum->Search(Photons_Spectrum_Tile4,10,"",0.05); //Search the peaks in the spectrum and displays them in the histogram
            peaks_position = spectrum->GetPositionX();
            std::cout<<"List of peak positions for Tile 4:"<<std::endl;
            std::cout<<PrintEnergyPeaks(peaks_position,peak_found)<<std::endl;
        }

        ESpectrum_Canvas->Modified();
        ESpectrum_Canvas->Update();
        gSystem->ProcessEvents();

        std::cout<<"Choose the position of the energy filter window for Tile 1 or 3. Press -1 if you don't want apply filter or if the Tile1 is the dead Tile"<<std::endl;
        std::cin>>energy_peak_Tile1;
        std::cout<<"For Tile 2 ?"<<std::endl;
        std::cin>>energy_peak_Tile2;
        std::cout<<"For Tile 4 ?"<<std::endl;
        std::cin>>energy_peak_Tile4;

        std::cout<<energy_peak_Tile1<<" "<<energy_peak_Tile2<<" "<<energy_peak_Tile4<<std::endl;
    }

    std::cout<< "Do you want to apply a temperature filter (Max peak value +/- 0.25"<<"\u00B0"<<"C)? (Y/N)"<<std::endl;
    std::cin>>temper_filter_answer;
    temper_filter_answer = towlower (temper_filter_answer);
    switch (temper_filter_answer)
    {
    case 'y': case 'yes':
        temper_filtering = true;
        break;

    case 'n': case 'no':
        temper_filtering = false;
        break;
    default: std::cout<<"Energy filter has not been correctly initialized. Please check that your choice is Y or N."<<std::endl;
        break;
    }

    if(temper_filtering == true){
        std::cout<<"Computing the reference temperature value"<<std::endl;

        for(int i=0; i<Tile1.size(); i++){
            Temper_Tile1->Fill(Tile1.at(i).temper);
        }
        for(int i=0; i<Tile2.size(); i++){
            Temper_Tile2->Fill(Tile2.at(i).temper);
        }
        for(int i=0; i<Tile4.size(); i++){
            Temper_Tile4->Fill(Tile4.at(i).temper);
        }

        if(Temper_Tile1->GetMaximumBin() == 1){
            temper_peak_Tile1 = -1;
        }
        else {temper_peak_Tile1 = Temper_Tile1->GetBinCenter(Temper_Tile1->GetMaximumBin());}

        if(Temper_Tile2->GetMaximumBin() == 1){
            temper_peak_Tile2 = -1;
        }
        else {temper_peak_Tile2 = Temper_Tile2->GetBinCenter(Temper_Tile2->GetMaximumBin());}

        if(Temper_Tile4->GetMaximumBin() == 1){
            temper_peak_Tile4 = -1;
        }
        else {temper_peak_Tile4 = Temper_Tile4->GetBinCenter(Temper_Tile4->GetMaximumBin());}
        std::cout<<"Temperature peaks position on T1 T2 and T4 : "<<temper_peak_Tile1<<" "<<temper_peak_Tile2<<" "<<temper_peak_Tile4<<std::endl;
    }

    std::cout << "Sorting the data of each tiles by the timestamp ..." <<std::endl;

    std::sort(Tile1.begin(), Tile1.end(), sort_function());
    std::sort(Tile2.begin(), Tile2.end(), sort_function());
    std::sort(Tile4.begin(), Tile4.end(), sort_function());

    int compteur =0;
    for(int i=1; i<Tile2.size(); i++){
        if(Tile2.at(i).htimestamp > Tile2.at(i-1).htimestamp){
            compteur++;
        }
    }

    std::cout<<"Nb d'evenement sans frame multiples sur T2 : "<<compteur<<std::endl;

    //////////////////////ONLY USED TO CHECK THE RAW DATA AND RAW 2D MAP ///////////////////
    std::cout<< "Generating raw data CRT between T2 and T4 ..."<<std::endl;
    if(Tile1.size() != 0 && Tile2.size() !=0){
        record_CRTA_B = Global_analysis_bis (Tile1, Tile2, Tile1_Img, Tile2_Img);
         //Tile 1 - 2 ----> If you want to see the raw data, you just have to write the histogram in the rootfile
        for(int i=0; i<record_CRTA_B[0].size(); i++){
            TA_TB->Fill(record_CRTA_B[0].at(i));
        }
        /*for(int i=0; i<record_CRTA_B[0].size(); i++){
            std::cout<<record_CRTA_B[0].at(i)<<std::endl;
        }*/
        if((Tile1.empty()!=true) && (Tile2.empty()!=true)){
            mapTAcoincTA_TB = map2D (record_CRTA_B[5], record_CRTA_B[7], DimXY);
            mapTBcoincTA_TB = map2D (record_CRTA_B[6], record_CRTA_B[8], DimXY);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "Research of coincidence events ..."<<std::endl;

    if((energy_filtering==false) && (temper_filtering==false)){
        /** Coincidence Tile 1 and Tile 4 */
        //if((Tile1.empty()!=true) && (Tile4.empty()!=true)){
        //record_CRT1_4 = CRT_filter (Tile1, Tile4, Tile1_Img, Tile4_Img);
        //}
        //else{std::cout<<"Tile 1 or 4 has not collected data. The reasearch of coincidence is then impossible !\n"<<std::endl;}

        /** Coincidence Tile 1 and Tile 2 */
        //if((Tile1.empty()!=true) && (Tile2.empty()!=true)){
        //record_CRT1_2 = CRT_filter (Tile1, Tile2, Tile1_Img, Tile2_Img);
        //}
        //else{std::cout<<"Tile 1 or 2 has not collected data. The reasearch of coincidence is then impossible !\n"<<std::endl;}

        /** Coincidence Tile 2 and Tile 4 */
        if((Tile2.empty()!=true) && (Tile4.empty()!=true)){
        record_CRT2_4 = CRT_filter (Tile2, Tile4, Tile2_Img, Tile4_Img);
        }
        else{std::cout<<"Tile 2 or 4 has not collected data. The reasearch of coincidence is then impossible !\n"<<std::endl;}

        /** Coincidence between Tile 2 and 4 with trigger signal provided by Tile 1 */
        //////////TODO///////////
    }
    else{
        /** Coincidence Tile 1 and Tile 4 */
        if((Tile1.empty()!=true) && (Tile4.empty()!=true)){
        record_CRT1_4 = Nrj_Temper_filter (energy_filtering, temper_filtering, energy_peak_Tile1, energy_peak_Tile4, temper_peak_Tile1, temper_peak_Tile4, Tile1, Tile4, Tile1_Img, Tile4_Img);
        }

        /** Coincidence Tile 1 and Tile 2 */
        if((Tile1.empty()!=true) && (Tile2.empty()!=true)){
        record_CRT1_2 = Nrj_Temper_filter (energy_filtering, temper_filtering, energy_peak_Tile1, energy_peak_Tile2, temper_peak_Tile1, temper_peak_Tile2, Tile1, Tile2, Tile1_Img, Tile2_Img);
        }

        /** Coincidence Tile 2 and Tile 4 */
        if((Tile2.empty()!=true) && (Tile4.empty()!=true)){
        record_CRT2_4 = Nrj_Temper_filter (energy_filtering, temper_filtering, energy_peak_Tile2, energy_peak_Tile4, temper_peak_Tile2, temper_peak_Tile4, Tile2, Tile4, Tile2_Img, Tile4_Img);
        }

        /** Coincidence between Tile 2 and 4 with trigger signal provided by Tile 1 */
        //////////TODO///////////
    }

    std::cout<< "Building the depth profiles ..."<<std::endl;
    std::cout<< "Choose the axis corresponding to the depth of the Tile2 detector (X,Y or Z) (see the orientation of the detector with respect to the source)"<<std::endl;
    std::cout<< "    y\n    |\n    |\n    |________z           top surface of the CeBr detectors : XY plane \n   /\n  /\n /x"<<std::endl;
    std::cin>>axis;
    if((Tile2.empty()!=true) && (Tile4.empty()!=true)){
    DepthProfileT2coincT2_T4 = DepthProfile(record_CRT2_4[5], record_CRT2_4[7], record_CRT2_4[9], axis);
    }
    std::cout<< "Choose the axis corresponding to the depth of the Tile4 detector"<<std::endl;
    std::cin>>axis;
    if((Tile2.empty()!=true) && (Tile4.empty()!=true)){
    DepthProfileT4coincT2_T4 = DepthProfile(record_CRT2_4[6], record_CRT2_4[8], record_CRT2_4[10], axis);
    }

    std::cout<< "Building the 2D map ..."<<std::endl;

    /** Tile4*/
    if((Tile1.empty()!=true) && (Tile4.empty()!=true)){
    mapT4coincT1_T4 = map2D (record_CRT1_4[6], record_CRT1_4[8], DimXY);
    }

    if((Tile2.empty()!=true) && (Tile4.empty()!=true)){
    mapT4coincT2_T4 = map2D (record_CRT2_4[6], record_CRT2_4[8], DimXY);
    }

    /** Tile2*/
    if((Tile1.empty()!=true) && (Tile2.empty()!=true)){
    mapT2coincT1_T2 = map2D (record_CRT1_2[6], record_CRT1_2[8], DimXY);
    }

    if((Tile2.empty()!=true) && (Tile4.empty()!=true)){
    mapT2coincT2_T4 = map2D (record_CRT2_4[5], record_CRT2_4[7], DimXY);
    }

    std::cout<< "Building the 3D map ..."<<std::endl;

    /** Tile4*/
    if((Tile1.empty()!=true) && (Tile4.empty()!=true)){
    map3DT4coincT1_T4 = map3D (record_CRT1_4[6], record_CRT1_4[8], record_CRT1_4[10], DimXY);
    }

    if((Tile2.empty()!=true) && (Tile4.empty()!=true)){
    map3DT4coincT2_T4 = map3D (record_CRT2_4[6], record_CRT2_4[8], record_CRT2_4[10], DimXY);
    }
    /** Tile2*/
    if((Tile1.empty()!=true) && (Tile2.empty()!=true)){
    map3DT2coincT1_T2 = map3D (record_CRT1_2[6], record_CRT1_2[8], record_CRT1_2[10], DimXY);
    }

    if((Tile2.empty()!=true) && (Tile4.empty()!=true)){
    map3DT2coincT2_T4 = map3D (record_CRT2_4[5], record_CRT2_4[7], record_CRT2_4[9], DimXY);
    }

    std::cout << "Creation of graphs and histograms ..." <<std::endl;

    if(record_CRT1_2[0].size()>0){
        mean_CRT_T1_T2 = std::accumulate(std::begin(record_CRT1_2[0]), std::end(record_CRT1_2[0]),0)/record_CRT1_2[0].size();
    }
    else{
        mean_CRT_T1_T2 = 0;
    }
    if(record_CRT1_4[0].size()>0){
        mean_CRT_T1_T4 = std::accumulate(std::begin(record_CRT1_4[0]), std::end(record_CRT1_4[0]),0)/record_CRT1_4[0].size();
    }
    else{
        mean_CRT_T1_T4 = 0;
    }

    if(record_CRT2_4[0].size()>0){
        mean_CRT_T2_T4 = (int)FindMedian(record_CRT2_4[0]);
    }
    else{
        mean_CRT_T2_T4 = 0;
    }

    std::cout<<mean_CRT_T1_T2<<" "<<mean_CRT_T1_T4<<" "<<mean_CRT_T2_T4<<std::endl;

    CRT_Tile1_Tile2= new TH1F("","", 400, mean_CRT_T1_T2-20, mean_CRT_T1_T2+20);
    CRT_Tile1_Tile2->SetTitle(";Time difference (ns);Counts");
    CRT_Tile1_Tile4= new TH1F("","", 400, mean_CRT_T1_T4-20, mean_CRT_T1_T4+20);
    CRT_Tile1_Tile4->SetTitle(";Time difference (ns);Counts");
    CRT_Tile2_Tile4= new TH1F("","", 400, mean_CRT_T2_T4-20, mean_CRT_T2_T4+20);
    CRT_Tile2_Tile4->SetTitle(";Time difference (ns);Counts");

    if((Tile2.empty()!=true) && (Tile4.empty()!=true)){
    PhotonSpectrumT2vsT4 = E_TileA_vs_E_TileB (record_CRT2_4[1], record_CRT2_4[2]);
    }

    for(int i=0; i<record_CRT1_4[0].size(); i++){
        CRT_Tile1_Tile4->Fill(record_CRT1_4[0].at(i));
    }

    for(int i=0; i<record_CRT1_4[1].size(); i++){
        Coinc_Photons_Spectrum_Tile14->Fill(record_CRT1_4[1].at(i));
    }

    for(int i=0; i<record_CRT1_4[2].size(); i++){
        Coinc_Photons_Spectrum_Tile41->Fill(record_CRT1_4[2].at(i));
    }

    for(int i=0; i<record_CRT1_2[0].size(); i++){
        CRT_Tile1_Tile2->Fill(record_CRT1_2[0].at(i));
    }

    for(int i=0; i<record_CRT1_2[1].size(); i++){
        Coinc_Photons_Spectrum_Tile12->Fill(record_CRT1_2[1].at(i));
    }

    for(int i=0; i<record_CRT1_2[2].size(); i++){
        Coinc_Photons_Spectrum_Tile21->Fill(record_CRT1_2[2].at(i));
    }

    for(int i=0; i<record_CRT2_4[0].size(); i++){
        CRT_Tile2_Tile4->Fill(record_CRT2_4[0].at(i));
    }

    for(int i=0; i<record_CRT2_4[1].size(); i++){
        Coinc_Photons_Spectrum_Tile24->Fill(record_CRT2_4[1].at(i));
    }

    for(int i=0; i<record_CRT2_4[2].size(); i++){
        Coinc_Photons_Spectrum_Tile42->Fill(record_CRT2_4[2].at(i));
    }

    /** This part of program is dedicated to recover the raw data positions to check the good match between the values and the data format structure*/
    /*for(int i=0; i<record_CRT2_4[8].size(); i++){
        X->Fill(record_CRT2_4[8].at(i));
    }
    for(int i=0; i<record_CRT2_4[9].size(); i++){
        Y->Fill(record_CRT2_4[9].at(i));
    }

    for(int i=0; i<record_CRT2_4[10].size(); i++){
        Z->Fill(record_CRT2_4[10].at(i));
    }*/
    /** End of addtionnal part of program **/

    //std::cout<< Form("%s.root", filename.erase(filename.end()-5, filename.end()))<<std::endl;*/
    std::cout<<"Press any key to save the data and close the program"<<std::endl;
    system("read");

    std::cout<<"Writing the output root file ..."<<std::endl;
    outroot.cd();
    CRT_Tile1_Tile4->Write("CRT_Diamond_CeBr1");
    CRT_Tile1_Tile2->Write("CRT_Diamond_CeBr2");
    CRT_Tile2_Tile4->Write("CRT_CeBr2_CeBr1");
    TA_TB->Write("Test");
    //mapTAcoincTA_TB->Write("Raw_Floodmap_T2");
    //mapTBcoincTA_TB->Write("Raw_Floodmap_T4");
    ESpectrum_Canvas->Write("EnergySpectrum");
    PhotonSpectrumT2vsT4->Write("E_T2_vs_E_T4");
    Coinc_Photons_Spectrum_Tile21->Write("ESpectrumT2whenT1");
    Coinc_Photons_Spectrum_Tile41->Write("ESpectrumT4whenT1");
    Coinc_Photons_Spectrum_Tile24->Write("ESpectrumT2whenT4");
    Coinc_Photons_Spectrum_Tile42->Write("ESpectrumT4whenT2");
    DepthProfileT2coincT2_T4->Write("Depth profile T2");
    DepthProfileT4coincT2_T4->Write("Depth profile T4");
    mapT2coincT2_T4->Write("2Dmap_T2");
    mapT4coincT2_T4->Write("2Dmap_T4");
    map3DT2coincT2_T4->Write("3Dmap_T2");
    map3DT4coincT2_T4->Write("3Dmap_T4");
    Temper_Tile1->Write("TemperT1");
    Temper_Tile2->Write("TemperT2");
    Temper_Tile4->Write("TemperT4");
    //X->Write("PosX");
    //Y->Write("PosY");
    //Z->Write("PosZ");

    return 0;
}

