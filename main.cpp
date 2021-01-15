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
#include <cmath>
#include <dirent.h> //class used to manipulate and traverse directories
#include <stdio.h>
#include <mutex>
#include <vector>
#include <array>
#include <algorithm> //class used for the max_element function

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
    TFile outroot("/Users/allegrini//Documents/11.CameraTemporal/Tests_DAMAVAN_10_2020/Analysis_output/test_2.root", "RECREATE");
    //TFile outroot(Form("%s.root", filename.erase(filename.end()-5, filename.end())), "RECREATE");

    /////////////////////// Histo/Graphs definition ////////////////////////

    /** Temporary graphs for raw data filtering **/
    TH1F* TA_TB = new TH1F("","",83886080,-41943040,41943040);

    /** Coincidence Time Resolution */
    TH1F* CRT_Tile1_Tile4 = new TH1F("CRT_Diamond_CeBr1","CRT_Diamond_CeBr1", 200, 100, 120);
    CRT_Tile1_Tile4->GetXaxis()->SetTitle("Time (ns)");
    CRT_Tile1_Tile4->GetYaxis()->SetTitle("Counts");
    TH1F* CRT_Tile1_Tile2 = new TH1F("CRT_Diamond_CeBr2","CRT_Diamond_CeBr2", 200, 100, 120);
    CRT_Tile1_Tile2->GetXaxis()->SetTitle("Time (ns)");
    CRT_Tile1_Tile2->GetYaxis()->SetTitle("Counts");
    TH1F* CRT_Tile2_Tile4 = new TH1F("CRT_CeBr2_CeBr1","CRT_CeBr2_CeBr1", 200, -10, 10);
    CRT_Tile2_Tile4->GetXaxis()->SetTitle("Time (ns)");
    CRT_Tile2_Tile4->GetYaxis()->SetTitle("Counts");

    /** Photons spectrum */
    TH1F* Photons_Spectrum_Tile1 = new TH1F("Photons_Spectrum_Tile1","Photons_Spectrum_Tile1", 1000, 0, 20000);
    TH1D* Photons_Spectrum_Tile2 = new TH1D("Photons_Spectrum_Tile2","Photons_Spectrum_Tile2", 1000, 0, 20000);
    TH1F* Photons_Spectrum_Tile4 = new TH1F("Photons_Spectrum_Tile4","Photons_Spectrum_Tile4", 1000, 0, 20000);
    TSpectrum *spectrum = new TSpectrum();

    TH1F* Coinc_Photons_Spectrum_Tile12 = new TH1F("Coinc_Photons_Spectrum_Tile1","Coinc_Photons_Spectrum_Tile1", 1000, 0, 20000);
    TH1F* Coinc_Photons_Spectrum_Tile14 = new TH1F("Coinc_Photons_Spectrum_Tile1","Coinc_Photons_Spectrum_Tile1", 1000, 0, 20000);
    TH1F* Coinc_Photons_Spectrum_Tile21 = new TH1F("Coinc_Photons_Spectrum_Tile2","Coinc_Photons_Spectrum_Tile2", 1000, 0, 20000);
    TH1F* Coinc_Photons_Spectrum_Tile24 = new TH1F("Coinc_Photons_Spectrum_Tile2","Coinc_Photons_Spectrum_Tile2", 1000, 0, 20000);
    TH1F* Coinc_Photons_Spectrum_Tile41 = new TH1F("Coinc_Photons_Spectrum_Tile4","Coinc_Photons_Spectrum_Tile4", 1000, 0, 20000);
    TH1F* Coinc_Photons_Spectrum_Tile42 = new TH1F("Coinc_Photons_Spectrum_Tile4","Coinc_Photons_Spectrum_Tile4", 1000, 0, 20000);

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

    if(reader.openReadFile(Qfilename) == FileReader::OPENNING_ERROR){
        std::cout << "Recherche du fichier: "<< filename <<" pas de fichier" << std::endl;
        return -1;
    }

    std::cout << "**** START ANALYSIS ****" << std::endl;

    std::cout << "Recovering data of each tiles ..." <<std::endl;
    while(reader.getNextEvent(&readTile) == FileReader::OK){

        if(readTile.tuile.tuile==1){
            Tile1.push_back(readTile.tuile);
            Tile1_Img.push_back(readTile.image);
        }

        else if(readTile.tuile.tuile==2){
            Tile2.push_back(readTile.tuile);
            Tile2_Img.push_back(readTile.image);
        }

        else if(readTile.tuile.tuile==4 || readTile.tuile.tuile==3){
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

        std::cout<<"Choose the position of the energy filter window for Tile 1. Press -1 if you don't want apply filter or if the Tile1 is the dead Tile"<<std::endl;
        std::cin>>energy_peak_Tile1;
        std::cout<<"For Tile 2 ?"<<std::endl;
        std::cin>>energy_peak_Tile2;
        std::cout<<"For Tile 3 or 4 ?"<<std::endl;
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
    }

    std::cout << "Sorting the data of each tiles by the timestamp ..." <<std::endl;

    std::sort(Tile1.begin(), Tile1.end(), sort_function());
    std::sort(Tile2.begin(), Tile2.end(), sort_function());
    std::sort(Tile4.begin(), Tile4.end(), sort_function());

    //////////////////////ONLY USED TO CHECK THE RAW DATA ///////////////////
    /*std::cout<< "Generating raw data CRT between T2 and T4 ..."<<std::endl;
    record_CRTA_B = Global_analysis (Tile2, Tile4);
     //Tile 1 - 2 ----> If you want to see the raw data, you just have to write the histogram in the rootfile
    for(int i=0; i<record_CRTA_B[0].size(); i++){
        TA_TB->Fill(record_CRTA_B[0].at(i));
    }*/

    /////////////////////////////////////////////////////////////////////////

    std::cout<< "Research of coincidence events ..."<<std::endl;

    if((energy_filtering==false) && (temper_filtering==false)){
        /** Coincidence Tile 1 and Tile 4 */
        record_CRT1_4 = CRT_filter (Tile1, Tile4);

        /** Coincidence Tile 1 and Tile 2 */
        record_CRT1_2 = CRT_filter (Tile1, Tile2);

        /** Coincidence Tile 2 and Tile 4 */
        record_CRT2_4 = CRT_filter (Tile2, Tile4);

        /** Coincidence between Tile 2 and 4 with trigger signal provided by Tile 1 */
        //////////TODO///////////
    }
    else{
        /** Coincidence Tile 1 and Tile 4 */
        record_CRT1_4 = Nrj_Temper_filter (energy_filtering, temper_filtering, energy_peak_Tile1, energy_peak_Tile4, temper_peak_Tile1, temper_peak_Tile4, Tile1, Tile4);

        /** Coincidence Tile 1 and Tile 2 */
        record_CRT1_2 = Nrj_Temper_filter (energy_filtering, temper_filtering, energy_peak_Tile1, energy_peak_Tile2, temper_peak_Tile1, temper_peak_Tile2, Tile1, Tile2);

        /** Coincidence Tile 2 and Tile 4 */
        record_CRT2_4 = Nrj_Temper_filter (energy_filtering, temper_filtering, energy_peak_Tile2, energy_peak_Tile4, temper_peak_Tile2, temper_peak_Tile4, Tile2, Tile4);

        /** Coincidence between Tile 2 and 4 with trigger signal provided by Tile 1 */
        //////////TODO///////////
    }

    std::cout<< "Building the depth profiles ..."<<std::endl;
    std::cout<< "Choose the axis corresponding to the depth of the Tile2 detector (X,Y or Z) (see the orientation of the detector with respect to the source)"<<std::endl;
    std::cout<< "    y\n    |\n    |\n    |________z           top surface of the CeBr detectors : XY plane \n   /\n  /\n /x"<<std::endl;
    std::cin>>axis;
    DepthProfileT2coincT2_T4 = DepthProfile(record_CRT2_4[5], record_CRT2_4[7], record_CRT2_4[9], axis);
    std::cout<< "Choose the axis corresponding to the depth of the Tile4 detector"<<std::endl;
    std::cin>>axis;
    DepthProfileT4coincT2_T4 = DepthProfile(record_CRT2_4[6], record_CRT2_4[8], record_CRT2_4[10], axis);

    std::cout<< "Building the 2D map ..."<<std::endl;

    /** Tile4*/
    mapT4coincT1_T4 = map2D (record_CRT1_4[6], record_CRT1_4[8], DimXY);

    mapT4coincT2_T4 = map2D (record_CRT2_4[6], record_CRT2_4[8], DimXY);
    /** Tile2*/
    mapT2coincT1_T2 = map2D (record_CRT1_2[6], record_CRT1_2[8], DimXY);

    mapT2coincT2_T4 = map2D (record_CRT2_4[5], record_CRT2_4[7], DimXY);

    std::cout<< "Building the 3D map ..."<<std::endl;

    /** Tile4*/
    map3DT4coincT1_T4 = map3D (record_CRT1_4[6], record_CRT1_4[8], record_CRT1_4[10], DimXY);

    map3DT4coincT2_T4 = map3D (record_CRT2_4[6], record_CRT2_4[8], record_CRT2_4[10], DimXY);
    /** Tile2*/
    map3DT2coincT1_T2 = map3D (record_CRT1_2[6], record_CRT1_2[8], record_CRT1_2[10], DimXY);

    map3DT2coincT2_T4 = map3D (record_CRT2_4[5], record_CRT2_4[7], record_CRT2_4[9], DimXY);

    std::cout << "Creation of graphs and histograms ..." <<std::endl;

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

    //std::cout<< Form("%s.root", filename.erase(filename.end()-5, filename.end()))<<std::endl;*/
    std::cout<<"Press any key to save the data and close the program"<<std::endl;
    system("read");

    std::cout<<"Writing the output root file ..."<<std::endl;
    outroot.cd();
    CRT_Tile1_Tile4->Write("CRT_Diamond_CeBr1");
    CRT_Tile1_Tile2->Write("CRT_Diamond_CeBr2");
    CRT_Tile2_Tile4->Write("CRT_CeBr2_CeBr1");
    //TA_TB->Write("Test");
    ESpectrum_Canvas->Write("EnergySpectrum");
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

    return 0;
}

