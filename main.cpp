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

    /** computation of temperature range filter*/
    char temper_filter_answer;
    bool temper_filtering;
    int temper_peak_Tile1, temper_peak_Tile2, temper_peak_Tile4;

    /** Output */
    TFile outroot("/Users/allegrini//Documents/11.CameraTemporal/Tests_DAMAVAN_10_2020/Analysis_output/test_2.root", "RECREATE");
    //TFile outroot(Form("%s.root", filename.erase(filename.end()-5, filename.end())), "RECREATE");

    /////////////////////// Histo/Graphs definition ////////////////////////

    /** Temporary graphs for raw data filtering **/
    //TH1F* T1_T2 = new TH1F("","",83886080,-41943040,41943040);

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
    TH1F* Photons_Spectrum_Tile1 = new TH1F("Photons_Spectrum_Tile1","Photons_Spectrum_Tile1", 5000, 0, 100000);
    TH1D* Photons_Spectrum_Tile2 = new TH1D("Photons_Spectrum_Tile2","Photons_Spectrum_Tile2", 5000, 0, 100000);
    TH1F* Photons_Spectrum_Tile4 = new TH1F("Photons_Spectrum_Tile4","Photons_Spectrum_Tile4", 5000, 0, 100000);
    TSpectrum *spectrum = new TSpectrum();

    TH1F* Coinc_Photons_Spectrum_Tile12 = new TH1F("Coinc_Photons_Spectrum_Tile1","Coinc_Photons_Spectrum_Tile1", 5000, 0, 100000);
    TH1F* Coinc_Photons_Spectrum_Tile14 = new TH1F("Coinc_Photons_Spectrum_Tile1","Coinc_Photons_Spectrum_Tile1", 5000, 0, 100000);
    TH1F* Coinc_Photons_Spectrum_Tile21 = new TH1F("Coinc_Photons_Spectrum_Tile2","Coinc_Photons_Spectrum_Tile2", 5000, 0, 100000);
    TH1F* Coinc_Photons_Spectrum_Tile24 = new TH1F("Coinc_Photons_Spectrum_Tile1","Coinc_Photons_Spectrum_Tile1", 5000, 0, 100000);
    TH1F* Coinc_Photons_Spectrum_Tile41 = new TH1F("Coinc_Photons_Spectrum_Tile2","Coinc_Photons_Spectrum_Tile4", 5000, 0, 100000);
    TH1F* Coinc_Photons_Spectrum_Tile42 = new TH1F("Coinc_Photons_Spectrum_Tile1","Coinc_Photons_Spectrum_Tile1", 5000, 0, 100000);

    /** Temperature spectrum */
    TH1F* Temper_Tile1 = new TH1F("Temper_Tile1","Temper_Tile1", 125, -2000, 500);
    TH1F* Temper_Tile2 = new TH1F("Temper_Tile2","Temper_Tile2", 125, -2000, 500);
    TH1F* Temper_Tile4 = new TH1F("Temper_Tile4","Temper_Tile4", 125, -2000, 500);

    /** 2D mapping of the position of the photons*/
    int DimXY=32;
    TH2D* mapT2coincT1_T2 = new TH2D("", "", DimXY, -0.5, 8.5, DimXY, 0.5, 8.5);
    TH2D* mapT4coincT1_T4 = new TH2D("", "", DimXY, -0.5, 8.5, DimXY, 0.5, 8.5);
    TH2D* mapT2coincT2_T4 = new TH2D("", "", DimXY, -0.5, 8.5, DimXY, 0.5, 8.5);
    TH2D* mapT4coincT2_T4 = new TH2D("", "", DimXY, -0.5, 8.5, DimXY, 0.5, 8.5);

    ///////////////////////   Canvas definition    ////////////////////////
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
        std::cout<<"Computing the reference energy peak position"<<std::endl;

        for(int i=0; i<Tile1.size(); i++){
            Photons_Spectrum_Tile1->Fill(Tile1.at(i).photons);
        }
        for(int i=0; i<Tile2.size(); i++){
            Photons_Spectrum_Tile2->Fill(Tile2.at(i).photons);
        }
        for(int i=0; i<Tile4.size(); i++){
            Photons_Spectrum_Tile4->Fill(Tile4.at(i).photons);
        }

        CRT_Canvas->cd();
        Photons_Spectrum_Tile2->Draw();

        //Int_t peak_found = spectrum->Search(Photons_Spectrum_Tile2,2,"",0.10);
        //std::cout<<peak_found<<std::endl;

        if((Photons_Spectrum_Tile1->GetBinCenter(Photons_Spectrum_Tile1->GetMaximumBin()))-((Photons_Spectrum_Tile1->GetBinWidth(Photons_Spectrum_Tile1->GetMaximumBin()))/2) == 0){
            energy_peak_Tile1 = -1;
        }
        else {energy_peak_Tile1 = Photons_Spectrum_Tile1->GetBinCenter(Photons_Spectrum_Tile1->GetMaximumBin());}

        if((Photons_Spectrum_Tile2->GetBinCenter(Photons_Spectrum_Tile2->GetMaximumBin()))-((Photons_Spectrum_Tile2->GetBinWidth(Photons_Spectrum_Tile2->GetMaximumBin()))/2) == 0){
            energy_peak_Tile2 = -1;
        }
        else {energy_peak_Tile2 = Photons_Spectrum_Tile2->GetBinCenter(Photons_Spectrum_Tile2->GetMaximumBin());}

        if((Photons_Spectrum_Tile4->GetBinCenter(Photons_Spectrum_Tile4->GetMaximumBin()))-((Photons_Spectrum_Tile4->GetBinWidth(Photons_Spectrum_Tile4->GetMaximumBin()))/2) == 0){
            energy_peak_Tile4 = -1;
        }
        else {energy_peak_Tile4 = -1/*Photons_Spectrum_Tile4->GetBinCenter(Photons_Spectrum_Tile4->GetMaximumBin())*/;}
        std::cout<<energy_peak_Tile1<<" "<<energy_peak_Tile2<<" "<<energy_peak_Tile4<<std::endl;
    }

    std::cout<< "Do you want to apply a temperature filter (Max peak value +/- 0.25 \370 C)? (Y/N)"<<std::endl;
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
        std::cout<<temper_peak_Tile1<<" "<<temper_peak_Tile2<<" "<<temper_peak_Tile4<<std::endl;
    }

    std::cout << "Sorting the data of each tiles by the timestamp ..." <<std::endl;

    std::sort(Tile1.begin(), Tile1.end(), sort_function());
    std::sort(Tile2.begin(), Tile2.end(), sort_function());
    std::sort(Tile4.begin(), Tile4.end(), sort_function());

    //////////////////////ONLY USED TO CHECK THE RAW DATA ///////////////////
    /*std::cout<< "Filtering of events under the time coincidence peak ..."<<std::endl;
    Coincidence Tile 1 and Tile 2
    record_CRT1_2 = Global_analysis (Tile1, Tile2);
     //Tile 1 - 2 ----> If you want to see the raw data, you just have to write the histogram in the rootfile
    for(int i=0; i<record_CRT1_2[0].size(); i++){
        T1_T2->Fill(record_CRT1_2[0].at(i));
    }*/

    /////////////////////////////////////////////////////////////////////////

    std::cout<< "Research of coincidence events ..."<<std::endl;
    std::cout<< "Peak location :"<<std::endl;

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



    std::cout<< "Building the 2D map ..."<<std::endl;

    /** Tile1*/
    mapT4coincT1_T4 = map2D (record_CRT1_4[6], record_CRT1_4[8], DimXY);
    /** Tile2*/
    mapT4coincT2_T4 = map2D (record_CRT2_4[6], record_CRT2_4[8], DimXY);
    /** Tile4*/
    mapT2coincT1_T2 = map2D (record_CRT1_2[6], record_CRT1_2[8], DimXY);

    mapT2coincT2_T4 = map2D (record_CRT2_4[5], record_CRT2_4[7], DimXY);



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

    std::cout<<"Writing the output root file ..."<<std::endl;
    outroot.cd();
    CRT_Tile1_Tile4->Write("CRT_Diamond_CeBr1");
    CRT_Tile1_Tile2->Write("CRT_Diamond_CeBr2");
    CRT_Tile2_Tile4->Write("CRT_CeBr2_CeBr1");
    //T1_T2->Write("Test");
    Photons_Spectrum_Tile1->Write("EnergySpectrumT1");
    Photons_Spectrum_Tile2->Write("EnergySpectrumT2");
    Photons_Spectrum_Tile4->Write("EnergySpectrumT4");
    Coinc_Photons_Spectrum_Tile21->Write("ESpectrumT2whenT1");
    Coinc_Photons_Spectrum_Tile41->Write("ESpectrumT4whenT1");
    Coinc_Photons_Spectrum_Tile24->Write("ESpectrumT2whenT4");
    Coinc_Photons_Spectrum_Tile42->Write("ESpectrumT4whenT2");
    mapT2coincT2_T4->Write("2Dmap_T2");
    mapT4coincT2_T4->Write("2Dmap_T4");
    Temper_Tile1->Write("TemperT1");
    Temper_Tile2->Write("TemperT2");
    Temper_Tile4->Write("TemperT4");

    return 0;
}

