/**
 * @author      : b.hjira@ip2i.in2p3.fr
 * @file        : main
 * @created     : Monday, February 15th, 2021 16:46:04 CEST
 * @comment     : Ce programme analyse les données enregistrées dans les détecteurs CeBr3 de la caméra temporelle classique et sortie une carte de correction
 */
// TEMPORAL CAMERA LIBRARY
#include "quadratorevt.h"
#include <filereader.h>

//FUNCTIONS
#include <functions_calib.h>

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
#include <dirent.h> //classe utilisée pour manipuler et parcourir les répertoires
#include <stdio.h>
#include <mutex>
#include <vector>
#include <array>
#include <algorithm> //classe utilisée pour la fonction max_element
//#include <bits/stdc++.h>
#include <chrono>

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

using namespace std;

int main(int argc, char **argv)
{
    /////////////////// Starting ROOT in the c++ program /////////////////
    TApplication *theApp = new TApplication("App", &argc, argv);
    if (gROOT->IsBatch())
    {
        fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
        return 1;
    }
    ///////////////////////////save the outputs////////////////////////
    std::string filename;
    string nomFichier("/home/hjira/calibration/output/calibration.root");
    TFile monFlux(nomFichier.c_str(), "RECREATE");
    ///////////////////////////variable declaration////////////////////
    int n, ki, kf;// use in for loop==> to choose all acquisitions or one
    int choice_tile = 2;
    int l = 0;// use in for loop==> to make the sum of number_events_pixel 
    int number_pixels_axis;
    double real_time_acquisition = 3600;
      

    /** computation of reference energy peak*/
    int peak_found = 0;
    double *peaks_position;
    int energy_peak_Tile;

    std::vector<TuileEvt> Tile;
    std::vector<TuileImg> Tile_Img;

    cout << "Enter the number of pixels per axis" << endl;
    cin >> number_pixels_axis;
    int number_pixels = pow(number_pixels_axis, 2);
    cout << "*** Number_pixels=***" << number_pixels << endl;
    cout << "*** Enter 0: for all acquisitions or a number (1--> 20) for each acquisition***" << endl;
    cin >> n;
    if (n == 0)
    {
        ki = 5, kf = 25;
    }
    else
    {
        ki = n, kf = n;
    }

    std::cout << "Choose number of tile : " <<std::endl;
    cin>> choice_tile;
    int acquisition_number = kf - ki + 1;
    vector<double> number_events_pixel_k(number_pixels * acquisition_number, 0.0);
    vector<double> total_event(number_pixels, 0.0);
    vector<double> number_events_pixel;
    double maxElement;
    //TCanvas *ESpectrum_Canvas = new TCanvas("ESpectrum", "ESpectrum", 1200, 1200);
    TCanvas* Map_tuile2_event= new TCanvas("Map_tuile2_photons","Map_tuile2_photons", 1600, 1000);
    TH2D *Map2D = new TH2D("Map2D", "Map2D", number_pixels_axis, 0.5, (number_pixels_axis + 0.5), number_pixels_axis, 0.5, (number_pixels_axis + 0.5));

    TH1F *Photons_Spectrum_Tile2 = new TH1F("Photons_Spectrum_Tile2", "Photons_Spectrum_Tile2", 1000, 0, 20000);
    Photons_Spectrum_Tile2->SetTitle(";Energy (a.u);Counts");
    Photons_Spectrum_Tile2->GetXaxis()->SetLabelFont(62);
    Photons_Spectrum_Tile2->GetXaxis()->SetTitleFont(62);
    Photons_Spectrum_Tile2->GetYaxis()->SetLabelFont(62);
    Photons_Spectrum_Tile2->GetYaxis()->SetTitleFont(62);
    Photons_Spectrum_Tile2->SetLineWidth(2);
    //Photons_Spectrum_Tile2->SetLineColor(kRed);
    TSpectrum *spectrum = new TSpectrum();

    int w = 0; //counter to count the number of acquisitions
    for (int k = ki; k <= kf; k++)
    {
        w += 1;
        cout << "******************************** n Acquisition : " << w << " *************************************" << endl;
        cout << "******************************** Acquisition : " << k << " *************************************" << endl;
        FileReader reader;
        TuileDam readTile;

        filename = Form("/home/hjira/calibration/data_lqua/Irradiation_XY_Homogene_%i_vect.lqua", k);

        cout << "path =  " << filename << "\n";
        QString Qfilename = QString::fromStdString(filename);

        if (reader.openReadFile(Qfilename) == FileReader::OPENNING_ERROR)
        {
            std::cout << "Recherche du fichier: " << filename << " pas de fichier" << std::endl;
            return -1;
        }

        while (reader.getNextEvent(&readTile) == FileReader::OK)
        {

            if(readTile.tuile.tuile == choice_tile)
            {
                Tile.push_back(readTile.tuile);
                Tile_Img.push_back(readTile.image);
            }
            else
            {
                std::cout << "ERREUR D'IDENTIFICATION DE LA TUILE: " << readTile.tuile.tuile << std::endl;
                return 0;
            }
        }

        std::cout << "Number of events in tile " << choice_tile << " : " << Tile.size() << std::endl;

        //===================================================================================
        //**************************called the functions********************************
        //===================================================================================

        for (int i = 0; i < Tile.size(); i++)
        {
            Photons_Spectrum_Tile2->Fill(Tile.at(i).photons);
        }
         ///////////////////////////choose  the position of the energy filter windows////////////////////////
        
        peak_found = spectrum->Search(Photons_Spectrum_Tile2,10,"",0.05); //Search the peaks in the spectrum and displays them in the histogram
        peaks_position = spectrum->GetPositionX();
        cout<<"peaks_position = "<<*peaks_position<<endl;   
        cout<<"peak_found = "<<peak_found <<endl;        
        std::cout<<PrintEnergyPeaks(peaks_position,peak_found)<<std::endl;
        std::cout<<"Choose the position of the energy filter window . Press 0 if you don't want apply filter"<<std::endl;
        std::cin>>energy_peak_Tile;                    
       
        //energy_peak_Tile = 9250; // without filter

        number_events_pixel = ray_number_events_pixel(Tile, number_pixels_axis, energy_peak_Tile);
        
        // print_event_vect(number_events_pixel);

        for (int i = 0; i < number_pixels; i++)
        {
            total_event[i] += number_events_pixel[i];
            number_events_pixel_k[l++] = number_events_pixel[i];
        }

    } //fine boucle k acquisition
    cout << "********************Total number of events per pixel*****************" << endl;
        print_event_vect(total_event);
    cout << "****************Normalized number of events per pixel*****************" << endl;
        vector<double> total_event_unit = ray_number_events_pixel_unit(total_event);
        print_event_vect(total_event_unit);
    cout << " ************************Profile 1d**********************" << endl;
        vector<double> profile_1d = profile_1d_vect(total_event);
    cout << " ********************Normalized 1d profile********************" << endl;
        vector<double> profile_1d_unit = ray_number_events_pixel_unit(profile_1d);
        print_profile_1d_vect(profile_1d_unit);
    cout << " ********************statistic********************" << endl;
    
    
    if((n==0)&&(energy_peak_Tile==0)){
        vector<double> statistic = ray_statistic(total_event, number_events_pixel_k, acquisition_number);
        print_event_vect_statistic(statistic);
    }else{
     cout << "standard deviation = 0 " << endl;
    }


    /***************************pdisplay by root********************************/
        Map2D = Map_number_events_pixel(total_event_unit, Map2D);
        Map2D->Draw("COLZ");
        monFlux.cd();
        Map2D->Write("Calib-Map");
    //------------------------------
        cout << "Photon spectrum exists" << endl;
        Photons_Spectrum_Tile2->Reset();
        Photons_Spectrum_Tile2 = Map_Photons_Spectrum(Tile, Photons_Spectrum_Tile2, energy_peak_Tile);
        Photons_Spectrum_Tile2->Draw(" ");
        monFlux.cd();
        Photons_Spectrum_Tile2->Write("Photons_Spectrum-Map");

    return 0;
}
