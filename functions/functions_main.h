/**
 * @author      : o.allegrini@ipnl.in2p3.fr
 * @file        : main
 * @created     : jeudi oct. 15, 2020 16:46:04 CEST
 * @comment     : This program analyses the data recorded in the CeBr3 detectors of the classic temporal camera and the data from the detector trigger (diamond hodoscope)
 */
// TEMPORAL CAMERA LIBRARY
#include "quadratorevt.h"

//C++ LIBRARIES
#include <iostream>
#include <filereader.h>
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
#include <algorithm> //class used for the max_element funtion

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
#include "TTimer.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TGWindow.h"
#include "TGClient.h"
#include "TPaveText.h"
#include "TMath.h"



/** sort_in_time function
 * @brief: sort the TuileEvt of a TuileDam from one of the arguments of the TuileEvt structure
 * @return: sorted vector of TuileEvt based on frame number and timestamp for events with the same frame number
 * @note:
 * */
struct sort_function
{
    inline bool operator() (const TuileEvt& struct1, const TuileEvt& struct2)
    {
        return ((struct1.htimestamp) < (struct2.htimestamp)) ||
               (((struct1.htimestamp) == (struct2.htimestamp)) && ((struct1.ltimestamp&0x00FFFFFF) > (struct2.ltimestamp&0x00FFFFFF)));
    }
};

/** Global analysis function
 * @brief: Apply an energy filter and/or temperature filters on the data, centered on the maximum peak +/- 10% for energy filter and +/- 0.25 Celcius degree for the temperature
 * @return: array of vectors containing: [0] CRT values [1]Energy data of Tile A [2] Energy data of Tile B
 * @note: Future improvements: Additionnal recording of the x,y,z positions after the different filtering options
 * */

std::array<std::vector<double>,3> Global_analysis (int energy_peak_TileA, int energy_peak_TileB, int temper_peak_TileA, int temper_peak_TileB, std::vector<TuileEvt> TileA, std::vector<TuileEvt> TileB, bool nrg_filtering, bool tmp_filtering){
    int it_TileA = 0, it_TileB = 0, CRT_true = 0, CRT = 0, energy_TileA, energy_TileB = 0, frameB = 0;
    bool real_value = false;
    std::array<std::vector<double>,3> record_data;
    while(((TileA.at(it_TileA+1).htimestamp) == (TileA.at(it_TileA).htimestamp)) || ((TileA.at(it_TileA).htimestamp) < (TileB.at(it_TileB).htimestamp))){
        it_TileA++;
    }
    while(((TileB.at(it_TileB).htimestamp) <= (TileA.at(it_TileA).htimestamp)) && ((it_TileA < (TileA.size()-1)) && (it_TileB < (TileB.size()-1)))){
        if ((TileB.at(it_TileB).htimestamp) == (TileA.at(it_TileA).htimestamp)){
            if (real_value == false){
                frameB = it_TileB;
                CRT = std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(it_TileB).ltimestamp&0x00FFFFFF)));
                energy_TileA = TileA.at(it_TileA).photons;
                energy_TileB = TileB.at(frameB).photons;
                CRT_true = (TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(it_TileB).ltimestamp&0x00FFFFFF);
                real_value = true;
                std::cout<<TileB.at(it_TileB).htimestamp<<" "<<energy_TileB<<" "<<CRT<<std::endl;
            }
            else if ((std::abs((int)(((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(it_TileB).ltimestamp&0x00FFFFFF))))) < CRT){
                frameB = it_TileB;
                CRT = std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(it_TileB).ltimestamp&0x00FFFFFF)));
                energy_TileA = TileA.at(it_TileA).photons;
                energy_TileB = TileB.at(frameB).photons;
                CRT_true = (TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF);
                std::cout<<"--"<<TileB.at(it_TileB).htimestamp<<" "<<energy_TileB<<" "<<CRT<<std::endl;
            }
        }
        //std::cout<<energy_TileB<<std::endl;
        it_TileB++;
        while ((TileB.at(it_TileB).htimestamp) > (TileA.at(it_TileA).htimestamp) && ((it_TileA < (TileA.size()-1)) && (it_TileB < (TileB.size()-1)))){
            if (((std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF)))) < CRT) && (TileA.at(it_TileA).htimestamp == TileB.at(frameB).htimestamp)){
                CRT = std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF)));
                energy_TileA = TileA.at(it_TileA).photons;
                energy_TileB = TileB.at(frameB).photons;
                CRT_true = (TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF);
                std::cout<<"---"<<TileB.at(frameB).htimestamp<<" "<<energy_TileB<<" "<<CRT<<std::endl;
            }
            //std::cout<<"-----"<<energy_TileB<<std::endl;
            if (((TileA.at(it_TileA+1).htimestamp) > (TileA.at(it_TileA).htimestamp)) && (TileA.at(it_TileA).htimestamp == TileB.at(frameB).htimestamp) && (real_value == true) && (it_TileA+1 <= TileA.size())){
                if(nrg_filtering == false && tmp_filtering == false){
                    record_data[0].push_back((double)CRT_true*5/256);
                    record_data[1].push_back(energy_TileA);
                    record_data[2].push_back(energy_TileB);
                }
                ///////////////// energy filtering only ///////////////////
                else if(nrg_filtering == true && tmp_filtering == false){
                    if((energy_peak_TileA != -1) && (energy_peak_TileB != -1)){
                        if((energy_TileA >= energy_peak_TileA*0.9) && (energy_TileA <= energy_peak_TileA*1.1) && (energy_TileB >= energy_peak_TileB*0.9) && (energy_TileB <= energy_peak_TileB*1.1)){
                            record_data[0].push_back((double)CRT_true*5/256);
                            record_data[1].push_back(energy_TileA);
                            record_data[2].push_back(energy_TileB);
                        }
                    }
                    else if(energy_peak_TileA == -1){
                        if((energy_TileB >= energy_peak_TileB*0.9) && (energy_TileB <= energy_peak_TileB*1.1)){
                            record_data[0].push_back((double)CRT_true*5/256);
                            record_data[1].push_back(energy_TileA);
                            record_data[2].push_back(energy_TileB);
                        }
                    }
                    else if(energy_peak_TileB == -1){
                        if((energy_TileA >= energy_peak_TileA*0.9) && (energy_TileA <= energy_peak_TileA*1.1)){
                            record_data[0].push_back((double)CRT_true*5/256);
                            record_data[1].push_back(energy_TileA);
                            record_data[2].push_back(energy_TileB);
                            std::cout<<"------------>"<<TileB.at(frameB).htimestamp<<" "<<energy_TileB<<" "<<CRT<<std::endl;
                        }
                    }
                }
                ////////////////// temperature filtering only ///////////////////
                else if(nrg_filtering == false && tmp_filtering == true){
                    if((temper_peak_TileA != -1) && (temper_peak_TileB != -1)){
                        if(((TileA.at(it_TileA).temper >= temper_peak_TileA-25) && (TileA.at(it_TileA).temper <= temper_peak_TileA+25)) && ((TileB.at(it_TileB-1).temper >= temper_peak_TileB-25) && (TileB.at(it_TileB-1).temper <= temper_peak_TileB*+25))){
                            record_data[0].push_back((double)CRT_true*5/256);
                            record_data[1].push_back(energy_TileA);
                            record_data[2].push_back(energy_TileB);
                        }
                    }
                    else if(temper_peak_TileA == -1){
                        if((TileB.at(it_TileB-1).temper >= temper_peak_TileB-25) && (TileB.at(it_TileB-1).temper <= temper_peak_TileB+25)){
                            record_data[0].push_back((double)CRT_true*5/256);
                            record_data[1].push_back(energy_TileA);
                            record_data[2].push_back(energy_TileB);
                        }
                    }
                    else if(temper_peak_TileB == -1){
                        if((TileA.at(it_TileA).temper >= temper_peak_TileA-25) && (TileA.at(it_TileA).temper >= temper_peak_TileA+25)){
                            record_data[0].push_back((double)CRT_true*5/256);
                            record_data[1].push_back(energy_TileA);
                            record_data[2].push_back(energy_TileB);
                        }
                    }
                }
                ////////////////// energy and temperature filtering ///////////////////
                else if(nrg_filtering == true && tmp_filtering == true){
                    if((energy_peak_TileA != -1) && (energy_peak_TileB != -1)){
                        if(((energy_TileA >= energy_peak_TileA*0.9) && (energy_TileA <= energy_peak_TileA*1.1)) && ((energy_TileB >= energy_peak_TileB*0.9) && (energy_TileB <= energy_peak_TileB*1.1)) && ((TileA.at(it_TileA).temper >= temper_peak_TileA-25) && (TileA.at(it_TileA).temper <= temper_peak_TileA+25)) && ((TileB.at(it_TileB-1).temper >= temper_peak_TileB-25) && (TileB.at(it_TileB-1).temper <= temper_peak_TileB+25))){
                            record_data[0].push_back((double)CRT_true*5/256);
                            record_data[1].push_back(energy_TileA);
                            record_data[2].push_back(energy_TileB);
                        }
                    }
                    else if(energy_peak_TileA == -1){
                        if(((energy_TileB >= energy_peak_TileB*0.9) && (energy_TileB <= energy_peak_TileB*1.1)) && ((TileB.at(it_TileB-1).temper >= temper_peak_TileB-25) && (TileB.at(it_TileB-1).temper <= temper_peak_TileB+25))){
                            record_data[0].push_back((double)CRT_true*5/256);
                            record_data[1].push_back(energy_TileA);
                            record_data[2].push_back(energy_TileB);
                        }
                    }
                    else if(energy_peak_TileB == -1){
                        if(((energy_TileA >= energy_peak_TileA*0.9) && (energy_TileA >= energy_peak_TileA*1.1)) && ((TileA.at(it_TileA).temper >= temper_peak_TileA-25) && (TileA.at(it_TileA).temper >= temper_peak_TileA+25))){
                            record_data[0].push_back((double)CRT_true*5/256);
                            record_data[1].push_back(energy_TileA);
                            record_data[2].push_back(energy_TileB);
                        }
                    }
                }
                CRT = 0;
                CRT_true = 0;
                energy_TileA = 0;
                energy_TileB = 0;
                real_value = false;

            }
            it_TileA++;
        }
    }
    return record_data;
}

/** Raw 2D map function
 * @brief: Reconstruction of the 2D with all events in a Tile without filter and/or CRT selection
 * @return: 2D array[8][8] containing the sum of the photons per pixel for the whole acquisition.
 * @note: Mainly used for developpement and debug. Recontsruction with CRT selection and filters is more appropriated.
 * @note: This function doesn't use the X Y Z coordinates obtained by vectorization. Potential 3D reconstruction is the impossible
 * */

 TH2I* Raw_map2D (std::vector<TuileImg> TileImg, int dim){
    //uint32_t[8][8]
    TH2I* map2D = new TH2I("2D_Map", "2D_Map", dim, -0.5, 8.5, dim, -0.5, 8.5);
    //uint32_t map2D [dim][dim];

    //loop on data
    for(int i=0; i<TileImg.size(); i++){
        //double loop on pixel X array and Y column
        for(int ii=0; ii<dim; ii++){
            for(int iii=0; iii<dim; iii++){
                //map2D[ii][iii] = map2D[ii][iii]+TileImg.at(i).g_pix[ii][iii];
                map2D->SetBinContent(ii,iii,map2D->GetBinContent(ii,iii)+TileImg.at(i).g_pix[ii][iii]);
            }
        }
    }
    return map2D;
}
