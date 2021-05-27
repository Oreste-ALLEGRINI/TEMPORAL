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
#include <cstddef>
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
#include "TPaveLabel.h"
#include "TMath.h"
#include "TH3.h"

struct TIME_PHOTON //object containing the time of a die and the sum of the photons of this die
{
    int time;
    int photon;  // Does not participate in comparisons
    bool operator() (TIME_PHOTON i,TIME_PHOTON j) { return (i.time<j.time); }
}time_photon;

/** Find_Median function
 * @brief: split the full path of filename and extension to add  the output repository in the path
 * @return: the full path of the output analysis file (without extension)
 * @note:
 * */
double FindMedian (std::vector<double> v)
{
  int median;
  std::sort(v.begin(),v.end());
  median = v.at(v.size()/2);
}

/** SplitFilename function
 * @brief: split the full path of filename and extension to add  the output repository in the path
 * @return: the full path of the output analysis file (without extension)
 * @note:
 * */
std::string SplitFilename (const std::string& str, std::string repository)
{
  std::string file;
  std::string interfile;
  std::size_t rm_path = str.find_last_of("/\\");
  file = str.substr(0,rm_path);
  file = (const std::string&)file;
  rm_path = file.find_last_not_of("/\\");
  file = file.substr(0,rm_path);
  rm_path = file.find_last_of("/\\");
  file = file.substr(0,rm_path);
  file+=repository;

  rm_path = str.find_last_of("/\\");
  interfile = str.substr(rm_path+1);
  interfile = (const std::string&)interfile;
  interfile = interfile.substr(0,rm_path);
  std::size_t rm_ext = interfile.find_last_of(".");
  interfile = interfile.substr(0,rm_ext);

  file+=interfile;
  return file;
}


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

/** CRT_calculation
 * @brief: This function uses the data contained in the TuileImage.g_die to compute the mean CRT of the 3 first die. First the hit time of each die is sorted and then the 3 first time are use to calculate the average
 * @return: The int value of the average of the time of the three first hitted die per Tile.
 * @note:
 * */
double CRT_calculation (int frame, std::vector<TuileImg> TileImg) {
    int entry_time[4][4];
    int entry_photons[4][4];
    double avr_CRT;
    int sum_die; // sum of the photons in the first three hitted dies
    std::vector<TIME_PHOTON> time_die;

    //Recovering the dies data of the events frameA and frameB
    for(int i=0; i<8; i+=2){
        for(int ii=0; ii<8; ii+=2){
            entry_time[i/2][ii/2] = TileImg.at(frame).g_die[i/2][ii/2];
            entry_photons[i/2][ii/2] = TileImg.at(frame).g_pix[i][ii]+TileImg.at(frame).g_pix[i][ii+1]+TileImg.at(frame).g_pix[i+1][ii]+TileImg.at(frame).g_pix[i+1][ii+1]; //Computing the sum of photons contained in each die
            time_die.push_back({entry_time[i/2][ii/2],entry_photons[i/2][ii/2]});
        }
    }

    std::sort(time_die.begin(), time_die.end(),time_photon);
    while(time_die.at(0).time == -1){
        time_die.emplace_back(time_die.at(0));
        time_die.erase(time_die.begin());
    }
    sum_die = time_die.at(0).photon + time_die.at(1).photon + time_die.at(2).photon;
    //std::sort((unsigned int*)&entry_time[0][0], (unsigned int *)&entry_time[3][3]+1);
    //avr_CRT=(entry_time[0][0]+entry_time[0][1]+entry_time[0][2])/3;
    avr_CRT=(time_die.at(0).time*time_die.at(0).photon+time_die.at(1).time*time_die.at(1).photon+time_die.at(2).time*time_die.at(2).photon)/sum_die;

    return avr_CRT;
}

/** Global analysis function
 * @brief: This function recovers the events which have a corresponding frame number (time coincidence)
 * @return: array of vectors containing: [0] CRT values [1]Energy data of Tile A [2] Energy data of Tile B [3] Temperature TileA [4] Temperature TileB [5] X position Tile A [6] Y position Tile A [7] Z position TileA [8][9][10] positions for tile B
 * @note: Future improvements: Additionnal recording of the x,y,z positions after the different filtering options
 * */

std::array<std::vector<double>,11> Global_analysis_bis (std::vector<TuileEvt> TileA, std::vector<TuileEvt> TileB, std::vector<TuileImg> TileImgA, std::vector<TuileImg> TileImgB, bool WeightedMeanTimeCorr){
    int counter = 0, it_TileA = 0, it_TileB = 0, it_TileB_eq = 0, CRT_true = 0, CRT = -1, energy_TileA =0, energy_TileB = 0, frameB = 0, temper_TileA = 0, temper_TileB = 0, posX_TileA = 0, posX_TileB = 0, posY_TileA = 0, posY_TileB = 0, posZ_TileA = 0, posZ_TileB = 0;
    bool real_value = false;
    float progress=0.0;
    double CRT_corrA= 0, CRT_corrB = 0;
    std::array<std::vector<double>,11> record_data;
    while(progress < 1.0) {
        int barWidth = 70;
        while(((TileA.at(it_TileA).htimestamp) < (TileB.at(it_TileB).htimestamp))){
            it_TileA++;
        }
        while(((TileB.at(it_TileB).htimestamp) <= (TileA.at(it_TileA).htimestamp)) && ((it_TileA < (TileA.size()-1)) && (it_TileB < (TileB.size()-1)))){
            if ((TileB.at(it_TileB).htimestamp) == (TileA.at(it_TileA).htimestamp)){
                it_TileB_eq = it_TileB;
                while((TileB.at(it_TileB_eq).htimestamp)==(TileB.at(it_TileB).htimestamp)){
                    if(CRT==-1 && real_value == false){
                        CRT = std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(it_TileB_eq).ltimestamp&0x00FFFFFF)));
                        CRT_true = (TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(it_TileB_eq).ltimestamp&0x00FFFFFF);
                        energy_TileA = TileA.at(it_TileA).photons;
                        energy_TileB = TileB.at(it_TileB_eq).photons;
                        temper_TileA = TileA.at(it_TileA).temper;
                        temper_TileB = TileB.at(it_TileB_eq).temper;
                        posX_TileA = TileA.at(it_TileA).mX;
                        posX_TileB = TileB.at(it_TileB_eq).mX;
                        posY_TileA = TileA.at(it_TileA).mY;
                        posY_TileB = TileB.at(it_TileB_eq).mY;
                        posZ_TileA = TileA.at(it_TileA).mZ;
                        posZ_TileB = TileB.at(it_TileB_eq).mZ;
                        real_value = true;
                    }
                    else if(real_value == true){
                        if(std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(it_TileB_eq).ltimestamp&0x00FFFFFF))) <= CRT){
                            CRT = std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(it_TileB_eq).ltimestamp&0x00FFFFFF)));
                            CRT_true = (TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(it_TileB_eq).ltimestamp&0x00FFFFFF);
                            energy_TileA = TileA.at(it_TileA).photons;
                            energy_TileB = TileB.at(it_TileB_eq).photons;
                            temper_TileA = TileA.at(it_TileA).temper;
                            temper_TileB = TileB.at(it_TileB_eq).temper;
                            posX_TileA = TileA.at(it_TileA).mX;
                            posX_TileB = TileB.at(it_TileB_eq).mX;
                            posY_TileA = TileA.at(it_TileA).mY;
                            posY_TileB = TileB.at(it_TileB_eq).mY;
                            posZ_TileA = TileA.at(it_TileA).mZ;
                            posZ_TileB = TileB.at(it_TileB_eq).mZ;
                            counter++;
                            frameB = it_TileB_eq;
                        }
                    }
                    ///////// MAYBE THIS PART CAN BE REMOVED. SHOULD BE TESTED WITH AN HIGH INTENSITY ACQUISITION WITH A LOT OF EVENTS PER FRAME AND AN IMPORTANT NOISE ////////
                    if((TileB.at(it_TileB_eq+1).htimestamp) > (TileB.at(it_TileB_eq).htimestamp) && (counter > 1) && (TileB.at(it_TileB_eq).htimestamp) != (TileB.at(frameB).htimestamp)){
                        if(WeightedMeanTimeCorr==1){
                            record_data[0].push_back((double)CRT_true*5/256 /*+ (CRT_calculation(it_TileA, TileImgA) - CRT_calculation(frameB, TileImgB))*5/256*/);
                        }
                        else {
                            record_data[0].push_back((double)CRT_true*5/256 /*+ CRT_calculation(frameB, TileImgB)*5/256*/);
                        }
                        record_data[1].push_back(energy_TileA);
                        record_data[2].push_back(energy_TileB);
                        record_data[3].push_back(temper_TileA);
                        record_data[4].push_back(temper_TileB);
                        record_data[5].push_back(posX_TileA);
                        record_data[6].push_back(posX_TileB);
                        record_data[7].push_back(posY_TileA);
                        record_data[8].push_back(posY_TileB);
                        record_data[9].push_back(posZ_TileA);
                        record_data[10].push_back(posZ_TileB);
                        CRT = -1;
                        CRT_true = 0;
                        energy_TileA = 0;
                        energy_TileB = 0;
                        temper_TileA = 0;
                        temper_TileB = 0;
                        posX_TileA = 0; posY_TileA = 0; posZ_TileA =0;
                        posX_TileB = 0; posY_TileB = 0; posZ_TileB =0;
                        real_value = false;
                        counter = 0;
                    }
                    //////////////////////////////////////////////////
                    else if ((TileB.at(it_TileB_eq+1).htimestamp) > (TileB.at(it_TileB_eq).htimestamp)){
                        if(WeightedMeanTimeCorr==1){
                            record_data[0].push_back((double)CRT_true*5/256 /*+ (CRT_calculation(it_TileA, TileImgA) - CRT_calculation(it_TileB_eq, TileImgB))*5/256*/);
                        }
                        else {
                            record_data[0].push_back((double)CRT_true*5/256 /*+ CRT_calculation(it_TileB_eq, TileImgB)*5/256*/);
                        }
                        record_data[1].push_back(energy_TileA);
                        record_data[2].push_back(energy_TileB);
                        record_data[3].push_back(temper_TileA);
                        record_data[4].push_back(temper_TileB);
                        record_data[5].push_back(posX_TileA);
                        record_data[6].push_back(posX_TileB);
                        record_data[7].push_back(posY_TileA);
                        record_data[8].push_back(posY_TileB);
                        record_data[9].push_back(posZ_TileA);
                        record_data[10].push_back(posZ_TileB);
                        CRT = -1;
                        CRT_true = 0;
                        energy_TileA = 0;
                        energy_TileB = 0;
                        temper_TileA = 0;
                        temper_TileB = 0;
                        posX_TileA = 0; posY_TileA = 0; posZ_TileA =0;
                        posX_TileB = 0; posY_TileB = 0; posZ_TileB =0;
                        real_value = false;
                        counter = 0;

                    }
                    it_TileB_eq++;
                }
                std::cout << "[";
                    int pos = barWidth * progress;
                    for (int i = 0; i < barWidth; ++i) {
                        if (i < pos) std::cout << "=";
                        else if (i == pos) std::cout << ">";
                        else std::cout << " ";
                    }
                std::cout << "] " << int(progress * 100.0) << " %\r";
                progress = (float)it_TileA/TileA.size();
                std::cout.flush();
                //CRT_calculation(it_TileA, TileImgA);
                //CRT_calculation(frameB, TileImgB);

            }
            it_TileB++;
            while ((TileB.at(it_TileB).htimestamp) > (TileA.at(it_TileA).htimestamp) && ((it_TileA < (TileA.size()-1)) && (it_TileB < (TileB.size()-1)))){
                it_TileA++;
            }
        }
        progress = 1;
        std::cout.flush();
    }
    std::cout << std::endl;
    return record_data;
}



/** Global analysis function
 * @brief: This function recovers the events which have a corresponding frame number (time coincidence)
 * @return: array of vectors containing: [0] CRT values [1]Energy data of Tile A [2] Energy data of Tile B [3] Temperature TileA [4] Temperature TileB [5] X position Tile A [6] Y position Tile A [7] Z position TileA [8][9][10] positions for tile B
 * @note: Future improvements: Additionnal recording of the x,y,z positions after the different filtering options
 * */

std::array<std::vector<double>,11> Global_analysis (std::vector<TuileEvt> TileA, std::vector<TuileEvt> TileB, std::vector<TuileImg> TileImgA, std::vector<TuileImg> TileImgB, bool WeightedMeanTimeCorr){
    int it_TileA = 0, it_TileB = 0, CRT_true = 0, CRT = 0, energy_TileA =0, energy_TileB = 0, frameB = 0, temper_TileA = 0, temper_TileB = 0, posX_TileA = 0, posX_TileB = 0, posY_TileA = 0, posY_TileB = 0, posZ_TileA = 0, posZ_TileB = 0;
    bool real_value = false;
    float progress=0.0;
    double CRT_corrA= 0, CRT_corrB = 0;
    std::array<std::vector<double>,11> record_data;
    while(progress < 1.0) {
        int barWidth = 70;
        while(((TileA.at(it_TileA+1).htimestamp) == (TileA.at(it_TileA).htimestamp)) || ((TileA.at(it_TileA).htimestamp) < (TileB.at(it_TileB).htimestamp))){
            it_TileA++;
        }
        while(((TileB.at(it_TileB).htimestamp) <= (TileA.at(it_TileA).htimestamp)) && ((it_TileA < (TileA.size()-1)) && (it_TileB < (TileB.size()-1)))){
            if ((TileB.at(it_TileB).htimestamp) == (TileA.at(it_TileA).htimestamp)){
                if(CRT==0 && real_value == false){
                    frameB = it_TileB;
                    CRT = std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF)));
                    CRT_true = (TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF);
                    energy_TileA = TileA.at(it_TileA).photons;
                    energy_TileB = TileB.at(frameB).photons;
                    temper_TileA = TileA.at(it_TileA).temper;
                    temper_TileB = TileB.at(frameB).temper;
                    posX_TileA = TileA.at(it_TileA).mX;
                    posX_TileB = TileB.at(frameB).mX;
                    posY_TileA = TileA.at(it_TileA).mY;
                    posY_TileB = TileB.at(frameB).mY;
                    posZ_TileA = TileA.at(it_TileA).mZ;
                    posZ_TileB = TileB.at(frameB).mZ;
                    real_value = true;
                }
                else if(real_value == true){
                    if(std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(it_TileB).ltimestamp&0x00FFFFFF))) <= CRT){
                        frameB = it_TileB;
                        CRT = std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF)));
                        CRT_true = (TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF);
                        energy_TileA = TileA.at(it_TileA).photons;
                        energy_TileB = TileB.at(frameB).photons;
                        temper_TileA = TileA.at(it_TileA).temper;
                        temper_TileB = TileB.at(frameB).temper;
                        posX_TileA = TileA.at(it_TileA).mX;
                        posX_TileB = TileB.at(frameB).mX;
                        posY_TileA = TileA.at(it_TileA).mY;
                        posY_TileB = TileB.at(frameB).mY;
                        posZ_TileA = TileA.at(it_TileA).mZ;
                        posZ_TileB = TileB.at(frameB).mZ;
                    }
                }
            }
            it_TileB++;
            while ((TileB.at(it_TileB).htimestamp) > (TileA.at(it_TileA).htimestamp) && ((it_TileA < (TileA.size()-1)) && (it_TileB < (TileB.size()-1)))){
                if (TileA.at(it_TileA).htimestamp == TileB.at(frameB).htimestamp){
                    if(std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF))) < CRT){
                        CRT = std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF)));
                        CRT_true = (TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF);
                        energy_TileA = TileA.at(it_TileA).photons;
                        energy_TileB = TileB.at(frameB).photons;
                        temper_TileA = TileA.at(it_TileA).temper;
                        temper_TileB = TileB.at(frameB).temper;
                        posX_TileA = TileA.at(it_TileA).mX;
                        posX_TileB = TileB.at(frameB).mX;
                        posY_TileA = TileA.at(it_TileA).mY;
                        posY_TileB = TileB.at(frameB).mY;
                        posZ_TileA = TileA.at(it_TileA).mZ;
                        posZ_TileB = TileB.at(frameB).mZ;
                    }
                    if (((TileA.at(it_TileA+1).htimestamp) > (TileA.at(it_TileA).htimestamp)) && (it_TileA+1 <= TileA.size())){
                        //CRT_corrA = CRT_calculation(it_TileA, TileImgA);
                        if(WeightedMeanTimeCorr==1){
                            record_data[0].push_back((double)CRT_true*5/256 /*+ (CRT_calculation(it_TileA, TileImgA) - CRT_calculation(frameB, TileImgB))*5/256*/);
                        }
                        else{
                            record_data[0].push_back((double)CRT_true*5/256 /*+ (CRT_calculation(it_TileA, TileImgA)*5/256*/);
                        }
                        record_data[1].push_back(energy_TileA);
                        record_data[2].push_back(energy_TileB);
                        record_data[3].push_back(temper_TileA);
                        record_data[4].push_back(temper_TileB);
                        record_data[5].push_back(posX_TileA);
                        record_data[6].push_back(posX_TileB);
                        record_data[7].push_back(posY_TileA);
                        record_data[8].push_back(posY_TileB);
                        record_data[9].push_back(posZ_TileA);
                        record_data[10].push_back(posZ_TileB);
                        CRT = 0;
                        CRT_true = 0;
                        energy_TileA = 0;
                        energy_TileB = 0;
                        temper_TileA = 0;
                        temper_TileB = 0;
                        posX_TileA = 0; posY_TileA = 0; posZ_TileA =0;
                        posX_TileB = 0; posY_TileB = 0; posZ_TileB =0;
                        real_value = false;

                        std::cout << "[";
                            int pos = barWidth * progress;
                            for (int i = 0; i < barWidth; ++i) {
                                if (i < pos) std::cout << "=";
                                else if (i == pos) std::cout << ">";
                                else std::cout << " ";
                            }
                        std::cout << "] " << int(progress * 100.0) << " %\r";
                            std::cout.flush();
                        //CRT_calculation(it_TileA, TileImgA);
                        //CRT_calculation(frameB, TileImgB);
                            progress = (float)it_TileA/TileA.size();

                    }
                }
                it_TileA++;
            }
        }
        progress = 1;
        std::cout.flush();
    }
    std::cout << std::endl;
    return record_data;
}

/** CRT filter function
 * @brief: This function computes the CRT peak and the baseline due to electronical noise in order to filter the data which correspond to physical events.
 * @return: array of vectors containing ([0] CRT [1] Energy TileA [2] Energy TileB) of physical events under the peak.
 * @note: The function Global_analysis is called to provide the data having a frame coincidence in the two tiles which are compared
 **/
std::array<std::vector<double>,11> CRT_filter (std::vector<TuileEvt> TileA, std::vector<TuileEvt> TileB, std::vector<TuileImg> TileImgA, std::vector<TuileImg> TileImgB, bool WeightedMeanTimeCorr){
    std::array<std::vector<double>,11> Raw_data;
    std::array<std::vector<double>,11> Clear_data;
    double CRT_peak;
    TH1F* Hist_Raw_data = new TH1F("","",83886080,-41943040,41943040);
    Raw_data = Global_analysis_bis(TileA, TileB, TileImgA, TileImgB, WeightedMeanTimeCorr);

    //Seek for the CRT peak
    for(int i=0; i<Raw_data[0].size(); i++){
        Hist_Raw_data->Fill(Raw_data[0].at(i));
    }
    CRT_peak = Hist_Raw_data->GetBinCenter(Hist_Raw_data->GetMaximumBin());
    //std::cout<<CRT_peak<<std::endl;  //Use it to find the CRT peak and redefine the histogram range

    //Data filtering (window +/- 5 ns which also corresponds to +/- 256 timestamp)
    for(int i=0; i<Raw_data[0].size(); i++){
        if((Raw_data[0].at(i) >= CRT_peak-5) && (Raw_data[0].at(i) <= CRT_peak+5)){
            Clear_data[0].push_back(Raw_data[0].at(i));
            Clear_data[1].push_back(Raw_data[1].at(i));
            Clear_data[2].push_back(Raw_data[2].at(i));
            Clear_data[3].push_back(Raw_data[3].at(i));
            Clear_data[4].push_back(Raw_data[4].at(i));
            Clear_data[5].push_back(Raw_data[5].at(i));
            Clear_data[6].push_back(Raw_data[6].at(i));
            Clear_data[7].push_back(Raw_data[7].at(i));
            Clear_data[8].push_back(Raw_data[8].at(i));
            Clear_data[9].push_back(Raw_data[9].at(i));
            Clear_data[10].push_back(Raw_data[10].at(i));
        }
    }
    return Clear_data;
}

/** Nrj Temper filter function
 * @brief: This function Apply an energy filter and/or temperature filters on the physical events, centered on the maximum peak +/- 10% for energy filter and +/- 0.25 Celcius degree for the temperature
 * @return: array of vectors containing ([0] CRT [1] Energy TileA [2] Energy TileB) of physical events under the peak.
 * @note: The function CRT_filter is called to provide the input data
 **/
std::array<std::vector<double>,11> Nrj_Temper_filter(bool nrg_filtering, bool tmp_filtering, int energy_peak_TileA, int energy_peak_TileB, int temper_peak_TileA, int temper_peak_TileB, std::vector<TuileEvt> TileA, std::vector<TuileEvt> TileB, std::vector<TuileImg> TileImgA, std::vector<TuileImg> TileImgB, bool WeightedMeanTimeCorr){
    std::array<std::vector<double>,11> input_data;
    std::array<std::vector<double>,11> record_data;
    input_data = CRT_filter(TileA, TileB, TileImgA, TileImgB, WeightedMeanTimeCorr);
    for(int i=0; i<input_data[0].size(); i++){
        ///////////////// energy filtering only ///////////////////
        if(nrg_filtering == true && tmp_filtering == false){
            if((energy_peak_TileA != -1) && (energy_peak_TileB != -1)){
                //if((input_data[1].at(i) + input_data[2].at(i) >= 0.9*(energy_peak_TileA+energy_peak_TileB)/2) && (input_data[1].at(i) + input_data[2].at(i) >= 1.2*(energy_peak_TileA+energy_peak_TileB)/2)){ //==> Filter to study the Compton events of 1275 keV gammas from 22Na source
                if((input_data[1].at(i) >= energy_peak_TileA*0.9) && (input_data[1].at(i) <= energy_peak_TileA*1.1) && (input_data[2].at(i) >= energy_peak_TileB*0.9) && (input_data[2].at(i) <= energy_peak_TileB*1.1)){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if((energy_peak_TileA == -1) && (energy_peak_TileB != -1)){
                if((input_data[2].at(i) >= energy_peak_TileB*0.9) && (input_data[2].at(i) <= energy_peak_TileB*1.1)){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if((energy_peak_TileB == -1) && (energy_peak_TileA != -1)){
                if((input_data[1].at(i) >= energy_peak_TileA*0.9) && (input_data[1].at(i) <= energy_peak_TileA*1.1)){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if((energy_peak_TileB == -1) && (energy_peak_TileA == -1)){
                record_data[0].push_back(input_data[0].at(i));
                record_data[1].push_back(input_data[1].at(i));
                record_data[2].push_back(input_data[2].at(i));
                record_data[3].push_back(input_data[3].at(i));
                record_data[4].push_back(input_data[4].at(i));
                record_data[5].push_back(input_data[5].at(i));
                record_data[6].push_back(input_data[6].at(i));
                record_data[7].push_back(input_data[7].at(i));
                record_data[8].push_back(input_data[8].at(i));
                record_data[9].push_back(input_data[9].at(i));
                record_data[10].push_back(input_data[10].at(i));
            }
        }
        ////////////////// temperature filtering only ///////////////////
        else if(nrg_filtering == false && tmp_filtering == true){
            if((temper_peak_TileA != -1) && (temper_peak_TileB != -1)){
                if(((input_data[3].at(i) >= temper_peak_TileA-25) && (input_data[3].at(i) <= temper_peak_TileA+25)) && ((input_data[4].at(i) >= temper_peak_TileB-25) && (input_data[4].at(i) <= temper_peak_TileB*+25))){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if(temper_peak_TileA == -1){
                if((input_data[4].at(i) >= temper_peak_TileB-25) && (input_data[4].at(i) <= temper_peak_TileB+25)){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if(temper_peak_TileB == -1){
                if((input_data[3].at(i) >= temper_peak_TileA-25) && (input_data[3].at(i) >= temper_peak_TileA+25)){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
        }
        ////////////////// energy and temperature filtering ///////////////////
        else if(nrg_filtering == true && tmp_filtering == true){
            if((energy_peak_TileA != -1) && (energy_peak_TileB != -1)){
                if(((input_data[1].at(i) >= energy_peak_TileA*0.9) && (input_data[1].at(i) <= energy_peak_TileA*1.1)) && ((input_data[2].at(i) >= energy_peak_TileB*0.9) && (input_data[2].at(i) <= energy_peak_TileB*1.1)) && ((input_data[3].at(i) >= temper_peak_TileA-25) && (input_data[3].at(i) <= temper_peak_TileA+25)) && ((input_data[4].at(i) >= temper_peak_TileB-25) && (input_data[4].at(i) <= temper_peak_TileB+25))){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if(energy_peak_TileA == -1){
                if(((input_data[2].at(i) >= energy_peak_TileB*0.9) && (input_data[2].at(i) <= energy_peak_TileB*1.1)) && ((input_data[4].at(i) >= temper_peak_TileB-25) && (input_data[4].at(i) <= temper_peak_TileB+25))){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if(energy_peak_TileB == -1){
                if(((input_data[1].at(i) >= energy_peak_TileA*0.9) && (input_data[1].at(i) <= energy_peak_TileA*1.1)) && ((input_data[3].at(i) >= temper_peak_TileA-25) && (input_data[3].at(i) <= temper_peak_TileA+25))){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
        }

    }
    return record_data;
}

/** 1D_DepthProfile function
 * @brief: Reconstruction of the 2D with all events in a Tile without filter and/or CRT selection
 * @return: TH1F of the interaction depth of the gamma in the detector.
 * @note: Mainly used for developpement and debug. Recontsruction with CRT selection and filters is more appropriated.
 * @note: This function doesn't use the X Y Z coordinates obtained by vectorization. Potential 3D reconstruction is the impossible
 * */

 TH1F* DepthProfile (std::vector<double> PosX, std::vector<double> PosY, std::vector<double> PosZ, char dimension){
    float dimX;
    int nbins;
    if((dimension == ('X' | 'x')) || (dimension == ('Y' | 'y'))){
        dimX = 32;
        nbins = 8;
    }
    else if(dimension == ('Z' | 'z')){
        dimX = 20;
        nbins = 32;
    }

    TH1F* histo_1d = new TH1F("profile_1D","profile_1D", nbins, 0, dimX);
    histo_1d->SetTitle(";Position (mm);Counts");
    histo_1d->GetXaxis()->SetLabelFont(22);
    histo_1d->GetXaxis()->SetTitleFont(22);
    histo_1d->GetYaxis()->SetLabelFont(22);
    histo_1d->GetYaxis()->SetTitleFont(22);
    histo_1d->SetLineWidth(2);
    if(dimension == ('X' | 'x')){
        for(int i=0; i<PosX.size(); i++){
            histo_1d->Fill((PosX.at(i)-500)/(8000/32));
        }
    }
    else if(dimension == ('Y' | 'y')){
        for(int i=0; i<PosY.size(); i++){
            histo_1d->Fill((PosY.at(i)-500)/(8000/32));
        }
    }
    else if(dimension == ('Z' | 'z')){
        for(int i=0; i<PosZ.size(); i++){
            histo_1d->Fill((PosZ.at(i)/1000));
        }
    }
    return histo_1d;
}

/** 2D map function
 * @brief: Reconstruction of the 2D with all events in a Tile without filter and/or CRT selection
 * @return: TH2D[8][8] containing the sum of the photons per pixel for the whole acquisition.
 * @note: Mainly used for developpement and debug. Reconstruction with CRT selection and filters is more appropriated.
 * @note: This function doesn't use the X Y Z coordinates obtained by vectorization. Potential 3D reconstruction is the impossible
 * */

 TH2D* map2D (std::vector<double> PosX, std::vector<double> PosY, int dim){
    TH2D* map2D = new TH2D("2D_Map", "2D_Map", dim, 0, 3.2, dim, 0, 3.2);
    map2D->SetTitle(";X (cm);Y (cm)");
    map2D->GetXaxis()->SetLabelFont(22);
    map2D->GetXaxis()->SetTitleFont(22);
    map2D->GetYaxis()->SetLabelFont(22);
    map2D->GetYaxis()->SetTitleFont(22);
    //uint32_t map2D [dim][dim];

    //loop on data
    for(int i =0; i< PosX.size(); i++){
        map2D->Fill((PosX.at(i)-500)/(8000/3.2), (PosY.at(i)-500)/(8000/3.2));
    }
    return map2D;
}

 /** 3D map function
  * @brief: Reconstruction of the 2D with all events in a Tile without filter and/or CRT selection
  * @return: TH3F [8][8][3] containing the sum of the photons per pixel for the whole acquisition.
  * @note: Mainly used for developpement and debug. Recontsruction with CRT selection and filters is more appropriated.
  * @note: This function doesn't use the X Y Z coordinates obtained by vectorization. Potential 3D reconstruction is the impossible
  * */

  TH3F* map3D (std::vector<double> PosX, std::vector<double> PosY, std::vector<double> PosZ, int dim){
     TH3F* map3D = new TH3F("3D_Map", "3D_Map", dim, 0, 3.2, dim, 0, 3.2, 18, 0, 2);
     map3D->SetTitle(";X (cm);Y (cm);Z (cm)");
     map3D->GetXaxis()->SetLabelFont(22);
     map3D->GetXaxis()->SetTitleFont(22);
     map3D->GetYaxis()->SetLabelFont(22);
     map3D->GetYaxis()->SetTitleFont(22);
     map3D->GetZaxis()->SetLabelFont(22);
     map3D->GetZaxis()->SetTitleFont(22);

     //loop on data
     for(int i =0; i< PosX.size(); i++){
         map3D->Fill((PosX.at(i)-500)/(8000/3.2), (PosY.at(i)-500)/(8000/3.2), PosZ.at(i)/(20000/2));
     }
     return map3D;
 }

  /** E_TileA_vs_E_TileB
   * @brief: Reconstruction of the 2D with all events in a Tile without filter and/or CRT selection
   * @return: TH2D[8][8] containing the sum of the photons per pixel for the whole acquisition.
   * @note: Mainly used for developpement and debug. Recontsruction with CRT selection and filters is more appropriated.
   * @note: This function doesn't use the X Y Z coordinates obtained by vectorization. Potential 3D reconstruction is the impossible
   * */

   TH2D* E_TileA_vs_E_TileB (std::vector<double> E_TileA, std::vector<double> E_TileB){
      TH2D* ETAvsETB = new TH2D("", "", 1200, 0, 12000, 1200, 0, 12000);
      ETAvsETB->SetTitle(";E_TileA (a.u);E_TileB (a.u)");
      ETAvsETB->GetXaxis()->SetLabelFont(22);
      ETAvsETB->GetXaxis()->SetTitleFont(22);
      ETAvsETB->GetYaxis()->SetLabelFont(22);
      ETAvsETB->GetYaxis()->SetTitleFont(22);
      //uint32_t map2D [dim][dim];

      //loop on data
      for(int i =0; i< E_TileA.size(); i++){
          ETAvsETB->Fill(E_TileA.at(i), E_TileB.at(i));
      }
      return ETAvsETB;
  }

 /** PrintEnergyPeaks function
  * @brief: Display the energy peak positions in a table and print it in the shell
  * @return: string for pretty printing
  * */
 inline std::string PrintEnergyPeaks(double* position, int nb_peak)
 {
   std::string result;
   int width   = 11;
   char corner = '+';
   char hline  = '-';
   char vline  = '|';
   char sep    = ' ';

   auto hLineSep = [&]() {
     for(int i=0; i<nb_peak+1;i++){
       result.push_back(corner);
       result.append(width,hline);
     }
     result.push_back(corner);
     result.push_back('\n');
   };
   auto hLineNumber = [&](){
       result.push_back(vline);
       std::string zero = std::to_string(-1);
       int zeroSize = zero.size();
       result.append((width-zeroSize)/2 + (zeroSize+1)%2,sep);
       result.append(zero);
       result.append((width-zeroSize)/2,sep);

       for(int i=0;i<nb_peak;i++){
           std::string num = std::to_string(i+1);
           int numSize = num.size();
           result.push_back(vline);
           result.append((width-numSize)/2 + (numSize+1)%2,sep);
           result.append(num);
           result.append((width-numSize)/2,sep);
       }
       result.push_back(vline);
       result.push_back('\n');
   };

   auto hLinePos = [&](){
       result.push_back(vline);
       std::string none = "none";
       int noneSize = none.size();
       result.append((width-noneSize)/2 + (noneSize+1)%2,sep);
       result.append(none);
       result.append((width-noneSize)/2,sep);

       for(int i=0;i<nb_peak;i++){
           std::string pos = std::to_string((int)position[i]);
           int posSize = pos.size();
           result.push_back(vline);
           result.append((width-posSize)/2 + (posSize+1)%2,sep);
           result.append(pos);
           result.append((width-posSize)/2,sep);
       }
       result.push_back(vline);
       result.push_back('\n');
   };

   //printing
     hLineSep();
     hLineNumber();
     hLineSep();
     hLinePos();
     hLineSep();
   return result;
 }

 /** DOI Time_distribution function
  * @brief: Reconstruction of the 1D time distribution of the proton hit in the Tile segmented into 5 mm slices
  * @return: TCanvas containing 4 or 6 TH1F corresponding to the depth ranges 0/5 5/10 10/15 15/20 (20/25 25/30) mm
  * @note: Mainly used for the DOI correction development.
  * @note: This function uses the raw ltimestamp data of the Tile.
  * */

  TCanvas* DOI_TimeDistribution (std::vector<TuileEvt> Tile, char dimension){
     TCanvas* time_z = new TCanvas("","",1600,1600);
     TH1F* time_0_5mm = new TH1F("0/5mm","0/5mm", 2000, -100, +100);
     time_0_5mm->SetTitle(";Time (ns);Counts");
     time_0_5mm->GetXaxis()->SetLabelFont(22);
     time_0_5mm->GetXaxis()->SetTitleFont(22);
     time_0_5mm->GetYaxis()->SetLabelFont(22);
     time_0_5mm->GetYaxis()->SetTitleFont(22);
     time_0_5mm->SetLineWidth(2);
     TH1F* time_5_10mm = new TH1F("5/10mm","5/10mm", 2000, -100, +100);
     time_5_10mm->SetTitle(";Time (ns);Counts");
     time_5_10mm->GetXaxis()->SetLabelFont(22);
     time_5_10mm->GetXaxis()->SetTitleFont(22);
     time_5_10mm->GetYaxis()->SetLabelFont(22);
     time_5_10mm->GetYaxis()->SetTitleFont(22);
     time_5_10mm->SetLineWidth(2);
     TH1F* time_10_15mm = new TH1F("10/15mm","10/15mm", 2000, -100, +100);
     time_10_15mm->SetTitle(";Time (ns);Counts");
     time_10_15mm->GetXaxis()->SetLabelFont(22);
     time_10_15mm->GetXaxis()->SetTitleFont(22);
     time_10_15mm->GetYaxis()->SetLabelFont(22);
     time_10_15mm->GetYaxis()->SetTitleFont(22);
     time_10_15mm->SetLineWidth(2);
     TH1F* time_15_20mm = new TH1F("15/20mm","15/20mm", 2000, -100, +100);
     time_15_20mm->SetTitle(";Time (ns);Counts");
     time_15_20mm->GetXaxis()->SetLabelFont(22);
     time_15_20mm->GetXaxis()->SetTitleFont(22);
     time_15_20mm->GetYaxis()->SetLabelFont(22);
     time_15_20mm->GetYaxis()->SetTitleFont(22);
     time_15_20mm->SetLineWidth(2);
     TH1F* time_20_25mm = new TH1F("20/25mm","20/25mm", 2000, -100, +100);
     time_20_25mm->SetTitle(";Time (ns);Counts");
     time_20_25mm->GetXaxis()->SetLabelFont(22);
     time_20_25mm->GetXaxis()->SetTitleFont(22);
     time_20_25mm->GetYaxis()->SetLabelFont(22);
     time_20_25mm->GetYaxis()->SetTitleFont(22);
     time_20_25mm->SetLineWidth(2);
     TH1F* time_25_30mm = new TH1F("25/30mm","25/30mm", 2000, -100, +100);
     time_25_30mm->SetTitle(";Time (ns);Counts");
     time_25_30mm->GetXaxis()->SetLabelFont(22);
     time_25_30mm->GetXaxis()->SetTitleFont(22);
     time_25_30mm->GetYaxis()->SetLabelFont(22);
     time_25_30mm->GetYaxis()->SetTitleFont(22);
     time_25_30mm->SetLineWidth(2);
     std::cout<<"Z : "<<Tile.at(15).mZ/2000;
     std::cout<<"X : "<<Tile.at(15).mX/8000;
     if(dimension == ('Z' | 'z')){
         for( int i=0; i<Tile.size(); i++){
             if (Tile.at(i).mZ/2000 <= 5){
                 time_0_5mm->Fill(Tile.at(i).ltimestamp&0x00FFFFFF);
             }
             else if ((Tile.at(i).mZ/2000 > 5) && (Tile.at(i).mZ/2000 <= 10)){
                 time_5_10mm->Fill(Tile.at(i).ltimestamp&0x00FFFFFF);
             }
             else if ((Tile.at(i).mZ/2000 > 10) && (Tile.at(i).mZ/2000 <= 15)){
                 time_10_15mm->Fill(Tile.at(i).ltimestamp&0x00FFFFFF);
             }
             else if (Tile.at(i).mZ/2000 > 15){
                 time_15_20mm->Fill(Tile.at(i).ltimestamp&0x00FFFFFF);
             }
         }
         time_z->Divide(2,2);
         time_z->cd(1);
         time_0_5mm->Draw();
         time_z->cd(2);
         time_5_10mm->Draw();
         time_z->cd(3);
         time_10_15mm->Draw();
         time_z->cd(4);
         time_15_20mm->Draw();
     }
     else if (dimension == ('X' | 'x')){
         for( int i=0; i<Tile.size(); i++){
             if (Tile.at(i).mX/8000 <= 5){
                 time_0_5mm->Fill(Tile.at(i).ltimestamp&0x00FFFFFF);
             }
             else if ((Tile.at(i).mX/8000 > 5) && (Tile.at(i).mX/8000 <= 10)){
                 time_5_10mm->Fill(Tile.at(i).ltimestamp&0x00FFFFFF);
             }
             else if ((Tile.at(i).mX/8000 > 10) && (Tile.at(i).mX/8000 <= 15)){
                 time_10_15mm->Fill(Tile.at(i).ltimestamp&0x00FFFFFF);
             }
             else if (Tile.at(i).mX/8000 > 15 && (Tile.at(i).mX/8000 <= 20)){
                 time_15_20mm->Fill(Tile.at(i).ltimestamp&0x00FFFFFF);
             }
             else if (Tile.at(i).mX/8000 > 20 && (Tile.at(i).mX/8000 <= 25)){
                 time_20_25mm->Fill(Tile.at(i).ltimestamp&0x00FFFFFF);
             }
             else if (Tile.at(i).mX/8000 > 25 && (Tile.at(i).mX/8000 <= 30)){
                 time_25_30mm->Fill(Tile.at(i).ltimestamp&0x00FFFFFF);
             }
         }
         time_z->Divide(2,3);
         time_z->cd(1);
         time_0_5mm->Draw();
         time_z->cd(2);
         time_5_10mm->Draw();
         time_z->cd(3);
         time_10_15mm->Draw();
         time_z->cd(4);
         time_15_20mm->Draw();
         time_z->cd(5);
         time_20_25mm->Draw();
         time_z->cd(6);
         time_25_30mm->Draw();
     }
     return time_z;
 }

