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

using namespace std;


/** ray_number_events_pixel function
 * @brief: Faire sortir les nombres des événements dans chaque pixel
 * @return: Résultat sous forme d’un vecteur à une dimension=nombre de pixel choisir 
 * @note:
 * */
//============================================================================
 vector <double> ray_number_events_pixel(std::vector<TuileEvt> TileA, int n_p_a,int p)
{

    double x_min=500;
	double x_max=8500;
	double x=0;
	double y=0;
	//double p=*p;
	double h=(x_max-x_min)/(n_p_a);
	double sum_event;
	vector <double> ray_numbre_event; 
    int r=0;
    for(int i=0;i<n_p_a;i++)
	{  
          for(int j=0;j<n_p_a;j++)
            {
                    sum_event=0;
                    
                    for(int k=0;k<TileA.size();k++)                    
                    {
                            if(p==0)                                    
                            {
                                            x=TileA.at(k).mX;
                                            y=TileA.at(k).mY;                                   
            
                                            if(((x>(x_min+h*i))&&(x<=(x_min+h*(i+1))))&&((y>(x_min+h*j))&&(y<=(x_min+h*(j+1)))))
                                            {				 
                                                    sum_event+=1;   					 							    
                                             }                          
         
                            }else{
                            
                                    if ((TileA.at(k).photons >= p*0.9)&&(TileA.at(k).photons <= p*1.1))// if ((k>=3850*0.9)&&(k<=3850*1.1))
                                    {                                       
                                            x=TileA.at(k).mX;
                                            y=TileA.at(k).mY;                                   
            
                                            if(((x>(x_min+h*i))&&(x<=(x_min+h*(i+1))))&&((y>(x_min+h*j))&&(y<=(x_min+h*(j+1)))))
                                             {				 
                                                   sum_event+=1;   					 							    
                                             }
                                    }                            
                            
                            }                                           
                     } 
                     
                
              ray_numbre_event.push_back(sum_event);                         
             //cout<<"Pixel_511{ "<<i+1<<" , "<<j+1<<" }  "<<sum_event<<endl;
                                           
            }     
    }   
  
return ray_numbre_event;     

}
//============================================================================
/** print_event_vect function
 * @brief: voir les nombres des événements dans chaque pixel 
 * @return: Résultat sous forme d’un vecteur à une dimension=nombre de pixel choisir 
 * @note:
 * */
void print_event_vect(vector<double> t)
{
     string const nomFichier("/home/hjira/calibration/output/number_events_pixel.txt");
     ofstream monFlux(nomFichier.c_str());
     if (monFlux){   
        
          for(int k=0;k<t.size();k++)
          {
            monFlux << t[k]<< endl; 
              //
              cout<<"Pixel { "<<k+1<<" }   "<<t[k]<<endl;           
          
          }                                                                                   
                                                                                                                               
    }
    monFlux.close();       
}    
//============================================================================
/** ray_number_events_pixel_unit function
 * @brief: Normaliser un vecteur de n dimension  
 * @return: Résultat sous forme d’un vecteur à une dimension=nombre de pixel choisir 
 * @note:
 * */
vector <double> ray_number_events_pixel_unit(vector<double> t)
{       
        double m=*max_element(t.begin(),t.end()); 
        vector<double> x;
        for(int i=0;i<t.size();i++)
        { 
          x.push_back(t[i]/m);                    
        }      
    return x; 
}
//============================================================================
/** profile_1d_vect function
 * @brief: Faire la somme des évènements pour chaque ligne d’une matrice 2D    
 * @return: Résultat sous forme d’un vecteur à une dimension= 1/2 nombre de pixel choisir 
 * @note:
 * */
vector <double> profile_1d_vect(vector<double> t)
{

         int size=sqrt(t.size());
            double sum=0;
            double tt[size][size];
            int k=0;
            int l=0;
            int m=0;
            
        vector<double> x;
            
         for (int i=0;i<size;i++)
            {   l+=1;
                sum=0;
                for (int j=0;j<size;j++)
                {
                    tt[i][j]=t[k++];
                    sum=sum+tt[i][j];
                          
                }
               x.push_back(sum);  
                
            }     
       
    return x; 
}

//============================================================================
/** print_profile_1d_vect function
 * @brief: voir la somme des événements dans chaque dans chaque ligne d’une matrice 2D    
 * @return: Résultat sous forme d’un vecteur à une dimension= 1/2 nombre de pixel choisir 
 * @note:
 * */
void print_profile_1d_vect(vector<double> t)
{
    string const nomFichier("/home/hjira/calibration/output/profile_1d.txt");
    ofstream monFlux(nomFichier.c_str());
    if (monFlux){ 
    
	                  for(int i=0;i<t.size();i++)
                        { 
                         monFlux<<t[i]<<endl;               
                         cout<<i+1<<" "<<t[i]<<endl;
                        }                                   
    }
     monFlux.close();       
}

//============================================================================
/** ray_standard_deviation function
 * @brief: Permettre de calculer la moyenne et l’écart type pour une carte de correction (la somme de n’acquisitions)   
 * @return: Sortie le rapport de l’écart type sur la valeur moyenne pour chaque pixel 
 * @note:
 * */
vector <double> ray_statistic(vector<double> t, vector<double> tt, double a_n)
{
        vector<double> mean_value;
        vector<double> standard_deviation;
        vector<double> standard_deviation_by_mean_value;
        int size=t.size();
        double mv;
        double ratio;
       // cout<<"a_n = "<<a_n<<endl;
        for(int i=0;i<t.size();i++)
        { 
                 
                 mv=t[i]/(a_n);
      
                 mean_value.push_back(mv);  
          
          
       //cout<<"i= "<<i+1<<" : valeur mean = "<<mv<<endl;  
        } 
         
        double sum=0;
        int m=0;    
        for(int j=0;j<t.size();j++)
        { 
                for(int l=0;l<a_n;l++)
                {         
                        
                sum+=(tt[j+(size*l+0)]-mean_value[j])*(tt[j+(size*l+0)]-mean_value[j]);
                
       //  cout<<"j= "<<j+1<<" : sum "<<sum<<endl; 
                }
             //cout<<"j= "<<j+1<<" : sum "<<sum<<endl; 
                sum=sum/a_n;
                sum=sqrt(sum);
                standard_deviation.push_back(sum);       
      // cout<<"j= "<<j+1<<" : standard_deviation = "<<sum<<endl; 
                if (mean_value[j]==0)
                {       
                            ratio=0;
                    }else{
                            ratio=standard_deviation[j]/mean_value[j];
                 }
            
         standard_deviation_by_mean_value.push_back(ratio);  
         
         //cout<<"j= "<<j+1<<" standard_deviation_by_mean_value = "<<standard_deviation_by_mean_value[j]<<endl;     
         //cout<<"j= "<<j+1<<" : standard_deviation_by_mean_value = "<<standard_deviation_by_mean_value[j]<<endl;                    
        }      
        
        
    return standard_deviation_by_mean_value; 
}

//============================================================================
/** print_event_vect function
 * @brief: voir le types de statistique indiquer dans la fonction ray_statistics dans chaque pixel 
 * @return: Résultat sous forme d’un vecteur à une dimension=nombre de pixel choisir 
 * @note:
 * */
void print_event_vect_statistic(vector<double> t)
{
     string const nomFichier1("/home/hjira/calibration/output/statistic.txt");
     ofstream monFlux1(nomFichier1.c_str());
     if (monFlux1){   
        
          for(int k=0;k<t.size();k++)
          {
            monFlux1 << t[k]<< endl; 
              //
              cout<<"Pixel { "<<k+1<<" }   "<<t[k]<<endl;           
          
          }                                                                                                                                                                                                                  
    }
    monFlux1.close();       
}   
//============================================================================
/**Map_number_events_pixe function
 * @brief: voir la carte des évenements en utilisant un histogramme 2D   
 * @return: 
 * @note:
 * */
TH2D* Map_number_events_pixel(vector<double> t, TH2D* Map2D)
{   
   int size=sqrt(t.size());
   int k=0;
   for(int i=0;i<size;++i)
    {
       for(int j=0;j<size;++j)
         {
            Map2D->Fill(i+1,j+1,t[k++]) ;
          } 
    }
     return Map2D;
}
//============================================================================
/**Map_Photons_Spectrum function
 * @brief: Voir le spectre total ou sélectionné pour une acquisition en utilisant un histogramme 1D     
 * @return: 
 * @note:
 * */
TH1F* Map_Photons_Spectrum(std::vector<TuileEvt> TileA, TH1F* Photons_Spectrum_Tile2, double p)// p=position
{   
    

   for(int k=0;k<TileA.size();++k)
    {
    
    
        if (p==0)
        {
             Photons_Spectrum_Tile2->Fill(TileA.at(k).photons);
            }else{
                if ((TileA.at(k).photons >= p*0.9)&&(TileA.at(k).photons <= p*1.1))// if ((k>=3850*0.9)&&(k<=3850*1.1))
                { 
                   Photons_Spectrum_Tile2->Fill(TileA.at(k).photons);
                }      
        }                
                        
    }
     return Photons_Spectrum_Tile2;
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

