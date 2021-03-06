///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//    Added 18->19 for reading in the luminosity table
//    Incoherent factor added to table --Joey
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>

#include "readinluminosity.h"
#include "starlightconstants.h"
#include "inputParameters.h"

using namespace std;


//______________________________________________________________________________
readLuminosity::readLuminosity()//:inputread(input)
  : _Warray(0), _Yarray(0), _Farray(0), _Farray1(0), _Farray2(0)
{
  //storing inputparameters into protected variables for the object to use them
  _ReadInputNPT=inputParametersInstance.nmbPtBinsInterference();
  _ReadInputnumy=inputParametersInstance.nmbRapidityBins();
  _ReadInputnumw=inputParametersInstance.nmbWBins();
  _ReadInputgg_or_gP=inputParametersInstance.productionMode();
  _ReadInputinterferencemode=inputParametersInstance.interferenceEnabled();
  
}


//______________________________________________________________________________
readLuminosity::~readLuminosity()
{ 
  if(_Warray) delete [] _Warray;
  if(_Yarray) delete [] _Yarray;
  if(_Farray) delete [] _Farray;
  if(_Farray1) delete [] _Farray1;
  if(_Farray2) delete [] _Farray2; 
}


//______________________________________________________________________________
void readLuminosity::read()
{
  
  if(!_Warray) _Warray = new double[_ReadInputnumw];
  if(!_Yarray) _Yarray = new double[_ReadInputnumy];
  if(!_Farray) 
  {
    _Farray = new double*[_ReadInputnumw];
    for(int i = 0; i < _ReadInputnumw; i++)
    {
      _Farray[i] = new double[_ReadInputnumy];
    }
  }
  if(!_Farray1) 
  {
    _Farray1 = new double*[_ReadInputnumw];
    for(int i = 0; i < _ReadInputnumw; i++)
    {
      _Farray1[i] = new double[_ReadInputnumy];
    }
  }
  if(!_Farray2) 
  {
    _Farray2 = new double*[_ReadInputnumw];
    for(int i = 0; i < _ReadInputnumw; i++)
    {
      _Farray2[i] = new double[_ReadInputnumy];
    }
  }

  double dummy[19]; //14//18

  //  double (*finterm)[starlightLimits::MAXWBINS]=new double[starlightLimits::MAXWBINS][starlightLimits::MAXYBINS];  
  //decreased from 1000*1000; too big! causes fault!

  double fpart =0.;
  double fptsum=0.;
  ifstream wylumfile;

  _f_max=0.0;
  _f_max1=0.0;
  _f_max2=0.0;

  wylumfile.open("slight.txt");
  for(int i=0;i < 19;i++){ // was 14; this is to account for sergei's additional parameters ie d-Au//was19
    wylumfile >> dummy[i];
  }
  int A_1 = dummy[1];
  int A_2 = dummy[3];

  for(int i=0;i<_ReadInputnumw;i++){
    wylumfile >> _Warray[i];
  }
  for(int i=0;i<_ReadInputnumy;i++){
    wylumfile >> _Yarray[i];
  }

  if( (_ReadInputgg_or_gP == 1) || (A_2 == 1 && A_1 != 1) || (A_1 ==1 && A_2 != 1) ){ 
    for(int i=0;i<_ReadInputnumw;i++){
      for(int j=0;j<_ReadInputnumy;j++){
        wylumfile >> _Farray[i][j];
        if( _Farray[i][j] > _f_max ) _f_max=_Farray[i][j];
      }
    }
    //Normalize farray 
    for(int i=0;i<_ReadInputnumw;i++){
      for(int j=0;j<_ReadInputnumy;j++){
        _Farray[i][j]=_Farray[i][j]/_f_max;
      }
    }
  } else {
    for(int i=0;i<_ReadInputnumw;i++){
      for(int j=0;j<_ReadInputnumy;j++){
        wylumfile >> _Farray1[i][j];
	//        if( _Farray1[i][j] > _f_max ) _f_max=_Farray1[i][j];
      }
    }
    for(int i=0;i<_ReadInputnumw;i++){
      for(int j=0;j<_ReadInputnumy;j++){
        wylumfile >> _Farray2[i][j];
        if( _Farray1[i][j] + _Farray2[i][j] > _f_max ) _f_max=(_Farray1[i][j] + _Farray2[i][j]);
      }
    }
    //Normalize farray, farray1, farray2 
    for(int i=0;i<_ReadInputnumw;i++){
      for(int j=0;j<_ReadInputnumy;j++){
        _Farray1[i][j]=_Farray1[i][j]/_f_max;
        _Farray2[i][j]=_Farray2[i][j]/_f_max;
        _Farray[i][j]=_Farray1[i][j]+_Farray2[i][j];
      }
    }
  }

  if (_ReadInputgg_or_gP != 1 && _ReadInputinterferencemode != 0) {
        // only numy/2 y bins here, from 0 (not -ymax) to ymax
        double **finterm  = new double*[starlightLimits::MAXWBINS];
        for (int i = 0; i < starlightLimits::MAXWBINS; i++) finterm[i] = new double[starlightLimits::MAXYBINS];
        for (int i=0;i<_ReadInputnumy/2;i++) {
            //fmax=0;
            //we want to convert _fptarray to an integral array where fpt(i,j) is near 0, and fpt(j,NPT) ~1. This will facilitate a simple table loookup
            fptsum=0.;
            for (int j=0;j<_ReadInputNPT;j++) {
                wylumfile >> fpart;
                finterm[i][j] = fpart;
                _fptarray[i][j]=0.;
                fptsum=fptsum+fpart;
            }
            //convert array to integral
            _fptarray[i][0]=finterm[i][0]/fptsum;
            for (int j=1;j<_ReadInputNPT;j++) {
                for (int k=0;k<=j;k++) {
                    _fptarray[i][j]=_fptarray[i][j]+finterm[i][k];
                }
                _fptarray[i][j]=_fptarray[i][j]/fptsum;
            }
        }
        delete [] finterm;

    }
  wylumfile >> _bwnormsave;
  wylumfile.close();
  //8delete[] finterm;	
  return;
}
