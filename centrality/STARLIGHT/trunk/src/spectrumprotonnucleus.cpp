/*
    <one line to give the library's name and an idea of what it does.>
    Copyright (C) 2011  Oystein Djuvsland <oystein.djuvsland@gmail.com>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <cmath>
#include "spectrumprotonnucleus.h"
#include "beambeamsystem.h"
#include "beam.h"
#include <iostream>

spectrumProtonNucleus::spectrumProtonNucleus(beamBeamSystem *b) : spectrum(b)
{
      _bMin = 4.0;
}

bool spectrumProtonNucleus::generateBreakupProbabilities()
{
    int nbbins = _nBbins;

    double b_min = _bMin;
    //double b_max = _bMax;

    //    double binc = (log(b_max/b_min))/(double)nbbins;
    double binc = exp((log(_bMax/_bMin))/(double)_nBbins);

    double b = b_min;

    _probOfBreakup.resize(nbbins);

    for (int i = 0; i < nbbins; i++)
    {
        //double bimp = b;
        //double rhad = 0;
	_beamBeamSystem->beam1().Z() > 1 ? _probOfBreakup[i] = exp(-getNucleonNucleonSigma()*_beamBeamSystem->beam1().thickness(b)) :
	_beamBeamSystem->beam2().Z() > 1 ? _probOfBreakup[i] = exp(-getNucleonNucleonSigma()*_beamBeamSystem->beam2().thickness(b)) :
	b < 7.76 ? _probOfBreakup[i] = 0 : _probOfBreakup[i] = 1; // Should always be true though
	//double hardCutoff = 0;
	
        b = b*binc;
    }

    return true;
}

double spectrumProtonNucleus::getSigma(double ) const
{
    return 0.11;
}



