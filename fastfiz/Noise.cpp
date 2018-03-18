/***************************************************************************
 *   Copyright (C) 2009 by Alon Altman   *
 *   epsalon@stanford.edu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Library General Public License as       *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "Noise.h"
#include <iomanip>
#include <iostream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace std;

Pool::Noise *Pool::Noise::Factory(istream &sourceStream) {
  int noiseType;
  sourceStream >> noiseType;
  NoiseType noise_type((NoiseType)noiseType); // Should cast things
  switch ((int)noise_type) {
  case NT_NONE: {
    NoNoise *newNoise = new NoNoise(sourceStream);
    return newNoise;
    break;
  }

  case NT_GAUSSIAN: {
    GaussianNoise *newNoise = new GaussianNoise(sourceStream);
    return newNoise;
    break;
  }

  default:
    cerr << "Unidentified noise." << endl;
    cerr << "NoiseType: " << noise_type << endl;
    break;
  }
  return NULL;
}

string Pool::Noise::toString() {
  ostringstream output;
  toStream(output);
  return output.str();
}

void Pool::Noise::toStream(ostream &out) const { out << noiseType(); }

void Pool::GaussianNoise::applyNoise(ShotParams &sp) const {
  sp.a += gsl_ran_gaussian(Utils::rng(), a);
  sp.b += gsl_ran_gaussian(Utils::rng(), b);
  sp.theta += gsl_ran_gaussian(Utils::rng(), theta);
  sp.phi += gsl_ran_gaussian(Utils::rng(), phi);
  sp.v += gsl_ran_gaussian(Utils::rng(), v);
}

void Pool::GaussianNoise::importFromStream(istream &sourceStream) {
  sourceStream >> a >> b >> theta >> phi >> v;
}

void Pool::GaussianNoise::toStream(ostream &out) const {
  Noise::toStream(out);
  out << " " << setprecision(20) << a << " " << setprecision(20) << b << " "
      << setprecision(20) << theta << " " << setprecision(20) << phi << " "
      << setprecision(20) << v;
}

ostream &Pool::operator<<(ostream &os, const Noise &obj) {
  obj.toStream(os);
  return os;
}
