#ifndef POOL_NOISE_H
#define POOL_NOISE_H

#ifdef SWIG
%{
#include "Noise.h"
%}
#endif /* SWIG */

#include "FastFiz.h"
#include <iostream>
#include <string>
using namespace std;


namespace Pool {

  /** Noise Types */
  enum NoiseType {NT_UNDEFINED,NT_NONE,NT_GAUSSIAN=2};
#ifndef SWIG
  inline istream & operator>> (istream & is, NoiseType nt)
  {
    int code;
    is >> code;
    nt = NoiseType(code);
    return is;
  }
#endif
  
  /** Base Noise class for implementing a generic type of noise */
  class Noise {
    public:
#ifndef SWIG
      /** Factory Function -- read noise info from a stream */
      static Noise* Factory(istream & sourceStream);
#endif
      /** Factory Function -- read noise info from a string */
      static Noise* Factory(string sourceString)
      {
              istringstream noiseSourceStream(sourceString);
              return Factory(noiseSourceStream);
      };
      
      
      ///////////////////////////////////////////////////
      // Input and output funzies.
      ///////////////////////////////////////////////////
      // Define friend operator functions
#ifndef SWIG
      friend ostream& operator<<(ostream& os, const Noise& obj);
      
      virtual void toStream(ostream & out) const ;
#endif
      
      virtual string toString();
      
      /** Return the type of noise used */
      virtual NoiseType noiseType() const  =0;
      
      /** Apply the noise to the given shot parameters.
       * @param sp Shot parameters -- will be changed.
       */
      virtual void applyNoise(ShotParams &sp) const = 0;

    protected:
      // Import from a stream
      virtual void importFromStream(istream & sourceStream) {};
      
  }; // Class Noise

#ifndef SWIG
  ostream& operator<<(ostream&os, const Noise& obj);
  
  /** Trivial class implementing null noise. */
  class NoNoise : public Noise
  {
    public:
      
      NoNoise () {
      }

      NoNoise (istream &sourceStream) {
        importFromStream(sourceStream);
      }
      
      // Parse common part of GameState from string
      NoNoise (string noiseString)
      {
        istringstream infoStream(noiseString);
        importFromStream(infoStream);
      };
      
      // Return type of game
      NoiseType noiseType() const  {return NT_NONE;};

      virtual void applyNoise(ShotParams &sp) const {};
  }; // Class NoNoise
  
  /** Class implementing Gaussian noise */
  class GaussianNoise : public Noise {
    public:
      
      GaussianNoise (double magnitude = 1.0) :
        a(0.5*magnitude), b(0.5*magnitude), theta(0.1*magnitude), phi(0.125*magnitude), v(0.075*magnitude)
      {};
      
      GaussianNoise (double _a, double _b, double _theta, double _phi, double _v) :
        a(_a), b(_b), theta(_theta), phi(_phi), v(_v) {};

      GaussianNoise (istream &sourceStream) {
        importFromStream(sourceStream);
      };
      
      // Parse common part of GameState from string
      GaussianNoise (string noiseString)
      {
        istringstream infoStream(noiseString);
        importFromStream(infoStream);
      };
      
      // Return type of game
      NoiseType noiseType() const  {return NT_GAUSSIAN;};

      virtual void toStream(ostream & out) const ;

      virtual void applyNoise(ShotParams &sp) const;
    protected:
      // Import from a stream
      virtual void importFromStream(istream & sourceStream);
    private:
      double a,b,theta,phi,v;
  };
#endif
  
} // Namespace pool

#endif
