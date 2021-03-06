#ifndef RecoDataProducts_CrvRecoPulse_hh
#define RecoDataProducts_CrvRecoPulse_hh
//
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include <vector>

namespace mu2e 
{
  class CrvRecoPulse
  {
    public:

    CrvRecoPulse() : _PEs(0), _PEsPulseHeight(0), _pulseTime(0), _pulseHeight(0), _pulseBeta(0), _pulseFitChi2(0), _LEtime(0), 
                     _scintillatorBarIndex(0), _SiPMNumber(0) {}

    CrvRecoPulse(int PEs, int PEsPulseHeight, double pulseTime, double pulseHeight, double pulseBeta, double pulseFitChi2, double LEtime, 
                 const std::vector<size_t> &waveformIndices, mu2e::CRSScintillatorBarIndex scintillatorBarIndex, int SiPMNumber) : 
                                                                            _PEs(PEs), 
                                                                            _PEsPulseHeight(PEsPulseHeight), 
                                                                            _pulseTime(pulseTime), 
                                                                            _pulseHeight(pulseHeight),
                                                                            _pulseBeta(pulseBeta),
                                                                            _pulseFitChi2(pulseFitChi2),
                                                                            _LEtime(LEtime),
                                                                            _waveformIndices(waveformIndices),
                                                                            _scintillatorBarIndex(scintillatorBarIndex),
                                                                            _SiPMNumber(SiPMNumber)
                                                                             {}

    int    GetPEs() const          {return _PEs;}
    int    GetPEsPulseHeight() const {return _PEsPulseHeight;}
    double GetPulseTime() const    {return _pulseTime;}
    double GetPulseHeight() const  {return _pulseHeight;}
    double GetPulseBeta() const    {return _pulseBeta;}
    double GetPulseFitChi2() const {return _pulseFitChi2;}
    double GetLEtime() const       {return _LEtime;}

    const std::vector<size_t>    &GetWaveformIndices() const      {return _waveformIndices;}
    mu2e::CRSScintillatorBarIndex GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    int                           GetSiPMNumber() const           {return _SiPMNumber;}

    private:

    int    _PEs;
    int    _PEsPulseHeight;  //used for PEs which were calculated using the pulse height and the pulse height calibration factor
    double _pulseTime;
    double _pulseHeight;
    double _pulseBeta;
    double _pulseFitChi2;
    double _LEtime;

    std::vector<size_t>            _waveformIndices;  //indices in the vector of the CrvDigiCollection (which is the same as the index in the CrvDigiMCCollection)
    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    int                            _SiPMNumber; 
  };
  typedef std::vector<mu2e::CrvRecoPulse> CrvRecoPulseCollection;
}

#endif /* RecoDataProducts_CrvRecoPulse_hh */
