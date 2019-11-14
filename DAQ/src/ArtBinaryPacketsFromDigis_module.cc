//
// EDProducer module for converting straw hit digis into DTC formatted
// packets (which are stored in DataBlock data products)
//
//
// Original author Tomonari Miyashita
//
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

#include <math.h>


// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Provenance.h"

// mu2e-artdaq-core includes
#include "mu2e-artdaq-core/Overlays/ArtFragment.hh"
#include "mu2e-artdaq-core/Overlays/ArtFragmentReader.hh"

//pci_linux_kernel_module includes
#include "dtcInterfaceLib/DTC.h"

// Mu2e includes.
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
//#include "DAQDataProducts/inc/DataBlockCollection.hh"

#include "SeedService/inc/SeedService.hh"

#include <fstream>
#include <stdexcept>


// Definitions and typedefs needed for compatibility with HLS codeblock
#define NUM_PRESAMPLES 4
#define START_SAMPLES 5 // 0 indexed
#define NUM_SAMPLES 15
#define LOWER_TDC 16000
#define UPPER_TDC 64000

typedef uint16_t tdc_type;
typedef uint8_t  tot_type;
typedef uint16_t adc_type;
typedef uint16_t calib_constant_type;
typedef uint8_t  flag_mask_type;
typedef uint64_t timestamp;

using namespace std;

using DataBlockHeader   = mu2e::ArtFragmentReader::DataBlockHeader;
using TrackerDataPacket = mu2e::ArtFragmentReader::TrackerDataPacket;
using  adc_t = mu2e::ArtFragment::adc_t;

 
namespace mu2e {
  constexpr int format_version = 6;
  enum  class PacketType : uint8_t {  DCSRequest = 0, Heartbeat = 1,  DataRequest = 2,   DCSReply = 3,   Dataheader = 5   };
 //--------------------------------------------------------------------
  //
  //
  class ArtBinaryPacketsFromDigis : public art::EDProducer {
  public:

    // using adc_t = mu2e::DataBlock::adc_t;
    // using dtc_id = mu2e::DataBlock::dtc_id;

    //data struct for the header, tracker, tacker and CRV
      
    typedef std::vector<TrackerDataPacket> TrackerDataPacketVector;
    typedef std::vector<DataBlockHeader>   DataBlockHeaderVector;

    explicit ArtBinaryPacketsFromDigis(fhicl::ParameterSet const& pset);

    virtual void beginJob() override;

    virtual void endJob();

    void produce( art::Event & ) override;

    flag_mask_type filter( TrackerDataPacket&   trkData, 
			   calib_constant_type clockstart, 
			   calib_constant_type panelTDCoffset, calib_constant_type hvoffset, calib_constant_type caloffset,
			   calib_constant_type energy_max_LSHIFT8, calib_constant_type energy_min_LSHIFT8,
			   calib_constant_type gain_RSHIFT15,
			   calib_constant_type inverse_ionization_energy_LSHIFT26);

  private:
    size_t _maxDMABlockSize;
    // Within each event (corresponding to a unique timestamp) the DataBlocks
    // are divided into Direct Memory Access (DMA) blocks with a max size
    // in bytes corresponding to _dmaBlockSize.
    // NOTE: THE DMA BLOCK SIZE INCLUDES THE DMA BLOCK HEADER !!!

    size_t _bufferSize;
    std::vector<adc_t> outputBuffer;

    size_t _generateTimestampTable;

    // Table used for mapping between DTC timestamp and art EventID
    std::string _tableFile;
    std::vector< std::pair<timestamp, timestamp> > tsTable;

    size_t  _timestampOffset;
    size_t  _numWordsWritten;
    size_t  _numEventsProcessed;
    int     _includeTracker;
    int     _includeCalorimeter;
    int     _includeCosmicRayVeto;
    
    int _includeDMAHeaders;

    // Set to 1 to save packet data to a binary file
    int _generateBinaryFile;

    string                _outputFile;
    ofstream              outputStream;
    
    //    const size_t number_of_rocs = 216;
    const size_t number_of_rocs = 240;
    const size_t number_of_straws_per_roc = 96; // Each panel in the tracker has 96 straws
    const size_t number_of_rocs_per_dtc = 6;
    const size_t numADCSamples = 15;

    // 96 straws per panel
    // 1 ROC per panel
    // 216 panels
    //
    // 6 ROCs per DTC
    // 36 DTCs

    int _generateTextFile;

    // Diagnostics level.
    int _diagLevel;

    int _enableFPGAEmulation;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Label of the module that made the hits.
    art::ProductToken<StrawDigiCollection> const  _sdtoken;
    art::ProductToken<CaloDigiCollection> const  _cdtoken;

    void   fillTrackerDataPacket(const StrawDigi& SD, TrackerDataPacket& TrkData);

    void   fillTrackerHeaderDataPacket(const StrawDigi& SD, DataBlockHeader& HeaderData, uint64_t &EventNum);

    void   fillEmptyTrackerDataPacket (TrackerDataPacket& TrkData);
    void   fillEmptyHeaderDataPacket  (DataBlockHeader& HeaderData, uint64_t& EventNum, uint8_t& ROCId, uint8_t& DTCId);
    //    void   fillCaloHeaderDataPacket(const CaloDigi&  SD, DataBlockHeader& HeaderData);

    void   processTrackerData(art::Event &evt, uint64_t& eventNum,
			      std::vector<TrackerDataPacketVector>& trkHitVectorByDTC, 
			      std::vector<DataBlockHeaderVector> & trkHeaderVectorByTDC,
			      size_t& curDMABlockSize,
			      size_t& numDataBlocksInCurDMABlock,
			      std::vector<size_t>& dataBlockPartition,
			      std::vector<size_t>& dataBlockPartitionSizes);
      
    
    void   fillTrackerDMABlocks(std::vector<TrackerDataPacketVector>& trkHitVectorByDTC, 
				std::vector<DataBlockHeaderVector>&   trkHeaderVectorByTDC,
				size_t& curDMABlockIdx, size_t& curDataBlockCount, size_t& curDMABlockByteCount,
				bool& atBeginningOfDMABlock, 
				std::vector<adc_t>& masterVector, std::vector<adc_t>& DMABlockVector,
				std::vector<adc_t>& headerlessMasterVector, std::vector<adc_t>& headerlessDMABlockVector,
				std::vector<size_t>& dataBlockPartition,
				std::vector<size_t>& dataBlockPartitionSizes);

    void   fillTrackerDataStream(std::vector<adc_t>&dataStream, TrackerDataPacket& curDataBlock, DataBlockHeader& headerDataBlock);


    void flushBuffer();
    
    std::vector<adc_t> generateDMABlockHeader (size_t theCount) const;
    std::vector<adc_t> generateEventByteHeader(size_t theCount) const;

    void   printTrackerHeader(DataBlockHeader& headerDataBlock);
    void   printTrackerData(TrackerDataPacket& curDataBlock);
  };
  
  void   ArtBinaryPacketsFromDigis::printTrackerHeader(DataBlockHeader& headerDataBlock){
    printf("[ArtBinaryPacketsFromDigis::printTrackerHEader] START header print  \n");
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] ByteCount      : %i \n", headerDataBlock.ByteCount    );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] Hopcount       : %i \n", headerDataBlock.Hopcount     );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] PacketType     : %i \n", headerDataBlock.PacketType   );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] ROCID 	   : %i \n", headerDataBlock.ROCID 	     );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] unused1	   : %i \n", headerDataBlock.unused1	     );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] SubsystemID    : %i \n", headerDataBlock.SubsystemID  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] Valid 	   : %i \n", headerDataBlock.Valid 	     );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] PacketCount    : %i \n", headerDataBlock.PacketCount  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] unused2        : %i \n", headerDataBlock.unused2      );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] TimestampLow   : %i \n", headerDataBlock.TimestampLow );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] TimestampMed   : %i \n", headerDataBlock.TimestampMed );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] TimestampHigh  : %i \n", headerDataBlock.TimestampHigh);
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] Status	   : %i \n", headerDataBlock.Status	     );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] FormatVersion  : %i \n", headerDataBlock.FormatVersion);
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] DTCID	   : %i \n", headerDataBlock.DTCID	     );
    printf("[ArtBinaryPacketsFromDigis::printTrackerHeader] EVBMode        : %i \n", headerDataBlock.EVBMode      );
  
  }
  
  void   ArtBinaryPacketsFromDigis::printTrackerData(TrackerDataPacket& trkData){
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] START tracker-data print \n");
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] StrawIndex    : %i \n", (int)trkData.StrawIndex	  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] TDC0		: %i \n", (int)trkData.TDC0		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] TDC1		: %i \n", (int)trkData.TDC1		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] TOT0		: %i \n", (int)trkData.TOT0		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] TOT1		: %i \n", (int)trkData.TOT1		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC00         : %i \n", (int)trkData.ADC00  		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC01  	: %i \n", (int)trkData.ADC01  		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC02  	: %i \n", (int)trkData.ADC02  		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC03  	: %i \n", (int)trkData.ADC03  		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC04  	: %i \n", (int)trkData.ADC04  		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC05A  	: %i \n", (int)trkData.ADC05A  	  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC05B 	: %i \n", (int)trkData.ADC05B 		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC06  	: %i \n", (int)trkData.ADC06  		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC07  	: %i \n", (int)trkData.ADC07  		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC08  	: %i \n", (int)trkData.ADC08  		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC09  	: %i \n", (int)trkData.ADC09  		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC10A 	: %i \n", (int)trkData.ADC10A 		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC10B  	: %i \n", (int)trkData.ADC10B  	  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC11  	: %i \n", (int)trkData.ADC11  		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC12  	: %i \n", (int)trkData.ADC12  		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC13 	: %i \n", (int)trkData.ADC13 		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] ADC14 	: %i \n", (int)trkData.ADC14 		  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] unused1 	: %i \n", (int)trkData.unused1 	  );
    printf("[ArtBinaryPacketsFromDigis::printTrackerData] PreprocessingFlags : %i \n", (int)trkData.PreprocessingFlags);

    
  }




  void   ArtBinaryPacketsFromDigis::fillTrackerDataStream(std::vector<adc_t>& dataStream, 
							  TrackerDataPacket&  trkDataBlock, 
							  DataBlockHeader&    headerDataBlock){
    
    auto sz = ceil((sizeof(TrackerDataPacket) + sizeof(DataBlockHeader))/sizeof(adc_t));
    //check that the trkDataBlock is not empty
    if (headerDataBlock.PacketCount == 0){
      sz = ceil( sizeof(DataBlockHeader)/sizeof(adc_t));
    }
    auto pos = dataStream.size();
    dataStream.resize(dataStream.size() + sz);

    memcpy(&dataStream[pos], &headerDataBlock, sizeof(DataBlockHeader));

    if (headerDataBlock.PacketCount == 2){
      pos += sizeof(DataBlockHeader) / sizeof(adc_t);
      memcpy(&dataStream[pos], &trkDataBlock, sizeof(TrackerDataPacket));
    }

    
  }

  

  void   ArtBinaryPacketsFromDigis::fillTrackerDMABlocks(std::vector<TrackerDataPacketVector>& trkHitVectorByDTC, 
							 std::vector<DataBlockHeaderVector>& trkHeaderVectorByTDC,
							 size_t& curDMABlockIdx, size_t& curDataBlockCount, size_t& curDMABlockByteCount,
							 bool& atBeginningOfDMABlock, 
							 std::vector<adc_t>& masterVector, std::vector<adc_t>& DMABlockVector,
							 std::vector<adc_t>& headerlessMasterVector, 
							 std::vector<adc_t>& headerlessDMABlockVector,
							 std::vector<size_t>& dataBlockPartition,
							 std::vector<size_t>& dataBlockPartitionSizes){
  
    for(size_t collectionIdx = 0; collectionIdx<trkHitVectorByDTC.size(); collectionIdx++) {
      TrackerDataPacketVector datablocks   = trkHitVectorByDTC   [collectionIdx];
      DataBlockHeaderVector   headerBlocks = trkHeaderVectorByTDC[collectionIdx];
      for(size_t dataBlockIdx = 0; dataBlockIdx<datablocks.size(); dataBlockIdx++) {
	TrackerDataPacket curDataBlock    = datablocks  [dataBlockIdx];
	DataBlockHeader   headerDataBlock = headerBlocks[dataBlockIdx];

	if(atBeginningOfDMABlock) {
	  atBeginningOfDMABlock = false;

	  DMABlockVector.clear();
	  headerlessDMABlockVector.clear();
	  
	  curDataBlockCount = 0;

	  if(_includeDMAHeaders>0) {
	    std::vector<adc_t> dma_header = generateDMABlockHeader( dataBlockPartitionSizes[curDMABlockIdx] );
	    DMABlockVector.insert(DMABlockVector.end(), dma_header.begin(), dma_header.end());

	    std::vector<adc_t> event_header = generateEventByteHeader( dataBlockPartitionSizes[curDMABlockIdx] - 8 - 8 );
	    DMABlockVector.insert(DMABlockVector.end(), event_header.begin(), event_header.end());	  
	    // The *exclusive* event byte count is equal to the *inclusive* DMA byte
	    // count minus the size of the DMA byte header (8 bytes)
	    // minus the size of the event byte count header (8 bytes)
	    // NOTE: Under the current specification, each DMA block contains only 1 event
	    
	    curDMABlockByteCount = 8 + 8;
	  }

	}

	// Add the current DataBlock to the current SuperBlock
	//curDataBlock.setTimestamp(ts); // Overwrite the timestamp
	
	//create a vector with all the adc_t values of the header and the trackerPacket
	std::vector<adc_t>  dataStream;

	fillTrackerDataStream(dataStream, curDataBlock, headerDataBlock);
	for(size_t adcNum = 0; adcNum < dataStream.size(); adcNum++) {
	  DMABlockVector.push_back(dataStream[adcNum]);
	  headerlessDMABlockVector.push_back(dataStream[adcNum]);
	}
	curDMABlockByteCount += dataStream.size() * sizeof(adc_t); 
	curDataBlockCount++;

	if ( _diagLevel > 1) {
	  if (dataBlockIdx==0){
	    std::cout << "================================================" << std::endl;
	    //std::cout << "\t\tTimestamp: " << ts << std::endl;
	    std::cout << "\t\tDTCID: " << (int)headerDataBlock.DTCID << std::endl;
	    std::cout << "\t\tSYSID: " << (int)headerDataBlock.SubsystemID << std::endl;
	  }
	  printTrackerHeader(headerDataBlock);
	  if (headerDataBlock.PacketCount>0) printTrackerData(curDataBlock);
	  
	}

	if(curDataBlockCount==dataBlockPartition[curDMABlockIdx]) {
	  // Reached end of current DMA block

	  if ( _diagLevel > 1 ) {
	    std::cout << "Number of bytes in DMABlock: " << curDMABlockByteCount << std::endl;
	  }

	  masterVector.insert(masterVector.end(), DMABlockVector.begin(), DMABlockVector.end());
	  headerlessMasterVector.insert(headerlessMasterVector.end(), headerlessDMABlockVector.begin(), headerlessDMABlockVector.end());
	  curDMABlockIdx++;

	  atBeginningOfDMABlock = true;
	}
	

      } // End loop over DataBlocks

    } // End loop over DTC collections
  }

  //--------------------------------------------------------------------------------
  //
  //--------------------------------------------------------------------------------
  void   ArtBinaryPacketsFromDigis::fillEmptyTrackerDataPacket (TrackerDataPacket& TrkData){
    bzero(&TrkData, sizeof(TrkData));
  }



  //--------------------------------------------------------------------------------  
  // 
  //--------------------------------------------------------------------------------
  void   ArtBinaryPacketsFromDigis::fillEmptyHeaderDataPacket(DataBlockHeader& headerData, uint64_t& EventNum, 
							      uint8_t& ROCId, uint8_t& DTCId){
  
    // Fill in the byte count field of the header packet
    // Word 0
    adc_t   nBytes = sizeof(DataBlockHeader);//ask Eric! //FIX ME!
    headerData.ByteCount = nBytes;
    // Word 1
    headerData.Hopcount = 0;//ask Eric!!!//FIX ME!
    headerData.PacketType = 5;//PacketType::Dataheader; 

    headerData.ROCID  = ROCId;
    headerData.unused1= 0;//ask Eric!
    
    headerData.SubsystemID = DTCLib::DTC_Subsystem_Tracker; //: 3;

    headerData.Valid    =  1;
    // Word 2
    headerData.PacketCount = 0;
    headerData.unused2 = 0;// : 5;
    // Word 3
    uint64_t timestamp = EventNum;
    headerData.TimestampLow  = static_cast<adc_t>(timestamp & 0xFFFF);
    // Word 4
    headerData.TimestampMed  = static_cast<adc_t>((timestamp >> 16) & 0xFFFF);
    // Word 5
    headerData.TimestampHigh = static_cast<adc_t>((timestamp >> 32) & 0xFFFF);
    // Word 6
    headerData.Status        = 0; // 0 corresponds to "TimeStamp had valid data"
    headerData.FormatVersion = format_version;
    // Word 7
    headerData.DTCID   = DTCId;
    uint8_t  evbMode = 0;//maybe off-spill vs on-spill?
    headerData.EVBMode = evbMode;
  }




  void   ArtBinaryPacketsFromDigis::fillTrackerHeaderDataPacket(const StrawDigi& SD, DataBlockHeader& HeaderData, 
								uint64_t& EventNum){
    // Word 0
    adc_t   nBytes = sizeof(DataBlockHeader) + sizeof(TrackerDataPacket);//ask Eric! //FIX ME!
    HeaderData.ByteCount = nBytes;
    // Word 1
    HeaderData.Hopcount = 0;//currently unused
    HeaderData.PacketType =  5;//PacketType::Dataheader;

    // 96 straws per ROC/panel
    // 6 panels / plane
    //      // 36 planes => 216 panels/ROCs
    // There are actually 40 planes in the input file => 240 panels/ROCs    
    int panel = SD.strawId().getPanel();
    int plane = SD.strawId().getPlane();

    // ROC ID, counting from 0 across all DTCs (for the tracker)
    //    uint8_t localROCID = panel;
    uint8_t globalROCID = (plane*6) + panel;

    // 240 ROCs total
    if(globalROCID >= number_of_rocs) {
      throw cet::exception("DATA") << " Global ROC ID " << globalROCID
				   << " exceeds limit of " << number_of_rocs;
    }
    HeaderData.ROCID  = panel;
    HeaderData.unused1= 0;
    HeaderData.SubsystemID = DTCLib::DTC_Subsystem_Tracker; 
    HeaderData.Valid    = 1 ;
    // Word 2
    HeaderData.PacketCount = 2;
    HeaderData.unused2     = 0;
    // Word 3
    uint64_t timestamp = EventNum;
    HeaderData.TimestampLow  = static_cast<adc_t>(timestamp & 0xFFFF);
    // Word 4
    HeaderData.TimestampMed  = static_cast<adc_t>((timestamp >> 16) & 0xFFFF);
    // Word 5
    HeaderData.TimestampHigh = static_cast<adc_t>((timestamp >> 32) & 0xFFFF);
    // Word 6
    HeaderData.Status        = 0; // 0 corresponds to "TimeStamp had valid data"
    HeaderData.FormatVersion = format_version;
    // Word 7
    HeaderData.DTCID   = static_cast<uint8_t>(globalROCID/number_of_rocs_per_dtc);
    uint8_t  evbMode = 0;//ask Eric
    HeaderData.EVBMode = evbMode;

  }

  void   ArtBinaryPacketsFromDigis::fillTrackerDataPacket(const StrawDigi& SD, TrackerDataPacket& TrkData){
    
    TrkData.StrawIndex = SD.strawId().asUint16();
    TrkData.TDC0       = SD.TDC(StrawEnd::cal);
    TrkData.TDC1       = SD.TDC(StrawEnd::hv);
    TrkData.TOT0       = SD.TOT(StrawEnd::cal);
    TrkData.TOT1       = SD.TOT(StrawEnd::hv);
    
    TrkTypes::ADCWaveform const& theWaveform = SD.adcWaveform();
    
    TrkData.ADC00   = (theWaveform[0] & 0xFFF);  // 12b
    TrkData.ADC01   = (theWaveform[1] & 0xFFF);  // 24b
    TrkData.ADC02   = (theWaveform[2] & 0xFFF);  // 36b
    TrkData.ADC03   = (theWaveform[3] & 0xFFF);  // 48b
    TrkData.ADC04   = (theWaveform[4] & 0xFFF);  // 60b
    TrkData.ADC05A  = (theWaveform[5] & 0xF);  // 64b
    TrkData.ADC05B  = ((theWaveform[5] >> 4) & 0xFF);  // 8b
    TrkData.ADC06   = (theWaveform[6] & 0xFFF);  // 20b
    TrkData.ADC07   = (theWaveform[7] & 0xFFF);  // 32b
    TrkData.ADC08   = (theWaveform[8] & 0xFFF);  // 44b
    TrkData.ADC09   = (theWaveform[9] & 0xFFF);  // 56b
    TrkData.ADC10A  = (theWaveform[10] & 0xFF);  // 64b
    TrkData.ADC10B  = ((theWaveform[10] >> 8) & 0xF);  //  4b
    TrkData.ADC11   = (theWaveform[11]& 0xFFF); // 16b
    TrkData.ADC12   = (theWaveform[12]& 0xFFF); // 28b
    TrkData.ADC13   = (theWaveform[13]& 0xFFF); // 40b
    TrkData.ADC14   = (theWaveform[14]& 0xFFF); // 52b
    
    TrkData.PreprocessingFlags  = 0;

    //    TrkData.unused1 = ;
    if(_enableFPGAEmulation) {
      // Note: Eventually, there will be a calibration database to provide these sorts of values
      // For now, these are just placeholders
      calib_constant_type clockstart = 320; //10 / .03125 (tdclsb)
      calib_constant_type maxenergyls8 = 583, minenergyls8 = 0;
      calib_constant_type gainrs15 = 1389;
      calib_constant_type inverseionizationenergyls26 = 633;
      calib_constant_type panelTDCoffset = 0,  hvoffset = 0,  caloffset = 0;


      flag_mask_type f = filter(TrkData,
				clockstart,panelTDCoffset,hvoffset,caloffset,
				maxenergyls8,minenergyls8,gainrs15,inverseionizationenergyls26);

      TrkData.PreprocessingFlags  =  f & 0xFF;
    }
  }


  ArtBinaryPacketsFromDigis::ArtBinaryPacketsFromDigis(fhicl::ParameterSet const& pset):
    art::EDProducer{pset},
    _maxDMABlockSize        (pset.get<size_t>("maxDMABlockSize",32000)), // Maximum size in bytes of a DMA block
    _bufferSize             (pset.get<size_t>("bufferSize",1000)),
    _generateTimestampTable (pset.get<size_t>("generateTimestampTable",0)),
    _tableFile              (pset.get<std::string>("tableFile","tsTable.bin")),
    _timestampOffset        (pset.get<size_t>("timestampOffset",0)),
    _numWordsWritten        (0),
    _numEventsProcessed     (0),
    _includeTracker         (pset.get<int>("includeTracker",1)),
    _includeCalorimeter     (pset.get<int>("includeCalorimeter",1)),
    _includeCosmicRayVeto   (pset.get<int>("includeCosmicRayVeto",0)),
    _includeDMAHeaders      (pset.get<int>("includeDMAHeaders",1)),
    _generateBinaryFile     (pset.get<int>("generateBinaryFile",1)),
    _outputFile             (pset.get<string>("outputFile","DTC_packets.bin")),
    _generateTextFile       (pset.get<int>("generateTextFile",0)),
    _diagLevel              (pset.get<int>("diagLevel",0)),
    _enableFPGAEmulation    (pset.get<int>("enableFPGAEmulation",0)),
    _maxFullPrint           (pset.get<int>("maxFullPrint",5)),
    _sdtoken{consumes<StrawDigiCollection>(pset.get<std::string>("StrawDigiCollection","makeSD"))},
    _cdtoken{consumes<CaloDigiCollection> (pset.get<std::string>("CaloDigiCollection","CaloDigiFromShower"))}{
      
      produces<timestamp>();
      produces< std::vector<adc_t> >();
      
      if(_generateBinaryFile == 1) {
	outputStream.open(_outputFile, std::ios::out | std::ios::binary);
      }
    }

  void ArtBinaryPacketsFromDigis::beginJob(){

    if ( _diagLevel > 0 ) {
      cout << "ArtBinaryPacketsFromDigis Diaglevel: "
           << _diagLevel << " "
           << _maxFullPrint
           << endl;

      if(_enableFPGAEmulation) {
	cout << "FPGA Emulation Enabled" << endl;
      }

    }

  }

  void ArtBinaryPacketsFromDigis::endJob(){
    if( _generateBinaryFile == 1 ) {
      flushBuffer();
      outputStream.close();
    }

    if(_generateTimestampTable) {
      std::ofstream tsTableStream;
      tsTableStream.open(_tableFile, std::ios::out | std::ios::binary);
      for(size_t idx=0; idx<tsTable.size(); idx++) {
	tsTableStream.write(reinterpret_cast<const char *>(&(tsTable[idx].first)), sizeof(timestamp));
	tsTableStream.write(reinterpret_cast<const char *>(&(tsTable[idx].second)), sizeof(timestamp));

	if (_diagLevel > 3) {
	  std::cout << "TIMESTAMP_MAPPING: timestamp: "
		    << tsTable[idx].first
		    << " uniqueid: "
		    << tsTable[idx].second
		    << std::endl;
	}

      }
      tsTableStream << std::flush;
      tsTableStream.close();
    }

    if (_diagLevel > 0) {
      std::cout << "BinaryPacketsFromDataBlocks: "
		<< "Finished writing "
		<< _numWordsWritten
		<< " words from "
		<< _numEventsProcessed
		<< " events to "
		<< _outputFile
		<< std::endl;
    }

  }

  
  void ArtBinaryPacketsFromDigis::flushBuffer() {
    
    for(size_t idx = 0; idx<outputBuffer.size(); idx++) {
      outputStream.write(reinterpret_cast<const char *>(&(outputBuffer[idx])), sizeof(adc_t));
    }
    outputStream << std::flush;

    _numWordsWritten+=outputBuffer.size()*2;
    outputBuffer.clear();
  }


  std::vector<adc_t> ArtBinaryPacketsFromDigis::generateDMABlockHeader(size_t theCount) const {

    uint64_t byteCount = theCount;

    std::vector<adc_t> header;
    header.push_back(static_cast<adc_t>( byteCount        & 0xFFFF));
    header.push_back(static_cast<adc_t>((byteCount >> 16) & 0xFFFF));
    header.push_back(static_cast<adc_t>((byteCount >> 32) & 0xFFFF));
    header.push_back(static_cast<adc_t>((byteCount >> 48) & 0xFFFF));
  
    return header;
  }

  std::vector<adc_t> ArtBinaryPacketsFromDigis::generateEventByteHeader(size_t theCount) const {

    uint64_t byteCount = theCount;

    std::vector<adc_t> header;
    header.push_back(static_cast<adc_t>( byteCount        & 0xFFFF));
    header.push_back(static_cast<adc_t>((byteCount >> 16) & 0xFFFF));
    header.push_back(static_cast<adc_t>((byteCount >> 32) & 0xFFFF));
    header.push_back(static_cast<adc_t>((byteCount >> 48) & 0xFFFF));
  
    return header;
  }


  void ArtBinaryPacketsFromDigis::produce(art::Event & evt) {

    // unique_ptr<DataBlockCollection> dtcPackets(new DataBlockCollection);

    uint64_t eventNum = evt.id().event();//is not unique! internal counter???//FIXME!
    uint64_t ts = _numEventsProcessed + _timestampOffset;

    if ( _diagLevel > 2 ) {
      cout << "ArtBinaryPacketsFromDigis: eventNum: " << eventNum << endl;
    }


    // Determine how to divide DataBlocks between DMABlocks within each timestamp
    std::vector<size_t> dataBlockPartition;
    std::vector<size_t> dataBlockPartitionSizes;

    size_t   curDMABlockSize(0);
    size_t   numDataBlocksInCurDMABlock(0);
 
    std::vector<TrackerDataPacketVector> trkHitVectorByDTC;   
    std::vector<DataBlockHeaderVector>   trkHeaderVectorByDTC;
 

    if(_includeTracker>0) {
      processTrackerData(evt, eventNum,
			 trkHitVectorByDTC, trkHeaderVectorByDTC, 
			 curDMABlockSize, numDataBlocksInCurDMABlock,
			 dataBlockPartition, dataBlockPartitionSizes);
    }
    
    // if (_includeCalorimeter>0){
    //   processCalorimterData();
    // }

    // Break the DataBlocks into DMABlocks and add DMABlock headers

    // Index of the current DMA block in the partition array
    size_t curDMABlockIdx = 0;

    // Number of DataBlocks added to the current DMA block
    size_t curDataBlockCount = 0; 
    // Number of DataBlocks added to the current DMA block
    size_t curDMABlockByteCount = 0; 

    bool atBeginningOfDMABlock = true;

    std::vector<adc_t> masterVector;
    std::vector<adc_t> DMABlockVector;

    std::vector<adc_t> headerlessMasterVector;
    std::vector<adc_t> headerlessDMABlockVector;

    if(_includeTracker>0) {
      
      fillTrackerDMABlocks(trkHitVectorByDTC, trkHeaderVectorByDTC,
			   curDMABlockIdx, curDataBlockCount, curDMABlockByteCount,
			   atBeginningOfDMABlock, 
			   masterVector, DMABlockVector,
			   headerlessMasterVector, headerlessDMABlockVector,
			   dataBlockPartition, dataBlockPartitionSizes);
    }

    // Write all values, including superblock header and DMA header values, to output buffer
    if( _generateBinaryFile == 1) {
      for ( size_t idx=0; idx<masterVector.size(); idx++ ) {
	if(outputBuffer.size()>= _bufferSize) {
	  flushBuffer();
	}
	outputBuffer.push_back(masterVector[idx]);
      }
    }

    _numEventsProcessed += 1;

    // Store the timestamp and DataBlockCollection in the event
    evt.put(std::unique_ptr<timestamp>(new timestamp( ts )));
    evt.put(std::make_unique< std::vector<adc_t> >(headerlessMasterVector));

  } // end of ::produce

  void ArtBinaryPacketsFromDigis::processTrackerData(art::Event &evt, uint64_t& eventNum,
						     std::vector<TrackerDataPacketVector>& trkHitVectorByDTC, 
						     std::vector<DataBlockHeaderVector> & trkHeaderVectorByDTC,
						     size_t& curDMABlockSize,
						     size_t& numDataBlocksInCurDMABlock,
						     std::vector<size_t>& dataBlockPartition,
						     std::vector<size_t>& dataBlockPartitionSizes){
    
    // bool gblresult;
    // art::Handle<StrawDigiCollection> strawDigisHandle;
    // gblresult = evt.getByLabel(_sdToken, strawDigisHandle);
    // if (!gblresult) throw cet::exception("DATA") << " Missing digi data";

    // StrawDigiCollection const& hits_SD =  *strawDigisHandle;
    
    auto  const &sdH = evt.getValidHandle(_sdtoken);
    const StrawDigiCollection& hits_SD(*sdH);

    // std::vector<TrackerDataPacket> trkHitVector;   // Vector of trk hit digi data
    // std::vector<DataBlockHeader>   trkHeaderVector;// Vector of trk hit header data
    std::vector<TrackerDataPacket> trkHitVector;   // Vector of trk hit digi data
    std::vector<DataBlockHeader>   trkHeaderVector;// Vector of trk hit header data
   
    for ( size_t i=0; i<hits_SD.size(); ++i ) {
      StrawDigi const& SD = hits_SD.at(i);

      // Fill struct with info for current hit
      TrackerDataPacket trkData;
      fillTrackerDataPacket(SD, trkData);
      DataBlockHeader   headerData;
      fillTrackerHeaderDataPacket(SD, headerData, eventNum);

      trkHitVector.push_back(trkData);
      trkHeaderVector.push_back(headerData);
    }
    
    if ( _diagLevel > 1 ) {
      std::cout << "[ArtBinaryPacketsFromDigis::processTrackerData ] Total number of tracker non-empty DataBlocks = " << 
	trkHitVector.size() << std::endl;
    }

    uint8_t max_dtc_id = number_of_rocs/number_of_rocs_per_dtc-1;
    if(number_of_rocs % number_of_rocs_per_dtc > 0) {
      max_dtc_id += 1;
    }
    
    // Loop over the DTC/ROC pairs and generate datablocks for each ROC
    for(uint8_t dtcID = 0; dtcID < max_dtc_id; dtcID++) {
      std::vector<TrackerDataPacket> currTrkHitVector;   // Vector of trk hit digi data
      std::vector<DataBlockHeader>   currTrkHeaderVector;// Vector of trk hit header data

      for(uint8_t rocID = 0; rocID < number_of_rocs_per_dtc; ++rocID) {
	// Find all hits for this event coming from the specified DTC/ROC combination
	for (size_t curHitIdx = 0; curHitIdx < trkHeaderVector.size(); curHitIdx++) {
	  if (trkHeaderVector[curHitIdx].DTCID == dtcID && 
	      trkHeaderVector[curHitIdx].ROCID == rocID ) {
	    currTrkHitVector.push_back(trkHitVector[curHitIdx]);
	    currTrkHeaderVector.push_back(trkHeaderVector[curHitIdx]);
	  }
	}

	if (currTrkHitVector.size() == 0) {
	  // No hits, so just fill a header packet and no data packets
	  DataBlockHeader   headerData;
	  TrackerDataPacket trkData;

	  fillEmptyHeaderDataPacket(headerData, eventNum, rocID, dtcID);
	  fillEmptyTrackerDataPacket(trkData);

	  currTrkHitVector.push_back(trkData);
	  currTrkHeaderVector.push_back(headerData);
	}
      } //Done looping over the ROCs in a given DTC
      
      trkHitVectorByDTC   .push_back(currTrkHitVector);
      trkHeaderVectorByDTC.push_back(currTrkHeaderVector);
    }


    // Loop again to evaluate the space of each block
    for(uint8_t dtcID = 0; dtcID < trkHeaderVectorByDTC.size(); dtcID++) {
      std::vector<TrackerDataPacket> currTrkHitVector    = trkHitVectorByDTC[dtcID];   // Vector of trk hit digi data
      std::vector<DataBlockHeader>   currTrkHeaderVector = trkHeaderVectorByDTC[dtcID];// Vector of trk hit header data

      //determine how to divide thedata blocks between DMVAblocks within each timestamp
      for (size_t dataBlockIdx=0; dataBlockIdx<currTrkHeaderVector.size(); ++dataBlockIdx){
	// TrackerDataPacket trkData    = currTrkHitVector[i];
	DataBlockHeader   headerData = currTrkHeaderVector[dataBlockIdx];

	if(numDataBlocksInCurDMABlock == 0) {
	  // Starting a new DMA Block, so allocate
	  // space for a new DMA block header (64 bits)
	  // and a new event byte count header (64 bits)
	  curDMABlockSize = 8 + 8;
	}

	numDataBlocksInCurDMABlock++; // Increment number of DataBlocks in the current DMA block
	curDMABlockSize += headerData.ByteCount * 2; // Size of current data block in 8bit words

	if(curDMABlockSize > _maxDMABlockSize) {
	  throw cet::exception("DATA") << "Current DMA Block size (" 
				       << curDMABlockSize << ") exceeds max DMA Block size ("
				       << _maxDMABlockSize << ")" << std::endl;
	}

	bool atEndOfDMABlock = false;

	if(dataBlockIdx == currTrkHeaderVector.size() - 1 &&
	   dtcID == trkHitVectorByDTC.size() - 1) {
	  // There are no more DataBlocks, so this marks the end of
	  // the current DMA block
	  atEndOfDMABlock = true;
	} 

	if( dataBlockIdx == currTrkHeaderVector.size() - 1 &&
	    dtcID < trkHitVectorByDTC.size() - 1 &&
	    trkHitVectorByDTC[dtcID+1].size()>0 &&
	    curDMABlockSize + (2 * trkHeaderVectorByDTC[dtcID+1][0].ByteCount) > _maxDMABlockSize) {
	  // Adding the DataBlock would go over the size limit
	  // so this marks the end of the current DMA block
	  atEndOfDMABlock = true;
	} 	 

	if(dataBlockIdx < currTrkHeaderVector.size() - 1 &&
	   curDMABlockSize + (2 * currTrkHeaderVector[dataBlockIdx + 1].ByteCount) > _maxDMABlockSize) {
	  // Adding the next DataBlock would put us over the limit so this is
	  // the end of the current DMA block
	  atEndOfDMABlock = true;
	}

	if(atEndOfDMABlock) {
	  dataBlockPartition.push_back(numDataBlocksInCurDMABlock);
	  dataBlockPartitionSizes.push_back(curDMABlockSize);
	  numDataBlocksInCurDMABlock = 0;
	}
	
      }
     
    } // Done looping of DTC/ROC pairs
  }



  flag_mask_type ArtBinaryPacketsFromDigis::filter( TrackerDataPacket& TrkData,
						    calib_constant_type clockstart, 
						    calib_constant_type panelTDCoffset, calib_constant_type hvoffset, calib_constant_type caloffset,
						    calib_constant_type energy_max_LSHIFT8, calib_constant_type energy_min_LSHIFT8,
						    calib_constant_type gain_RSHIFT15,
						    calib_constant_type inverse_ionization_energy_LSHIFT26) {

    int failed_time = 0;
    int failed_energy = 0;
    int thistdc = TrkData.TDC0 < TrkData.TDC1 ? TrkData.TDC0 + hvoffset : TrkData.TDC1 + caloffset;
    thistdc = thistdc + clockstart + panelTDCoffset;
    if (thistdc < LOWER_TDC || thistdc > UPPER_TDC) {
      failed_time = 1;
    }
    std::array<adc_t, 15> const waveform = TrkData.Waveform();
    //sum up presamples
    int pedsum = 0;
    for (int i = 0; i < NUM_PRESAMPLES; i++) {
      pedsum += waveform[i];
    }
    int pedestal = pedsum / NUM_PRESAMPLES;
    
    int peak = 0;
    for (int i = START_SAMPLES; i < NUM_SAMPLES; i++) {
      if (TrkData.Waveform().at(i) > peak) {
	peak = waveform[i];
      } else {
	break;
      }
    }
    
    int energy = peak - pedestal;
    
    int energy_max_adjusted = ((((energy_max_LSHIFT8 * gain_RSHIFT15) >> 9) * inverse_ionization_energy_LSHIFT26) >> 10);
    int energy_min_adjusted = ((((energy_min_LSHIFT8 * gain_RSHIFT15) >> 9) * inverse_ionization_energy_LSHIFT26) >> 10);
    if (energy > energy_max_adjusted || energy < energy_min_adjusted) {
      failed_energy = 1;
    }
    
    return (failed_energy<<1) | failed_time;
    
  }
  
}


using mu2e::ArtBinaryPacketsFromDigis;
DEFINE_ART_MODULE(ArtBinaryPacketsFromDigis);
