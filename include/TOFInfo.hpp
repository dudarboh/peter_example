#ifndef TOFInfo_h
#define TOFInfo_h 1

#include <vector>
#include "marlin/Processor.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "UTIL/LCRelationNavigator.h"
#include "DDRec/Vector3D.h"
#include "DD4hep/Detector.h"

class TOFInfo : public marlin::Processor {
    public:
        marlin::Processor* newProcessor() { return new TOFInfo; }
        TOFInfo();
        void init();
        void processEvent(EVENT::LCEvent* evt);        

        void printOutputFCN(EVENT::ReconstructedParticle* pfo);
        void printOutputGNN(EVENT::ReconstructedParticle* pfo);
        void printOutputCNN(EVENT::ReconstructedParticle*);
        void drawPFO(EVENT::LCEvent* event, EVENT::ReconstructedParticle* pfo);

        dd4hep::rec::Vector3D getHelixMomAtTrackState(const EVENT::TrackState& ts);
        EVENT::MCParticle* getMaxTrackWeightMC(EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator nav);
        std::vector<EVENT::CalorimeterHit*> selectFrankEcalHits( EVENT::Cluster* cluster, EVENT::Track* track, int maxEcalLayer);
        double getTofClosest( EVENT::Cluster* cluster, EVENT::Track* track, double timeResolution);

    private:
        dd4hep::Detector& _detector = dd4hep::Detector::getInstance();
        int _nEvent{};
        double _bField{};
};

#endif
