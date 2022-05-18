#include "TOFInfo.hpp"

#include <limits>
#include "EVENT/LCCollection.h"
#include "UTIL/PIDHandler.h"

#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"
#include "marlin/VerbosityLevels.h"
#include "marlinutil/GeometryUtil.h"
#include "marlinutil/DDMarlinCED.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/Randomize.h"
#include "HelixClass.h"
#include "marlinutil/CalorimeterHitType.h"

using namespace EVENT;
using dd4hep::rec::Vector3D;
using std::cout, std::endl, std::vector, std::string;
using UTIL::LCRelationNavigator;
TOFInfo aTOFInfo ;

TOFInfo::TOFInfo() : marlin::Processor("TOFInfo") {}


void TOFInfo::init(){
    marlin::Global::EVENTSEEDER->registerProcessor(this);
    DDMarlinCED::init(this);
    _bField = MarlinUtil::getBzAtOrigin();
}


void TOFInfo::processEvent(EVENT::LCEvent * event){
    CLHEP::RandGauss::setTheSeed( marlin::Global::EVENTSEEDER->getSeed(this) );
    cout<<endl;
    cout<<"***********************************************"<<endl;
    cout<<"****************Event "<<(++_nEvent)<<"**************************"<<endl;
    cout<<"***********************************************"<<endl;

    LCCollection* pfos = event->getCollection("PandoraPFOs");
    LCRelationNavigator pfoToMc( event->getCollection("RecoMCTruthLink") );
    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );
        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();
        
        // Analyze only simple pfos. Otherwise ignore. It covers ~99.9% pfos
        if( nClusters != 1 || nTracks != 1) continue;
        
        // we can also filter only pions/kaons/protons
        MCParticle* mc = getMaxTrackWeightMC(pfo, pfoToMc);
        bool isPion = std::abs( mc->getPDG() ) == 211;
        bool isKaon = std::abs( mc->getPDG() ) == 321;
        bool isProton = std::abs( mc->getPDG() ) == 2212;
        if ( not (isPion || isKaon || isProton) ) continue;

        cout<<endl<<"        ======PFO "<<i+1<<"===== PDG: "<<mc->getPDG()<<" ======"<<endl;
        //the main fun in these functions
        printOutputFCN(pfo);
        // printOutputGNN(pfo);
        // printOutputCCN(pfo);
        
        drawPFO(event, pfo);
    }
}

void TOFInfo::printOutputFCN(ReconstructedParticle* pfo){
    /* For this one we agreed to have 10 Frank hits*/
    Cluster* cluster = pfo->getClusters()[0];
    Track* track = pfo->getTracks()[0];
    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D trackPosAtEcal ( tsEcal->getReferencePoint() );
    Vector3D trackMomAtEcal = getHelixMomAtTrackState(*tsEcal);
    //consider this to be "TRUE" reference value for the training
    double tof = getTofClosest(cluster, track, 0.);
    cout<<std::setprecision(4)<<"Track p: "<<trackMomAtEcal.r()<<"    pT: "<<trackMomAtEcal.trans()<<"    pz: "<<trackMomAtEcal.z()<<" (GeV)    theta: "<<trackMomAtEcal.theta()*180./M_PI<<" (deg)"<<"    \"TRUE\" TOF: "<<tof<<" ns"<<endl;


    int nLayers = 10;
    std::vector<CalorimeterHit*> frankHits = selectFrankEcalHits(cluster, track, nLayers);
    cout<<"Hit #"<<"    time true (ns)    time 50ps (ns)    energy (GeV)    d (mm)    d|| (mm)    dT (mm)"<<endl;

    // loop over Frank hits and print their info
    for(size_t i=0; i < frankHits.size(); ++i){
        CalorimeterHit* hit  = frankHits[i];
        Vector3D hitPos( hit->getPosition() );
        double dToTrack = (hitPos - trackPosAtEcal).r(); // you can take this
        double dPerp = (hitPos - trackPosAtEcal).cross(trackMomAtEcal.unit()).r();
        double dTau = std::sqrt(dToTrack*dToTrack - dPerp*dPerp);
        double hitTime = hit->getTime(); // in ns
        double timeResolution = 0.05; // in ns consider this input parameter, you can change that to whatever you want...
        double hitTimeSmeared = CLHEP::RandGauss::shoot(hitTime, timeResolution);
        cout<<"    "<<i+1;
        cout<<std::setprecision(4)<<"         "<<hitTime;
        cout<<std::setprecision(4)<<"          "<<hitTimeSmeared;
        cout<<std::setprecision(4)<<"             "<<hit->getEnergy();
        cout<<std::setprecision(4)<<"       "<<dToTrack;
        cout<<std::setprecision(4)<<"       "<<dTau;
        cout<<std::setprecision(4)<<"       "<<dPerp<<endl;        
    }
}


void TOFInfo::printOutputGNN(ReconstructedParticle* pfo){
    /* It is complete copy-paste of the FCN output, but instead of "Frank hits" I will use all ECAL hits from the cluster.*/
    Cluster* cluster = pfo->getClusters()[0];
    Track* track = pfo->getTracks()[0];
    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D trackPosAtEcal ( tsEcal->getReferencePoint() );
    Vector3D trackMomAtEcal = getHelixMomAtTrackState(*tsEcal);
    //consider this to be "TRUE" reference value for the training
    double tof = getTofClosest(cluster, track, 0.);
    cout<<std::setprecision(4)<<"Track p: "<<trackMomAtEcal.r()<<"    pT: "<<trackMomAtEcal.trans()<<"    pz: "<<trackMomAtEcal.z()<<" (GeV)    theta: "<<trackMomAtEcal.theta()*180./M_PI<<" (deg)"<<"    \"TRUE\" TOF: "<<tof<<" ns"<<endl;

    cout<<"Hit #"<<"    time true (ns)    time 50ps (ns)    energy (GeV)    d (mm)    d|| (mm)    dT (mm)"<<endl;

    // loop over cluster hits and print their info
    vector<CalorimeterHit*> hits = cluster->getCalorimeterHits();
    for(size_t i=0; i < hits.size(); ++i){
        CalorimeterHit* hit  = hits[i];
        //ignore non-ECAL hits
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if ( not isECALHit ) continue;
        //you can also add a limit how many layers we go in depth e.g.:
        // int layer = hitType.layer();
        // if ( layer >= 10 ) continue;

        Vector3D hitPos( hit->getPosition() );
        double dToTrack = (hitPos - trackPosAtEcal).r(); // you can take this
        double dPerp = (hitPos - trackPosAtEcal).cross(trackMomAtEcal.unit()).r();
        double dTau = std::sqrt(dToTrack*dToTrack - dPerp*dPerp);
        double hitTime = hit->getTime(); // in ns
        double timeResolution = 0.05; // in ns consider this input parameter, you can change that to whatever you want...
        double hitTimeSmeared = CLHEP::RandGauss::shoot(hitTime, timeResolution);
        cout<<"    "<<i+1;
        cout<<std::setprecision(4)<<"         "<<hitTime;
        cout<<std::setprecision(4)<<"          "<<hitTimeSmeared;
        cout<<std::setprecision(4)<<"             "<<hit->getEnergy();
        cout<<std::setprecision(4)<<"       "<<dToTrack;
        cout<<std::setprecision(4)<<"       "<<dTau;
        cout<<std::setprecision(4)<<"       "<<dPerp<<endl;        
    }
}

void TOFInfo::printOutputCNN(EVENT::ReconstructedParticle*){
    /* Hell, no, good luck with that one. See examples about how to extract individual hit information from the cluster*/
    return;
}

void TOFInfo::drawPFO(EVENT::LCEvent* event, EVENT::ReconstructedParticle* pfo){
    /* This one draws PFO in the event display for debug purposes and adds a pause after each PFO*/
    DDMarlinCED::newEvent(this, event);
    DDMarlinCED::drawDD4hepDetector(_detector, false, vector<string>{""});
    DDCEDPickingHandler& pHandler= DDCEDPickingHandler::getInstance();
    pHandler.update(event);
    
    //draw PFO using momentum at the ECAL
    Track* track = pfo->getTracks()[0];
    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D pos ( tsEcal->getReferencePoint() );
    Vector3D mom = getHelixMomAtTrackState(*tsEcal);

    int marker = 1;
    int size = 1;
    unsigned long color = 0x4242f5;
    DDMarlinCED::drawHelix(-_bField, -pfo->getCharge(), pos.x(), pos.y(), pos.z(), mom.x(), mom.y(), mom.z(), marker, size, color, 0., 3000., 3000., 0.);
    DDMarlinCED::drawHelix(_bField, pfo->getCharge(), pos.x(), pos.y(), pos.z(), mom.x(), mom.y(), mom.z(), marker, size, color, 0., 3000., 3000., 0.);

    Cluster* cluster = pfo->getClusters()[0];
    int type = 0; // point
    int pointSize = 5;
    for( auto hit : cluster->getCalorimeterHits() ){
        Vector3D hitPos (hit->getPosition());

        unsigned long hitColor;
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if ( isECALHit ) hitColor = 0x09ed1c;
        else  hitColor = 0x7a1d1d;
        
        ced_hit_ID(hitPos.x(), hitPos.y(), hitPos.z(), type, 0, pointSize, hitColor, 0);
    }

    DDMarlinCED::draw(this);

    return;
}


EVENT::MCParticle* TOFInfo::getMaxTrackWeightMC(EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator nav){
    const std::vector<LCObject*>& mcs = nav.getRelatedToObjects(pfo);
    const std::vector<float>& weights = nav.getRelatedToWeights(pfo);
    if ( mcs.size() == 0 ) return nullptr;
    //get index of the highest TRACK weight MC particle
    auto trackWeightComparator = [](float a, float b){return (int(a)%10000)/1000. < (int(b)%10000)/1000.;};
    int i = std::max_element(weights.begin(), weights.end(), trackWeightComparator) - weights.begin();
    float trackWeight = (int(weights[i])%10000)/1000.;

    // track doesn't exist
    if ( trackWeight == 0) return nullptr;
    return static_cast<MCParticle*> ( mcs[i] );
}


std::vector<EVENT::CalorimeterHit*> TOFInfo::selectFrankEcalHits( EVENT::Cluster* cluster, EVENT::Track* track, int maxEcalLayer){
    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D trackPosAtEcal ( tsEcal->getReferencePoint() );
    Vector3D trackMomAtEcal = getHelixMomAtTrackState(*tsEcal);

    vector<CalorimeterHit*> selectedHits(maxEcalLayer, nullptr);
    vector<double> minDistances(maxEcalLayer, std::numeric_limits<double>::max());

    for ( auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        int layer = hitType.layer();
        if ( (!isECALHit) || ( layer >= maxEcalLayer) ) continue;

        Vector3D hitPos( hit->getPosition() );
        double dToLine = (hitPos - trackPosAtEcal).cross(trackMomAtEcal.unit()).r();
        if ( dToLine < minDistances[layer] ){
            minDistances[layer] = dToLine;
            selectedHits[layer] = hit;
        }
    }
    // remove empty layers (you can comment it, to get always array of length 10 hits including nullptrs)
    selectedHits.erase( std::remove_if( selectedHits.begin(), selectedHits.end(), [](CalorimeterHit* h) { return h == nullptr; } ), selectedHits.end() );

    return selectedHits;
}


dd4hep::rec::Vector3D TOFInfo::getHelixMomAtTrackState(const EVENT::TrackState& ts){
    double phi = ts.getPhi();
    double d0 = ts.getD0();
    double z0 = ts.getZ0();
    double omega = ts.getOmega();
    double tanL = ts.getTanLambda();

    HelixClass helix;
    helix.Initialize_Canonical(phi, d0, z0, omega, tanL, _bField);
    return helix.getMomentum();
}

double TOFInfo::getTofClosest( EVENT::Cluster* cluster, EVENT::Track* track, double timeResolution){
    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D trackPosAtEcal ( tsEcal->getReferencePoint() );

    double hitTime = std::numeric_limits<double>::max();
    double closestDistance = std::numeric_limits<double>::max();
    for( auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if (! isECALHit) continue;

        Vector3D hitPos( hit->getPosition() );
        double dToTrack = (hitPos - trackPosAtEcal).r();
        if( dToTrack < closestDistance ){
            closestDistance = dToTrack;
            hitTime = hit->getTime();
        }
    }

    if ( hitTime == std::numeric_limits<double>::max() ) return 0.;
    return CLHEP::RandGauss::shoot(hitTime, timeResolution) - closestDistance/CLHEP::c_light;
}
