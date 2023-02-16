#include "trackqa_processor.h"

#include <JANA/JApplication.h>
#include <JANA/JEvent.h>

#include <Math/GenVector/PxPyPzM4D.h>

#include <spdlog/spdlog.h>

#include <algorithms/tracking/JugTrack/TrackingResultTrajectory.hpp>
#include <algorithms/tracking/ParticlesFromTrackFitResult.h>

#include <services/rootfile/RootFile_service.h>
#include <services/geometry/acts/ACTSGeo_service.h>

#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>

#include <edm4hep/MCParticle.h>


//------------------
// (Constructor)
//------------------
trackqa_processor::trackqa_processor(JApplication *app) :
	JEventProcessor(app)
{
}

//------------------
// Init
//------------------
void trackqa_processor::Init()
{
    std::string plugin_name=("track_qa");

    // Get JANA application
    auto app = GetApplication();

    // Ask service locator a file to write histograms to
    auto root_file_service = app->GetService<RootFile_service>();

    // Get TDirectory for histograms root file
    auto globalRootLock = app->GetService<JGlobalRootLock>();
    globalRootLock->acquire_write_lock();
    auto file = root_file_service->GetHistFile();
    globalRootLock->release_lock();

    // Create a directory for this plugin. And subdirectories for series of histograms
    m_dir_main = file->mkdir(plugin_name.c_str());
    m_dir_sub = m_dir_main->mkdir("eta_bins");

    //Define histograms
    h1a = new TH2D("h1a","",100,0,25,100,0,25);
    h1a->GetXaxis()->SetTitle("True Momentum [GeV/c]");h1a->GetXaxis()->CenterTitle();
    h1a->GetYaxis()->SetTitle("Rec. Track Momentum [GeV/c]");h1a->GetYaxis()->CenterTitle();
    h1a->SetDirectory(m_dir_main);

    hchi2 = new TH1D("hchi2","",50,0,50);
    hchi2->GetXaxis()->SetTitle("Track #Chi^{2} Sum");hchi2->GetXaxis()->CenterTitle();
    hchi2->GetYaxis()->SetTitle("Counts");hchi2->GetYaxis()->CenterTitle();
    hchi2->SetLineWidth(2);hchi2->SetLineColor(kBlue);
    hchi2->SetDirectory(m_dir_main);

    heta = new TH1D("heta","",50,-4,4);
    heta->GetXaxis()->SetTitle("#eta (Generated)");heta->GetXaxis()->CenterTitle();
    heta->GetYaxis()->SetTitle("Counts");heta->GetYaxis()->CenterTitle();
    heta->SetLineWidth(2);heta->SetLineColor(kBlue);
    heta->SetDirectory(m_dir_main);

    hp = new TH1D("hp","",44,0,11);
    hp->GetXaxis()->SetTitle("Momentum [GeV]");hp->GetXaxis()->CenterTitle();
    hp->GetYaxis()->SetTitle("Counts");hp->GetYaxis()->CenterTitle();
    hp->SetLineWidth(2);hp->SetLineColor(kBlue);
    hp->SetDirectory(m_dir_main);

    hpt = new TH1D("hpt","",50,0,10);
    hpt->GetXaxis()->SetTitle("Transverse Momentum [GeV]");hpt->GetXaxis()->CenterTitle();
    hpt->GetYaxis()->SetTitle("Counts");hpt->GetYaxis()->CenterTitle();
    hpt->SetLineWidth(2);hpt->SetLineColor(kBlue);
    hpt->SetDirectory(m_dir_main);

    hhits = new TH1D("hhits","",50,0,50);
    hhits->GetXaxis()->SetTitle("Number of Hits Found in Tracker");hhits->GetXaxis()->CenterTitle();
    hhits->GetYaxis()->SetTitle("Counts");hhits->GetYaxis()->CenterTitle();
    hhits->SetLineWidth(2);hhits->SetLineColor(kBlue);
    hhits->SetDirectory(m_dir_main);

    hNDF = new TH1D("hNDF","",25,0,25);
    hNDF->GetXaxis()->SetTitle("NDF");hNDF->GetXaxis()->CenterTitle();
    hNDF->GetYaxis()->SetTitle("Counts");hNDF->GetYaxis()->CenterTitle();
    hNDF->SetLineWidth(2);hNDF->SetLineColor(kBlue);
    hNDF->SetDirectory(m_dir_main);

    hchi2_by_hits = new TH1D("hchi2_by_hits","",50,0,10);
    hchi2_by_hits->GetXaxis()->SetTitle("Track #Chi^{2} Sum/Number of Hits");hchi2_by_hits->GetXaxis()->CenterTitle();
    hchi2_by_hits->GetYaxis()->SetTitle("Counts");hchi2_by_hits->GetYaxis()->CenterTitle();
    hchi2_by_hits->SetLineWidth(2);hchi2_by_hits->SetLineColor(kBlue);
    hchi2_by_hits->SetDirectory(m_dir_main);

    hchi2_by_NDF = new TH1D("hchi2_by_NDF","",50,0,10);
    hchi2_by_NDF->GetXaxis()->SetTitle("Track #Chi^{2} Sum/NDF");hchi2_by_NDF->GetXaxis()->CenterTitle();
    hchi2_by_NDF->GetYaxis()->SetTitle("Counts");hchi2_by_NDF->GetYaxis()->CenterTitle();
    hchi2_by_NDF->SetLineWidth(2);hchi2_by_NDF->SetLineColor(kBlue);
    hchi2_by_NDF->SetDirectory(m_dir_main);

    hchi2_by_meas = new TH1D("hchi2_by_meas","",50,0,10);
    hchi2_by_meas->GetXaxis()->SetTitle("Track #Chi^{2} Sum/meas");hchi2_by_meas->GetXaxis()->CenterTitle();
    hchi2_by_meas->GetYaxis()->SetTitle("Counts");hchi2_by_meas->GetYaxis()->CenterTitle();
    hchi2_by_meas->SetLineWidth(2);hchi2_by_meas->SetLineColor(kBlue);
    hchi2_by_meas->SetDirectory(m_dir_main);

    //chi^2 and number of hits
    hchi2_vs_eta = new TH2D("hchi2_vs_eta","",50,-4,4,50,0,50);
    hchi2_vs_eta->GetXaxis()->SetTitle("#eta (Generated)");hchi2_vs_eta->GetXaxis()->CenterTitle();
    hchi2_vs_eta->GetYaxis()->SetTitle("Track #Chi^{2} Sum");hchi2_vs_eta->GetYaxis()->CenterTitle();
    hchi2_vs_eta->SetDirectory(m_dir_main);

    hchi2_vs_hits = new TH2D("hchi2_vs_hits","",50,0,50,50,0,50);
    hchi2_vs_hits->GetXaxis()->SetTitle("Number of Hits");hchi2_vs_hits->GetXaxis()->CenterTitle();
    hchi2_vs_hits->GetYaxis()->SetTitle("Track #Chi^{2} Sum");hchi2_vs_hits->GetYaxis()->CenterTitle();
    hchi2_vs_hits->SetDirectory(m_dir_main);

    hchi2_vs_hits_zoomed = new TH2D("hchi2_vs_hits_zoomed","",10,0,10,50,0,5);
    hchi2_vs_hits_zoomed->GetXaxis()->SetTitle("Number of Hits");hchi2_vs_hits_zoomed->GetXaxis()->CenterTitle();
    hchi2_vs_hits_zoomed->GetYaxis()->SetTitle("Track #Chi^{2} Sum");hchi2_vs_hits_zoomed->GetYaxis()->CenterTitle();
    hchi2_vs_hits_zoomed->SetDirectory(m_dir_main);

    hhits_vs_eta = new TH2D("hhits_vs_eta","",50,-4,4,50,0,50);
    hhits_vs_eta->GetXaxis()->SetTitle("#eta (Generated)");hhits_vs_eta->GetXaxis()->CenterTitle();
    hhits_vs_eta->GetYaxis()->SetTitle("Number of Hits");hhits_vs_eta->GetYaxis()->CenterTitle();
    hhits_vs_eta->SetDirectory(m_dir_main);

    hhits_vs_eta_1 = new TH2D("hhits_vs_eta_1","At least one track reconstructed",50,-4,4,50,0,50);
    hhits_vs_eta_1->GetXaxis()->SetTitle("#eta (Generated)");hhits_vs_eta_1->GetXaxis()->CenterTitle();
    hhits_vs_eta_1->GetYaxis()->SetTitle("Number of Hits");hhits_vs_eta_1->GetYaxis()->CenterTitle();
    hhits_vs_eta_1->SetDirectory(m_dir_main);

    htracks_vs_eta = new TH2D("htracks_vs_eta","",50,-4,4,10,0,10);
    htracks_vs_eta->GetXaxis()->SetTitle("#eta (Generated)");htracks_vs_eta->GetXaxis()->CenterTitle();
    htracks_vs_eta->GetYaxis()->SetTitle("Number of Tracks");htracks_vs_eta->GetYaxis()->CenterTitle();
    htracks_vs_eta->SetDirectory(m_dir_main);

    heta_vs_p_vs_chi2 = new TH3D("heta_vs_p_vs_chi2","",50,-4,4,44,0,11,50,0,50); //x: eta, y: p, z: chi^2
    heta_vs_p_vs_chi2->GetXaxis()->SetTitle("#eta (Generated)"); heta_vs_p_vs_chi2->GetXaxis()->CenterTitle();
    heta_vs_p_vs_chi2->GetYaxis()->SetTitle("Momentum [GeV]"); heta_vs_p_vs_chi2->GetYaxis()->CenterTitle();
    heta_vs_p_vs_chi2->GetZaxis()->SetTitle("Track #Chi^{2} Sum"); heta_vs_p_vs_chi2->GetZaxis()->CenterTitle();
    heta_vs_p_vs_chi2->SetDirectory(m_dir_main);

    hNDF_states = new TH2D("hNDF_states","",25,0,25,25,0,25);
    hNDF_states->GetXaxis()->SetTitle("NDF");hNDF_states->GetXaxis()->CenterTitle();
    hNDF_states->GetYaxis()->SetTitle("States");hNDF_states->GetYaxis()->CenterTitle();
    hNDF_states->SetDirectory(m_dir_main);

    hmeasptrack_vs_eta = new TH2D("hmeasptrack_vs_eta","",50,-4,4,10,0,10);
    hmeasptrack_vs_eta->GetXaxis()->SetTitle("#eta (Generated)");hmeasptrack_vs_eta->GetXaxis()->CenterTitle();
    hmeasptrack_vs_eta->GetYaxis()->SetTitle("Number of Measurements per Track");hmeasptrack_vs_eta->GetYaxis()->CenterTitle();
    hmeasptrack_vs_eta->SetDirectory(m_dir_main);

    hmeasptrack_vs_hits = new TH2D("hmeasptrack_vs_hits","",50,0,50,10,0,10);
    hmeasptrack_vs_hits->GetXaxis()->SetTitle("Number of Hits");hmeasptrack_vs_hits->GetXaxis()->CenterTitle();
    hmeasptrack_vs_hits->GetYaxis()->SetTitle("Number of Measurements per Track");hmeasptrack_vs_hits->GetYaxis()->CenterTitle();
    hmeasptrack_vs_hits->SetDirectory(m_dir_main);

    hmeasptrack_vs_chi2perNDF = new TH2D("hmeasptrack_vs_chi2perNDF","",50,0,50,10,0,10);
    hmeasptrack_vs_chi2perNDF->GetXaxis()->SetTitle("Track #Chi^{2} Sum/NDF");hmeasptrack_vs_chi2perNDF->GetXaxis()->CenterTitle();
    hmeasptrack_vs_chi2perNDF->GetYaxis()->SetTitle("Number of Measurements per Track");hmeasptrack_vs_chi2perNDF->GetYaxis()->CenterTitle();
    hmeasptrack_vs_chi2perNDF->SetDirectory(m_dir_main);

    hmeasptrack_vs_calstates = new TH2D("hmeasptrack_vs_calstates","",15,0,15,15,0,15);
    hmeasptrack_vs_calstates->GetXaxis()->SetTitle("Number of Measurements per Track");hmeasptrack_vs_calstates->GetXaxis()->CenterTitle();
    hmeasptrack_vs_calstates->GetYaxis()->SetTitle("Number of Calibrated States");hmeasptrack_vs_calstates->GetYaxis()->CenterTitle();
    hmeasptrack_vs_calstates->SetDirectory(m_dir_main);

    hmeaschi2_vs_chi2 = new TH2D("hmeaschi2_vs_chi2","",50,0,50,50,0,20);
    hmeaschi2_vs_chi2->GetXaxis()->SetTitle("Track #Chi^{2} Sum");hmeaschi2_vs_chi2->GetXaxis()->CenterTitle();
    hmeaschi2_vs_chi2->GetYaxis()->SetTitle("#Chi^{2} Individual Measurements");hmeaschi2_vs_chi2->GetYaxis()->CenterTitle();
    hmeaschi2_vs_chi2->SetDirectory(m_dir_main);

    hmeaschi2_vs_eta = new TH2D("hmeaschi2_vs_eta","",50,-4,4,50,0,20);
    hmeaschi2_vs_eta->GetXaxis()->SetTitle("#eta");hmeaschi2_vs_eta->GetXaxis()->CenterTitle();
    hmeaschi2_vs_eta->GetYaxis()->SetTitle("#Chi^{2} Individual Measurements");hmeaschi2_vs_eta->GetYaxis()->CenterTitle();
    hmeaschi2_vs_eta->SetDirectory(m_dir_main);

    hmeaschi2_vs_hits = new TH2D("hmeaschi2_vs_hits","",50,0,50,50,0,20);
    hmeaschi2_vs_hits->GetXaxis()->SetTitle("Number of Hits");hmeaschi2_vs_hits->GetXaxis()->CenterTitle();
    hmeaschi2_vs_hits->GetYaxis()->SetTitle("#Chi^{2} Individual Measurements");hmeaschi2_vs_hits->GetYaxis()->CenterTitle();
    hmeaschi2_vs_hits->SetDirectory(m_dir_main);

        //use 0.5i-4 to get lowerbound and 0.5i-3.5 to get upper bound
    for (int i=0; i<16; i++){
        TH2 *htemp = new TH2D(TString::Format("hchi2_vs_hits_eta_%.1f_%.1f", 0.5*i-4,0.5*i-3.5),
                    "",50,0,50,50,0,50);
        htemp->GetXaxis()->SetTitle("Number of Hits");htemp->GetXaxis()->CenterTitle();
        htemp->GetYaxis()->SetTitle("Track #Chi^{2} Sum");htemp->GetYaxis()->CenterTitle();
        hchi2_vs_hits_etabins.push_back(htemp);
        hchi2_vs_hits_etabins[i]->SetDirectory(m_dir_sub);

        TH2 *htemp1 = new TH2D(TString::Format("hmeasptrack_vs_hits_eta_%.1f_%.1f", 0.5*i-4,0.5*i-3.5),
                    "",50,0,50,10,0,10);
        htemp1->GetXaxis()->SetTitle("Number of Hits");htemp1->GetXaxis()->CenterTitle();
        htemp1->GetYaxis()->SetTitle("Number of Measurements per Track");htemp1->GetYaxis()->CenterTitle();
        hmeasptrack_vs_hits_etabins.push_back(htemp1);
        hmeasptrack_vs_hits_etabins[i]->SetDirectory(m_dir_sub);

        TH2 *htemp1b = new TH2D(TString::Format("hmeasptrack_vs_hits_eta_zoomed_%.1f_%.1f", 0.5*i-4,0.5*i-3.5),
                    "",10,0,10,10,0,10);
        htemp1b->GetXaxis()->SetTitle("Number of Hits");htemp1b->GetXaxis()->CenterTitle();
        htemp1b->GetYaxis()->SetTitle("Number of Measurements per Track");htemp1b->GetYaxis()->CenterTitle();
        hmeasptrack_vs_hits_etabins_zoomed.push_back(htemp1b);
        hmeasptrack_vs_hits_etabins_zoomed[i]->SetDirectory(m_dir_sub);

        TH2 *htemp2 = new TH2D(TString::Format("hmeasptrack_vs_chi2perNDF_eta_%.1f_%.1f", 0.5*i-4,0.5*i-3.5),
                    "",50,0,50,10,0,10);
        htemp2->GetXaxis()->SetTitle("Track #Chi^{2} Sum/NDF");htemp2->GetXaxis()->CenterTitle();
        htemp2->GetYaxis()->SetTitle("Number of Measurements per Track");htemp2->GetYaxis()->CenterTitle();
        hmeasptrack_vs_chi2perNDF_etabins.push_back(htemp2);
        hmeasptrack_vs_chi2perNDF_etabins[i]->SetDirectory(m_dir_sub);
    }
    // hchi2_vs_hits_etabins->SetDirectory(m_dir_main);

    hholes_vs_hits = new TH2D("hholes_vs_hits","",50,0,50,5,0,5);
    hholes_vs_hits->GetXaxis()->SetTitle("Number of Hits");hholes_vs_hits->GetXaxis()->CenterTitle();
    hholes_vs_hits->GetYaxis()->SetTitle("Number of Holes");hholes_vs_hits->GetYaxis()->CenterTitle();
    hholes_vs_hits->SetDirectory(m_dir_main);

    houtliers_vs_hits = new TH2D("houtliers_vs_hits","",50,0,50,5,0,5);
    houtliers_vs_hits->GetXaxis()->SetTitle("Number of Hits");houtliers_vs_hits->GetXaxis()->CenterTitle();
    houtliers_vs_hits->GetYaxis()->SetTitle("Number of Outliers");houtliers_vs_hits->GetYaxis()->CenterTitle();
    houtliers_vs_hits->SetDirectory(m_dir_main);

    hsummation = new TH2D("hsummation","",15,0,15,15,0,15);
    hsummation->GetXaxis()->SetTitle("Number of Meas per Track + Number of Outliers");hsummation->GetXaxis()->CenterTitle();
    hsummation->GetYaxis()->SetTitle("Number of Calibrated States");hsummation->GetYaxis()->CenterTitle();
    hsummation->SetDirectory(m_dir_main);

    hsummation2 = new TH2D("hsummation2","",50,0,50,50,0,50);
    hsummation2->GetXaxis()->SetTitle("Number of Meas per Track + Number of Outliers");hsummation2->GetXaxis()->CenterTitle();
    hsummation2->GetYaxis()->SetTitle("Number of Hits");hsummation2->GetYaxis()->CenterTitle();
    hsummation2->SetDirectory(m_dir_main);

    hsummation3 = new TH2D("hsummation3","",15,0,15,15,0,15);
    hsummation3->GetXaxis()->SetTitle("Number of Meas per Track + Number of Outliers + Number of Holes");hsummation3->GetXaxis()->CenterTitle();
    hsummation3->GetYaxis()->SetTitle("Number of States");hsummation3->GetYaxis()->CenterTitle();
    hsummation3->SetDirectory(m_dir_main);


    // Get log level from user parameter or default
    InitLogger(plugin_name);

    auto acts_service = app->GetService<ACTSGeo_service>();
    m_geo_provider = acts_service->actsGeoProvider();

}


//------------------
// Process
//------------------
// This function is called every event
void trackqa_processor::Process(const std::shared_ptr<const JEvent>& event)
{
    m_log->trace("");
    m_log->trace("trackqa_processor event");

    //Generated particles (only one particle generated)
    double mceta = 0; //particle eta
    double mcphi = 0; //particle phi
    double mcp = 0; //total momentum
    double mcpt = 0; //transverse momentum
    double mce = 0; //energy
    int mcid = 0; //PDG id
    int num_primary = 0; //Number of primary particles

    auto mcParticles = event->Get<edm4hep::MCParticle>("MCParticles");

    for( size_t iParticle=0;iParticle<mcParticles.size(); iParticle++ ){
	
	    auto mcparticle = mcParticles[iParticle];
        if(mcparticle->getGeneratorStatus() != 1) continue;
        auto& mom = mcparticle->getMomentum();
        
        mceta = -log(tan(atan2(sqrt(mom.x*mom.x+mom.y*mom.y),mom.z)/2.));
        mcphi = atan2(mom.y, mom.x);
	    mcp = sqrt(mom.x*mom.x+mom.y*mom.y+mom.z*mom.z);
        mcpt = sqrt(mom.x*mom.x+mom.y*mom.y);
        mce = sqrt(mcp*mcp + mcparticle->getMass()*mcparticle->getMass());	

        mcid = mcparticle->getPDG();

	    num_primary++;

    }

    //Print generated info to log
    m_log->trace("-------------------------");
    m_log->trace("Number of primary generated particles:");
    m_log->trace("{:>10}",num_primary);
    m_log->trace("Generated particle id, eta, p, E:");
    m_log->trace("{:>10} {:>10.2f} {:>10.2f} {:>10.2f}", mcid, mceta, mcp, mce);

    //Tracker hits
    std::vector<std::string> m_data_names = {
        "SiBarrelTrackerRecHits",         // Barrel Tracker
        "SiBarrelVertexRecHits",          // Vertex
        "SiEndcapTrackerRecHits",         // End Cap tracker
        "MPGDBarrelRecHits",              // MPGD
        "MPGDDIRCRecHits",                // MPGD DIRC
        "TOFBarrelRecHit",                // Barrel TOF
        "TOFEndcapRecHits"                // End Cap TOF
    };

    m_log->trace("-------------------------");
    int nHitsallTrackers = 0;
    
    for(size_t name_index = 0; name_index < m_data_names.size(); name_index++ ) {
        auto data_name = m_data_names[name_index];
        auto hits = event->Get<edm4eic::TrackerHit>(data_name);
        m_log->trace("Detector {} has {} digitized hits.",data_name,hits.size());
        
        int nHits_detector = 0;
        for(auto hit: hits) {
            auto cell_id = hit->getCellID(); //FIXME: convert to volume id and layer id
            auto x = hit->getPosition().x;
            auto y = hit->getPosition().y;
            auto z = hit->getPosition().z;
            auto r = sqrt(x*x + y*y);
            auto etahit = -log(tan(atan2(r,z)/2.));

            nHits_detector++;
            m_log->trace("For digitized hit number {}:",nHits_detector);
            m_log->trace("Cell Id is {}",cell_id);
            m_log->trace("Hit x, y, z, r, eta:");
            m_log->trace("{:>10.2f} {:>10.2f} {:>10.2f} {:>10.2f} {:>10.2f}",x,y,z,r,etahit);
            nHitsallTrackers++;
        }

        m_log->trace("");
    }

    m_log->trace("Total number of tracker hits is {}",nHitsallTrackers);
    m_log->trace("-------------------------");

    //ACTS Trajectories
    auto trajectories = event->Get<eicrecon::TrackingResultTrajectory>("CentralCKFTrajectories");
    //auto trajectories = event->Get<eicrecon::TrackingResultTrajectory>("CentralCKFSeededTrajectories"); // for realistic seeded trajectories

    m_log->trace("Number of ACTS Trajectories: {}", trajectories.size());
    m_log->trace("");

    // Loop over the trajectories
    int num_traj = 0;

    for (const auto& traj : trajectories) {

        // Get the entry index for the single trajectory
        // The trajectory entry indices and the multiTrajectory
        const auto &mj = traj->multiTrajectory();
        const auto &trackTips = traj->tips();

        // Skip empty
        if (trackTips.empty()) {
            m_log->trace("Empty multiTrajectory.");
            continue;
        }
        auto &trackTip = trackTips.front();

        // Collect the trajectory summary info
        auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
        int m_nStates = trajState.nStates;
        int m_nMeasurements = trajState.nMeasurements;
        int m_nOutliers = trajState.nOutliers;
        auto m_chi2Sum = trajState.chi2Sum;
        int m_NDF = trajState.NDF;
        auto m_measurementChi2 = trajState.measurementChi2;
        int m_nHoles = trajState.nHoles;

        //General Trajectory Information
        m_log->trace("Number of elements in trackTips {}", trackTips.size());
        m_log->trace("Number of states in trajectory     : {}", m_nStates);
        m_log->trace("Number of measurements in trajectory: {}", m_nMeasurements);
        m_log->trace("Number of outliers in trajectory     : {}", m_nOutliers);
        m_log->trace("Total chi-square of trajectory    : {:>8.2f}",m_chi2Sum);
   
        //Initial Trajectory Parameters
        const auto &initial_bound_parameters = traj->trackParameters(trackTip);
        const auto &parameter  = initial_bound_parameters.parameters();
        const auto &covariance = *initial_bound_parameters.covariance();

        m_log->trace("{:>8} {:>8} {:>8} {:>8} {:>8} {:>10} {:>10} {:>10}",
                    "[loc 0]","[loc 1]", "[phi]", "[theta]", "[q/p]", "[err phi]", "[err th]", "[err q/p]");
        m_log->trace("{:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>10.2g} {:>10.2g} {:>10.2g}",
                    parameter[Acts::eBoundLoc0],
                    parameter[Acts::eBoundLoc1],
                    parameter[Acts::eBoundPhi],
                    parameter[Acts::eBoundTheta],
                    parameter[Acts::eBoundQOverP],
                    sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)),
                    sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)),
                    sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP))
                    );

        auto p_traj = fabs(1. / parameter[Acts::eBoundQOverP]);
        auto eta_traj = -log( tan(parameter[Acts::eBoundTheta]/2.) );
        m_log->trace("Trajectory p, eta:");
        m_log->trace("{:>10.2f} {:>10.2f}",p_traj,eta_traj);
        m_log->trace("");

        //Information at tracking layers
        int m_nCalibrated = 0;
        int state_counter = 0;
        // visit the track points
        mj.visitBackwards(trackTip, [&](auto &&trackstate) {

            state_counter++;
            m_log->trace("Now at State number {}",state_counter);

            // get volume info
            auto geoID = trackstate.referenceSurface().geometryId();
            auto volume = geoID.volume();
            auto layer = geoID.layer();
            m_log->trace("Volume id is {}, layer id is {}",volume,layer);

            if (trackstate.hasCalibrated()) {
                m_nCalibrated++;
                m_log->trace("This is a calibrated state.");
            }
            else{
                m_log->trace("This is NOT a calibrated state.");
            }

            // get track state parameters and their covariances
            const auto &state_params = trackstate.predicted();
            const auto &state_covar = trackstate.predictedCovariance();

            //First print same information as for initial track parameters
            m_log->trace("{:>8} {:>8} {:>8} {:>8} {:>8} {:>10} {:>10} {:>10}",
                    "[loc 0]","[loc 1]", "[phi]", "[theta]", "[q/p]", "[err phi]", "[err th]", "[err q/p]");
            m_log->trace("{:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>10.2g} {:>10.2g} {:>10.2g}",
                    state_params[Acts::eBoundLoc0],
                    state_params[Acts::eBoundLoc1],
                    state_params[Acts::eBoundPhi],
                    state_params[Acts::eBoundTheta],
                    state_params[Acts::eBoundQOverP],
                    sqrt(state_covar(Acts::eBoundPhi, Acts::eBoundPhi)),
                    sqrt(state_covar(Acts::eBoundTheta, Acts::eBoundTheta)),
                    sqrt(state_covar(Acts::eBoundQOverP, Acts::eBoundQOverP))
                    );

            // convert local to global
            auto global = trackstate.referenceSurface().localToGlobal(
                    m_geo_provider->getActsGeometryContext(),
                    {parameter[Acts::eBoundLoc0], parameter[Acts::eBoundLoc1]},
                    {0, 0, 0} );
            auto global_r = sqrt(global.x()*global.x()+global.y()*global.y());
            m_log->trace("State global x, y, z, r and pathlength:");
            m_log->trace("{:>10.2f} {:>10.2f} {:>10.2f} {:>10.2f} {:>10.2f}",global.x(),global.y(),global.z(),global_r,trackstate.pathLength());
            
            m_log->trace("");
        }); //End visiting track points
     
        m_log->trace("Number of calibrated states: {}",m_nCalibrated);
        num_traj++;
        
        //Fill histograms
        if(num_primary==1){
            h1a->Fill(mcp,p_traj);
            hchi2->Fill(m_chi2Sum);
            hNDF->Fill(m_NDF);
            hchi2_by_hits->Fill(m_chi2Sum/nHitsallTrackers);
            hchi2_by_NDF->Fill(m_chi2Sum/m_NDF);
            hchi2_by_meas->Fill(m_chi2Sum/m_nMeasurements);
            
            hchi2_vs_eta->Fill(mceta, m_chi2Sum);
            hchi2_vs_hits->Fill(nHitsallTrackers, m_chi2Sum);
            hchi2_vs_hits_zoomed->Fill(nHitsallTrackers, m_chi2Sum);
            heta_vs_p_vs_chi2->Fill(mceta, mcp, m_chi2Sum);

            hNDF_states->Fill(m_NDF,state_counter);
            
            hmeasptrack_vs_eta->Fill(mceta, m_nMeasurements);
            hmeasptrack_vs_hits->Fill(nHitsallTrackers, m_nMeasurements);
            hmeasptrack_vs_chi2perNDF->Fill(m_chi2Sum/m_NDF, m_nMeasurements);
            hmeasptrack_vs_calstates->Fill(m_nCalibrated, m_nMeasurements);
            
            //floor(2*eta + 8) should give the index where bounds are [beg,end)
            int index = floor(2*mceta+8);
            hchi2_vs_hits_etabins[index]->Fill(nHitsallTrackers, m_chi2Sum);
            hmeasptrack_vs_hits_etabins[index]->Fill(nHitsallTrackers, m_nMeasurements);
            hmeasptrack_vs_hits_etabins_zoomed[index]->Fill(nHitsallTrackers, m_nMeasurements);
            hmeasptrack_vs_chi2perNDF_etabins[index]->Fill(m_chi2Sum/m_NDF, m_nMeasurements);
            

            for (int j=0; j<m_measurementChi2.size(); j++){
                hmeaschi2_vs_chi2->Fill(m_chi2Sum, m_measurementChi2[j]);
                hmeaschi2_vs_eta->Fill(mceta, m_measurementChi2[j]);
                hmeaschi2_vs_hits->Fill(nHitsallTrackers, m_measurementChi2[j]);
            }

            hholes_vs_hits->Fill(nHitsallTrackers, m_nHoles);
            houtliers_vs_hits->Fill(nHitsallTrackers, m_nOutliers);
            hsummation->Fill(m_nMeasurements + m_nOutliers, m_nCalibrated);
            hsummation2->Fill(m_nOutliers + m_nMeasurements, nHitsallTrackers);
            hsummation3->Fill(m_nMeasurements + m_nOutliers + m_nHoles,m_nStates);
            
        }

    } //End loop over trajectories

    if(num_primary==1){
        heta->Fill(mceta);
        hp->Fill(mcp);
        hpt->Fill(mcpt);
        hhits->Fill(nHitsallTrackers);
        htracks_vs_eta->Fill(mceta, num_traj);

        hhits_vs_eta->Fill(mceta, nHitsallTrackers);
        if(num_traj>0) hhits_vs_eta_1->Fill(mceta, nHitsallTrackers);
    }

    m_log->trace("-------------------------");

}

//------------------
// Finish
//------------------
void trackqa_processor::Finish()
{
	m_log->trace("trackqa_processor finished\n");
}

