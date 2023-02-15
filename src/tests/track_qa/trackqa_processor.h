#ifndef EICRECON_TRK_QA_PROCESSOR_H
#define EICRECON_TRK_QA_PROCESSOR_H

#include <JANA/JEventProcessor.h>
#include <JANA/JEventProcessorSequentialRoot.h>

#include <services/log/Log_service.h>
#include <extensions/spdlog/SpdlogMixin.h>
#include "algorithms/tracking/ActsGeometryProvider.h"

#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <iostream>
// #include <fstream>
// #include <string>
#include <vector>
// #include <sstream>
// #include <cstdlib>
using namespace std;

class JEvent;
class JApplication;

class trackqa_processor:
        public JEventProcessor,
        public eicrecon::SpdlogMixin<trackqa_processor>   // this automates proper log initialization
{
public:
    explicit trackqa_processor(JApplication *);
    ~trackqa_processor() override = default;

    //----------------------------
    // Init
    //
    // This is called once before the first call to the Process method
    // below. You may, for example, want to open an output file here.
    // Only one thread will call this.
    void Init() override;

    //----------------------------
    // Process
    //
    // This is called for every event. Multiple threads may call this
    // simultaneously. If you write something to an output file here
    // then make sure to protect it with a mutex or similar mechanism.
    // Minimize what is done while locked since that directly affects
    // the multi-threaded performance.
    void Process(const std::shared_ptr<const JEvent>& event) override;

    //----------------------------
    // Finish
    //
    // This is called once after all events have been processed. You may,
    // for example, want to close an output file here.
    // Only one thread will call this.
    void Finish() override;

    //Histograms
    TH2 *h1a; //Rec track momentum vs. true momentum
    TH1 *hchi2; //chi^2 histogram
    TH1 *heta; //eta histogram
    TH1 *hp; //p histogram
    TH1 *hpt; //pt histogram
    TH1 *hhits; //number of hits
    TH1 *hNDF; //number of degrees of freedom
    TH1 *hchi2_by_hits; //chi^2 divided by the # of hits
    TH1 *hchi2_by_NDF; //chi^2 divided by the # of degrees of freedom

        //Looking at chi^2 and # of hits
    TH2 *hchi2_vs_eta; //chi^2 vs eta
    TH2 *hchi2_vs_hits; //chi^2 vs # of hits (chi^2 is in y-axis, # of hits is in x-axis)
    TH2 *hchi2_vs_hits_zoomed; //chi^2 vs # of hits (chi^2 is in y-axis, # of hits is in x-axis) zoomed in
    vector<TH2*> hchi2_vs_hits_etabins; //chi^2 vs # of hits in 16 bins of eta
    TH2 *hhits_vs_eta; //# of hits vs eta
    TH2 *htracks_vs_eta; //# of tracks vs eta
    TH3 *heta_vs_p_vs_chi2; //#eta vs p vs chi^2

        //now looking at the number of measurements per track
    TH2 *hmeasptrack_vs_eta; //# of measurements per track vs eta
    TH2 *hmeasptrack_vs_hits; //# of measurements per track vs # of hits
    TH2 *hmeasptrack_vs_chi2perNDF; //# of measurements per track vs chi^2
    vector<TH2*> hmeasptrack_vs_hits_etabins; //# of measurements per track vs # of hits in 16 bins of eta
    vector<TH2*> hmeasptrack_vs_hits_etabins_zoomed; //# of measurements per track vs # of hits in 16 bins of eta zoomed in (smaller range)
    vector<TH2*> hmeasptrack_vs_chi2perNDF_etabins; //# of measurements per track vs chi^2 in 16 bins of eta
    TH2 *hmeasptrack_vs_calstates; //# of measurements per tracks vs # of calibrated states

    TH2 *hholes_vs_hits; //# of holes per # of hits
    TH2 *houtliers_vs_hits; //# of outliers [er # of hits
    TH2 *hsummation; //confirm that # calibrated states = # of trajectories + # outliers
    TH2 *hsummation2; //confirm that # hits = # of outliers + # meas per track

        //look at the individual chi^2 per measurement
    TH2 *hmeaschi2_vs_chi2;
    TH2 *hmeaschi2_vs_eta;
    TH2 *hmeaschi2_vs_hits;

private:

    std::shared_ptr<const ActsGeometryProvider> m_geo_provider;

    /// Directory to store histograms to
    TDirectory *m_dir_main{};
    TDirectory *m_dir_sub{};

};

#endif //EICRECON_TRK_QA_PROCESSOR_H
