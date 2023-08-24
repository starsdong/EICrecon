// Created by Dmitry Romanov
// Subject to the terms in the LICENSE file found in the top-level directory.
//

#include <edm4eic/RawTrackerHit.h>
#include <edm4eic/TrackerHitCollection.h>
#include <JANA/JEvent.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "TrackerHitReconstruction_factory.h"
#include "services/geometry/dd4hep/JDD4hep_service.h"
#include "services/log/Log_service.h"
#include "extensions/spdlog/SpdlogExtensions.h"

void TrackerHitReconstruction_factory::Init() {

    auto app = GetApplication();
    auto param_prefix = GetDefaultParameterPrefix();

    // Set time resolution
    auto cfg = GetDefaultConfig();
    app->SetDefaultParameter(param_prefix + ":TimeResolution", cfg.time_resolution, "threshold");
    m_reco_algo.applyConfig(cfg);

    // Init logger from default or user parameters
    InitLogger(app, param_prefix);

    // Init input collections tags and read from user parameters
    InitDataTags(param_prefix);

    // get geometry service
    auto geo_service = app->GetService<JDD4hep_service>();

    // Initialize reconstruction algorithm
    m_reco_algo.init(geo_service->detector(), m_log);
}

void TrackerHitReconstruction_factory::ChangeRun(const std::shared_ptr<const JEvent> &event) {
    // nothing here
}

void TrackerHitReconstruction_factory::Process(const std::shared_ptr<const JEvent> &event) {

    // Get RawTrackerHit-s with the proper tag
    auto raw_hits = event->Get<edm4eic::RawTrackerHit>(GetInputTags()[0]);

    // Output array
    std::vector<edm4eic::TrackerHit*> hits;

    try {
        // Create output hits using TrackerHitReconstruction algorithm
        for(auto raw_hit: raw_hits){
            hits.push_back(m_reco_algo.produce(raw_hit));
        }
        Set(hits);
        m_log->debug("End of process. Hits count: {}", hits.size());
    }
    catch(std::exception &e) {
        throw JException(e.what());
    }
}
