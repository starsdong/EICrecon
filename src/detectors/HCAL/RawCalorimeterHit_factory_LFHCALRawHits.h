// Copyright 2022, David Lawrence
// Subject to the terms in the LICENSE file found in the top-level directory.
//

#ifndef _RawCalorimeterHit_factory_LFHCALRawHits_h_
#define _RawCalorimeterHit_factory_LFHCALRawHits_h_

#include <random>

#include <JANA/JEvent.h>
#include <JANA/JFactoryT.h>
#include <services/geometry/dd4hep/JDD4hep_service.h>
#include <algorithms/calorimetry/CalorimeterHitDigi.h>
#include <edm4hep/SimCalorimeterHit.h>
#include <edm4hep/RawCalorimeterHit.h>
#include <Evaluator/DD4hepUnits.h>
#include <services/log/Log_service.h>
#include <extensions/spdlog/SpdlogExtensions.h>



class RawCalorimeterHit_factory_LFHCALRawHits : public JFactoryT<edm4hep::RawCalorimeterHit>, CalorimeterHitDigi {

public:

    //------------------------------------------
    // Constructor
    RawCalorimeterHit_factory_LFHCALRawHits() {
        SetTag("LFHCALRawHits");
        m_log = japp->GetService<Log_service>()->logger(GetTag());
    }

    //------------------------------------------
    // Init
    void Init() override {
        auto app = GetApplication();

        // Set default values for all config. parameters in CalorimeterHitDigi algorithm
        m_input_tag = "LFHCALHits";
        u_eRes = {};
        m_tRes = 0.0; // in ns
        m_capTime = 100 ; // given in ns
        m_capADC=65536;//2^16
        m_dyRangeADC=1 * dd4hep::GeV;
        m_pedMeanADC=20;
        m_pedSigmaADC = 0.8;
        m_resolutionTDC = 10 * dd4hep::picosecond;
        m_corrMeanScale = 1.0;
//         u_fields={"rlayerz","layerz"};
//         u_refs={1, 0};
        u_fields={"layerz"};
        u_refs={0};
        m_geoSvcName = "geoServiceName";
        m_readout = "LFHCALHits";
        m_geoSvc = app->GetService<JDD4hep_service>(); // TODO: implement named geometry service?

        // This is another option for exposing the data members as JANA configuration parameters.
//        app->SetDefaultParameter("HCAL:tag",              m_input_tag);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:energyResolutions",u_eRes);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:maxTimeInt",       m_capTime);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:timeResolution",   m_tRes);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:capacityADC",      m_capADC);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:dynamicRangeADC",  m_dyRangeADC);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:pedestalMean",     m_pedMeanADC);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:pedestalSigma",    m_pedSigmaADC);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:resolutionTDC",    m_resolutionTDC);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:scaleResponse",    m_corrMeanScale);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:signalSumFields",  u_fields);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:fieldRefNumbers",  u_refs);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:geoServiceName",   m_geoSvcName);
        app->SetDefaultParameter("HCAL:LFHCALRawHits:readoutClass",     m_readout);

        // Call Init for generic algorithm
        AlgorithmInit(m_log);
    }

    //------------------------------------------
    // ChangeRun
    void ChangeRun(const std::shared_ptr<const JEvent> &event) override {
        AlgorithmChangeRun();
    }

    //------------------------------------------
    // Process
    void Process(const std::shared_ptr<const JEvent> &event) override {
        // Prefill inputs
        simhits = event->Get<edm4hep::SimCalorimeterHit>(m_input_tag);

        // Call Process for generic algorithm
        AlgorithmProcess();

        // Hand owner of algorithm objects over to JANA
        Set(rawhits);
        rawhits.clear(); // not really needed, but better to not leave dangling pointers around
    }

};

#endif // _RawCalorimeterHit_factory_LFHCALRawHits_h_