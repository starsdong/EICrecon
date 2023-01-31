// Copyright 2022, Christopher Dilks, adapted from Alexander Kiselev's Juggler implementation `IRTAlgorithm`
// Subject to the terms in the LICENSE file found in the top-level directory.

#include "IrtCherenkovParticleID.h"

// AlgorithmInit
//---------------------------------------------------------------------------
void eicrecon::IrtCherenkovParticleID::AlgorithmInit(
    CherenkovDetectorCollection*     irt_det_coll,
    std::shared_ptr<spdlog::logger>& logger
    )
{
  // members
  m_irt_det_coll = irt_det_coll;
  m_log          = logger;
  m_init_failed  = false;

  // print the configuration parameters
  m_cfg.Print(m_log, spdlog::level::debug);

  // extract the the relevant `CherenkovDetector`, set to `m_irt_det`
  auto &detectors = m_irt_det_coll->GetDetectors();
  if(detectors.size() == 0) {
    m_log->error("No CherenkovDetectors found in input collection `irt_det_coll`");
    m_init_failed = true;
  }
  if(detectors.size() > 1)
    m_log->warn("IrtCherenkovParticleID currently only supports 1 CherenkovDetector at a time; taking the first");
  auto this_detector = *detectors.begin();
  m_det_name = this_detector.first;
  m_irt_det = this_detector.second;
  m_log->debug("Initializing IrtCherenkovParticleID algorithm for CherenkovDetector '{}'", m_det_name);

  // readout decoding
  m_cell_mask = m_irt_det->GetReadoutCellMask();
  m_log->debug("readout cellMask = {:#X}", m_cell_mask);

  // rebin refractive index tables to have `m_cfg.numRIndexBins` bins
  m_log->trace("Rebinning refractive index tables to have {} bins",m_cfg.numRIndexBins);
  for(auto [rad_name,irt_rad] : m_irt_det->Radiators()) {
    auto ri_lookup_table_orig = irt_rad->m_ri_lookup_table;
    irt_rad->m_ri_lookup_table.clear();
    irt_rad->m_ri_lookup_table = Tools::ApplyFineBinning( ri_lookup_table_orig, m_cfg.numRIndexBins );
    // m_log->trace("- {}", rad_name);
    // for(auto [energy,rindex] : irt_rad->m_ri_lookup_table) m_log->trace("  {:>5} eV   {:<}", energy, rindex);
  }

  // build `m_pid_radiators`, the list of radiators to use for PID
  for(auto [rad_name,irt_rad] : m_irt_det->Radiators())
    if(rad_name!="Filter")
      m_pid_radiators.insert({ std::string(rad_name), irt_rad });

  // check radiators' configuration, and pass it to `m_irt_det`'s radiators
  for(auto [rad_name,irt_rad] : m_pid_radiators) {
    // find `cfg_rad`, the associated `IrtCherenkovParticleIDConfig` radiator
    auto cfg_rad_it = m_cfg.radiators.find(rad_name);
    if(cfg_rad_it != m_cfg.radiators.end()) {
      auto cfg_rad = cfg_rad_it->second;
      // pass `cfg_rad` params to `irt_rad`, the IRT radiator
      irt_rad->m_ID = Tools::GetRadiatorID(std::string(rad_name));
      irt_rad->m_AverageRefractiveIndex = cfg_rad.referenceRIndex;
      irt_rad->SetReferenceRefractiveIndex(cfg_rad.referenceRIndex);
      if(cfg_rad.attenuation>0)
        irt_rad->SetReferenceAttenuationLength(cfg_rad.attenuation);
      if(cfg_rad.smearing>0) {
        if(cfg_rad.smearingMode=="uniform")
          irt_rad->SetUniformSmearing(cfg_rad.smearing);
        else if(cfg_rad.smearingMode=="gaussian")
          irt_rad->SetGaussianSmearing(cfg_rad.smearing);
        else
          m_log->error("Unknown smearing mode '{}' for {} radiator",cfg_rad.smearingMode,rad_name);
      }
      irt_rad->SetTrajectoryBinCount(cfg_rad.zbins - 1);
    } else m_log->error("Cannot find radiator '{}' in IrtCherenkovParticleIDConfig instance",rad_name);
  }

  // get PDG info for the particles we want to identify in PID
  // FIXME: cannot use `TDatabasePDG` since it is not thread safe; until we
  // have a proper PDG database service, we hard-code the masses in Tools.h
  m_log->debug("List of particles for PID:");
  for(auto pdg : m_cfg.pdgList) {
    auto mass = Tools::GetPDGMass(pdg);
    m_pdg_mass.insert({ pdg, mass });
    m_log->debug("  {:>8}  M={} GeV", pdg, mass);
  }

}


// AlgorithmChangeRun
//---------------------------------------------------------------------------
void eicrecon::IrtCherenkovParticleID::AlgorithmChangeRun() {
}


// AlgorithmProcess
//---------------------------------------------------------------------------
std::vector<edm4eic::CherenkovParticleID*> eicrecon::IrtCherenkovParticleID::AlgorithmProcess(
    std::vector<const edm4eic::RawPMTHit*>& in_raw_hits,
    std::map<std::string,std::vector<const edm4eic::TrackSegment*>>& in_charged_particles
    )
{
  // logging
  m_log->trace("{:=^70}"," call IrtCherenkovParticleID::AlgorithmProcess ");
  m_log->trace("number of raw sensor hits: {}", in_raw_hits.size());

  // annoy the user, if a cheat mode is enabled
  m_cfg.PrintCheats(m_log);

  // start output collections
  std::vector<edm4eic::CherenkovParticleID*> out_cherenkov_pids;
  if(m_init_failed) return out_cherenkov_pids;

  // check `in_charged_particles`: each radiator should have the same number of TrackSegments
  long num_charged_particles = -1;
  for(const auto& [rad_name,charged_particle_list] : in_charged_particles) {
    if(num_charged_particles<0) {
      num_charged_particles = charged_particle_list.size();
      m_log->trace("number of reconstructed charged particles: {}", charged_particle_list.size());
    }
    else if(num_charged_particles != charged_particle_list.size()) {
      m_log->error("radiators have differing numbers of TrackSegments");
      return out_cherenkov_pids;
    }
  }


  // loop over raw hits ***************************************************
  // - make `irt_photons`: a list of `OpticalPhoton`s for the IRT algorithm
  std::vector<std::shared_ptr<OpticalPhoton>> irt_photons;
  m_log->trace("{:#<70}","### SENSOR HITS ");
  for(const auto& raw_hit : in_raw_hits) {

    // get sensor and pixel info
    auto cell_id       = raw_hit->getCellID();
    uint64_t sensor_id = cell_id & m_cell_mask;
    TVector3 pixel_pos = m_irt_det->m_ReadoutIDToPosition(cell_id);
    m_log->trace("cell_id={:#X}  copy={}  pixel_pos=( {:>10.2f} {:>10.2f} {:>10.2f} )",
        cell_id, sensor_id, pixel_pos.x(), pixel_pos.y(), pixel_pos.z());
    // FIXME: `pixel_pos` is slightly different from juggler (but who is right?)

    // start new IRT photon
    auto mc_photon  = raw_hit->getPhoton(); // get MC photon (typically only used by cheat modes)
    auto irt_sensor = m_irt_det->m_PhotonDetectors[0]; // FIXME: assumes one sensor type
    auto irt_photon = std::make_shared<OpticalPhoton>();
    irt_photon->SetVolumeCopy(sensor_id);
    irt_photon->SetDetectionPosition(pixel_pos);
    irt_photon->SetPhotonDetector(irt_sensor);
    irt_photon->SetDetected(true);

    // cheat mode: get photon vertex info from MC truth
    if(m_cfg.cheatPhotonVertex || m_cfg.cheatTrueRadiator) {
      irt_photon->SetVertexPosition(Tools::PodioVector3_to_TVector3(mc_photon.getVertex()));
      irt_photon->SetVertexMomentum(Tools::PodioVector3_to_TVector3(mc_photon.getMomentum()));
    }

    // add to `irt_photons`
    irt_photons.push_back(std::move(irt_photon));

  } // end `in_raw_hits` loop


  // loop over charged particles ********************************************
  m_log->trace("{:#<70}","### CHARGED PARTICLES ");
  for(long i=0; i<num_charged_particles; i++) {
    m_log->trace("{:-<70}",fmt::format("--- charged particle #{} ",i));

    // start an `irt_particle` and a local map of `irt_rad` -> charged particle object (for output)
    auto irt_particle = std::make_unique<ChargedParticle>();
    std::unordered_map<CherenkovRadiator*, const edm4eic::TrackSegment*> out_charged_particles;

    // loop over radiators
    for(auto [rad_name,irt_rad] : m_pid_radiators) {

      // start a new IRT `RadiatorHistory` (it will be destroyed when `irt_particle` is destroyed)
      auto irt_rad_history = new RadiatorHistory();
      irt_particle->StartRadiatorHistory({ irt_rad, irt_rad_history });

      // get the `charged_particle` for this radiator
      auto charged_particle_list_it = in_charged_particles.find(rad_name);
      if(charged_particle_list_it == in_charged_particles.end()) {
        m_log->error("Cannot find radiator '{}' in `in_charged_particles`",rad_name);
        continue;
      }
      auto charged_particle_list = charged_particle_list_it->second;
      auto charged_particle      = charged_particle_list.at(i);
      out_charged_particles.insert({ irt_rad, charged_particle });

      // loop over `TrackPoint`s of this `charged_particle`, adding each to the IRT radiator
      irt_rad->ResetLocations();
      m_log->trace("TrackPoints in '{}' radiator:",rad_name);
      for(const auto& point : charged_particle->getPoints()) {
        TVector3 position = Tools::PodioVector3_to_TVector3(point.position);
        TVector3 momentum = Tools::PodioVector3_to_TVector3(point.momentum);
        // irt_rad_history->AddStep( new ChargedParticleStep(position,momentum) ); // FIXME: not needed, if we can just use AddLocation?
        irt_rad->AddLocation(position,momentum);
        m_log->trace(" point: x=( {:>10.2f} {:>10.2f} {:>10.2f} )", position.x(), position.y(), position.z());
        m_log->trace("        p=( {:>10.2f} {:>10.2f} {:>10.2f} )", momentum.x(), momentum.y(), momentum.z());
      }

      // add each `irt_photon` to the radiator history
      // - we don't know which photon goes with which radiator, thus we add them all to each radiator;
      //   the radiators' photons are mixed in ChargedParticle::PIDReconstruction
      for(auto irt_photon : irt_photons) {

        // cheat mode, for testing only: use MC to get the actual radiator
        if(m_cfg.cheatTrueRadiator) {
          auto mc_rad = m_irt_det->GuessRadiator(
                irt_photon->GetVertexPosition(),
                irt_photon->GetVertexPosition() // FIXME: assumes IP is at (0,0,0)
                );
          if(mc_rad != irt_rad) continue;
        }

        // cheat mode: retrieve a refractive index estimate; it is not exactly the one, which 
        // was used in GEANT, but should be very close
        if(m_cfg.cheatPhotonVertex) {
          double ri;
          auto ri_set = Tools::GetFinelyBinnedTableEntry(
              irt_rad->m_ri_lookup_table,
              1e9 * irt_photon->GetVertexMomentum().Mag(),
              &ri
              );
          if(ri_set) irt_photon->SetVertexRefractiveIndex(ri);
        }

        // add copy of `irt_photon` to the history, since ~RadiatorHistory() will destroy it
        irt_rad_history->AddOpticalPhoton(new OpticalPhoton(*irt_photon));
        /* FIXME: this considers ALL of the `irt_photons`... we can limit this
         * once we add the ability to get a fiducial volume for each track, i.e.,
         * a region of sensors where we expect to see this `irt_particle`'s
         * Cherenkov photons; this should also combat sensor noise
         */
      }

    } // end radiator loop


    // particle identification +++++++++++++++++++++++++++++++++++++++++++++++++++++

    // define a mass hypothesis for each particle we want to check
    m_log->trace("{:+^70}"," PARTICLE IDENTIFICATION ");
    CherenkovPID irt_pid;
    std::unordered_map<int,MassHypothesis*> pdg_to_hyp; // `pdg` -> hypothesis
    for(auto [pdg,mass] : m_pdg_mass) {
      irt_pid.AddMassHypothesis(mass);
      pdg_to_hyp.insert({ pdg, irt_pid.GetHypothesis(irt_pid.GetHypothesesCount()-1) });
    }

    // run IRT PID
    irt_particle->PIDReconstruction(irt_pid);
    m_log->trace("{:-^70}"," IRT RESULTS ");

    // loop over radiators
    for(auto [rad_name,irt_rad] : m_pid_radiators) {

      // Cherenkov angle (theta) estimate
      unsigned npe            = 0;
      double   theta_sum      = 0.0;
      double   rindex_sum     = 0.0;
      double   wavelength_sum = 0.0;
      double   weight_sum     = 0.0;
      std::vector<std::pair<double,double>> phot_theta_phi;

      // loop over this radiator's photons, and decide which to include in the theta estimate
      auto irt_rad_history = irt_particle->FindRadiatorHistory(irt_rad);
      for(auto irt_photon : irt_rad_history->Photons()) {

        // check whether this photon was selected by at least one mass hypothesis
        bool selected = false;
        for(auto irt_photon_sel : irt_photon->_m_Selected)
          if(irt_photon_sel.second==irt_rad) { selected=true; break; }
        if(!selected) continue;

        // get this photon's theta and phi estimates
        auto phot_theta = irt_photon->_m_PDF[irt_rad].GetAverage();
        auto phot_phi   = irt_photon->m_Phi[irt_rad];

        // set weight
        auto phot_weight = 1.0;
        // auto phot_weight = fabs(sin(phi));

        // add to the total
        npe++;
        weight_sum += phot_weight;
        theta_sum  += phot_weight * phot_theta;
        phot_theta_phi.push_back({ phot_theta, phot_phi });
        if(m_cfg.cheatPhotonVertex) {
          rindex_sum     += irt_photon->GetVertexRefractiveIndex();
          wavelength_sum += 1239.84 / (1e9 * irt_photon->GetVertexMomentum().Mag()); // [GeV]->[eV] conversion
        }

      } // end loop over this radiator's photons

      // compute averages
      auto ComputeAve     = [] (auto sum, auto count) { return count>0 ? sum / count : 0.0; };
      auto theta_ave      = ComputeAve(theta_sum,      weight_sum);
      auto rindex_ave     = ComputeAve(rindex_sum,     npe);
      auto wavelength_ave = ComputeAve(wavelength_sum, npe);


      // fill output collections -----------------------------------------------

      // fill Cherenkov angle estimate
      edm4eic::MutableCherenkovParticleID out_cherenkov_pid;
      out_cherenkov_pid.setRadiator(   decltype(edm4eic::CherenkovParticleIDData::radiator)   (irt_rad->m_ID)  );
      out_cherenkov_pid.setNpe(        decltype(edm4eic::CherenkovParticleIDData::npe)        (npe)            );
      out_cherenkov_pid.setTheta(      decltype(edm4eic::CherenkovParticleIDData::theta)      (theta_ave)      );
      out_cherenkov_pid.setRindex(     decltype(edm4eic::CherenkovParticleIDData::rindex)     (rindex_ave)     );
      out_cherenkov_pid.setWavelength( decltype(edm4eic::CherenkovParticleIDData::wavelength) (wavelength_ave) );
      for(auto [phot_theta,phot_phi] : phot_theta_phi) {
        edm4hep::Vector2f theta_phi{ float(phot_theta), float(phot_phi) };
        out_cherenkov_pid.addToThetaPhiPhotons(theta_phi);
      }
      m_log->trace("-> {} Radiator (ID={}):", rad_name, irt_rad->m_ID);
      m_log->trace("  Cherenkov Angle Estimate:");
      m_log->trace("    {:>16}:  {:<10}",     "NPE",          npe);
      m_log->trace("    {:>16}:  {:<10.3}",   "<theta>",      theta_ave);
      m_log->trace("    {:>16}:  {:<10.3}",   "<rindex>",     rindex_ave);
      m_log->trace("    {:>16}:  {:<10.3}",   "<wavelength>", wavelength_ave);
      m_log->trace("    {:>16}:  {:<10.3}",   "weight_sum",   weight_sum);

      // relate mass hypotheses
      m_log->trace("  Mass Hypotheses:");
      m_log->trace("    {:>6}  {:>10}  {:>10}", "PDG", "Weight", "NPE");
      for(auto [pdg,mass] : m_pdg_mass) {
        
        // get hypothesis results
        auto irt_hypothesis = pdg_to_hyp.at(pdg);
        auto hyp_weight     = irt_hypothesis->GetWeight(irt_rad);
        auto hyp_npe        = irt_hypothesis->GetNpe(irt_rad);
        m_log->trace("    {:>6}  {:>10.8}  {:>10.8}", pdg, hyp_weight, hyp_npe);

        // fill `ParticleID` output collection
        edm4eic::CherenkovPdgHypothesis out_hypothesis;
        out_hypothesis.pdg    = decltype(edm4eic::CherenkovPdgHypothesis::pdg)    (pdg);
        out_hypothesis.weight = decltype(edm4eic::CherenkovPdgHypothesis::weight) (hyp_weight);
        out_hypothesis.npe    = decltype(edm4eic::CherenkovPdgHypothesis::npe)    (hyp_npe);

        // relate
        out_cherenkov_pid.addToHypotheses(out_hypothesis);

      } // end hypothesis loop

      // relate charged particle
      auto out_charged_particle = *out_charged_particles.at(irt_rad);
      out_cherenkov_pid.setChargedParticle(out_charged_particle);

      // append
      out_cherenkov_pids.push_back(new edm4eic::CherenkovParticleID(out_cherenkov_pid)); // force immutable

    } // end radiator loop

  } // end `in_charged_particles` loop

  return out_cherenkov_pids;
}