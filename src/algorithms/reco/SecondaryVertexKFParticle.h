// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Xin Dong

#pragma once

#include <DD4hep/Detector.h>
#include <algorithms/algorithm.h>
#include <algorithms/geo.h>
#include <edm4eic/ReconstructedParticleCollection.h>
#include <edm4eic/VertexCollection.h>
#include <string>      // for basic_string
#include <string_view> // for string_view

#include "algorithms/interfaces/WithPodConfig.h"
#include "algorithms/reco/SecondaryVertexKFParticleConfig.h"

namespace eicrecon {

using SecondaryVertexKFParticleAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4eic::VertexCollection, edm4eic::ReconstructedParticleCollection>,
    algorithms::Output<edm4eic::VertexCollection>>;

/**
 * @brief Two-track secondary (V0) vertex finder using the KFParticle package
 *
 * Same inputs, outputs and selection logic as SecondaryVerticesHelix, but the
 * daughters are built as KFParticle objects from the full perigee track
 * covariance and the V0 vertex is obtained from a KFParticle fit.
 */
class SecondaryVertexKFParticle : public SecondaryVertexKFParticleAlgorithm,
                                  public WithPodConfig<SecondaryVertexKFParticleConfig> {

public:
  SecondaryVertexKFParticle(std::string_view name)
      : SecondaryVertexKFParticleAlgorithm{name,
                                           {"inputVertices", "inputParticles"},
                                           {"outputSecondaryVertices"},
                                           "Reconstruct secondary vertices with KFParticle"} {}

  void init() final;
  void process(const Input&, const Output&) const final;

private:
  const dd4hep::Detector* m_det{algorithms::GeoSvc::instance().detector()};
};

} // namespace eicrecon
