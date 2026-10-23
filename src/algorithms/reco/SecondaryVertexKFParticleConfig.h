// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Xin Dong

#pragma once

#include <edm4eic/unit_system.h>
#include <vector>

namespace eicrecon {

struct SecondaryVertexKFParticleConfig {

  bool unlikesign   = true;
  float minDca      = 0.03 * edm4eic::unit::mm; // mm, daughter to pVtx
  float maxDca12    = 1. * edm4eic::unit::mm;   // mm, dca between daughter 1 and 2
  float maxDca      = 1. * edm4eic::unit::mm;   // mm, dca of V0 to pVtx
  float minCostheta = 0.8; // costheta, theta: angle of V0 decay direction and momentum
  float maxChi2Ndf  = 10.; // chi2/ndf of the KFParticle V0 vertex fit

  // V0 mass hypotheses refitted for every vertex: K0S, Lambda, anti-Lambda
  std::vector<int> v0Hypotheses = {310, 3122, -3122};
};

} // namespace eicrecon
