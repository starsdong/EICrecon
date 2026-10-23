// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Xin Dong

#include <DD4hep/Fields.h>
#include <Evaluator/DD4hepUnits.h>
#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>
#include <KFVertex.h>
#include <edm4eic/Cov4f.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/TrackCollection.h>
#include <edm4eic/TrajectoryCollection.h>
#include <edm4eic/unit_system.h>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

#include "algorithms/reco/SecondaryVertexKFParticle.h"
#include "algorithms/reco/SecondaryVertexKFParticleConfig.h"

namespace eicrecon {

namespace {

  // KFParticle works in cm, GeV and kGauss; EDM4eic uses mm and GeV
  constexpr double kMmToCm = edm4eic::unit::mm / edm4eic::unit::cm;
  constexpr double kCmToMm = edm4eic::unit::cm / edm4eic::unit::mm;

  /**
   * Build a KFPTrack from the perigee parameters of the first track attached to
   * a reconstructed particle. The ACTS bound covariance [l0,l1,phi,theta,q/p]
   * is propagated to the Cartesian basis [x,y,z,px,py,pz] that KFParticle uses.
   * Track::positionMomentumCovariance is not filled by ActsToTracks, so the
   * bound covariance is the only full covariance available.
   */
  std::optional<KFPTrack> make_kfp_track(const edm4eic::ReconstructedParticle& p) {
    if (p.getTracks().empty()) {
      return std::nullopt;
    }
    const auto track      = p.getTracks()[0];
    const auto trajectory = track.getTrajectory();
    if (!trajectory.isAvailable() || trajectory.getTrackParameters().empty()) {
      return std::nullopt;
    }
    const auto pars = trajectory.getTrackParameters()[0];

    const double qop = pars.getQOverP();
    if (!std::isfinite(qop) || qop == 0.) {
      return std::nullopt;
    }
    const double q     = std::copysign(1., qop);
    const double pAbs  = 1. / std::abs(qop);
    const double l0    = pars.getLoc().a; // mm
    const double phi   = pars.getPhi();
    const double theta = pars.getTheta();
    const double sphi = std::sin(phi), cphi = std::cos(phi);
    const double sth = std::sin(theta), cth = std::cos(theta);

    // Jacobian d(x,y,z,px,py,pz)/d(l0,l1,phi,theta,q/p) for an ACTS perigee
    // surface along z: x = -l0 sin(phi), y = l0 cos(phi), z = l1 (+ center)
    const double dpdqop = -q * pAbs * pAbs;
    const std::array<std::array<double, 5>, 6> J{{
        {-sphi, 0., -l0 * cphi, 0., 0.},
        {cphi, 0., -l0 * sphi, 0., 0.},
        {0., 1., 0., 0., 0.},
        {0., 0., -pAbs * sth * sphi, pAbs * cth * cphi, dpdqop * sth * cphi},
        {0., 0., pAbs * sth * cphi, pAbs * cth * sphi, dpdqop * sth * sphi},
        {0., 0., 0., -pAbs * sth, dpdqop * cth},
    }};
    const auto& cb = pars.getCovariance();

    // C = J Cb J^T, packed lower-triangular as KFParticle expects, with
    // positions converted from mm to cm
    std::array<float, 21> C{};
    for (std::size_t i = 0, k = 0; i < 6; ++i) {
      for (std::size_t j = 0; j <= i; ++j, ++k) {
        double c = 0.;
        for (std::size_t a = 0; a < 5; ++a) {
          for (std::size_t b = 0; b < 5; ++b) {
            c += J[i][a] * cb(a, b) * J[j][b];
          }
        }
        if (i < 3) {
          c *= kMmToCm;
        }
        if (j < 3) {
          c *= kMmToCm;
        }
        C[k] = static_cast<float>(c);
      }
    }

    const auto pos = track.getPosition(); // mm, at the perigee
    KFPTrack kt;
    kt.SetParameters(pos.x * kMmToCm, pos.y * kMmToCm, pos.z * kMmToCm, pAbs * sth * cphi,
                     pAbs * sth * sphi, pAbs * cth);
    kt.SetCovarianceMatrix(C.data());
    kt.SetCharge(static_cast<int>(q));
    kt.SetNDF(static_cast<int>(track.getNdf()));
    kt.SetChi2(track.getChi2());
    return kt;
  }

} // namespace

/**
   * @brief Initialize the SecondaryVertexKFParticle Algorithm
   *
   * KFParticle is built with HomogeneousField: the field is a single static Bz
   * value, taken here from the detector field map at the nominal IP.
   */
void SecondaryVertexKFParticle::init() {
  auto field           = m_det->field().magneticField({0., 0., 0.});
  const float bz_tesla = field.z() / dd4hep::tesla;
  KFParticle::SetField(bz_tesla * 10.); // kGauss
  debug("KFParticle homogeneous field Bz = {} tesla", bz_tesla);
}

/**
   * @brief Produce a list of secondary vertex candidates
   *
   * @param rcvtx  - input collection of all vertex candidates
   * @return edm4eic::VertexCollection
   */
void SecondaryVertexKFParticle::process(const SecondaryVertexKFParticle::Input& input,
                                        const SecondaryVertexKFParticle::Output& output) const {
  const auto [rcvtx, rcparts]   = input;
  auto [out_secondary_vertices] = output;

  if (rcvtx->empty()) {
    debug("No primary vertex in this event");
    return;
  }

  // Primary vertex as a KFVertex, in cm
  const auto pv     = (*rcvtx)[0];
  const auto& pvCov = pv.getPositionError();
  KFPVertex kfpPv;
  kfpPv.SetXYZ(pv.getPosition().x * kMmToCm, pv.getPosition().y * kMmToCm,
               pv.getPosition().z * kMmToCm);
  kfpPv.SetCovarianceMatrix(pvCov(0, 0) * kMmToCm * kMmToCm, pvCov(1, 0) * kMmToCm * kMmToCm,
                            pvCov(1, 1) * kMmToCm * kMmToCm, pvCov(2, 0) * kMmToCm * kMmToCm,
                            pvCov(2, 1) * kMmToCm * kMmToCm, pvCov(2, 2) * kMmToCm * kMmToCm);
  kfpPv.SetChi2(pv.getChi2());
  kfpPv.SetNDF(pv.getNdf());
  const KFVertex pVtx(kfpPv);

  debug("Primary vertex = ({},{},{})cm", pVtx.GetX(), pVtx.GetY(), pVtx.GetZ());

  std::vector<KFParticle> kVec;
  std::vector<unsigned int> indexVec;
  for (unsigned int i = 0; i < rcparts->size(); ++i) {
    const auto& p = (*rcparts)[i];
    if (p.getCharge() == 0) {
      continue;
    }
    const auto kt = make_kfp_track(p);
    if (!kt) {
      continue;
    }
    // Mass hypothesis from the particle PDG; KFParticle falls back to the
    // pion mass for codes it does not know (including 0)
    KFParticle kp(*kt, p.getPDG());
    double dca = kp.GetDistanceFromVertex(pVtx) * edm4eic::unit::cm;
    if (dca < m_cfg.minDca) {
      continue;
    }

    kVec.push_back(kp);
    indexVec.push_back(i);
  }

  debug("\tVector size {}, {}", kVec.size(), indexVec.size());

  for (unsigned int i1 = 0; i1 < kVec.size(); ++i1) {
    for (unsigned int i2 = i1 + 1; i2 < kVec.size(); ++i2) {
      const auto& p1 = (*rcparts)[indexVec[i1]];
      const auto& p2 = (*rcparts)[indexVec[i2]];

      if (!(m_cfg.unlikesign && p1.getCharge() + p2.getCharge() == 0))
        continue;

      const auto& k1 = kVec[i1];
      const auto& k2 = kVec[i2];

      double dca12 = k1.GetDistanceFromParticle(k2) * edm4eic::unit::cm;
      if (std::isnan(dca12))
        continue;
      if (dca12 > m_cfg.maxDca12)
        continue;

      // KFParticle V0 vertex fit
      KFParticle v0;
      const KFParticle* daughters[2] = {&k1, &k2};
      v0.Construct(daughters, 2);
      if (v0.GetNDF() <= 0 || !std::isfinite(v0.GetChi2()))
        continue;
      if (v0.GetChi2() / v0.GetNDF() > m_cfg.maxChi2Ndf)
        continue;

      // Pointing: angle between V0 momentum and PV -> V0 decay direction
      const double dx     = v0.GetX() - pVtx.GetX();
      const double dy     = v0.GetY() - pVtx.GetY();
      const double dz     = v0.GetZ() - pVtx.GetZ();
      const double decayL = std::sqrt(dx * dx + dy * dy + dz * dz); // cm
      const double pv0    = v0.GetP();
      if (decayL <= 0. || pv0 <= 0.)
        continue;
      const double costheta = (dx * v0.GetPx() + dy * v0.GetPy() + dz * v0.GetPz()) / decayL / pv0;
      if (costheta < m_cfg.minCostheta)
        continue;

      double dca2pv = v0.GetDistanceFromVertex(pVtx) * edm4eic::unit::cm;
      if (dca2pv > m_cfg.maxDca)
        continue;

      const double beta = pv0 / v0.GetE();
      const double time = decayL * dd4hep::cm / (beta * dd4hep::c_light) /
                          dd4hep::ns; // ns, as for the helix finder

      constexpr double kCm2ToMm2 = kCmToMm * kCmToMm;
      const edm4eic::Cov4f cov(
          v0.GetCovariance(0, 0) * kCm2ToMm2, v0.GetCovariance(1, 1) * kCm2ToMm2,
          v0.GetCovariance(2, 2) * kCm2ToMm2, 0., v0.GetCovariance(0, 1) * kCm2ToMm2,
          v0.GetCovariance(0, 2) * kCm2ToMm2, 0., v0.GetCovariance(1, 2) * kCm2ToMm2, 0., 0.);

      auto vtx = out_secondary_vertices->create();
      vtx.setType(2); // 2 for secondary
      vtx.setChi2(v0.GetChi2());
      vtx.setNdf(v0.GetNDF());
      vtx.setPosition({(float)(v0.GetX() * kCmToMm), (float)(v0.GetY() * kCmToMm),
                       (float)(v0.GetZ() * kCmToMm), (float)time});
      vtx.setPositionError(cov);
      vtx.addToAssociatedParticles(p1);
      vtx.addToAssociatedParticles(p2);

      debug("One secondary vertex found at (x,y,z) = ({}, {}, {}) mm, chi2/ndf = {}/{}.",
            v0.GetX() * kCmToMm, v0.GetY() * kCmToMm, v0.GetZ() * kCmToMm, v0.GetChi2(),
            v0.GetNDF());

    } // end i2
  } // end i1

} // end process

} // namespace eicrecon
