// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Xin Dong

#pragma once

#include <optional>
#include <utility>

namespace eicrecon {

/// Daughter PDG codes (positive, negative) for a supported V0 mass hypothesis
inline std::optional<std::pair<int, int>> v0_daughter_pdgs(int v0_pdg) {
  switch (v0_pdg) {
  case 310:
    return std::pair{211, -211}; // K0S -> pi+ pi-
  case 3122:
    return std::pair{2212, -211}; // Lambda -> p pi-
  case -3122:
    return std::pair{211, -2212}; // anti-Lambda -> pi+ pbar
  default:
    return std::nullopt;
  }
}

} // namespace eicrecon
