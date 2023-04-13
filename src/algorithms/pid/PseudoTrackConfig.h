// Copyright 2023, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.

#pragma once

#include <spdlog/spdlog.h>

namespace eicrecon {

  class PseudoTrackConfig {
    public:

      /////////////////////////////////////////////////////
      // CONFIGURATION PARAMETERS
      //   NOTE: some defaults are hard-coded here; override externally

      unsigned numPoints = 10; // number of pseudo-track points; set to 0 to use ALL of them

      //
      /////////////////////////////////////////////////////

      // print all parameters
      void Print(
          std::shared_ptr<spdlog::logger> m_log,
          spdlog::level::level_enum lvl=spdlog::level::debug
          )
      {
        m_log->log(lvl, "{:=^60}"," PseudoTrackConfig Settings ");
        auto print_param = [&m_log, &lvl] (auto name, auto val) {
          m_log->log(lvl, "  {:>20} = {:<}", name, val);
        };
        print_param("numPoints",numPoints);
        m_log->log(lvl, "{:=^60}","");
      }

  };
}
