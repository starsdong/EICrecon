// Copyright 2023, Alexander Kiselev, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.
//
// Ported from Juggler's JugPID `IRTAlgorithmServices`
//

#pragma once

// general
#include <map>
#include <math.h>

// ROOT
#include <TVector3.h>

// IRT
#include <IRT/ParametricSurface.h>


namespace eicrecon {


  // local PDG mass database
  // FIXME: cannot use `TDatabasePDG` since it is not thread safe; until we
  // have a proper PDG database service, we hard-code the masses we need;
  // use Tools::GetPDGMass for access
  const std::unordered_map<int,double> g_pdg_db_for_pid {
    { -11,  0.000510999 },
    { 211,  0.13957     },
    { 321,  0.493677    },
    { 2212, 0.938272    }
  };

  const std::unordered_map<int,std::string> g_radiator_ids {
    { 0, "Aerogel" },
    { 1, "Gas" }
  };


  // Tools class, filled with miscellanous helper functions
  class Tools {
    public:

      // -------------------------------------------------------------------------------------
      // Radiator IDs
      static std::string GetRadiatorName(int id) {
        std::string name;
        try { name = g_radiator_ids.at(id); }
        catch(const std::out_of_range& e) {
          throw std::runtime_error(fmt::format("RUNTIME ERROR: unknown radiator ID={} in algorithms/pid/Tools::GetRadiatorName",id));
        }
        return name;
      }
      static int GetRadiatorID(std::string name) {
        for(auto& [id,name_tmp] : g_radiator_ids)
          if(name==name_tmp) return id;
        throw std::runtime_error(fmt::format("RUNTIME ERROR: unknown radiator '{}' in algorithms/pid/Tools::GetRadiatorID",name));
        return -1;
      }
      static std::unordered_map<int,std::string> GetRadiatorIDs() { return g_radiator_ids; };


      // -------------------------------------------------------------------------------------
      // PDG mass lookup
      static double GetPDGMass(int pdg) {
        double mass;
        try { mass = g_pdg_db_for_pid.at(pdg); }
        catch(const std::out_of_range& e) {
          throw std::runtime_error(fmt::format("RUNTIME ERROR: unknown PDG={} in algorithms/pid/Tools::GetPDGMass",pdg));
        }
        return mass;
      }
      static int GetNumPDGs() { return g_pdg_db_for_pid.size(); };


      // -------------------------------------------------------------------------------------
      static std::vector<std::pair<double, double>> ApplyFineBinning(
          const std::vector<std::pair<double,double>> &input,
          unsigned nbins
          )
      {
        std::vector<std::pair<double, double>> ret;

        // Well, could have probably just reordered the initial vector;
        std::map<double, double> buffer;

        for(auto entry: input)
          buffer[entry.first] = entry.second;

        // Sanity checks; return empty map in case do not pass them;
        if (buffer.size() < 2 || nbins < 2) return ret;

        double from = (*buffer.begin()).first;
        double to   = (*buffer.rbegin()).first;
        // Will be "nbins+1" equidistant entries;
        double step = (to - from) / nbins;

        for(auto entry: buffer) {
          double e1 = entry.first;
          double qe1 = entry.second;

          if (!ret.size())
            ret.push_back(std::make_pair(e1, qe1));
          else {
            const auto &prev = ret[ret.size()-1];

            double e0 = prev.first;
            double qe0 = prev.second;
            double a = (qe1 - qe0) / (e1 - e0);
            double b = qe0 - a*e0;
            // FIXME: check floating point accuracy when moving to a next point; do we actually
            // care whether the overall number of bins will be "nbins+1" or more?;
            for(double e = e0+step; e<e1; e+=step)
              ret.push_back(std::make_pair(e, a*e + b));
          } //if
        } //for entry

        return ret;
      }


      // -------------------------------------------------------------------------------------
      static bool GetFinelyBinnedTableEntry(
          const std::vector<std::pair<double, double>> &table,
          double argument,
          double *entry
          )
      {
        // Get the tabulated table reference; perform sanity checks;
        //const std::vector<std::pair<double, double>> &qe = u_quantumEfficiency.value();
        unsigned dim = table.size(); if (dim < 2) return false;

        // Find a proper bin; no tricks, they are all equidistant;
        auto const &from = table[0];
        auto const &to = table[dim-1];
        double emin = from.first;
        double emax = to.first;
        double step = (emax - emin) / (dim - 1);
        int ibin = (int)floor((argument - emin) / step);

        //printf("%f vs %f, %f -> %d\n", ev, from.first, to. first, ibin);

        // Out of range check;
        if (ibin < 0 || ibin >= int(dim)) return false;

        *entry = table[ibin].second;
        return true;
      }

      // -------------------------------------------------------------------------------------
      // convert PODIO vector datatype to ROOT TVector3
      template<class PodioVector3>
        static TVector3 PodioVector3_to_TVector3(const PodioVector3 v) {
          return TVector3(v.x, v.y, v.z);
        }
      // convert ROOT::Math::Vector to ROOT TVector3
      template<class MathVector3>
        static TVector3 MathVector3_to_TVector3(MathVector3 v) {
          return TVector3(v.x(), v.y(), v.z());
        }

      // -------------------------------------------------------------------------------------
      // printing: convert objects to strings
      static std::string TVector3_to_string(std::string name, TVector3 v) {
        return fmt::format("{:>30} = ( {:>10.2f} {:>10.2f} {:>10.2f} )", name, v.x(), v.y(), v.z());
      }

  }; // class Tools
} // namespace eicrecon