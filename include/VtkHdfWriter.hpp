#ifndef _VTKHDFWRITER_HPP_
#define _VTKHDFWRITER_HPP_

/**
 * VtkHdfWriter.hpp
 *
 * Write a whole time series of an unstructured grid into ONE VTKHDF (.vtkhdf)
 * file, instead of a directory full of per-step .vtu/.pvtu files plus a .pvd
 * collection. VTKHDF is a single HDF5 container that ParaView (modern versions)
 * opens natively: all steps, mesh topology, and point-data arrays live in one
 * binary file, so there is exactly one artifact per run to download.
 *
 * This implements the *temporal* VTKHDF "UnstructuredGrid" layout (VTK File
 * Formats, VTKHDF v2.0). Every step is stored as a self-contained partition:
 * its points/cells/connectivity are concatenated into the global datasets and a
 * "Steps" group records, per step, the offsets into those globals plus the time
 * value. Because each step is self-contained, a step's geometry may change from
 * one step to the next (moving walls, re-meshed tubes) -- exactly our use case.
 *
 * Usage (accumulate in memory, then write once):
 *   glymphbie::VtkHdfTimeSeries ts;
 *   for (each step) {
 *     glymphbie::VtkHdfStep s;
 *     s.AppendVTU(vtu_a, "BC", 0);   // merge one or more sctl::VTUData blocks,
 *     s.AppendVTU(vtu_b, "BC", 1);   //   each tagged with an integer ObjectId
 *     ts.AddStep(time, std::move(s));
 *   }
 *   ts.Write("vis/walls");           // -> vis/walls.vtkhdf
 *
 * A VTUData carries a single flat "value" array whose component count is
 * inferred as value.Dim()/(coord.Dim()/3) (the same rule sctl::VTUData uses),
 * stored as the named point-data array; an "ObjectId" scalar is added per node
 * so merged blocks can be split/colored apart in ParaView.
 *
 * Header-only, but requires linking the HDF5 C library (see CMakeLists).
 */

#include <sctl.hpp>
#include <hdf5.h>
#include <cstdint>
#include <cstring>
#include <map>
#include <string>
#include <vector>

namespace glymphbie {

// One time step's worth of unstructured-grid data, assembled by merging one or
// more sctl::VTUData blocks. Offsets are kept in the VTK "offsets" convention
// used by sctl::VTUData -- one END offset per cell, with NO leading zero; the
// writer prepends the single leading 0 the VTKHDF spec requires.
struct VtkHdfStep {
  std::vector<float>   points;        // AoS xyz, length 3*NP
  std::vector<int64_t> connectivity;  // point ids, length NConn
  std::vector<int64_t> offsets;       // per-cell END offset (no leading 0), length NC
  std::vector<uint8_t> types;         // VTK cell type per cell, length NC
  // Named point-data arrays: name -> (ncomp, data[NP*ncomp]). Kept aligned so
  // every array covers every point in the step.
  std::map<std::string, std::pair<int, std::vector<float>>> pdata;

  int64_t NumPoints()   const { return (int64_t)points.size() / 3; }
  int64_t NumCells()    const { return (int64_t)types.size(); }
  int64_t NumConnIds()  const { return (int64_t)connectivity.size(); }

  // Merge one sctl::VTUData block into this step. Its coords/cells are appended
  // (connectivity re-based onto the running point count, offsets re-based onto
  // the running connectivity length), its flat `value` array is stored under
  // `value_name`, and a constant `object_id` is written into an "ObjectId"
  // scalar for the block's nodes.
  void AppendVTU(const sctl::VTUData& vtu, const std::string& value_name, int object_id) {
    const int64_t np0        = NumPoints();               // point-id base for this block
    const int64_t conn_base  = (int64_t)connectivity.size();  // connectivity base for offsets
    const int64_t np         = (int64_t)vtu.coord.Dim() / 3;

    for (sctl::Long i = 0; i < vtu.coord.Dim(); i++) points.push_back((float)vtu.coord[i]);

    for (sctl::Long i = 0; i < vtu.connect.Dim(); i++)
      connectivity.push_back((int64_t)vtu.connect[i] + np0);
    for (sctl::Long i = 0; i < vtu.offset.Dim(); i++)
      offsets.push_back((int64_t)vtu.offset[i] + conn_base);
    for (sctl::Long i = 0; i < vtu.types.Dim(); i++)
      types.push_back((uint8_t)vtu.types[i]);

    // Named value array (component count inferred exactly as sctl::VTUData does).
    const int ncomp = (np > 0) ? (int)(vtu.value.Dim() / np) : 0;
    if (np > 0 && ncomp > 0) {
      auto& arr = pdata[value_name];
      arr.first = ncomp;
      for (sctl::Long i = 0; i < vtu.value.Dim(); i++) arr.second.push_back((float)vtu.value[i]);
    }
    // ObjectId scalar per node of this block.
    auto& oid = pdata["ObjectId"];
    oid.first = 1;
    for (int64_t i = 0; i < np; i++) oid.second.push_back((float)object_id);
  }
};

// Accumulates VtkHdfStep snapshots and writes them as one temporal VTKHDF
// UnstructuredGrid file.
class VtkHdfTimeSeries {
 public:
  void AddStep(double time, VtkHdfStep step) {
    times_.push_back(time);
    steps_.push_back(std::move(step));
  }
  sctl::Long NumSteps() const { return (sctl::Long)steps_.size(); }

  // Write "<fname_prefix>.vtkhdf". Returns true on success. Serial (single
  // partition per step); rank 0 should be the only caller.
  bool Write(const std::string& fname_prefix) const {
    const std::string fname = fname_prefix + ".vtkhdf";
    const size_t NSteps = steps_.size();
    if (NSteps == 0) return false;

    // ---- Concatenate the per-step geometry into the global VTKHDF datasets ----
    std::vector<int64_t> nPoints, nCells, nConnIds;   // per-step counts
    std::vector<float>   Points;
    std::vector<int64_t> Connectivity;
    std::vector<int64_t> Offsets;                     // [0, ends...] per step, concatenated
    std::vector<uint8_t> Types;
    std::vector<double>  Values(times_.begin(), times_.end());
    std::vector<int64_t> partOffsets, numberOfParts, pointOffsets;
    std::vector<int64_t> cellOffsets, connIdOffsets;  // (NSteps,1) row-major

    // Point-data arrays are stored concatenated across steps too; collect the
    // union of names (with component counts) so every array spans every step.
    std::map<std::string, int> pd_ncomp;
    for (const auto& s : steps_)
      for (const auto& kv : s.pdata) pd_ncomp[kv.first] = kv.second.first;
    std::map<std::string, std::vector<float>>   pd_data;
    std::map<std::string, std::vector<int64_t>> pd_offsets;  // Steps/PointDataOffsets/<name>

    int64_t cumPoints = 0, cumCells = 0, cumConn = 0;
    for (size_t t = 0; t < NSteps; t++) {
      const VtkHdfStep& s = steps_[t];
      const int64_t NP = s.NumPoints(), NC = s.NumCells(), NCI = s.NumConnIds();

      nPoints.push_back(NP);
      nCells.push_back(NC);
      nConnIds.push_back(NCI);

      Points.insert(Points.end(), s.points.begin(), s.points.end());
      Connectivity.insert(Connectivity.end(), s.connectivity.begin(), s.connectivity.end());
      Types.insert(Types.end(), s.types.begin(), s.types.end());
      // Per-step offsets: the VTKHDF spec wants NC+1 entries starting at 0.
      Offsets.push_back(0);
      Offsets.insert(Offsets.end(), s.offsets.begin(), s.offsets.end());

      // Steps metadata (single partition per step).
      partOffsets.push_back((int64_t)t);
      numberOfParts.push_back(1);
      pointOffsets.push_back(cumPoints);
      cellOffsets.push_back(cumCells);
      connIdOffsets.push_back(cumConn);

      // Point data (pad any array a step happens to be missing with zeros so all
      // arrays stay length NP for this step).
      for (const auto& nc : pd_ncomp) {
        const std::string& name = nc.first;
        const int ncomp = nc.second;
        pd_offsets[name].push_back(cumPoints);
        std::vector<float>& dst = pd_data[name];
        auto it = s.pdata.find(name);
        if (it != s.pdata.end() && (int64_t)it->second.second.size() == NP * ncomp) {
          dst.insert(dst.end(), it->second.second.begin(), it->second.second.end());
        } else {
          dst.insert(dst.end(), (size_t)(NP * ncomp), 0.0f);
        }
      }

      cumPoints += NP;
      cumCells  += NC;
      cumConn   += NCI;
    }

    // ---- Write the HDF5 container ----
    hid_t file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) return false;
    hid_t root = H5Gcreate2(file, "VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    { // Version attribute = [2, 0]
      const int64_t ver[2] = {2, 0};
      const hsize_t d = 2;
      hid_t sp = H5Screate_simple(1, &d, nullptr);
      hid_t at = H5Acreate2(root, "Version", H5T_STD_I64LE, sp, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(at, H5T_NATIVE_INT64, ver);
      H5Aclose(at); H5Sclose(sp);
    }
    { // Type attribute = "UnstructuredGrid" (fixed-length ASCII, matches VTK's writer)
      const char* type = "UnstructuredGrid";
      hid_t st = H5Tcopy(H5T_C_S1);
      H5Tset_size(st, std::strlen(type));
      H5Tset_strpad(st, H5T_STR_NULLPAD);
      H5Tset_cset(st, H5T_CSET_ASCII);
      hid_t sp = H5Screate(H5S_SCALAR);
      hid_t at = H5Acreate2(root, "Type", st, sp, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(at, st, type);
      H5Aclose(at); H5Sclose(sp); H5Tclose(st);
    }

    write_1d(root, "NumberOfPoints",          H5T_STD_I64LE, H5T_NATIVE_INT64, nPoints.data(),  nPoints.size());
    write_1d(root, "NumberOfCells",           H5T_STD_I64LE, H5T_NATIVE_INT64, nCells.data(),   nCells.size());
    write_1d(root, "NumberOfConnectivityIds", H5T_STD_I64LE, H5T_NATIVE_INT64, nConnIds.data(), nConnIds.size());
    write_2d(root, "Points", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, Points.data(), Points.size() / 3, 3);
    write_1d(root, "Connectivity", H5T_STD_I64LE, H5T_NATIVE_INT64, Connectivity.data(), Connectivity.size());
    write_1d(root, "Offsets",      H5T_STD_I64LE, H5T_NATIVE_INT64, Offsets.data(),      Offsets.size());
    write_1d(root, "Types",        H5T_STD_U8LE,  H5T_NATIVE_UINT8, Types.data(),        Types.size());

    // Point-data arrays.
    hid_t pd_group = H5Gcreate2(root, "PointData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (const auto& nc : pd_ncomp) {
      const std::string& name = nc.first;
      const int ncomp = nc.second;
      const std::vector<float>& d = pd_data[name];
      if (ncomp == 1) write_1d(pd_group, name.c_str(), H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, d.data(), d.size());
      else            write_2d(pd_group, name.c_str(), H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, d.data(), d.size() / ncomp, ncomp);
    }
    H5Gclose(pd_group);

    // Steps group (temporal metadata).
    hid_t steps = H5Gcreate2(root, "Steps", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    { const int64_t ns = (int64_t)NSteps;
      hid_t sp = H5Screate(H5S_SCALAR);
      hid_t at = H5Acreate2(steps, "NSteps", H5T_STD_I64LE, sp, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(at, H5T_NATIVE_INT64, &ns);
      H5Aclose(at); H5Sclose(sp);
    }
    write_1d(steps, "Values",        H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, Values.data(),        Values.size());
    write_1d(steps, "PartOffsets",   H5T_STD_I64LE,  H5T_NATIVE_INT64,  partOffsets.data(),   partOffsets.size());
    write_1d(steps, "NumberOfParts", H5T_STD_I64LE,  H5T_NATIVE_INT64,  numberOfParts.data(), numberOfParts.size());
    write_1d(steps, "PointOffsets",  H5T_STD_I64LE,  H5T_NATIVE_INT64,  pointOffsets.data(),  pointOffsets.size());
    write_2d(steps, "CellOffsets",          H5T_STD_I64LE, H5T_NATIVE_INT64, cellOffsets.data(),   cellOffsets.size(),   1);
    write_2d(steps, "ConnectivityIdOffsets", H5T_STD_I64LE, H5T_NATIVE_INT64, connIdOffsets.data(), connIdOffsets.size(), 1);

    hid_t pdo = H5Gcreate2(steps, "PointDataOffsets", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (const auto& kv : pd_offsets)
      write_1d(pdo, kv.first.c_str(), H5T_STD_I64LE, H5T_NATIVE_INT64, kv.second.data(), kv.second.size());
    H5Gclose(pdo);

    H5Gclose(steps);
    H5Gclose(root);
    H5Fclose(file);
    return true;
  }

 private:
  // Create + write a 1D dataset of length n (no-op if n == 0, so ParaView never
  // sees a zero-length required array).
  static void write_1d(hid_t loc, const char* name, hid_t ftype, hid_t mtype,
                        const void* data, size_t n) {
    if (n == 0) return;
    const hsize_t dim = (hsize_t)n;
    hid_t sp = H5Screate_simple(1, &dim, nullptr);
    hid_t ds = H5Dcreate2(loc, name, ftype, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(ds); H5Sclose(sp);
  }
  // Create + write a 2D dataset of shape (n, m).
  static void write_2d(hid_t loc, const char* name, hid_t ftype, hid_t mtype,
                        const void* data, size_t n, size_t m) {
    if (n == 0) return;
    const hsize_t dims[2] = {(hsize_t)n, (hsize_t)m};
    hid_t sp = H5Screate_simple(2, dims, nullptr);
    hid_t ds = H5Dcreate2(loc, name, ftype, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(ds); H5Sclose(sp);
  }

  std::vector<double>     times_;
  std::vector<VtkHdfStep> steps_;
};

}  // namespace glymphbie

#endif  // _VTKHDFWRITER_HPP_
