#ifndef STARFIELD_CATALOG_H
#define STARFIELD_CATALOG_H

#include "types.h"
#include "constants.h"

#include <cstdint>
#include <string>
#include <vector>

namespace starfield {

// ---------------------------------------------------------------------------
// Delta-encoded Hipparcos Catalog
// ---------------------------------------------------------------------------

/// Encode a sorted star catalog to delta-encoded binary format.
/// Stars MUST be sorted by RA before encoding.
/// @param stars  sorted star catalog
/// @returns binary blob (header + delta records)
std::vector<uint8_t> encodeCatalog(const std::vector<Star>& stars);

/// Decode a delta-encoded binary catalog back to full Star records.
/// @param data  binary blob from encodeCatalog()
/// @returns decoded star catalog
std::vector<Star> decodeCatalog(const uint8_t* data, size_t size);
std::vector<Star> decodeCatalog(const std::vector<uint8_t>& data);

/// Sort stars by RA for optimal delta encoding.
void sortByRA(std::vector<Star>& stars);

/// Compute compression statistics.
struct CompressionStats {
    size_t rawSize       = 0;  // bytes if stored uncompressed
    size_t encodedSize   = 0;  // bytes delta-encoded
    double ratio         = 0;  // compression ratio
    double maxDeltaRA_arcsec  = 0;  // max RA delta [arcsec]
    double maxDeltaDec_arcsec = 0;  // max Dec delta [arcsec]
    int    overflows     = 0;  // deltas that exceeded int16 range
};

CompressionStats computeStats(const std::vector<Star>& stars);

// ---------------------------------------------------------------------------
// Star Catalog Queries
// ---------------------------------------------------------------------------

/// Find all stars within a cone (for FOV queries).
/// @param catalog  decoded star catalog
/// @param center   cone center unit vector (ICRS J2000)
/// @param halfAngle  cone half-angle [rad]
/// @param magLimit  faintest magnitude to include
/// @returns matching stars, sorted by magnitude (brightest first)
std::vector<Star> queryCone(
    const std::vector<Star>& catalog,
    const Vector3& center,
    double halfAngle,
    double magLimit = 6.5);

/// Find stars in a rectangular FOV (for star tracker).
/// @param catalog  decoded star catalog
/// @param boresight  boresight unit vector (ICRS J2000)
/// @param up  up direction in sensor frame
/// @param fovX  horizontal FOV [rad]
/// @param fovY  vertical FOV [rad]
/// @param magLimit  faintest magnitude
std::vector<Star> queryRectFOV(
    const std::vector<Star>& catalog,
    const Vector3& boresight,
    const Vector3& up,
    double fovX, double fovY,
    double magLimit = 6.5);

/// Apply proper motion to a star direction vector.
/// @param star  star record
/// @param epochYears  years since J2000.0
/// @returns corrected unit direction vector
Vector3 applyProperMotion(const Star& star, double epochYears);

// ---------------------------------------------------------------------------
// GPU Buffer Generation
// ---------------------------------------------------------------------------

/// Convert star catalog to GPU-ready vertex buffer for Cesium rendering.
/// @param stars  catalog stars
/// @param magLimit  faintest star to include in vertex buffer
/// @returns vertex buffer (StarVertex[])
std::vector<StarVertex> toVertexBuffer(
    const std::vector<Star>& stars,
    float magLimit = 8.0f);

}  // namespace starfield

#endif  // STARFIELD_CATALOG_H
