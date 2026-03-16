#ifndef STARFIELD_TYPES_H
#define STARFIELD_TYPES_H

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace starfield {

using Vector3 = std::array<double, 3>;
using Quaternion = std::array<double, 4>;  // [w, x, y, z]
using Matrix3x3 = std::array<std::array<double, 3>, 3>;

// ---------------------------------------------------------------------------
// Star record (full precision, in-memory)
// ---------------------------------------------------------------------------

struct Star {
    uint32_t hipId       = 0;       // Hipparcos catalog ID
    double   ra          = 0.0;     // Right ascension J2000 [rad]
    double   dec         = 0.0;     // Declination J2000 [rad]
    float    vmag        = 0.0f;    // Visual magnitude (Hp)
    float    bv          = 0.0f;    // B-V color index
    float    pmRA        = 0.0f;    // Proper motion RA*cos(dec) [mas/yr]
    float    pmDec       = 0.0f;    // Proper motion Dec [mas/yr]
    float    parallax    = 0.0f;    // Parallax [mas]
    uint8_t  spectralIdx = 0;       // Spectral type index (O=0..M=6)

    // Derived: J2000 unit direction vector (ICRS)
    Vector3 direction() const;
};

// ---------------------------------------------------------------------------
// Delta-encoded catalog (compact binary for WASM)
// ---------------------------------------------------------------------------
// Layout:
//   Header (16 bytes):
//     magic[4] = "HIPS"
//     version  = 1
//     count    = number of stars
//     flags    = encoding flags
//
//   Reference star (full precision):
//     ra_ref   (float64)
//     dec_ref  (float64)
//
//   Per-star records (sorted by RA for spatial coherence):
//     delta_ra   (int16)  — offset from previous star [0.01 arcsec units]
//     delta_dec  (int16)  — offset from previous star [0.01 arcsec units]
//     vmag       (int16)  — magnitude * 1000 (milli-mag)
//     bv         (int16)  — color index * 1000
//     pmRA       (int16)  — proper motion * 10 [0.1 mas/yr]
//     pmDec      (int16)  — proper motion * 10 [0.1 mas/yr]
//     parallax   (uint16) — parallax * 100 [0.01 mas]
//     hipId      (uint32) — catalog ID
//   = 20 bytes per star
//
//   Total for 118,218 stars: ~2.25 MB (vs ~8 MB uncompressed)

struct CatalogHeader {
    char     magic[4]   = {'H','I','P','S'};
    uint16_t version    = 1;
    uint32_t count      = 0;
    uint16_t flags      = 0;
    double   ra_ref     = 0.0;
    double   dec_ref    = 0.0;
};

struct DeltaStarRecord {
    int16_t  delta_ra;     // 0.01 arcsec units from previous
    int16_t  delta_dec;    // 0.01 arcsec units from previous
    int16_t  vmag;         // milli-magnitudes
    int16_t  bv;           // milli-color-index
    int16_t  pmRA;         // 0.1 mas/yr
    int16_t  pmDec;        // 0.1 mas/yr
    uint16_t parallax;     // 0.01 mas
    uint32_t hipId;
};

static_assert(sizeof(DeltaStarRecord) == 20, "DeltaStarRecord must be 20 bytes");

// ---------------------------------------------------------------------------
// Star tracker types
// ---------------------------------------------------------------------------

/// Star tracker field of view configuration
struct StarTrackerFOV {
    double fovX       = 8.0;    // horizontal FOV [deg]
    double fovY       = 8.0;    // vertical FOV [deg]
    double magLimit   = 6.5;    // faintest detectable star
    int    maxStars   = 50;     // max stars to report
    double exclusionAngleSun  = 30.0;  // sun exclusion [deg]
    double exclusionAngleMoon = 10.0;  // moon exclusion [deg]
    double exclusionAngleEarth = 20.0; // earth limb exclusion [deg]
};

/// Observed star in sensor frame
struct ObservedStar {
    uint32_t hipId;
    double   sensorX;     // x pixel coordinate (normalized -1..1)
    double   sensorY;     // y pixel coordinate (normalized -1..1)
    double   magnitude;
    Vector3  bodyVector;  // unit vector in body frame
};

/// Star tracker observation result
struct StarTrackerObservation {
    std::vector<ObservedStar> stars;
    int    starsInFOV     = 0;
    int    starsDetected  = 0;
    double epoch          = 0.0;  // Julian Date
};

/// Attitude determination result
struct AttitudeResult {
    Quaternion quaternion = {1,0,0,0};  // ECI-to-body quaternion
    Matrix3x3  dcm       = {};          // Direction Cosine Matrix
    double     error_arcsec = 0.0;      // attitude error estimate [arcsec]
    bool       valid      = false;
    std::string method;                 // "TRIAD", "QUEST", "ESOQ2"
    int        starsUsed  = 0;
};

/// Camera pointing for test validation
struct CameraPointing {
    Vector3    positionECEF = {};  // camera position [m] in ECEF
    Quaternion attitudeECI  = {1,0,0,0};  // ECI-to-camera quaternion
    double     epoch        = 0.0;  // Julian Date (TDB)
    StarTrackerFOV fov;
};

// ---------------------------------------------------------------------------
// Rendering types (for Cesium integration)
// ---------------------------------------------------------------------------

/// GPU-ready star vertex data (matches shader attributes)
struct StarVertex {
    float posX, posY, posZ;   // J2000 unit direction
    float magnitude;
    float colorIndex;
    float pmRA, pmDec;
    float parallax;
    float hipId;
};

static_assert(sizeof(StarVertex) == 36, "StarVertex must be 36 bytes");

}  // namespace starfield

#endif  // STARFIELD_TYPES_H
