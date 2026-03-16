#ifndef STARFIELD_STAR_TRACKER_H
#define STARFIELD_STAR_TRACKER_H

#include "types.h"
#include "catalog.h"

#include <vector>

namespace starfield {

// ---------------------------------------------------------------------------
// Coordinate Transforms
// ---------------------------------------------------------------------------

/// Convert ECEF position to geodetic (lat, lon, alt).
struct Geodetic { double lat, lon, alt; };  // rad, rad, m
Geodetic ecefToGeodetic(const Vector3& ecef);

/// Convert ECEF position + velocity to ECI at given Julian Date.
/// Uses simple sidereal rotation (adequate for star tracker).
Vector3 ecefToECI(const Vector3& ecef, double jd);

/// Compute Greenwich Mean Sidereal Time [rad] from Julian Date.
double gmst(double jd);

/// Rotate a vector from ECI to ECEF.
Vector3 eciToECEF(const Vector3& eci, double jd);

/// Convert ECI direction to RA/Dec [rad].
void eciToRADec(const Vector3& dir, double& ra, double& dec);

/// Convert RA/Dec [rad] to ECI unit vector.
Vector3 raDecToECI(double ra, double dec);

/// Quaternion to DCM.
Matrix3x3 quaternionToDCM(const Quaternion& q);

/// DCM to quaternion.
Quaternion dcmToQuaternion(const Matrix3x3& dcm);

/// Rotate vector by quaternion.
Vector3 rotateByQuaternion(const Vector3& v, const Quaternion& q);

/// Quaternion inverse.
Quaternion quaternionInverse(const Quaternion& q);

/// Quaternion multiply.
Quaternion quaternionMultiply(const Quaternion& a, const Quaternion& b);

// ---------------------------------------------------------------------------
// Star Tracker Simulation
// ---------------------------------------------------------------------------

/// Simulate a star tracker observation from a given camera pointing.
/// 1. Determines which catalog stars fall in the FOV
/// 2. Projects them to sensor coordinates
/// 3. Returns body-frame unit vectors for each detected star
StarTrackerObservation simulateObservation(
    const std::vector<Star>& catalog,
    const CameraPointing& pointing);

// ---------------------------------------------------------------------------
// Attitude Determination Algorithms
// ---------------------------------------------------------------------------

/// TRIAD method: determines attitude from exactly 2 vector observations.
/// Uses the two brightest stars. Fast, simple, O(1).
/// @param obs  observation with at least 2 stars
/// @param catalog  full catalog for reference vectors
AttitudeResult solveTriad(
    const StarTrackerObservation& obs,
    const std::vector<Star>& catalog);

/// QUEST method: optimal single-point attitude from N vector observations.
/// Minimizes Wahba's loss function via quaternion eigenvector.
/// @param obs  observation with at least 2 stars
/// @param catalog  full catalog for reference vectors
/// @param weights  optional per-star weights (default: magnitude-based)
AttitudeResult solveQUEST(
    const StarTrackerObservation& obs,
    const std::vector<Star>& catalog,
    const std::vector<double>& weights = {});

/// ESOQ2 method: fast optimal attitude from N vector observations.
/// Mortari's efficient formulation — O(N) after setup.
AttitudeResult solveESOQ2(
    const StarTrackerObservation& obs,
    const std::vector<Star>& catalog,
    const std::vector<double>& weights = {});

/// Lost-in-space: identify stars without prior attitude knowledge.
/// Uses triangle pattern matching against the catalog.
/// @param obs  observation (stars have unknown identity)
/// @param catalog  full catalog for matching
/// @param angularTolerance  matching tolerance [arcsec]
AttitudeResult solveLostInSpace(
    const StarTrackerObservation& obs,
    const std::vector<Star>& catalog,
    double angularTolerance = 10.0);

// ---------------------------------------------------------------------------
// Cesium World-Space Integration
// ---------------------------------------------------------------------------

/// Compute star tracker FOV cone in Cesium world coordinates (ECEF).
/// Returns 4 corner vectors + boresight for rendering the FOV overlay.
struct FOVCone {
    Vector3 boresightECEF;
    Vector3 corners[4];    // TL, TR, BR, BL in ECEF
    double  halfAngle;     // cone half-angle [rad]
};

FOVCone computeFOVCone(const CameraPointing& pointing);

/// Project catalog stars to Cesium world-space for overlay rendering.
/// Converts observed stars from body frame to ECEF for display.
struct WorldSpaceStar {
    uint32_t hipId;
    Vector3  positionECEF;    // direction * large_distance in ECEF
    double   magnitude;
    double   sensorX, sensorY;  // sensor frame coords
    bool     identified;       // matched to catalog
};

std::vector<WorldSpaceStar> projectToWorldSpace(
    const StarTrackerObservation& obs,
    const CameraPointing& pointing);

}  // namespace starfield

#endif  // STARFIELD_STAR_TRACKER_H
