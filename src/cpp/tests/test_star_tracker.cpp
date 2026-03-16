#include "starfield/star_tracker.h"
#include "starfield/catalog.h"
#include "starfield/constants.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

namespace {

using namespace starfield;

void assertNear(double actual, double expected, double tol,
                const char* label) {
    if (std::abs(actual - expected) > tol) {
        std::cerr << "FAIL [" << label << "]: expected " << expected
                  << ", got " << actual
                  << " (diff=" << std::abs(actual - expected) << ")\n";
        assert(false);
    }
}

// Create a small catalog of well-known bright stars for testing
std::vector<Star> createBrightStarCatalog() {
    std::vector<Star> stars;

    // Sirius (HIP 32349) — brightest star
    stars.push_back({32349, 101.287*DEG_TO_RAD, -16.716*DEG_TO_RAD,
                     -1.46f, 0.00f, -546.0f, -1223.0f, 379.2f, 0});

    // Canopus (HIP 30438)
    stars.push_back({30438, 95.988*DEG_TO_RAD, -52.696*DEG_TO_RAD,
                     -0.74f, 0.15f, 19.9f, 23.2f, 10.6f, 0});

    // Arcturus (HIP 69673)
    stars.push_back({69673, 213.915*DEG_TO_RAD, 19.182*DEG_TO_RAD,
                     -0.05f, 1.23f, -1093.4f, -1999.4f, 88.8f, 0});

    // Vega (HIP 91262)
    stars.push_back({91262, 279.235*DEG_TO_RAD, 38.784*DEG_TO_RAD,
                     0.03f, 0.00f, 200.9f, 286.2f, 130.2f, 0});

    // Capella (HIP 24608)
    stars.push_back({24608, 79.172*DEG_TO_RAD, 45.998*DEG_TO_RAD,
                     0.08f, 0.80f, 75.5f, -427.1f, 77.3f, 0});

    // Rigel (HIP 24436)
    stars.push_back({24436, 78.634*DEG_TO_RAD, -8.202*DEG_TO_RAD,
                     0.13f, -0.03f, 1.3f, -0.5f, 3.8f, 0});

    // Procyon (HIP 37279)
    stars.push_back({37279, 114.826*DEG_TO_RAD, 5.225*DEG_TO_RAD,
                     0.34f, 0.42f, -714.6f, -1036.8f, 284.6f, 0});

    // Betelgeuse (HIP 27989)
    stars.push_back({27989, 88.793*DEG_TO_RAD, 7.407*DEG_TO_RAD,
                     0.42f, 1.85f, 27.3f, 10.9f, 6.6f, 0});

    // Aldebaran (HIP 21421)
    stars.push_back({21421, 68.980*DEG_TO_RAD, 16.509*DEG_TO_RAD,
                     0.87f, 1.54f, 62.8f, -189.4f, 48.9f, 0});

    // Polaris (HIP 11767)
    stars.push_back({11767, 37.954*DEG_TO_RAD, 89.264*DEG_TO_RAD,
                     1.98f, 0.60f, 44.2f, -11.7f, 7.5f, 0});

    // Add more stars for denser coverage (fill with synthetic)
    std::mt19937 rng(123);
    std::uniform_real_distribution<double> raDist(0, TWO_PI);
    std::uniform_real_distribution<double> decDist(-PI/2, PI/2);
    std::uniform_real_distribution<float> magDist(1.0f, 6.5f);
    std::uniform_real_distribution<float> bvDist(-0.3f, 1.8f);

    for (int i = 0; i < 500; ++i) {
        Star s;
        s.hipId = 100000 + i;
        s.ra = raDist(rng);
        s.dec = decDist(rng);
        s.vmag = magDist(rng);
        s.bv = bvDist(rng);
        s.pmRA = 0.0f;
        s.pmDec = 0.0f;
        s.parallax = 10.0f;
        stars.push_back(s);
    }

    return stars;
}

// ===== GMST =====

void testGMST() {
    // J2000.0 epoch: GMST should be ~280.46 degrees
    double g = gmst(J2000_JD);
    double gDeg = g * RAD_TO_DEG;
    assertNear(gDeg, 280.46, 1.0, "gmst_j2000");
    std::cout << "  GMST at J2000.0: " << gDeg << " deg ✓\n";
}

// ===== ECEF ↔ ECI =====

void testECEFtoECI() {
    // At J2000.0 GMST~280 deg, so ECEF and ECI are rotated by ~280 deg
    Vector3 ecef = {R_EARTH, 0, 0};  // on equator at 0° lon
    Vector3 eci = ecefToECI(ecef, J2000_JD);

    // Should be rotated by GMST
    double g = gmst(J2000_JD);
    assertNear(eci[0], R_EARTH * std::cos(g), 1e3, "eci_x");
    assertNear(eci[1], R_EARTH * std::sin(g), 1e3, "eci_y");
    assertNear(eci[2], 0.0, 1.0, "eci_z");

    // Round-trip
    Vector3 ecef2 = eciToECEF(eci, J2000_JD);
    assertNear(ecef2[0], ecef[0], 1.0, "ecef_roundtrip_x");
    assertNear(ecef2[1], ecef[1], 1.0, "ecef_roundtrip_y");

    std::cout << "  ECEF↔ECI round-trip ✓\n";
}

// ===== Camera Pointing + Star Tracker Observation =====
// Test: ISS-like orbit, camera pointing at known star field

void testCameraPointingISS() {
    auto catalog = createBrightStarCatalog();

    // ISS position in ECEF (approx: over the equator at 0° lon, 400 km alt)
    CameraPointing pointing;
    pointing.positionECEF = {R_EARTH + 400e3, 0, 0};
    pointing.epoch = J2000_JD + 365.25 * 26;  // ~2026

    // Camera pointing toward Orion (RA~5.5h, Dec~0°)
    // In ECI: direction toward RA=82.5°, Dec=0°
    double targetRA = 82.5 * DEG_TO_RAD;
    double targetDec = 0.0;
    Vector3 boresightECI = raDecToECI(targetRA, targetDec);

    // Build quaternion: rotate body Z to point at boresight
    // Simple: use boresight as Z, north as Y hint
    Vector3 z = boresightECI;
    Vector3 yHint = {0, 0, 1};  // north
    Vector3 x = {z[1]*yHint[2]-z[2]*yHint[1],
                 z[2]*yHint[0]-z[0]*yHint[2],
                 z[0]*yHint[1]-z[1]*yHint[0]};
    double xn = std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    x = {x[0]/xn, x[1]/xn, x[2]/xn};
    Vector3 y = {z[1]*x[2]-z[2]*x[1], z[2]*x[0]-z[0]*x[2], z[0]*x[1]-z[1]*x[0]};

    // DCM columns are x, y, z (ECI-to-body)
    Matrix3x3 dcm = {{{x[0], y[0], z[0]},
                       {x[1], y[1], z[1]},
                       {x[2], y[2], z[2]}}};
    pointing.attitudeECI = dcmToQuaternion(dcm);

    pointing.fov.fovX = 20.0;  // wide FOV for more stars
    pointing.fov.fovY = 20.0;
    pointing.fov.magLimit = 6.5;
    pointing.fov.maxStars = 50;

    auto obs = simulateObservation(catalog, pointing);

    std::cout << "  ISS camera pointing toward Orion:\n"
              << "    Stars in FOV: " << obs.starsInFOV << "\n"
              << "    Stars detected: " << obs.starsDetected << "\n";

    assert(obs.starsDetected > 0);

    // All detected stars should have valid body vectors
    for (const auto& s : obs.stars) {
        double n = std::sqrt(s.bodyVector[0]*s.bodyVector[0] +
                            s.bodyVector[1]*s.bodyVector[1] +
                            s.bodyVector[2]*s.bodyVector[2]);
        assertNear(n, 1.0, 1e-10, "body_vec_unit");

        // Body vector Z should be positive (in front of camera)
        assert(s.bodyVector[2] > 0);

        // Sensor coords should be within [-1, 1]
        assert(std::abs(s.sensorX) <= 1.1);
        assert(std::abs(s.sensorY) <= 1.1);
    }

    std::cout << "  Star tracker observation ✓\n";
}

// ===== Water Ion Visibility in ECEF =====
// Test: simulate observing water ions (OH emission at 309nm) at known
// ECEF coordinates and verify they appear in Cesium world-space correctly

void testWaterIonECEF() {
    // Water ion source in ECEF: ~100km altitude, mesospheric
    // Location: above equator at 90°E longitude
    double ionAlt = 100e3;  // mesosphere
    double ionLat = 0.0;
    double ionLon = 90.0 * DEG_TO_RAD;
    Vector3 ionECEF = {(R_EARTH + ionAlt) * std::cos(ionLat) * std::cos(ionLon),
                        (R_EARTH + ionAlt) * std::cos(ionLat) * std::sin(ionLon),
                        (R_EARTH + ionAlt) * std::sin(ionLat)};

    // Camera in ECEF: ISS above same longitude at 400km
    Vector3 camECEF = {(R_EARTH + 400e3) * std::cos(ionLon),
                        (R_EARTH + 400e3) * std::sin(ionLon),
                        0};

    // Direction from camera to ion source
    Vector3 toIon = {ionECEF[0] - camECEF[0],
                     ionECEF[1] - camECEF[1],
                     ionECEF[2] - camECEF[2]};
    double dist = std::sqrt(toIon[0]*toIon[0]+toIon[1]*toIon[1]+toIon[2]*toIon[2]);
    Vector3 ionDir = {toIon[0]/dist, toIon[1]/dist, toIon[2]/dist};

    // Verify the ion source is below the camera (looking down at the atmosphere)
    double ionR = std::sqrt(ionECEF[0]*ionECEF[0]+ionECEF[1]*ionECEF[1]+ionECEF[2]*ionECEF[2]);
    double camR = std::sqrt(camECEF[0]*camECEF[0]+camECEF[1]*camECEF[1]+camECEF[2]*camECEF[2]);
    assert(ionR < camR);

    // Range should be ~300 km (400km - 100km roughly)
    assert(dist > 200e3 && dist < 500e3);

    // Convert ion ECEF to ECI at current epoch
    double jd = J2000_JD + 365.25 * 26;
    Vector3 ionECI = ecefToECI(ionECEF, jd);

    // Convert camera ECEF to ECI
    Vector3 camECI = ecefToECI(camECEF, jd);

    // Direction in ECI
    Vector3 toIonECI = {ionECI[0]-camECI[0], ionECI[1]-camECI[1], ionECI[2]-camECI[2]};
    double distECI = std::sqrt(toIonECI[0]*toIonECI[0]+toIonECI[1]*toIonECI[1]+toIonECI[2]*toIonECI[2]);

    // Distance should be preserved (rotation doesn't change distances)
    assertNear(distECI, dist, 1.0, "ecef_eci_distance");

    // Round-trip: ECI back to ECEF
    Vector3 ionECEF2 = eciToECEF(ionECI, jd);
    assertNear(ionECEF2[0], ionECEF[0], 1.0, "ion_ecef_roundtrip_x");
    assertNear(ionECEF2[1], ionECEF[1], 1.0, "ion_ecef_roundtrip_y");
    assertNear(ionECEF2[2], ionECEF[2], 1.0, "ion_ecef_roundtrip_z");

    // Get RA/Dec of ion direction from camera
    double ra, dec;
    Vector3 ionDirECI = {toIonECI[0]/distECI, toIonECI[1]/distECI, toIonECI[2]/distECI};
    eciToRADec(ionDirECI, ra, dec);

    std::cout << "  Water ion test:\n"
              << "    Ion ECEF: (" << ionECEF[0]/1e3 << ", " << ionECEF[1]/1e3
              << ", " << ionECEF[2]/1e3 << ") km\n"
              << "    Camera ECEF: (" << camECEF[0]/1e3 << ", " << camECEF[1]/1e3
              << ", " << camECEF[2]/1e3 << ") km\n"
              << "    Range: " << dist/1e3 << " km\n"
              << "    Ion direction (RA,Dec): (" << ra*RAD_TO_DEG << "°, "
              << dec*RAD_TO_DEG << "°)\n"
              << "    ECEF↔ECI round-trip ✓\n";
}

// ===== TRIAD Attitude Determination =====

void testTRIAD() {
    auto catalog = createBrightStarCatalog();

    // Create known attitude (identity = body aligned with ECI)
    CameraPointing pointing;
    pointing.positionECEF = {R_EARTH + 400e3, 0, 0};
    pointing.epoch = J2000_JD;
    pointing.attitudeECI = {1, 0, 0, 0};  // identity
    pointing.fov.fovX = 60.0;  // very wide to guarantee stars
    pointing.fov.fovY = 60.0;
    pointing.fov.magLimit = 6.5;

    auto obs = simulateObservation(catalog, pointing);
    if (obs.starsDetected < 2) {
        std::cout << "  TRIAD: not enough stars detected, skipping\n";
        return;
    }

    auto result = solveTriad(obs, catalog);

    std::cout << "  TRIAD: valid=" << result.valid
              << " stars=" << result.starsUsed
              << " error=" << result.error_arcsec << " arcsec\n";

    assert(result.valid);

    // The solved quaternion should be close to identity
    assertNear(std::abs(result.quaternion[0]), 1.0, 0.1, "triad_qw");
    std::cout << "  TRIAD quaternion: w=" << result.quaternion[0]
              << " x=" << result.quaternion[1]
              << " y=" << result.quaternion[2]
              << " z=" << result.quaternion[3] << " ✓\n";
}

// ===== QUEST Attitude Determination =====

void testQUEST() {
    auto catalog = createBrightStarCatalog();

    CameraPointing pointing;
    pointing.positionECEF = {R_EARTH + 400e3, 0, 0};
    pointing.epoch = J2000_JD;
    pointing.attitudeECI = {1, 0, 0, 0};
    pointing.fov.fovX = 60.0;
    pointing.fov.fovY = 60.0;
    pointing.fov.magLimit = 6.5;

    auto obs = simulateObservation(catalog, pointing);
    if (obs.starsDetected < 3) {
        std::cout << "  QUEST: not enough stars, skipping\n";
        return;
    }

    auto result = solveQUEST(obs, catalog);

    std::cout << "  QUEST: valid=" << result.valid
              << " stars=" << result.starsUsed
              << " error=" << result.error_arcsec << " arcsec\n";

    assert(result.valid);
    std::cout << "  QUEST quaternion: w=" << result.quaternion[0]
              << " x=" << result.quaternion[1]
              << " y=" << result.quaternion[2]
              << " z=" << result.quaternion[3] << " ✓\n";
}

// ===== FOV Cone in Cesium World-Space =====

void testFOVCone() {
    CameraPointing pointing;
    pointing.positionECEF = {R_EARTH + 400e3, 0, 0};
    pointing.epoch = J2000_JD;
    pointing.attitudeECI = {1, 0, 0, 0};
    pointing.fov.fovX = 8.0;
    pointing.fov.fovY = 8.0;

    auto cone = computeFOVCone(pointing);

    // Boresight should be a unit vector
    double bn = std::sqrt(cone.boresightECEF[0]*cone.boresightECEF[0] +
                          cone.boresightECEF[1]*cone.boresightECEF[1] +
                          cone.boresightECEF[2]*cone.boresightECEF[2]);
    assertNear(bn, 1.0, 1e-6, "fov_boresight_unit");

    // Corners should be unit vectors
    for (int i = 0; i < 4; ++i) {
        double cn = std::sqrt(cone.corners[i][0]*cone.corners[i][0] +
                              cone.corners[i][1]*cone.corners[i][1] +
                              cone.corners[i][2]*cone.corners[i][2]);
        assertNear(cn, 1.0, 1e-6, "fov_corner_unit");
    }

    // Half angle should match the larger FOV dimension
    assertNear(cone.halfAngle, 4.0 * DEG_TO_RAD, 1e-6, "fov_half_angle");

    std::cout << "  FOV cone in ECEF ✓\n";
}

// ===== World-Space Star Projection =====

void testWorldSpaceProjection() {
    auto catalog = createBrightStarCatalog();

    CameraPointing pointing;
    pointing.positionECEF = {R_EARTH + 400e3, 0, 0};
    pointing.epoch = J2000_JD;
    pointing.attitudeECI = {1, 0, 0, 0};
    pointing.fov.fovX = 60.0;
    pointing.fov.fovY = 60.0;
    pointing.fov.magLimit = 6.5;

    auto obs = simulateObservation(catalog, pointing);
    auto wsStars = projectToWorldSpace(obs, pointing);

    assert(wsStars.size() == obs.stars.size());

    for (const auto& ws : wsStars) {
        // ECEF position should be at large distance
        double dist = std::sqrt(ws.positionECEF[0]*ws.positionECEF[0] +
                                ws.positionECEF[1]*ws.positionECEF[1] +
                                ws.positionECEF[2]*ws.positionECEF[2]);
        assert(dist > 1e8);  // > 100,000 km
    }

    std::cout << "  World-space projection: " << wsStars.size() << " stars ✓\n";
}

}  // namespace

int main() {
    std::cout << "=== test_star_tracker ===\n";
    testGMST();
    testECEFtoECI();
    testCameraPointingISS();
    testWaterIonECEF();
    testTRIAD();
    testQUEST();
    testFOVCone();
    testWorldSpaceProjection();
    std::cout << "All star tracker tests passed.\n";
    return 0;
}
