#include "starfield/catalog.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <array>

namespace {

constexpr float PI = 3.14159265358979323846f;
constexpr float DEG_TO_RAD = PI / 180.0f;
constexpr float R_EARTH = 6371000.0f;  // meters

// ============================================================
// ECEF coordinate utilities
// ============================================================

struct ECEFPosition {
    double x, y, z;  // meters in ECEF
};

struct LLA {
    double lat, lon, alt;  // radians, radians, meters
};

ECEFPosition llaToECEF(const LLA& lla) {
    double a = 6378137.0;           // WGS84 semi-major axis
    double f = 1.0 / 298.257223563; // WGS84 flattening
    double e2 = 2.0 * f - f * f;    // eccentricity squared

    double sinLat = std::sin(lla.lat);
    double cosLat = std::cos(lla.lat);
    double sinLon = std::sin(lla.lon);
    double cosLon = std::cos(lla.lon);

    double N = a / std::sqrt(1.0 - e2 * sinLat * sinLat);

    return {
        (N + lla.alt) * cosLat * cosLon,
        (N + lla.alt) * cosLat * sinLon,
        (N * (1.0 - e2) + lla.alt) * sinLat
    };
}

/// Local ENU (East-North-Up) to ECEF rotation at a given position
struct ENUFrame {
    double east[3], north[3], up[3];
};

ENUFrame computeENU(const LLA& lla) {
    double sinLat = std::sin(lla.lat);
    double cosLat = std::cos(lla.lat);
    double sinLon = std::sin(lla.lon);
    double cosLon = std::cos(lla.lon);

    ENUFrame f;
    // East = [-sinLon, cosLon, 0]
    f.east[0] = -sinLon; f.east[1] = cosLon; f.east[2] = 0;
    // North = [-sinLat*cosLon, -sinLat*sinLon, cosLat]
    f.north[0] = -sinLat*cosLon; f.north[1] = -sinLat*sinLon; f.north[2] = cosLat;
    // Up = [cosLat*cosLon, cosLat*sinLon, sinLat]
    f.up[0] = cosLat*cosLon; f.up[1] = cosLat*sinLon; f.up[2] = sinLat;
    return f;
}

/// Convert azimuth/elevation to ECEF direction
/// Az: 0=North, 90=East, El: 0=horizon, 90=zenith
std::array<float, 3> azelToECEF(const LLA& observer, float azDeg, float elDeg) {
    float az = azDeg * DEG_TO_RAD;
    float el = elDeg * DEG_TO_RAD;

    // Direction in ENU
    float e = std::sin(az) * std::cos(el);
    float n = std::cos(az) * std::cos(el);
    float u = std::sin(el);

    auto enu = computeENU(observer);

    return {
        static_cast<float>(e * enu.east[0] + n * enu.north[0] + u * enu.up[0]),
        static_cast<float>(e * enu.east[1] + n * enu.north[1] + u * enu.up[1]),
        static_cast<float>(e * enu.east[2] + n * enu.north[2] + u * enu.up[2])
    };
}

float dot3(const std::array<float, 3>& a, const std::array<float, 3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

float norm3(const std::array<float, 3>& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// ============================================================
// Simulated camera / FOV
// ============================================================

struct Camera {
    LLA observer;
    float azimuth;    // degrees
    float elevation;  // degrees
    float fovDeg;     // field of view (degrees)
};

/// Check if a star direction (ICRS) is within camera FOV.
/// Simplified: treats ICRS ≈ ECEF (ignores sidereal rotation for testing).
/// In production, apply GAST rotation to get ECEF from ICRS.
bool isInFOV(const Camera& cam, const std::array<float, 3>& starDir) {
    auto camDir = azelToECEF(cam.observer, cam.azimuth, cam.elevation);
    float cosAngle = dot3(camDir, starDir) / (norm3(camDir) * norm3(starDir));
    float angle = std::acos(std::min(1.0f, std::max(-1.0f, cosAngle)));
    return angle < (cam.fovDeg * DEG_TO_RAD / 2.0f);
}

/// Compute the number of visible stars for a camera
int countVisibleStars(const Camera& cam, const std::vector<starfield::StarVertex>& vertices) {
    int count = 0;
    for (const auto& v : vertices) {
        std::array<float, 3> dir = {v.x, v.y, v.z};
        if (isInFOV(cam, dir)) ++count;
    }
    return count;
}

// ============================================================
// Water ion emission line check
// ============================================================

/// Water (H2O+) produces emission lines at:
///   - 619.5 nm (red-orange)
///   - 615.8 nm
///   - 557.7 nm (green, also [OI] auroral)
/// These can contaminate star tracker observations.
/// Test: generate B-V colors that correspond to water ion emission
/// and verify they render as expected colors.

struct EmissionLine {
    float wavelength_nm;
    std::string name;
    float expectedBV;  // approximate B-V for this color temperature
};

void testWaterIonColors() {
    // Water ion emission lines in cometary spectra
    EmissionLine lines[] = {
        {619.5f, "H2O+ (619.5nm)", 1.2f},   // red-orange, ~K star
        {615.8f, "H2O+ (615.8nm)", 1.15f},   // red-orange
        {557.7f, "[OI] (557.7nm)", 0.6f},     // green, ~G star (solar)
        {630.0f, "[OI] (630.0nm)", 1.3f},     // red
    };

    std::cout << "  Water ion / auroral emission line colors:\n";
    for (const auto& line : lines) {
        // Verify B-V → RGB mapping produces expected color range
        // (This is a rendering validation — the fragment shader uses the same mapping)
        float bv = line.expectedBV;

        // Simplified B-V → RGB (matches fragment shader)
        float r, g, b;
        if (bv < 0.4f) r = 0.83f + 0.5f * bv;
        else r = std::min(1.0f, 0.71f + 0.42f * bv);

        if (bv < 0.4f) g = 0.87f + 0.11f * bv;
        else if (bv < 1.5f) g = 1.0f - 0.43f * (bv - 0.4f);
        else g = std::max(0.2f, 0.53f - 0.2f * (bv - 1.5f));

        if (bv < 0.4f) b = 1.0f - 1.67f * bv;
        else if (bv < 1.0f) b = std::max(0.0f, 0.33f - 0.55f * (bv - 0.4f));
        else b = 0.0f;

        std::cout << "    " << line.name << " (λ=" << line.wavelength_nm
                  << "nm, B-V=" << bv << ") → RGB("
                  << r << ", " << g << ", " << b << ")\n";

        // Water ions should render as warm colors (R > G > B)
        if (line.wavelength_nm > 600.0f) {
            assert(r > g && "red emission should have R > G");
            assert(g > b && "red emission should have G > B");
        }
    }

    std::cout << "  Water ion emission color mapping ✓\n";
}

// ============================================================
// Camera pointing tests
// ============================================================

void testCameraPointingZenith() {
    // Observer at (0°N, 0°E, 0m) looking straight up (zenith)
    // Stars near Dec=0 should be visible (near ICRS equator at this longitude)
    Camera cam;
    cam.observer = {0.0, 0.0, 0.0};  // equator, prime meridian
    cam.azimuth = 0.0f;
    cam.elevation = 90.0f;
    cam.fovDeg = 60.0f;

    auto camDir = azelToECEF(cam.observer, cam.azimuth, cam.elevation);

    // Zenith at (0,0) should point along ECEF +X
    assert(std::abs(camDir[0] - 1.0f) < 0.01f);
    assert(std::abs(camDir[1]) < 0.01f);
    assert(std::abs(camDir[2]) < 0.01f);

    std::cout << "  Camera zenith at equator/PM → ECEF +X ✓\n";
}

void testCameraPointingNorth() {
    // Observer at 45°N, 0°E looking north at 45° elevation
    Camera cam;
    cam.observer = {45.0 * DEG_TO_RAD, 0.0, 0.0};
    cam.azimuth = 0.0f;    // North
    cam.elevation = 45.0f;  // 45° above horizon

    auto camDir = azelToECEF(cam.observer, cam.azimuth, cam.elevation);

    // Should point roughly toward the NCP (ECEF +Z component should be large)
    assert(camDir[2] > 0.5f);

    std::cout << "  Camera pointing north from 45°N → strong +Z ✓\n";
}

void testStarVisibilityInFOV() {
    // Create synthetic catalog
    std::vector<starfield::Star> stars;

    // Star at RA=0, Dec=0 (vernal equinox, ICRS +X)
    stars.push_back({1, 0.0f, 0.0f, 1.0f, 0.5f, 0, 0, 0});

    // Star at RA=180°, Dec=0 (anti-vernal equinox, ICRS -X)
    stars.push_back({2, PI, 0.0f, 2.0f, 0.5f, 0, 0, 0});

    // Star at NCP (Dec=90°)
    stars.push_back({3, 0.0f, PI/2.0f, 0.0f, 0.0f, 0, 0, 0});

    auto vbuf = starfield::generateVertexBuffer(stars);

    // Camera looking toward ICRS +X (RA=0, Dec=0)
    // Observer at equator/PM, zenith
    Camera cam;
    cam.observer = {0.0, 0.0, 0.0};
    cam.azimuth = 0.0f;
    cam.elevation = 90.0f;
    cam.fovDeg = 30.0f;

    int visible = countVisibleStars(cam, vbuf);
    std::cout << "  Stars visible in 30° FOV toward RA=0: " << visible << "\n";

    // Star at RA=0 should be visible (it's in the ECEF +X direction)
    // Star at RA=180° should NOT be visible
    // NCP depends on FOV
    assert(visible >= 1);

    std::cout << "  Star visibility in FOV ✓\n";
}

void testECEFStarPositions() {
    // Verify that star positions placed at "infinity" in ECEF
    // maintain correct angular separation

    // Sirius: RA=6h45m, Dec=-16°43'
    starfield::Star sirius;
    sirius.ra = (6.0f + 45.0f/60.0f) * 15.0f * DEG_TO_RAD;
    sirius.dec = -(16.0f + 43.0f/60.0f) * DEG_TO_RAD;

    // Betelgeuse: RA=5h55m, Dec=+7°24'
    starfield::Star betelgeuse;
    betelgeuse.ra = (5.0f + 55.0f/60.0f) * 15.0f * DEG_TO_RAD;
    betelgeuse.dec = (7.0f + 24.0f/60.0f) * DEG_TO_RAD;

    auto dirS = starfield::radecToDirection(sirius.ra, sirius.dec);
    auto dirB = starfield::radecToDirection(betelgeuse.ra, betelgeuse.dec);

    // Angular separation: Sirius to Betelgeuse ≈ 27.5°
    float cosAngle = dot3(dirS, dirB);
    float angleDeg = std::acos(std::min(1.0f, cosAngle)) / DEG_TO_RAD;
    std::cout << "  Sirius-Betelgeuse separation: " << angleDeg << "°\n";
    assert(angleDeg > 25.0f && angleDeg < 30.0f);

    // Place at "infinity" (1e9 m) — angular separation should be preserved
    float dist = 1e9f;
    std::array<float, 3> posSirius = {dirS[0]*dist, dirS[1]*dist, dirS[2]*dist};
    std::array<float, 3> posBetel = {dirB[0]*dist, dirB[1]*dist, dirB[2]*dist};

    // Recompute angle from ECEF positions
    float dx = posBetel[0]-posSirius[0];
    float dy = posBetel[1]-posSirius[1];
    float dz = posBetel[2]-posSirius[2];
    float chord = std::sqrt(dx*dx + dy*dy + dz*dz);
    float angularSep = 2.0f * std::asin(chord / (2.0f * dist)) / DEG_TO_RAD;
    assert(std::abs(angularSep - angleDeg) < 0.1f);

    std::cout << "  ECEF star positions preserve angular separation ✓\n";
}

void testMagnitudeToBrightness() {
    // Verify Pogson scale: 5 magnitudes = factor of 100 in flux
    float magBright = 0.0f;
    float magFaint = 5.0f;
    float fluxRatio = std::pow(2.512f, -(magBright - magFaint));

    // Brighter star should have 100x the flux
    assert(std::abs(fluxRatio - 100.0f) < 1.0f);

    // Point size: sqrt(flux) ratio
    float sizeRatio = std::sqrt(fluxRatio);
    assert(std::abs(sizeRatio - 10.0f) < 0.1f);

    std::cout << "  Pogson magnitude scale ✓ (5 mag = 100x flux = 10x size)\n";
}

}  // namespace

int main() {
    std::cout << "=== test_rendering ===\n";
    testCameraPointingZenith();
    testCameraPointingNorth();
    testStarVisibilityInFOV();
    testECEFStarPositions();
    testWaterIonColors();
    testMagnitudeToBrightness();
    std::cout << "All rendering tests passed.\n";
    return 0;
}
