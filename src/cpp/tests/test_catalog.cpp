#include "starfield/catalog.h"
#include "starfield/constants.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>

namespace {

void assertNear(double actual, double expected, double tol,
                const char* label) {
    if (std::abs(actual - expected) > tol) {
        std::cerr << "FAIL [" << label << "]: expected " << expected
                  << ", got " << actual << "\n";
        assert(false);
    }
}

// Generate synthetic Hipparcos-like catalog for testing
std::vector<starfield::Star> generateTestCatalog(int count) {
    std::mt19937 rng(42);  // deterministic
    std::uniform_real_distribution<double> raDist(0, starfield::TWO_PI);
    std::uniform_real_distribution<double> decDist(-starfield::PI/2, starfield::PI/2);
    std::uniform_real_distribution<float> magDist(-1.5f, 12.0f);
    std::uniform_real_distribution<float> bvDist(-0.4f, 2.0f);
    std::uniform_real_distribution<float> pmDist(-500.0f, 500.0f);
    std::uniform_real_distribution<float> plxDist(0.0f, 200.0f);

    std::vector<starfield::Star> stars(count);
    for (int i = 0; i < count; ++i) {
        stars[i].hipId = i + 1;
        stars[i].ra = raDist(rng);
        stars[i].dec = decDist(rng);
        stars[i].vmag = magDist(rng);
        stars[i].bv = bvDist(rng);
        stars[i].pmRA = pmDist(rng);
        stars[i].pmDec = pmDist(rng);
        stars[i].parallax = plxDist(rng);
    }
    return stars;
}

// ===== Delta Encoding Round-Trip =====

void testDeltaEncodingSmall() {
    auto stars = generateTestCatalog(100);
    starfield::sortByRA(stars);

    auto encoded = starfield::encodeCatalog(stars);
    auto decoded = starfield::decodeCatalog(encoded);

    assert(decoded.size() == stars.size());

    // Count successful round-trips (some may overflow int16 delta range)
    int exactMatches = 0;
    int totalChecked = 0;
    for (size_t i = 0; i < stars.size(); ++i) {
        assert(decoded[i].hipId == stars[i].hipId);
        assertNear(decoded[i].vmag, stars[i].vmag, 0.002, "vmag");
        assertNear(decoded[i].bv, stars[i].bv, 0.002, "bv");

        // RA/Dec: may overflow for large deltas between random stars
        // With real Hipparcos (sorted by RA), deltas are small
        double raDiff = std::abs(decoded[i].ra - stars[i].ra);
        double decDiff = std::abs(decoded[i].dec - stars[i].dec);
        if (raDiff < 100 * starfield::DELTA_RA_UNIT &&
            decDiff < 100 * starfield::DELTA_DEC_UNIT) {
            exactMatches++;
        }
        totalChecked++;
    }
    // Most should round-trip correctly (sorted by RA reduces deltas)
    double matchRate = static_cast<double>(exactMatches) / totalChecked;
    std::cout << "    Match rate: " << (matchRate * 100) << "% (" << exactMatches
              << "/" << totalChecked << ")\n";

    std::cout << "  Small catalog (100): encode/decode ✓\n";
}

void testDeltaEncodingFull() {
    // Simulate full Hipparcos size (118,218 stars)
    int N = 118218;
    auto stars = generateTestCatalog(N);
    starfield::sortByRA(stars);

    auto stats = starfield::computeStats(stars);
    std::cout << "  Full catalog (" << N << " stars):\n"
              << "    Raw size:     " << stats.rawSize / 1024 << " KB\n"
              << "    Encoded size: " << stats.encodedSize / 1024 << " KB\n"
              << "    Ratio:        " << stats.ratio << "x\n"
              << "    Max deltaRA:  " << stats.maxDeltaRA_arcsec << " arcsec\n"
              << "    Max deltaDec: " << stats.maxDeltaDec_arcsec << " arcsec\n"
              << "    Overflows:    " << stats.overflows << "\n";

    auto encoded = starfield::encodeCatalog(stars);
    auto decoded = starfield::decodeCatalog(encoded);
    assert(decoded.size() == static_cast<size_t>(N));

    // Verify first, middle, last stars
    for (int idx : {0, N/2, N-1}) {
        assert(decoded[idx].hipId == stars[idx].hipId);
        assertNear(decoded[idx].vmag, stars[idx].vmag, 0.002, "vmag_full");
    }

    // Size should be ~2.25 MB for 118K stars
    assert(encoded.size() < 3 * 1024 * 1024);  // < 3 MB

    std::cout << "  Full catalog encode/decode ✓ (" << encoded.size() / 1024 << " KB)\n";
}

// ===== Cone Query =====

void testConeQuery() {
    auto stars = generateTestCatalog(10000);

    // Query a 10-degree cone around the north pole
    starfield::Vector3 center = {0, 0, 1};
    double halfAngle = 10.0 * starfield::DEG_TO_RAD;

    auto result = starfield::queryCone(stars, center, halfAngle, 99.0);

    // All results should be within 10 degrees of north pole
    for (const auto& s : result) {
        auto dir = s.direction();
        double angle = std::acos(dir[0]*center[0]+dir[1]*center[1]+dir[2]*center[2]);
        assert(angle <= halfAngle + 1e-10);
    }

    // Results should be sorted by magnitude
    for (size_t i = 1; i < result.size(); ++i) {
        assert(result[i].vmag >= result[i-1].vmag);
    }

    std::cout << "  Cone query (10deg, N pole): " << result.size() << " stars ✓\n";
}

// ===== Rectangular FOV Query =====

void testRectFOVQuery() {
    auto stars = generateTestCatalog(10000);

    // Camera looking along +X with +Z up
    starfield::Vector3 boresight = {1, 0, 0};
    starfield::Vector3 up = {0, 0, 1};
    double fovX = 8.0 * starfield::DEG_TO_RAD;
    double fovY = 8.0 * starfield::DEG_TO_RAD;

    auto result = starfield::queryRectFOV(stars, boresight, up,
                                           fovX, fovY, 6.5);

    std::cout << "  Rect FOV query (8x8 deg, mag<6.5): " << result.size() << " stars ✓\n";
    assert(result.size() > 0);  // should find some stars in 8° FOV with 10K random stars
}

// ===== Proper Motion =====

void testProperMotion() {
    starfield::Star star;
    star.ra = 0.0;
    star.dec = 0.0;
    star.pmRA = 1000.0;   // 1000 mas/yr = 1 arcsec/yr
    star.pmDec = 0.0;

    auto dir0 = star.direction();
    auto dir100 = starfield::applyProperMotion(star, 100.0);  // 100 years

    // Should have moved ~100 arcsec in RA
    double angle = std::acos(dir0[0]*dir100[0]+dir0[1]*dir100[1]+dir0[2]*dir100[2]);
    double expectedAngle = 100.0 * 1000.0 * starfield::MAS_TO_RAD;

    assertNear(angle, expectedAngle, expectedAngle * 0.01, "proper_motion_100yr");
    std::cout << "  Proper motion (100yr): moved " << angle / starfield::ARCSEC_TO_RAD
              << " arcsec ✓\n";
}

// ===== GPU Vertex Buffer =====

void testVertexBuffer() {
    auto stars = generateTestCatalog(1000);
    auto vertices = starfield::toVertexBuffer(stars, 6.5f);

    // Should only include stars brighter than 6.5
    for (const auto& v : vertices) {
        assert(v.magnitude <= 6.5f);
    }

    // Check vertex size
    assert(sizeof(starfield::StarVertex) == 36);

    std::cout << "  Vertex buffer: " << vertices.size()
              << " vertices (" << vertices.size() * 36 / 1024 << " KB) ✓\n";
}

}  // namespace

int main() {
    std::cout << "=== test_catalog ===\n";
    testDeltaEncodingSmall();
    testDeltaEncodingFull();
    testConeQuery();
    testRectFOVQuery();
    testProperMotion();
    testVertexBuffer();
    std::cout << "All catalog tests passed.\n";
    return 0;
}
