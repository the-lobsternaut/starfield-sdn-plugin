/**
 * starfield-sdn-plugin WASM API
 *
 * Star catalog, FOV queries, attitude determination, GPU vertex buffers.
 * Binary catalog decoding + JSON-in/JSON-out for queries and attitude.
 */

#include "starfield/types.h"
#include "starfield/catalog.h"
#include "starfield/star_tracker.h"
#include "starfield/constants.h"

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#include <emscripten/bind.h>
using namespace emscripten;
#endif

#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>

// ============================================================================
// Version
// ============================================================================

static std::string version() {
    return "0.1.0";
}

// ============================================================================
// Catalog state (held in WASM memory for query reuse)
// ============================================================================

static std::vector<starfield::Star> g_catalog;

/// Load delta-encoded catalog from binary data
static int load_catalog(uintptr_t data_ptr, size_t size) {
    const uint8_t* data = reinterpret_cast<const uint8_t*>(data_ptr);
    g_catalog = starfield::decodeCatalog(data, size);
    return static_cast<int>(g_catalog.size());
}

/// Get loaded catalog size
static int catalog_size() {
    return static_cast<int>(g_catalog.size());
}

// ============================================================================
// Catalog queries (JSON output)
// ============================================================================

/// Query stars in a cone — returns JSON array
static std::string query_cone_json(double cx, double cy, double cz,
                                    double half_angle_deg, double mag_limit) {
    starfield::Vector3 center = {cx, cy, cz};
    double half_angle_rad = half_angle_deg * starfield::DEG_TO_RAD;
    auto stars = starfield::queryCone(g_catalog, center, half_angle_rad, mag_limit);

    std::ostringstream oss;
    oss.precision(10);
    oss << "{\"count\":" << stars.size() << ",\"stars\":[";
    for (size_t i = 0; i < stars.size(); ++i) {
        if (i > 0) oss << ",";
        auto dir = stars[i].direction();
        oss << "{\"hipId\":" << stars[i].hipId
            << ",\"ra\":" << stars[i].ra * starfield::RAD_TO_DEG
            << ",\"dec\":" << stars[i].dec * starfield::RAD_TO_DEG
            << ",\"vmag\":" << stars[i].vmag
            << ",\"bv\":" << stars[i].bv
            << ",\"dir\":[" << dir[0] << "," << dir[1] << "," << dir[2] << "]}";
    }
    oss << "]}";
    return oss.str();
}

/// Get compression stats for the loaded catalog
static std::string compression_stats_json() {
    auto stats = starfield::computeStats(g_catalog);
    std::ostringstream oss;
    oss << "{\"rawSize\":" << stats.rawSize
        << ",\"encodedSize\":" << stats.encodedSize
        << ",\"ratio\":" << stats.ratio
        << ",\"maxDeltaRA_arcsec\":" << stats.maxDeltaRA_arcsec
        << ",\"maxDeltaDec_arcsec\":" << stats.maxDeltaDec_arcsec
        << ",\"overflows\":" << stats.overflows << "}";
    return oss.str();
}

// ============================================================================
// GPU vertex buffer generation (binary output)
// ============================================================================

/// Generate GPU-ready vertex buffer, returns byte count written to output_ptr
static int to_vertex_buffer(uintptr_t output_ptr, size_t output_capacity,
                             float mag_limit) {
    auto verts = starfield::toVertexBuffer(g_catalog, mag_limit);
    size_t bytes = verts.size() * sizeof(starfield::StarVertex);
    if (bytes > output_capacity) return -2;
    std::memcpy(reinterpret_cast<void*>(output_ptr), verts.data(), bytes);
    return static_cast<int>(bytes);
}

// ============================================================================
// Star tracker / attitude determination (JSON output)
// ============================================================================

/// Simulate star tracker observation → JSON
static std::string simulate_observation_json(
    double posX, double posY, double posZ,       // ECEF position [m]
    double qw, double qx, double qy, double qz, // ECI-to-camera quat
    double epoch_jd,
    double fovX_deg, double fovY_deg, double mag_limit) {

    starfield::CameraPointing cp;
    cp.positionECEF = {posX, posY, posZ};
    cp.attitudeECI = {qw, qx, qy, qz};
    cp.epoch = epoch_jd;
    cp.fov.fovX = fovX_deg;
    cp.fov.fovY = fovY_deg;
    cp.fov.magLimit = mag_limit;

    auto obs = starfield::simulateObservation(g_catalog, cp);

    std::ostringstream oss;
    oss.precision(10);
    oss << "{\"starsInFOV\":" << obs.starsInFOV
        << ",\"starsDetected\":" << obs.starsDetected
        << ",\"epoch\":" << obs.epoch
        << ",\"stars\":[";
    for (size_t i = 0; i < obs.stars.size(); ++i) {
        if (i > 0) oss << ",";
        oss << "{\"hipId\":" << obs.stars[i].hipId
            << ",\"sensorX\":" << obs.stars[i].sensorX
            << ",\"sensorY\":" << obs.stars[i].sensorY
            << ",\"magnitude\":" << obs.stars[i].magnitude << "}";
    }
    oss << "]}";
    return oss.str();
}

// ============================================================================
// Coordinate transforms
// ============================================================================

static std::string ecef_to_geodetic_json(double x, double y, double z) {
    starfield::Vector3 ecef = {x, y, z};
    auto geo = starfield::ecefToGeodetic(ecef);
    std::ostringstream oss;
    oss.precision(12);
    oss << "{\"lat\":" << geo.lat * starfield::RAD_TO_DEG
        << ",\"lon\":" << geo.lon * starfield::RAD_TO_DEG
        << ",\"alt\":" << geo.alt << "}";
    return oss.str();
}

static double compute_gmst(double jd) {
    return starfield::gmst(jd);
}

// ============================================================================
// SDN Plugin ABI (C exports)
// ============================================================================

#ifdef __EMSCRIPTEN__

extern "C" {

EMSCRIPTEN_KEEPALIVE
void* sdn_malloc(size_t size) {
    return malloc(size);
}

EMSCRIPTEN_KEEPALIVE
void sdn_free(void* ptr) {
    free(ptr);
}

} // extern "C"

// ============================================================================
// Embind Exports
// ============================================================================

EMSCRIPTEN_BINDINGS(sdn_starfield) {
    function("version", &version);

    // Catalog management
    function("loadCatalog", &load_catalog, allow_raw_pointers());
    function("catalogSize", &catalog_size);
    function("compressionStats", &compression_stats_json);

    // Queries
    function("queryCone", &query_cone_json);
    function("toVertexBuffer", &to_vertex_buffer, allow_raw_pointers());

    // Star tracker
    function("simulateObservation", &simulate_observation_json);

    // Coordinate transforms
    function("ecefToGeodetic", &ecef_to_geodetic_json);
    function("computeGMST", &compute_gmst);
}

#endif
