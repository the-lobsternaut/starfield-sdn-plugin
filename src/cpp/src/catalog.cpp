#include "starfield/catalog.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace starfield {

// ---------------------------------------------------------------------------
// Star methods
// ---------------------------------------------------------------------------

Vector3 Star::direction() const {
    return {std::cos(dec) * std::cos(ra),
            std::cos(dec) * std::sin(ra),
            std::sin(dec)};
}

// ---------------------------------------------------------------------------
// Delta Encoding
// ---------------------------------------------------------------------------

void sortByRA(std::vector<Star>& stars) {
    std::sort(stars.begin(), stars.end(),
              [](const Star& a, const Star& b) { return a.ra < b.ra; });
}

std::vector<uint8_t> encodeCatalog(const std::vector<Star>& stars) {
    if (stars.empty()) return {};

    size_t headerSize = sizeof(CatalogHeader);
    size_t recordSize = sizeof(DeltaStarRecord);
    size_t totalSize = headerSize + stars.size() * recordSize;

    std::vector<uint8_t> data(totalSize);

    // Write header
    CatalogHeader header;
    header.count = static_cast<uint32_t>(stars.size());
    header.ra_ref = stars[0].ra;
    header.dec_ref = stars[0].dec;
    std::memcpy(data.data(), &header, headerSize);

    // Write delta-encoded records
    double prevRA = header.ra_ref;
    double prevDec = header.dec_ref;

    for (size_t i = 0; i < stars.size(); ++i) {
        const Star& s = stars[i];

        DeltaStarRecord rec;

        // Delta from previous star
        double deltaRA = (s.ra - prevRA) / DELTA_RA_UNIT;
        double deltaDec = (s.dec - prevDec) / DELTA_DEC_UNIT;

        // Clamp to int16 range
        rec.delta_ra = static_cast<int16_t>(
            std::max(-32767.0, std::min(32767.0, deltaRA)));
        rec.delta_dec = static_cast<int16_t>(
            std::max(-32767.0, std::min(32767.0, deltaDec)));

        rec.vmag = static_cast<int16_t>(s.vmag / VMAG_UNIT);
        rec.bv = static_cast<int16_t>(s.bv / BV_UNIT);
        rec.pmRA = static_cast<int16_t>(s.pmRA / PM_UNIT);
        rec.pmDec = static_cast<int16_t>(s.pmDec / PM_UNIT);
        rec.parallax = static_cast<uint16_t>(
            std::max(0.0f, s.parallax / static_cast<float>(PARALLAX_UNIT)));
        rec.hipId = s.hipId;

        std::memcpy(data.data() + headerSize + i * recordSize,
                     &rec, recordSize);

        prevRA = s.ra;
        prevDec = s.dec;
    }

    return data;
}

std::vector<Star> decodeCatalog(const uint8_t* data, size_t size) {
    size_t headerSize = sizeof(CatalogHeader);
    if (size < headerSize) {
        throw std::runtime_error("Invalid catalog: too small for header");
    }

    CatalogHeader header;
    std::memcpy(&header, data, headerSize);

    if (header.magic[0] != 'H' || header.magic[1] != 'I' ||
        header.magic[2] != 'P' || header.magic[3] != 'S') {
        throw std::runtime_error("Invalid catalog: bad magic");
    }

    size_t recordSize = sizeof(DeltaStarRecord);
    size_t expectedSize = headerSize + header.count * recordSize;
    if (size < expectedSize) {
        throw std::runtime_error("Invalid catalog: truncated data");
    }

    std::vector<Star> stars(header.count);
    double prevRA = header.ra_ref;
    double prevDec = header.dec_ref;

    for (uint32_t i = 0; i < header.count; ++i) {
        DeltaStarRecord rec;
        std::memcpy(&rec, data + headerSize + i * recordSize, recordSize);

        Star& s = stars[i];
        s.ra = prevRA + rec.delta_ra * DELTA_RA_UNIT;
        s.dec = prevDec + rec.delta_dec * DELTA_DEC_UNIT;
        s.vmag = rec.vmag * static_cast<float>(VMAG_UNIT);
        s.bv = rec.bv * static_cast<float>(BV_UNIT);
        s.pmRA = rec.pmRA * static_cast<float>(PM_UNIT);
        s.pmDec = rec.pmDec * static_cast<float>(PM_UNIT);
        s.parallax = rec.parallax * static_cast<float>(PARALLAX_UNIT);
        s.hipId = rec.hipId;

        prevRA = s.ra;
        prevDec = s.dec;
    }

    return stars;
}

std::vector<Star> decodeCatalog(const std::vector<uint8_t>& data) {
    return decodeCatalog(data.data(), data.size());
}

CompressionStats computeStats(const std::vector<Star>& stars) {
    CompressionStats stats;
    stats.rawSize = stars.size() * sizeof(Star);
    stats.encodedSize = sizeof(CatalogHeader) +
                        stars.size() * sizeof(DeltaStarRecord);
    stats.ratio = static_cast<double>(stats.rawSize) / stats.encodedSize;

    double prevRA = stars.empty() ? 0 : stars[0].ra;
    double prevDec = stars.empty() ? 0 : stars[0].dec;

    for (size_t i = 1; i < stars.size(); ++i) {
        double dra = std::abs(stars[i].ra - prevRA) / ARCSEC_TO_RAD;
        double ddec = std::abs(stars[i].dec - prevDec) / ARCSEC_TO_RAD;
        stats.maxDeltaRA_arcsec = std::max(stats.maxDeltaRA_arcsec, dra);
        stats.maxDeltaDec_arcsec = std::max(stats.maxDeltaDec_arcsec, ddec);

        double draUnits = std::abs(stars[i].ra - prevRA) / DELTA_RA_UNIT;
        double ddecUnits = std::abs(stars[i].dec - prevDec) / DELTA_DEC_UNIT;
        if (draUnits > 32767 || ddecUnits > 32767) stats.overflows++;

        prevRA = stars[i].ra;
        prevDec = stars[i].dec;
    }

    return stats;
}

// ---------------------------------------------------------------------------
// Catalog Queries
// ---------------------------------------------------------------------------

static double dot3(const Vector3& a, const Vector3& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static double norm3(const Vector3& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static Vector3 normalize3(const Vector3& v) {
    double n = norm3(v);
    if (n < 1e-15) return {0,0,0};
    return {v[0]/n, v[1]/n, v[2]/n};
}

static Vector3 cross3(const Vector3& a, const Vector3& b) {
    return {a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]};
}

std::vector<Star> queryCone(
    const std::vector<Star>& catalog,
    const Vector3& center, double halfAngle, double magLimit) {

    double cosHalf = std::cos(halfAngle);
    std::vector<Star> result;

    for (const auto& star : catalog) {
        if (star.vmag > magLimit) continue;
        Vector3 dir = star.direction();
        if (dot3(dir, center) >= cosHalf) {
            result.push_back(star);
        }
    }

    // Sort by magnitude (brightest first)
    std::sort(result.begin(), result.end(),
              [](const Star& a, const Star& b) { return a.vmag < b.vmag; });

    return result;
}

std::vector<Star> queryRectFOV(
    const std::vector<Star>& catalog,
    const Vector3& boresight, const Vector3& up,
    double fovX, double fovY, double magLimit) {

    // Build sensor frame: boresight = Z, up projected = Y, X = Y cross Z
    Vector3 z = normalize3(boresight);
    Vector3 y_raw = normalize3(up);
    Vector3 x = normalize3(cross3(y_raw, z));
    Vector3 y = cross3(z, x);

    double tanHalfX = std::tan(fovX / 2.0);
    double tanHalfY = std::tan(fovY / 2.0);

    std::vector<Star> result;

    for (const auto& star : catalog) {
        if (star.vmag > magLimit) continue;
        Vector3 dir = star.direction();

        // Project onto sensor frame
        double dz = dot3(dir, z);
        if (dz <= 0.0) continue;  // behind camera

        double dx = dot3(dir, x) / dz;
        double dy = dot3(dir, y) / dz;

        if (std::abs(dx) <= tanHalfX && std::abs(dy) <= tanHalfY) {
            result.push_back(star);
        }
    }

    std::sort(result.begin(), result.end(),
              [](const Star& a, const Star& b) { return a.vmag < b.vmag; });

    return result;
}

Vector3 applyProperMotion(const Star& star, double epochYears) {
    double dra = star.pmRA * MAS_TO_RAD * epochYears;
    double ddec = star.pmDec * MAS_TO_RAD * epochYears;

    double ra = star.ra + dra / std::cos(star.dec);
    double dec = star.dec + ddec;

    return {std::cos(dec) * std::cos(ra),
            std::cos(dec) * std::sin(ra),
            std::sin(dec)};
}

// ---------------------------------------------------------------------------
// GPU Buffer
// ---------------------------------------------------------------------------

std::vector<StarVertex> toVertexBuffer(
    const std::vector<Star>& stars, float magLimit) {

    std::vector<StarVertex> vertices;
    vertices.reserve(stars.size());

    for (const auto& s : stars) {
        if (s.vmag > magLimit) continue;
        Vector3 dir = s.direction();
        StarVertex v;
        v.posX = static_cast<float>(dir[0]);
        v.posY = static_cast<float>(dir[1]);
        v.posZ = static_cast<float>(dir[2]);
        v.magnitude = s.vmag;
        v.colorIndex = s.bv;
        v.pmRA = s.pmRA;
        v.pmDec = s.pmDec;
        v.parallax = s.parallax;
        v.hipId = static_cast<float>(s.hipId);
        vertices.push_back(v);
    }

    return vertices;
}

}  // namespace starfield
