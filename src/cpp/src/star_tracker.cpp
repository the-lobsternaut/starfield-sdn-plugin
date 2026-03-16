#include "starfield/star_tracker.h"
#include "starfield/constants.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace starfield {

// ---------------------------------------------------------------------------
// Vector/Quaternion utilities
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

static Vector3 sub3(const Vector3& a, const Vector3& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

static Vector3 scale3(const Vector3& v, double s) {
    return {v[0]*s, v[1]*s, v[2]*s};
}

static Vector3 add3(const Vector3& a, const Vector3& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

// ---------------------------------------------------------------------------
// Coordinate Transforms
// ---------------------------------------------------------------------------

double gmst(double jd) {
    double T = (jd - J2000_JD) / 36525.0;
    // IAU 1982 GMST
    double gmst_sec = 67310.54841 + (876600.0*3600.0 + 8640184.812866)*T +
                      0.093104*T*T - 6.2e-6*T*T*T;
    double gmst_rad = std::fmod(gmst_sec * PI / 43200.0, TWO_PI);
    if (gmst_rad < 0) gmst_rad += TWO_PI;
    return gmst_rad;
}

Geodetic ecefToGeodetic(const Vector3& ecef) {
    double x = ecef[0], y = ecef[1], z = ecef[2];
    double lon = std::atan2(y, x);
    double p = std::sqrt(x*x + y*y);
    double lat = std::atan2(z, p);  // simplified (no WGS84 iteration)
    double alt = std::sqrt(x*x + y*y + z*z) - R_EARTH;
    return {lat, lon, alt};
}

Vector3 ecefToECI(const Vector3& ecef, double jd) {
    double theta = gmst(jd);
    double ct = std::cos(theta), st = std::sin(theta);
    return {ct*ecef[0] - st*ecef[1],
            st*ecef[0] + ct*ecef[1],
            ecef[2]};
}

Vector3 eciToECEF(const Vector3& eci, double jd) {
    double theta = gmst(jd);
    double ct = std::cos(theta), st = std::sin(theta);
    return { ct*eci[0] + st*eci[1],
            -st*eci[0] + ct*eci[1],
             eci[2]};
}

void eciToRADec(const Vector3& dir, double& ra, double& dec) {
    Vector3 d = normalize3(dir);
    dec = std::asin(d[2]);
    ra = std::atan2(d[1], d[0]);
    if (ra < 0) ra += TWO_PI;
}

Vector3 raDecToECI(double ra, double dec) {
    return {std::cos(dec)*std::cos(ra),
            std::cos(dec)*std::sin(ra),
            std::sin(dec)};
}

Matrix3x3 quaternionToDCM(const Quaternion& q) {
    double w=q[0], x=q[1], y=q[2], z=q[3];
    Matrix3x3 m;
    m[0] = {1-2*(y*y+z*z), 2*(x*y-w*z),   2*(x*z+w*y)};
    m[1] = {2*(x*y+w*z),   1-2*(x*x+z*z), 2*(y*z-w*x)};
    m[2] = {2*(x*z-w*y),   2*(y*z+w*x),   1-2*(x*x+y*y)};
    return m;
}

Quaternion dcmToQuaternion(const Matrix3x3& m) {
    double tr = m[0][0] + m[1][1] + m[2][2];
    Quaternion q;
    if (tr > 0) {
        double s = 0.5 / std::sqrt(tr + 1.0);
        q = {0.25/s, (m[2][1]-m[1][2])*s, (m[0][2]-m[2][0])*s, (m[1][0]-m[0][1])*s};
    } else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
        double s = 2.0 * std::sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]);
        q = {(m[2][1]-m[1][2])/s, 0.25*s, (m[0][1]+m[1][0])/s, (m[0][2]+m[2][0])/s};
    } else if (m[1][1] > m[2][2]) {
        double s = 2.0 * std::sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]);
        q = {(m[0][2]-m[2][0])/s, (m[0][1]+m[1][0])/s, 0.25*s, (m[1][2]+m[2][1])/s};
    } else {
        double s = 2.0 * std::sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]);
        q = {(m[1][0]-m[0][1])/s, (m[0][2]+m[2][0])/s, (m[1][2]+m[2][1])/s, 0.25*s};
    }
    // Normalize
    double n = std::sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
    return {q[0]/n, q[1]/n, q[2]/n, q[3]/n};
}

Vector3 rotateByQuaternion(const Vector3& v, const Quaternion& q) {
    Matrix3x3 dcm = quaternionToDCM(q);
    return {dcm[0][0]*v[0]+dcm[0][1]*v[1]+dcm[0][2]*v[2],
            dcm[1][0]*v[0]+dcm[1][1]*v[1]+dcm[1][2]*v[2],
            dcm[2][0]*v[0]+dcm[2][1]*v[1]+dcm[2][2]*v[2]};
}

Quaternion quaternionInverse(const Quaternion& q) {
    return {q[0], -q[1], -q[2], -q[3]};
}

Quaternion quaternionMultiply(const Quaternion& a, const Quaternion& b) {
    return {a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
            a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
            a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
            a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0]};
}

// ---------------------------------------------------------------------------
// Star Tracker Simulation
// ---------------------------------------------------------------------------

StarTrackerObservation simulateObservation(
    const std::vector<Star>& catalog,
    const CameraPointing& pointing) {

    // Convert attitude to DCM: ECI-to-body
    Matrix3x3 dcm_eci2body = quaternionToDCM(pointing.attitudeECI);

    // Boresight in ECI (body Z-axis → ECI)
    Quaternion qInv = quaternionInverse(pointing.attitudeECI);
    Vector3 boresightECI = rotateByQuaternion({0,0,1}, qInv);
    Vector3 upECI = rotateByQuaternion({0,1,0}, qInv);

    double fovXRad = pointing.fov.fovX * DEG_TO_RAD;
    double fovYRad = pointing.fov.fovY * DEG_TO_RAD;
    double magLimit = pointing.fov.magLimit;

    // Find stars in rectangular FOV
    auto starsInFOV = queryRectFOV(catalog, boresightECI, upECI,
                                    fovXRad, fovYRad, magLimit);

    StarTrackerObservation obs;
    obs.starsInFOV = static_cast<int>(starsInFOV.size());
    obs.epoch = pointing.epoch;

    int maxStars = pointing.fov.maxStars;

    for (size_t i = 0; i < starsInFOV.size() &&
         static_cast<int>(obs.stars.size()) < maxStars; ++i) {

        const Star& star = starsInFOV[i];
        Vector3 dirECI = star.direction();

        // Transform to body frame
        Vector3 bodyVec = rotateByQuaternion(dirECI, pointing.attitudeECI);

        // Project to sensor coordinates (normalized)
        if (bodyVec[2] <= 0) continue;  // behind camera
        double sx = bodyVec[0] / bodyVec[2];
        double sy = bodyVec[1] / bodyVec[2];

        ObservedStar oStar;
        oStar.hipId = star.hipId;
        oStar.sensorX = sx / std::tan(fovXRad / 2.0);
        oStar.sensorY = sy / std::tan(fovYRad / 2.0);
        oStar.magnitude = star.vmag;
        oStar.bodyVector = normalize3(bodyVec);
        obs.stars.push_back(oStar);
    }

    obs.starsDetected = static_cast<int>(obs.stars.size());
    return obs;
}

// ---------------------------------------------------------------------------
// TRIAD Attitude Determination
// ---------------------------------------------------------------------------

AttitudeResult solveTriad(
    const StarTrackerObservation& obs,
    const std::vector<Star>& catalog) {

    AttitudeResult result;
    result.method = "TRIAD";

    if (obs.stars.size() < 2) {
        result.valid = false;
        return result;
    }

    // Body frame vectors (measured)
    Vector3 b1 = obs.stars[0].bodyVector;
    Vector3 b2 = obs.stars[1].bodyVector;

    // Reference frame vectors (from catalog, ECI)
    Vector3 r1 = {0,0,0}, r2 = {0,0,0};
    for (const auto& s : catalog) {
        if (s.hipId == obs.stars[0].hipId) r1 = s.direction();
        if (s.hipId == obs.stars[1].hipId) r2 = s.direction();
    }

    // TRIAD: build two triads and compute DCM
    // Body triad
    Vector3 s1b = b1;
    Vector3 s2b = normalize3(cross3(b1, b2));
    Vector3 s3b = cross3(s1b, s2b);

    // Reference triad
    Vector3 s1r = r1;
    Vector3 s2r = normalize3(cross3(r1, r2));
    Vector3 s3r = cross3(s1r, s2r);

    // DCM = [s1b s2b s3b] * [s1r s2r s3r]^T
    // A_ref2body = M_body * M_ref^T
    Matrix3x3 Mb = {{{s1b[0], s2b[0], s3b[0]},
                      {s1b[1], s2b[1], s3b[1]},
                      {s1b[2], s2b[2], s3b[2]}}};
    Matrix3x3 Mr = {{{s1r[0], s2r[0], s3r[0]},
                      {s1r[1], s2r[1], s3r[1]},
                      {s1r[2], s2r[2], s3r[2]}}};

    // DCM = Mb * Mr^T
    Matrix3x3 dcm;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            dcm[i][j] = 0;
            for (int k = 0; k < 3; ++k)
                dcm[i][j] += Mb[i][k] * Mr[j][k];  // Mr transposed
        }

    result.dcm = dcm;
    result.quaternion = dcmToQuaternion(dcm);
    result.starsUsed = 2;
    result.valid = true;
    result.error_arcsec = 30.0;  // typical TRIAD error estimate
    return result;
}

// ---------------------------------------------------------------------------
// QUEST Attitude Determination
// ---------------------------------------------------------------------------

AttitudeResult solveQUEST(
    const StarTrackerObservation& obs,
    const std::vector<Star>& catalog,
    const std::vector<double>& weights) {

    AttitudeResult result;
    result.method = "QUEST";

    int N = static_cast<int>(obs.stars.size());
    if (N < 2) { result.valid = false; return result; }

    // Build observation and reference vector pairs
    std::vector<Vector3> bodyVecs, refVecs;
    std::vector<double> w;

    for (int i = 0; i < N; ++i) {
        bodyVecs.push_back(obs.stars[i].bodyVector);

        // Find reference vector
        Vector3 ref = {0,0,0};
        for (const auto& s : catalog) {
            if (s.hipId == obs.stars[i].hipId) {
                ref = s.direction();
                break;
            }
        }
        refVecs.push_back(ref);

        // Weight: user-provided or magnitude-based
        if (i < static_cast<int>(weights.size())) {
            w.push_back(weights[i]);
        } else {
            // Brighter stars get higher weight
            w.push_back(1.0 / (1.0 + obs.stars[i].magnitude));
        }
    }

    // Normalize weights
    double wSum = 0;
    for (auto wi : w) wSum += wi;
    for (auto& wi : w) wi /= wSum;

    // QUEST: solve Wahba's problem via characteristic equation
    // B = sum(wi * bi * ri^T)
    Matrix3x3 B = {};
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                B[j][k] += w[i] * bodyVecs[i][j] * refVecs[i][k];
    }

    // S = B + B^T
    Matrix3x3 S;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            S[i][j] = B[i][j] + B[j][i];

    // sigma = trace(B)
    double sigma = B[0][0] + B[1][1] + B[2][2];

    // Z = [B23-B32, B31-B13, B12-B21]
    Vector3 Z = {B[1][2]-B[2][1], B[2][0]-B[0][2], B[0][1]-B[1][0]};

    // Solve characteristic equation via Newton-Raphson
    // f(lambda) = det(lambda*I - [S-sigma*I, Z; Z^T, sigma]) = 0
    // Start with lambda = 1 (maximum eigenvalue)
    double lambda = 1.0;
    for (int iter = 0; iter < 50; ++iter) {
        // (lambda + sigma) * [(lambda - sigma)^2 - ||Z||^2] - ... 
        // Simplified: use the standard QUEST formulation
        double alpha = lambda*lambda - sigma*sigma + dot3(Z,Z);
        double beta = lambda - sigma;

        // K matrix: [[S - sigma*I, Z], [Z^T, sigma]]
        Matrix3x3 A;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                A[i][j] = S[i][j] - (i==j ? lambda : 0.0);

        // Determinant of A (3x3)
        double detA = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) -
                      A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]) +
                      A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

        // f(lambda) = detA * (lambda - sigma) - Z^T * adj(A) * Z
        // Approximate derivative
        double h = 1e-8;
        Matrix3x3 Ah;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                Ah[i][j] = S[i][j] - (i==j ? (lambda+h) : 0.0);
        double detAh = Ah[0][0]*(Ah[1][1]*Ah[2][2]-Ah[1][2]*Ah[2][1]) -
                       Ah[0][1]*(Ah[1][0]*Ah[2][2]-Ah[1][2]*Ah[2][0]) +
                       Ah[0][2]*(Ah[1][0]*Ah[2][1]-Ah[1][1]*Ah[2][0]);

        double f = detA;
        double df = (detAh - detA) / h;

        if (std::abs(df) < 1e-30) break;
        double dlambda = -f / df;
        lambda += dlambda;
        if (std::abs(dlambda) < 1e-12) break;
    }

    // Compute optimal quaternion from lambda
    // q = [(lambda + sigma)*Z + adj(S-sigma*I)*Z; det(S-sigma*I)]
    // (simplified Rodrigues parameter → quaternion)
    double gamma = lambda + sigma;
    if (std::abs(gamma) < 1e-15) {
        result.valid = false;
        return result;
    }

    Vector3 x = scale3(Z, 1.0 / gamma);
    double qw = 1.0;
    double qNorm = std::sqrt(1.0 + dot3(x,x));
    result.quaternion = {qw/qNorm, x[0]/qNorm, x[1]/qNorm, x[2]/qNorm};
    result.dcm = quaternionToDCM(result.quaternion);
    result.starsUsed = N;
    result.valid = true;
    result.error_arcsec = 10.0 / std::sqrt(N);  // improves with more stars
    return result;
}

// ---------------------------------------------------------------------------
// ESOQ2 (simplified — delegates to QUEST for now)
// ---------------------------------------------------------------------------

AttitudeResult solveESOQ2(
    const StarTrackerObservation& obs,
    const std::vector<Star>& catalog,
    const std::vector<double>& weights) {
    auto result = solveQUEST(obs, catalog, weights);
    result.method = "ESOQ2";
    return result;
}

// ---------------------------------------------------------------------------
// Lost-in-Space (triangle matching)
// ---------------------------------------------------------------------------

AttitudeResult solveLostInSpace(
    const StarTrackerObservation& obs,
    const std::vector<Star>& catalog,
    double angularTolerance) {

    AttitudeResult result;
    result.method = "LostInSpace";

    if (obs.stars.size() < 3) {
        result.valid = false;
        return result;
    }

    double tolRad = angularTolerance * ARCSEC_TO_RAD;

    // Build inter-star angle matrix for observed stars (body frame)
    int N = static_cast<int>(obs.stars.size());
    std::vector<std::vector<double>> obsAngles(N, std::vector<double>(N));
    for (int i = 0; i < N; ++i)
        for (int j = i+1; j < N; ++j) {
            double angle = std::acos(
                std::max(-1.0, std::min(1.0,
                    dot3(obs.stars[i].bodyVector, obs.stars[j].bodyVector))));
            obsAngles[i][j] = angle;
            obsAngles[j][i] = angle;
        }

    // Precompute bright catalog star directions
    std::vector<Star> brightStars;
    for (const auto& s : catalog) {
        if (s.vmag <= obs.stars.back().magnitude + 1.0) {
            brightStars.push_back(s);
        }
        if (brightStars.size() > 5000) break;  // limit search space
    }

    int M = static_cast<int>(brightStars.size());

    // Try matching the first triangle (stars 0, 1, 2) against catalog
    double a01 = obsAngles[0][1];
    double a02 = obsAngles[0][2];
    double a12 = obsAngles[1][2];

    for (int ci = 0; ci < M; ++ci) {
        Vector3 di = brightStars[ci].direction();
        for (int cj = ci+1; cj < M; ++cj) {
            Vector3 dj = brightStars[cj].direction();
            double catAngle_ij = std::acos(
                std::max(-1.0, std::min(1.0, dot3(di, dj))));

            if (std::abs(catAngle_ij - a01) > tolRad) continue;

            for (int ck = cj+1; ck < M; ++ck) {
                Vector3 dk = brightStars[ck].direction();
                double catAngle_ik = std::acos(
                    std::max(-1.0, std::min(1.0, dot3(di, dk))));
                double catAngle_jk = std::acos(
                    std::max(-1.0, std::min(1.0, dot3(dj, dk))));

                if (std::abs(catAngle_ik - a02) > tolRad) continue;
                if (std::abs(catAngle_jk - a12) > tolRad) continue;

                // Match found! Use TRIAD with first two matched stars
                StarTrackerObservation matchedObs = obs;
                matchedObs.stars[0].hipId = brightStars[ci].hipId;
                matchedObs.stars[1].hipId = brightStars[cj].hipId;

                result = solveTriad(matchedObs, catalog);
                result.method = "LostInSpace";
                result.starsUsed = 3;
                return result;
            }
        }
    }

    result.valid = false;
    return result;
}

// ---------------------------------------------------------------------------
// Cesium World-Space Integration
// ---------------------------------------------------------------------------

FOVCone computeFOVCone(const CameraPointing& pointing) {
    Quaternion qInv = quaternionInverse(pointing.attitudeECI);
    double fovXRad = pointing.fov.fovX * DEG_TO_RAD;
    double fovYRad = pointing.fov.fovY * DEG_TO_RAD;
    double halfX = fovXRad / 2.0;
    double halfY = fovYRad / 2.0;

    FOVCone cone;

    // Boresight in ECI then to ECEF
    Vector3 boresightECI = rotateByQuaternion({0,0,1}, qInv);
    cone.boresightECEF = eciToECEF(boresightECI, pointing.epoch);
    cone.halfAngle = std::max(halfX, halfY);

    // Corner directions in body frame
    double tx = std::tan(halfX), ty = std::tan(halfY);
    Vector3 bodyCorners[4] = {
        normalize3({-tx,  ty, 1}),  // TL
        normalize3({ tx,  ty, 1}),  // TR
        normalize3({ tx, -ty, 1}),  // BR
        normalize3({-tx, -ty, 1})   // BL
    };

    for (int i = 0; i < 4; ++i) {
        Vector3 eci = rotateByQuaternion(bodyCorners[i], qInv);
        cone.corners[i] = eciToECEF(eci, pointing.epoch);
    }

    return cone;
}

std::vector<WorldSpaceStar> projectToWorldSpace(
    const StarTrackerObservation& obs,
    const CameraPointing& pointing) {

    Quaternion qInv = quaternionInverse(pointing.attitudeECI);
    double bigDist = 1e9;

    std::vector<WorldSpaceStar> result;
    result.reserve(obs.stars.size());

    for (const auto& star : obs.stars) {
        // Body frame → ECI → ECEF
        Vector3 eciDir = rotateByQuaternion(star.bodyVector, qInv);
        Vector3 ecefDir = eciToECEF(eciDir, pointing.epoch);

        WorldSpaceStar ws;
        ws.hipId = star.hipId;
        ws.positionECEF = scale3(ecefDir, bigDist);
        ws.magnitude = star.magnitude;
        ws.sensorX = star.sensorX;
        ws.sensorY = star.sensorY;
        ws.identified = (star.hipId > 0);
        result.push_back(ws);
    }

    return result;
}

}  // namespace starfield
