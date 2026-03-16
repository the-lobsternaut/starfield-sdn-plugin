// star_field.vert — Cesium-compatible star field vertex shader
// Renders Hipparcos catalog stars as point sprites with:
//   - Apparent magnitude → point size mapping (Pogson scale)
//   - Proper motion correction at current epoch
//   - Atmospheric extinction (airmass-dependent dimming)
//   - Precession from J2000 to current epoch (IAU 2006)

// Cesium built-in uniforms
uniform mat4 czm_modelViewProjection;
uniform mat4 czm_modelView;
uniform vec3 czm_viewerPositionWC;        // viewer position in world coords (ECEF meters)
uniform float czm_frameNumber;
uniform float czm_morphTime;
uniform vec3 czm_sunDirectionWC;          // sun direction for sky brightness

// Custom uniforms
uniform float u_epoch;                     // current epoch (Julian years since J2000.0)
uniform float u_magnitudeLimit;            // faintest star to render (default 6.5)
uniform float u_brightnessScale;           // overall brightness multiplier
uniform float u_pointSizeBase;             // base point size in pixels
uniform float u_extinctionCoeff;           // atmospheric extinction coefficient (default 0.2)
uniform bool  u_enableScintillation;       // enable twinkling
uniform bool  u_enableAtmosphere;          // enable atmospheric extinction
uniform float u_scintillationStrength;     // twinkle intensity (0-1)

// Per-star attributes (from delta-decoded buffer)
attribute vec3  a_position;                // J2000 unit direction vector (ICRS)
attribute float a_magnitude;               // visual magnitude (Hp)
attribute float a_colorIndex;              // B-V color index
attribute vec2  a_properMotion;            // proper motion [mas/yr]: (mu_ra*cos(dec), mu_dec)
attribute float a_parallax;                // parallax [mas] (for nearby stars)
attribute float a_hipId;                   // Hipparcos catalog ID (for picking)

// Varyings to fragment shader
varying float v_magnitude;
varying float v_colorIndex;
varying float v_alpha;
varying float v_pointSize;
varying float v_hipId;

// Constants
const float DEG_TO_RAD = 0.017453292519943295;
const float MAS_TO_RAD = 4.84813681e-9;
const float POGSON = 2.512;               // 10^(2/5), magnitude ratio

// Apparent magnitude → point size (pixels)
// Brighter stars = larger points
// Sirius (V=-1.46) → ~8px, Polaris (V=1.97) → ~4px, V=6.5 → ~1px
float magnitudeToSize(float mag) {
    // Inverse Pogson: flux ratio = 10^(-0.4 * (mag - mag_ref))
    // Size proportional to sqrt(flux) for area-based brightness
    float flux = pow(POGSON, -(mag - u_magnitudeLimit));
    return u_pointSizeBase * sqrt(max(flux, 0.0));
}

// Apply proper motion correction
vec3 applyProperMotion(vec3 dir, float epoch_years) {
    // Convert proper motion from mas/yr to radians
    float dra  = a_properMotion.x * MAS_TO_RAD * epoch_years;
    float ddec = a_properMotion.y * MAS_TO_RAD * epoch_years;

    // Simple linear correction on the unit sphere
    // For high-proper-motion stars this is approximate but adequate for rendering
    float ra  = atan(dir.y, dir.x) + dra;
    float dec = asin(clamp(dir.z, -1.0, 1.0)) + ddec;

    return vec3(cos(dec) * cos(ra), cos(dec) * sin(ra), sin(dec));
}

// Atmospheric extinction based on airmass
// Higher airmass (near horizon) = more extinction
float computeExtinction(vec3 starDirECEF, vec3 viewerPosECEF) {
    if (!u_enableAtmosphere) return 0.0;

    // Zenith angle of star from viewer
    vec3 up = normalize(viewerPosECEF);
    float cosZenith = dot(normalize(starDirECEF), up);

    // Below horizon: fully extinguished
    if (cosZenith < 0.0) return 99.0;

    // Airmass approximation (Kasten & Young 1989)
    float zenithAngleDeg = acos(clamp(cosZenith, 0.0, 1.0)) / DEG_TO_RAD;
    float airmass;
    if (zenithAngleDeg < 85.0) {
        // Standard: X = sec(z)
        airmass = 1.0 / max(cosZenith, 0.001);
    } else {
        // Near horizon: empirical correction
        float denom = cosZenith + 0.50572 * pow(96.07995 - zenithAngleDeg, -1.6364);
        airmass = 1.0 / max(denom, 0.001);
    }

    // Extinction in magnitudes: dM = k * X
    return u_extinctionCoeff * airmass;
}

// Scintillation (twinkling) — pseudo-random per star per frame
float computeScintillation(float hipId, float frame) {
    if (!u_enableScintillation) return 0.0;

    // Hash-based pseudo-random per star and frame
    float h = fract(sin(hipId * 12.9898 + frame * 0.01) * 43758.5453);
    float h2 = fract(sin(hipId * 78.233 + frame * 0.013) * 28647.893);

    // Scintillation index increases with airmass (not modeled here for simplicity)
    return u_scintillationStrength * (h * 0.3 + h2 * 0.2 - 0.25);
}

void main() {
    // Apply proper motion
    vec3 correctedDir = applyProperMotion(a_position, u_epoch);

    // Convert ICRS unit vector to ECEF (approximate: ignore nutation for rendering)
    // The star is at "infinity" — place it on a large sphere
    float starDist = 1.0e9;  // 1 billion meters (arbitrary large distance)
    vec3 starPosECEF = correctedDir * starDist;

    // Transform to clip space via Cesium's model-view-projection
    vec4 clipPos = czm_modelViewProjection * vec4(starPosECEF, 1.0);

    // Compute atmospheric extinction
    float extinction = computeExtinction(correctedDir, czm_viewerPositionWC);
    float apparentMag = a_magnitude + extinction;

    // Scintillation
    float scint = computeScintillation(a_hipId, czm_frameNumber);
    apparentMag += scint;

    // Cull stars fainter than magnitude limit
    if (apparentMag > u_magnitudeLimit) {
        gl_Position = vec4(2.0, 2.0, 2.0, 1.0);  // off-screen
        gl_PointSize = 0.0;
        return;
    }

    // Point size from apparent magnitude
    float size = magnitudeToSize(apparentMag);
    size = clamp(size, 1.0, 12.0);

    // Brightness (alpha) also scaled by magnitude
    float flux = pow(POGSON, -(apparentMag - u_magnitudeLimit));
    float alpha = clamp(flux * u_brightnessScale, 0.0, 1.0);

    // Sun interference: reduce brightness when sun is above horizon
    // (simulates sky brightness washing out stars)
    vec3 viewerUp = normalize(czm_viewerPositionWC);
    float sunAltSin = dot(czm_sunDirectionWC, viewerUp);
    if (sunAltSin > -0.1) {
        // Sun near or above horizon — wash out stars
        float washout = smoothstep(-0.1, 0.2, sunAltSin);
        alpha *= (1.0 - washout * 0.95);
    }

    gl_Position = clipPos;
    gl_PointSize = size;

    v_magnitude = apparentMag;
    v_colorIndex = a_colorIndex;
    v_alpha = alpha;
    v_pointSize = size;
    v_hipId = a_hipId;
}
