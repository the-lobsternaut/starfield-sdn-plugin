// star_field.frag — Cesium-compatible star field fragment shader
// Renders each star as a point sprite with:
//   - B-V color index → realistic RGB color mapping
//   - Airy disk point spread function (PSF) for natural star appearance
//   - Diffraction spikes (optional, for bright stars)
//   - Smooth alpha falloff at edges

#ifdef GL_FRAGMENT_PRECISION_HIGH
    precision highp float;
#else
    precision mediump float;
#endif

// Varyings from vertex shader
varying float v_magnitude;
varying float v_colorIndex;
varying float v_alpha;
varying float v_pointSize;
varying float v_hipId;

// Custom uniforms
uniform bool  u_enableDiffraction;         // show diffraction spikes on bright stars
uniform float u_diffractionThreshold;      // magnitude threshold for spikes (default 2.0)
uniform float u_gammaCorrection;           // display gamma (default 2.2)

// ============================================================
// B-V Color Index → RGB Color
// Based on Ballesteros (2012) formula + empirical correction
// Maps the full range of stellar spectral types:
//   O stars (B-V < -0.3): blue-white
//   B stars (B-V ~ -0.2): blue
//   A stars (B-V ~ 0.0):  white
//   F stars (B-V ~ 0.3):  yellow-white
//   G stars (B-V ~ 0.6):  yellow (Sun)
//   K stars (B-V ~ 1.0):  orange
//   M stars (B-V > 1.4):  red
// ============================================================

vec3 bvToRGB(float bv) {
    // Clamp to valid range
    bv = clamp(bv, -0.4, 2.0);

    float r, g, b;

    // Red channel
    if (bv < -0.2) {
        r = 0.61 + 0.11 * bv + 0.1 * bv * bv;
    } else if (bv < 0.0) {
        r = 0.83 + 0.7 * bv;
    } else if (bv < 0.4) {
        r = 0.83 + 0.5 * bv;
    } else {
        r = min(1.0, 0.71 + 0.42 * bv);
    }

    // Green channel
    if (bv < 0.0) {
        g = 0.70 + 0.07 * bv + 0.1 * bv * bv;
    } else if (bv < 0.4) {
        g = 0.87 + 0.11 * bv;
    } else if (bv < 1.5) {
        g = 1.0 - 0.43 * (bv - 0.4);
    } else {
        g = max(0.2, 0.53 - 0.2 * (bv - 1.5));
    }

    // Blue channel
    if (bv < -0.2) {
        b = 1.0;
    } else if (bv < 0.4) {
        b = 1.0 - 1.67 * bv;
    } else if (bv < 1.0) {
        b = max(0.0, 0.33 - 0.55 * (bv - 0.4));
    } else {
        b = 0.0;
    }

    return vec3(clamp(r, 0.0, 1.0), clamp(g, 0.0, 1.0), clamp(b, 0.0, 1.0));
}

// ============================================================
// Point Spread Function (PSF)
// Approximation of an Airy disk for natural star appearance
// Core: bright Gaussian center
// Halo: wider Gaussian wings
// ============================================================

float starPSF(vec2 uv, float magnitude) {
    float dist = length(uv);

    // Core Gaussian (tight, bright center)
    float core = exp(-dist * dist * 8.0);

    // Halo (wider, dimmer wings — more visible for bright stars)
    float haloStrength = max(0.0, (3.0 - magnitude) * 0.1);
    float halo = haloStrength * exp(-dist * dist * 2.0);

    return core + halo;
}

// ============================================================
// Diffraction Spikes (for bright stars)
// 4-pointed cross pattern
// ============================================================

float diffractionSpikes(vec2 uv, float magnitude) {
    if (!u_enableDiffraction) return 0.0;
    if (magnitude > u_diffractionThreshold) return 0.0;

    float brightness = max(0.0, (u_diffractionThreshold - magnitude) * 0.15);
    float dist = length(uv);

    // 4-pointed cross
    float spike1 = exp(-abs(uv.x) / 0.02) * exp(-uv.y * uv.y * 200.0);
    float spike2 = exp(-abs(uv.y) / 0.02) * exp(-uv.x * uv.x * 200.0);

    // Diagonal spikes (45°)
    vec2 rotUV = vec2(uv.x * 0.707 + uv.y * 0.707,
                      -uv.x * 0.707 + uv.y * 0.707);
    float spike3 = exp(-abs(rotUV.x) / 0.02) * exp(-rotUV.y * rotUV.y * 200.0);
    float spike4 = exp(-abs(rotUV.y) / 0.02) * exp(-rotUV.x * rotUV.x * 200.0);

    return brightness * (spike1 + spike2 + spike3 * 0.5 + spike4 * 0.5) * exp(-dist * 3.0);
}

void main() {
    // Point sprite UV: (0,0) at center, range [-0.5, 0.5]
    vec2 uv = gl_PointCoord - vec2(0.5);
    float dist = length(uv);

    // Discard outside circle
    if (dist > 0.5) discard;

    // Get star color from B-V index
    vec3 starColor = bvToRGB(v_colorIndex);

    // PSF (brightness profile)
    float psf = starPSF(uv, v_magnitude);

    // Diffraction spikes
    float spikes = diffractionSpikes(uv, v_magnitude);

    // Combine
    float brightness = psf + spikes;
    brightness *= v_alpha;

    // Apply gamma correction
    vec3 color = starColor * brightness;
    color = pow(color, vec3(1.0 / u_gammaCorrection));

    // Final alpha: smooth falloff
    float alpha = brightness * smoothstep(0.5, 0.3, dist);

    if (alpha < 0.001) discard;

    gl_FragColor = vec4(color, alpha);
}
