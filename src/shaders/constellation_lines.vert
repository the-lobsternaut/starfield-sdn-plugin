// constellation_lines.vert — Cesium-compatible constellation line shader
// Renders IAU constellation boundaries and asterism stick figures

uniform mat4 czm_modelViewProjection;
uniform vec3 czm_viewerPositionWC;
uniform vec3 czm_sunDirectionWC;
uniform float u_epoch;               // for proper motion of endpoint stars
uniform float u_lineOpacity;         // base opacity (default 0.3)

attribute vec3 a_position;           // J2000 unit direction vector (ICRS)
attribute float a_constellation;     // constellation index (for coloring)

varying float v_opacity;
varying vec3  v_color;

// Constellation colors: 12 distinct hues cycling through the zodiac
vec3 constellationColor(float idx) {
    float hue = mod(idx * 0.618033988749895, 1.0);  // golden ratio spacing
    // HSV to RGB (simplified)
    float h = hue * 6.0;
    float c = 0.6;
    float x = c * (1.0 - abs(mod(h, 2.0) - 1.0));
    vec3 rgb;
    if      (h < 1.0) rgb = vec3(c, x, 0.0);
    else if (h < 2.0) rgb = vec3(x, c, 0.0);
    else if (h < 3.0) rgb = vec3(0.0, c, x);
    else if (h < 4.0) rgb = vec3(0.0, x, c);
    else if (h < 5.0) rgb = vec3(x, 0.0, c);
    else              rgb = vec3(c, 0.0, x);
    return rgb + vec3(0.3);  // lighten
}

void main() {
    float starDist = 1.0e9;
    vec3 starPosECEF = a_position * starDist;
    gl_Position = czm_modelViewProjection * vec4(starPosECEF, 1.0);

    // Fade when sun is up
    vec3 viewerUp = normalize(czm_viewerPositionWC);
    float sunAltSin = dot(czm_sunDirectionWC, viewerUp);
    float dayFade = 1.0 - smoothstep(-0.1, 0.1, sunAltSin);

    v_opacity = u_lineOpacity * dayFade;
    v_color = constellationColor(a_constellation);
}
