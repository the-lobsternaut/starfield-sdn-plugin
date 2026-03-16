// constellation_lines.frag — Constellation line fragment shader

#ifdef GL_FRAGMENT_PRECISION_HIGH
    precision highp float;
#else
    precision mediump float;
#endif

varying float v_opacity;
varying vec3  v_color;

void main() {
    if (v_opacity < 0.01) discard;
    gl_FragColor = vec4(v_color, v_opacity);
}
