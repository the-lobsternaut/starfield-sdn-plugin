#ifndef STARFIELD_CONSTANTS_H
#define STARFIELD_CONSTANTS_H

namespace starfield {

constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 6.28318530717958647692;
constexpr double DEG_TO_RAD = 0.017453292519943295;
constexpr double RAD_TO_DEG = 57.29577951308232;
constexpr double ARCSEC_TO_RAD = 4.84813681109536e-6;
constexpr double MAS_TO_RAD = 4.84813681109536e-9;

// J2000.0 epoch as Julian Date
constexpr double J2000_JD = 2451545.0;
constexpr double JULIAN_YEAR = 365.25;  // days
constexpr double SECONDS_PER_DAY = 86400.0;

// Earth
constexpr double R_EARTH = 6.3781e6;   // equatorial radius [m]
constexpr double OMEGA_EARTH = 7.2921159e-5;  // rotation rate [rad/s]

// Delta encoding units
constexpr double DELTA_RA_UNIT = 0.01 * ARCSEC_TO_RAD;   // 0.01 arcsec in radians
constexpr double DELTA_DEC_UNIT = 0.01 * ARCSEC_TO_RAD;

// Magnitude encoding
constexpr double VMAG_UNIT = 0.001;   // milli-magnitudes
constexpr double BV_UNIT = 0.001;     // milli-color-index
constexpr double PM_UNIT = 0.1;       // 0.1 mas/yr
constexpr double PARALLAX_UNIT = 0.01; // 0.01 mas

}  // namespace starfield

#endif  // STARFIELD_CONSTANTS_H
