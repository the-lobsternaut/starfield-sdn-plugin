/**
 * Starfield Plugin — 6DOF Integration Tests
 */

#include "starfield/sixdof_core.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace sixdof;

void testAttitudePropConsistency() {
    State s;
    s.quat = qidentity();
    s.omega = {0.02, -0.01, 0.015};
    s.mass = 50;

    InertiaTensor I = inertiaDiag(5, 8, 6);
    auto coastFn = [](const State&, double) -> ForcesTorques { return {}; };
    double dt = 0.01, t = 0;

    for (int i = 0; i < 3000; i++) {
        s = rk4Step(s, I, dt, t, coastFn); t += dt;
        assert(std::abs(qnorm(s.quat) - 1.0) < 1e-6);
    }

    Vec3 body_x = {1, 0, 0};
    Vec3 inertial = qrotate(s.quat, body_x);
    assert(std::abs(v3norm(inertial) - 1.0) < 1e-10);

    std::cout << "  Attitude propagation ✓ (30s)\n";
}

void testReactionWheelDesat() {
    State s;
    s.quat = qidentity();
    s.omega = {0.1, 0, 0};
    s.mass = 20;

    InertiaTensor I = inertiaDiag(2, 3, 2.5);
    double dt = 0.01, t = 0;

    auto desatFn = [](const State& st, double) -> ForcesTorques {
        ForcesTorques ft;
        if (st.omega[0] > 0.001) ft.torque_body = {-0.01, 0, 0};
        return ft;
    };

    for (int i = 0; i < 2000; i++) { s = rk4Step(s, I, dt, t, desatFn); t += dt; }

    assert(std::abs(s.omega[0]) < 0.01);
    std::cout << "  Reaction wheel desat ✓ (final_omega=" << s.omega[0] << ")\n";
}

void testBodyToInertialRoundtrip() {
    Quat q = {0.5, 0.5, 0.5, 0.5}; // 120° about (1,1,1)
    assert(std::abs(qnorm(q) - 1.0) < 1e-10);

    Vec3 v_body = {1.5, -2.3, 0.7};
    Vec3 v_inertial = qrotate(q, v_body);
    Vec3 v_back = qrotateInv(q, v_inertial);

    for (int i = 0; i < 3; i++) assert(std::abs(v_body[i] - v_back[i]) < 1e-10);

    std::cout << "  Body↔Inertial roundtrip ✓\n";
}

int main() {
    std::cout << "=== starfield 6DOF tests ===\n";
    testAttitudePropConsistency();
    testReactionWheelDesat();
    testBodyToInertialRoundtrip();
    std::cout << "All starfield 6DOF tests passed.\n";
    return 0;
}
