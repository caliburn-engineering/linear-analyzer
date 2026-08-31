// tests/cascade_fixture.h
#pragma once

// The three facts every closed-loop test on the Ball-Balancer Cascade needs:
// which preset it is, what gain the application opens on, and how big the ball
// is.  They live here rather than in each test file because a second
// transcription is how two tests come to disagree about the plant they are
// both claiming to measure — the same reason `cascadeMechanism` and friends
// are the single source for the application itself.
//
// The simulation LOOPS are deliberately not shared.  `test_auto_balance`
// measures settling and peak against a setpoint; `test_attract_mode` records a
// whole trace and needs the disturbance schedule inside the loop.  One harness
// serving both would take a parameter per difference and answer neither
// question clearly.

#include "analysis/lqr.h"
#include "analysis/model_library.h"
#include "auto_balance.h"
#include "ball_sim.h"
#include "rolling_dynamics.h"
#include "table_kinematics.h"
#include "test_helpers.h"

#include <Eigen/Core>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace caliburn {

/// The simulator's ball, not a copy of it: the friction dead band that
/// `test_auto_balance` asserts about is a claim on this exact rolling
/// resistance.
constexpr double kFixtureBallRadius = kPlateBall.radius;
constexpr double kFixtureBallFriction = kPlateBall.rolling_friction;

inline const ModelEntry& cascadeModel(const std::vector<ModelEntry>& models) {
    for (const auto& m : models)
        if (isCascadeModel(m)) return m;
    std::fprintf(stderr, "FAIL: no cascade model in the library\n");
    std::exit(1);
}

/// The plate itself.  `legCommand` needs it because the servo travel limits
/// are a box and the workspace is not — see issue #22 — so every test that
/// evaluates the loop needs a mechanism to check the command against.
inline const TableKinematics& cascadeKinematics() {
    static const TableKinematics tk(
        cascadeMechanism(cascadeModel(getBuiltinModels()).params));
    return tk;
}

inline Eigen::MatrixXd gainFor(const ModelEntry& e,
                               const Eigen::VectorXd& q,
                               const Eigen::VectorXd& r) {
    const LqrResult res = computeLQR(e.system,
                                     q.asDiagonal().toDenseMatrix(),
                                     r.asDiagonal().toDenseMatrix());
    ASSERT_TRUE(res.success);
    return res.K;
}

/// The gain the application actually opens on — the weights from
/// `defaultLqrStateWeights` / `defaultLqrInputWeights`, solved against this
/// preset.  Every "does the demo work" claim is a claim about this matrix.
inline Eigen::MatrixXd defaultGain(const ModelEntry& e) {
    return gainFor(e, defaultLqrStateWeights(7), defaultLqrInputWeights(3));
}

}  // namespace caliburn
