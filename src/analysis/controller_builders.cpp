// src/analysis/controller_builders.cpp
#include "controller_builders.h"
#include "system_connect.h"

namespace caliburn {

LinearSystem buildPID(const PIDParams& p) {
    // Emit a state only when its output coefficient is nonzero.  A state with a
    // zero entry in C is unobservable, survives seriesConnect and
    // feedbackConnect untouched, and shows up as a phantom pole in the
    // closed-loop pole plot: PI with Ki = 0 puts a marginally-stable pole at
    // the origin, PD with Kd = 0 puts one at -1/tau_f.  Emitting the state does
    // not waste space, it lies to the plots.
    //
    // The threshold is exact zero, no epsilon.  The integrator's A block is [0]
    // regardless of Ki (Ki enters only C), so a tiny Ki is numerically harmless
    // and genuinely is a near-integrator; a dead band would silently contradict
    // a slider reading 1e-9.  The left slider stop writes exact 0.0f.
    //
    // No defensive guard on tau_f <= 0: it is a documented precondition, the
    // sliders enforce the range, and a silent clamp here would hide a
    // programming error at the one place it is cheap to catch.
    const bool has_i = (p.Ki != 0.0f);
    const bool has_d = (p.Kd != 0.0f);
    const int n = (has_i ? 1 : 0) + (has_d ? 1 : 0);

    const double Kp = p.Kp, Ki = p.Ki, Kd = p.Kd, tau_f = p.tau_f;

    LinearSystem c;
    c.A = Eigen::MatrixXd::Zero(n, n);
    c.B = Eigen::MatrixXd::Zero(n, 1);
    c.C = Eigen::MatrixXd::Zero(1, n);
    c.D = Eigen::MatrixXd::Constant(1, 1, Kp);

    int k = 0;
    if (has_i) {  // integrator: Ki / s
        c.A(k, k) = 0.0;
        c.B(k, 0) = 1.0;
        c.C(0, k) = Ki;
        ++k;
    }
    if (has_d) {  // filtered derivative: Kd*s / (tau_f*s + 1)
        c.A(k, k) = -1.0 / tau_f;
        c.B(k, 0) = 1.0 / tau_f;
        c.C(0, k) = -Kd / tau_f;
        c.D(0, 0) += Kd / tau_f;
    }
    return c;
}

namespace {

// One first-order section (s + z) / (s + p), z = 1/T, p = 1/(alpha*T).
// Same observability rule as buildPID: the output coefficient is z - p, which
// vanishes at alpha == 1, so that case drops to a pure gain with no state.
// Unreachable through the sliders (alpha_lead is clamped to [0.01, 0.99],
// alpha_lag to [1.01, 100]) — handled anyway, because a builder that is total
// in alpha does not depend on slider clamps staying where they are.
LinearSystem leadLagSection(double alpha, double T) {
    const double z = 1.0 / T;
    const double pole = 1.0 / (alpha * T);
    const double coeff = z - pole;
    const int n = (coeff != 0.0) ? 1 : 0;

    LinearSystem s;
    s.A = Eigen::MatrixXd::Zero(n, n);
    s.B = Eigen::MatrixXd::Zero(n, 1);
    s.C = Eigen::MatrixXd::Zero(1, n);
    s.D = Eigen::MatrixXd::Constant(1, 1, 1.0);
    if (n == 1) {
        s.A(0, 0) = -pole;
        s.B(0, 0) = 1.0;
        s.C(0, 0) = coeff;
    }
    return s;
}

}  // anonymous namespace

LinearSystem buildLeadLag(const LeadLagParams& p) {
    // Kc == 0 zeroes every output coefficient at once, so the same rule says
    // emit no states.  Also unreachable through the sliders (Kc >= 0.1).
    if (p.Kc == 0.0f) {
        LinearSystem z;
        z.A = Eigen::MatrixXd::Zero(0, 0);
        z.B = Eigen::MatrixXd::Zero(0, 1);
        z.C = Eigen::MatrixXd::Zero(1, 0);
        z.D = Eigen::MatrixXd::Zero(1, 1);
        return z;
    }

    LinearSystem c;
    switch (p.mode) {
        case CompensatorMode::Lead:
            c = leadLagSection(p.alpha_lead, p.T_lead);
            break;
        case CompensatorMode::Lag:
            c = leadLagSection(p.alpha_lag, p.T_lag);
            break;
        case CompensatorMode::LeadLag:
            // Each section drops independently: alpha_lead == 1 gives 1 state,
            // both == 1 gives 0.  seriesConnect handles n1 == 0 exactly (#5).
            c = seriesConnect(leadLagSection(p.alpha_lead, p.T_lead),
                              leadLagSection(p.alpha_lag, p.T_lag));
            break;
    }
    c.C *= static_cast<double>(p.Kc);
    c.D *= static_cast<double>(p.Kc);
    return c;
}

LinearSystem buildLoopController(const std::vector<Loop>& loops,
                                 LoopKind kind,
                                 const LinearSystem& plant) {
    const int p = plant.outputs();
    const int m = plant.inputs();

    // Build each loop's SISO block first — the total state count is not known
    // until every block has decided how minimal it is.  The assembler calls the
    // SISO builders rather than inlining the realization tables a second time,
    // so the observability rule lives in exactly one place.
    std::vector<LinearSystem> blocks;
    std::vector<const Loop*> kept;
    blocks.reserve(loops.size());
    kept.reserve(loops.size());
    int n_total = 0;

    for (const auto& l : loops) {
        // Out-of-range index: skip this loop and build the rest.  Unreachable
        // in principle — the (p, m) dimension guard reseeds wholesale before
        // any recompute — so this only decides what happens when it occurs
        // anyway.  Asserting would crash a UI that must keep rendering while
        // shipped builds scattered out of bounds; clamping would silently
        // control a different physical channel while carrying the old
        // channel's gains; bailing to all-zero would discard every other
        // loop's tuning.  Skipping degrades locally and behaves identically in
        // both build configs.
        if (l.out < 0 || l.out >= p || l.in < 0 || l.in >= m) continue;
        blocks.push_back(kind == LoopKind::PID ? buildPID(l.pid)
                                               : buildLeadLag(l.leadlag));
        kept.push_back(&l);
        n_total += blocks.back().states();
    }

    // Block-diagonal A; B row-block k scattered into column out_k; C
    // column-block k into row in_k; D accumulated at (in_k, out_k).  Only D
    // accumulates — the B and C scatters are unambiguous because each state
    // block belongs to exactly one loop.  An empty (or entirely skipped) list
    // yields the all-zero m x p system with states() == 0, which is issue #2's
    // "empty loop list is equivalent to ControllerType::None".
    LinearSystem c;
    c.A = Eigen::MatrixXd::Zero(n_total, n_total);
    c.B = Eigen::MatrixXd::Zero(n_total, p);
    c.C = Eigen::MatrixXd::Zero(m, n_total);
    c.D = Eigen::MatrixXd::Zero(m, p);

    int off = 0;
    for (std::size_t k = 0; k < blocks.size(); ++k) {
        const LinearSystem& b = blocks[k];
        const int nk = b.states();
        const int out = kept[k]->out;
        const int in = kept[k]->in;

        if (nk > 0) {
            c.A.block(off, off, nk, nk) = b.A;
            c.B.block(off, out, nk, 1) = b.B;
            c.C.block(in, off, 1, nk) = b.C;
        }
        c.D(in, out) += b.D(0, 0);
        off += nk;
    }
    return c;
}

}  // namespace caliburn
