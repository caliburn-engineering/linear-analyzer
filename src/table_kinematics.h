#pragma once

#include <Eigen/Dense>
#include <array>
#include <optional>
#include <cmath>

/// Parameters defining the 3-RRS parallel mechanism geometry
struct TableParams {
    double R_ground;    // Ground circle radius [m]
    double R_table;     // Table circle radius [m]
    double L1;          // Lower leg length [m]
    double L2;          // Upper leg length [m]
    double alpha_min;   // Minimum servo angle [rad] (≈ 10°)
    double alpha_max;   // Maximum servo angle [rad] (≈ 80°)
    double gamma_offset = 0.0; // Angular offset between ground and table triads [rad]
};

/// Table pose: roll (phi), pitch (theta), heave (z_c)
struct TablePose {
    double phi;     // Roll about X-axis [rad]
    double theta;   // Pitch about Y-axis [rad]
    double z_c;     // Table center height [m]
};

/// Result of inverse kinematics
struct IKResult {
    std::array<double, 3> alpha;  // Servo angles [rad]
    bool feasible;                // All solutions exist and within limits
};

/// Result of forward kinematics
struct FKResult {
    TablePose pose;
    int iterations;
    double residual_norm;
    bool converged;
};

class TableKinematics {
public:
    explicit TableKinematics(const TableParams& params);

    // --- Geometry ---

    /// Ground attachment point for leg i (i = 0, 1, 2)
    Eigen::Vector3d ground_point(int i) const;

    /// Radial unit vector at leg i (points outward from center)
    Eigen::Vector3d radial_dir(int i) const;

    /// Tangent unit vector at leg i (servo axis, tangent to ground circle)
    Eigen::Vector3d tangent_dir(int i) const;

    /// Angular position of leg i on the ground circle [rad]
    double beta(int i) const;

    /// Angular position of leg i on the table circle [rad]
    double gamma(int i) const;

    // --- Knee Positions ---

    /// Knee position given servo angle (simplified: L1 in radial-vertical plane)
    Eigen::Vector3d knee_position(int i, double alpha_i) const;

    // --- Table Attachment Points ---

    /// Table attachment point in table-local frame
    Eigen::Vector3d table_point_local(int i) const;

    /// Table rotation matrix from local to world: Ry(theta) * Rx(phi)
    Eigen::Matrix3d table_rotation(double phi, double theta) const;

    /// Table attachment point in world frame
    Eigen::Vector3d table_point_world(int i, const TablePose& pose) const;

    // --- Inverse Kinematics (closed-form) ---

    /// Solve for servo angles given desired table pose.
    /// Returns nullopt if any leg has no geometric solution.
    IKResult inverse_kinematics(const TablePose& pose) const;

    /// Solve IK for a single leg. Returns up to 2 solutions.
    /// Picks the one within [alpha_min, alpha_max], or nearest.
    std::optional<double> inverse_kinematics_leg(int i, const TablePose& pose) const;

    /// Check if a table pose is reachable (all IK solutions exist and within limits)
    bool is_feasible(const TablePose& pose) const;

    // --- Forward Kinematics (Newton-Raphson) ---

    /// Solve for table pose given servo angles.
    /// x0 is the initial guess. Uses Newton-Raphson iteration.
    FKResult forward_kinematics(const std::array<double, 3>& alpha,
                                 const TablePose& x0,
                                 int max_iter = 20,
                                 double tol = 1e-10) const;

    /// The level pose for a leg triple: `home_pose` of their mean angle.
    ///
    /// The seed every forward solve falls back to, and the safe end of every
    /// retreat.  It was written out at six call sites before it had a name.
    TablePose level_pose(const std::array<double, 3>& alpha) const;

    /// Can the mechanism be assembled at this leg triple at all?
    ///
    /// The servo travel limits are a box and the workspace is not — barely
    /// half the box has an assembly — so a per-leg clamp can ask for a
    /// configuration that does not exist.
    bool can_assemble(const std::array<double, 3>& alpha) const;

    /// Forward kinematics that stays on the assembly the machine is built in.
    ///
    /// `forward_kinematics` is a Newton solver, and the constraint equations
    /// have two roots.  Beside the built assembly — table above its knees —
    /// there is a FOLDED one with the table lying flat on the base, and for
    /// `R_ground == R_table` with `L1 == L2` that folded root is exact:
    /// `z_c = 0` satisfies every leg's length constraint to machine precision
    /// at every servo angle.  A residual norm cannot tell the two apart, so a
    /// converged solve is not by itself an answer about a mechanism.
    ///
    /// This warm-starts from `seed` for the two-iteration convergence a 60 Hz
    /// loop wants, then checks which root it landed on.  If it is the folded
    /// one — because the seed was already folded, or because a hard manoeuvre
    /// carried the pose down through the knee plane — it re-solves from the
    /// analytic level pose, which is on the built assembly by construction.
    ///
    /// Returns `converged = false` when no built assembly exists for `alpha`
    /// at all.  That is a real case, not a solver failure: the servo travel
    /// limits are a box and the workspace is not, so a per-leg clamp can ask
    /// for a configuration the mechanism cannot make.  Callers should keep the
    /// pose they had — a mechanism driven into a singularity binds and stops,
    /// it does not lie flat.
    ///
    /// See [#22](https://github.com/caliburn-engineering/caliburn/issues/22).
    FKResult solve_pose(const std::array<double, 3>& alpha,
                        const TablePose& seed) const;

    /// The heave a pose has to clear to be the built assembly rather than the
    /// folded one: a tenth of the mean height of the three knees.
    ///
    /// A tenth looks arbitrary and is not.  The folded assembly is not a
    /// family of poses — it is exactly one, `(phi, theta, z_c) = (0, 0, 0)`,
    /// at every servo angle, because a table sitting concentric on the base is
    /// exactly `L2` from each knee when `R_ground == R_table` and
    /// `L1 == L2`.  So the floor is separating a population from a single
    /// point at the origin, and the only way to get it wrong is to set it high
    /// enough to reject something real.
    ///
    /// Measured, over 300k samples of the servo box, against assemblies
    /// verified by round-tripping the pose back through closed-form IK: the
    /// built population's median sits at 1.98 mean knee heights and its
    /// minimum at 0.163, with 35 of 142301 below 0.5.  Those are near-singular
    /// configurations with the table close to its own knee plane, and they are
    /// real.  A floor at the full mean knee height would have thrown away
    /// every one of them; at a tenth, the lowest clears by 1.6x and the folded
    /// root still has to climb from nothing.
    ///
    /// Both ends are pinned by `test_assembly_mode`, so a parameter change
    /// that moves either population fails loudly rather than quietly eroding
    /// the margin.
    double assembly_floor(const std::array<double, 3>& alpha) const;

    /// FK residual: f_i = ||P_i - K_i||^2 - L2^2
    Eigen::Vector3d fk_residual(const std::array<double, 3>& alpha,
                                 const TablePose& pose) const;

    /// FK Jacobian: J_ij = df_i / dx_j, x = [phi, theta, z_c]
    Eigen::Matrix3d fk_jacobian(const std::array<double, 3>& alpha,
                                 const TablePose& pose) const;

    // --- Velocity Jacobian ---

    /// Constraint Jacobian w.r.t. servo angles: J_alpha(i,i) = df_i/d(alpha_i)
    /// Diagonal matrix (each constraint depends on only one servo)
    Eigen::Matrix3d constraint_jacobian_alpha(const std::array<double, 3>& alpha,
                                               const TablePose& pose) const;

    /// Derivative of knee position w.r.t. servo angle
    Eigen::Vector3d dKnee_dAlpha(int i, double alpha_i) const;

    /// Velocity Jacobian: d(pose)/dt = J_v * d(alpha)/dt
    /// Maps servo velocities to table pose velocities.
    /// J_v = -J_pose^{-1} * J_alpha
    Eigen::Matrix3d velocity_jacobian(const std::array<double, 3>& alpha,
                                       const TablePose& pose) const;

    // --- Singularity Analysis ---

    /// Manipulability index: sqrt(det(J_v * J_v^T))
    /// Zero at singularity, higher = better omnidirectional control
    double manipulability(const std::array<double, 3>& alpha,
                          const TablePose& pose) const;

    /// Condition number of velocity Jacobian (sigma_max / sigma_min)
    /// 1 = isotropic (ideal), infinity = singular
    double condition_number(const std::array<double, 3>& alpha,
                            const TablePose& pose) const;

    /// Determinant of the constraint Jacobian w.r.t. pose (J_pose)
    /// Zero = type-1 (forward) singularity — table has uncontrollable motion
    double det_J_pose(const std::array<double, 3>& alpha,
                      const TablePose& pose) const;

    /// Determinant of the constraint Jacobian w.r.t. alpha (J_alpha)
    /// Zero = type-2 (inverse) singularity — servo motion doesn't affect constraint
    double det_J_alpha(const std::array<double, 3>& alpha,
                       const TablePose& pose) const;

    // --- Home position ---

    /// Default pose when all servos at given angle (typically 45°)
    TablePose home_pose(double alpha_home = M_PI / 4.0) const;

    const TableParams& params() const { return params_; }

private:
    TableParams params_;

    /// Partial derivative of table rotation w.r.t. phi
    Eigen::Matrix3d dR_dphi(double phi, double theta) const;

    /// Partial derivative of table rotation w.r.t. theta
    Eigen::Matrix3d dR_dtheta(double phi, double theta) const;
};
