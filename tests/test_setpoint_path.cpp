// tests/test_setpoint_path.cpp
//
// The paths the ball is asked to trace.  Nothing here touches the controller —
// these are claims about the shapes themselves, and they are worth making
// separately because a tracking failure is much easier to read when the thing
// being tracked is known to be right.
#include "setpoint_path.h"

#include "test_helpers.h"

#include <cmath>

using namespace caliburn;

namespace {

SetpointPath shape(PathShape s, double r = 0.12, double period = 10.0) {
    SetpointPath p;
    p.shape = s;
    p.radius_m = r;
    p.period_s = period;
    return p;
}

// Every shape closes: a lap returns exactly where it started, at every size
// and period.  A path that drifts would send the ball off the plate eventually
// and the drift would be invisible for the first few laps.
void test_every_path_closes() {
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        for (double r : {0.04, 0.12, 0.20}) {
            for (double T : {3.0, 10.0, 30.0}) {
                const SetpointPath p = shape(s, r, T);
                for (int lap = 1; lap <= 4; ++lap) {
                    const Eigen::Vector2d a = pathPoint(p, 0.7);
                    const Eigen::Vector2d b = pathPoint(p, 0.7 + lap * T);
                    ASSERT_NEAR((a - b).norm(), 0.0, 1e-12);
                }
            }
        }
    }
}

// Nothing ever goes outside the radius asked for.  This is what lets the plate
// bound be set once, against the radius slider, rather than per shape.
void test_no_path_leaves_its_radius() {
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        const SetpointPath p = shape(s, 0.12);
        for (int i = 0; i < 2000; ++i) {
            const double t = 10.0 * i / 2000.0;
            ASSERT_TRUE(pathPoint(p, t).norm() <= 0.12 + 1e-12);
        }
    }
}

// And each one actually reaches it, at its corners.  Sizing from the
// circumradius is what makes the shapes comparable at the same setting.
void test_every_shape_reaches_its_radius() {
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        const SetpointPath p = shape(s, 0.12);
        double furthest = 0.0;
        for (int i = 0; i < 2000; ++i)
            furthest = std::max(furthest, pathPoint(p, 10.0 * i / 2000.0).norm());
        ASSERT_NEAR(furthest, 0.12, 1e-6);
    }
}

// Constant speed, not constant angle.  A polygon traversed by sweeping an
// angle would crawl through its corners, which is precisely where the
// interesting behaviour is.
void test_the_polygons_are_walked_at_constant_speed() {
    for (PathShape s : {PathShape::Square, PathShape::Triangle}) {
        const SetpointPath p = shape(s, 0.12, 6.0);
        const double expect = pathLength(p) / p.period_s;
        for (int i = 0; i < 600; ++i) {
            const double t = 6.0 * (i + 0.5) / 600.0;   // off the corners
            ASSERT_NEAR(pathVelocity(p, t).norm(), expect, 1e-9);
        }
    }
}

// The velocity is the position's derivative, which is the claim a feedforward
// term would rest on.  Checked against a difference of the position, away from
// the corners where the derivative genuinely does not exist.
void test_velocity_is_the_derivative_of_position() {
    const double h = 1e-6;
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        const SetpointPath p = shape(s, 0.12, 7.0);
        for (int i = 0; i < 200; ++i) {
            const double t = 7.0 * (i + 0.5) / 200.0;
            const Eigen::Vector2d fd = (pathPoint(p, t + h) - pathPoint(p, t - h)) / (2 * h);
            ASSERT_TRUE((fd - pathVelocity(p, t)).norm() < 1e-4);
        }
    }
}

// A corner IS a step in the reference velocity — that is what a corner means,
// and it is why the ball rounds one.  Worth asserting rather than assuming:
// if the path ever stopped having corners the demo would lose its point.
void test_a_corner_is_a_step_in_the_reference_velocity() {
    const SetpointPath p = shape(PathShape::Square, 0.12, 8.0);
    const double corner_t = p.period_s / 4.0;      // end of the first edge
    const Eigen::Vector2d before = pathVelocity(p, corner_t - 1e-4);
    const Eigen::Vector2d after = pathVelocity(p, corner_t + 1e-4);
    // A square's corner turns the velocity through a right angle.
    ASSERT_NEAR(before.normalized().dot(after.normalized()), 0.0, 1e-6);
    ASSERT_NEAR(before.norm(), after.norm(), 1e-9);
}

// A degenerate period must not divide, and `Fixed` must stay put — it is the
// setpoint the balance loop has always had, and the path layer has to leave it
// alone rather than dragging it to the origin every frame.
void test_the_degenerate_cases_behave() {
    SetpointPath p = shape(PathShape::Circle, 0.12, 0.0);
    ASSERT_TRUE(std::isfinite(pathPoint(p, 3.0).norm()));
    ASSERT_NEAR(pathVelocity(p, 3.0).norm(), 0.0, 1e-15);

    p = shape(PathShape::Fixed);
    ASSERT_NEAR(pathPoint(p, 3.0).norm(), 0.0, 1e-15);
    ASSERT_NEAR(pathLength(p), 0.0, 1e-15);
}

// Drawing: a polygon is drawn by its corners, a circle by samples.  A square
// drawn as a chord of samples is a square drawn wrong.
void test_the_outline_is_drawable() {
    Eigen::Matrix<double, 2, Eigen::Dynamic> pts;
    pathOutline(shape(PathShape::Square), 64, &pts);
    ASSERT_EQ((int)pts.cols(), 5);                     // four corners, closed
    ASSERT_NEAR((pts.col(0) - pts.col(4)).norm(), 0.0, 1e-12);

    pathOutline(shape(PathShape::Triangle), 64, &pts);
    ASSERT_EQ((int)pts.cols(), 4);

    pathOutline(shape(PathShape::Circle), 64, &pts);
    ASSERT_EQ((int)pts.cols(), 65);
    ASSERT_NEAR(pts.col(0).norm(), 0.12, 1e-12);

    pathOutline(shape(PathShape::Fixed), 64, &pts);
    ASSERT_EQ((int)pts.cols(), 0);
}

}  // namespace

int main() {
    test_every_path_closes();
    test_no_path_leaves_its_radius();
    test_every_shape_reaches_its_radius();
    test_the_polygons_are_walked_at_constant_speed();
    test_velocity_is_the_derivative_of_position();
    test_a_corner_is_a_step_in_the_reference_velocity();
    test_the_degenerate_cases_behave();
    test_the_outline_is_drawable();
    std::printf("test_setpoint_path: all passed\n");
    return 0;
}
