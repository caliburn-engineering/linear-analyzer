// tests/test_helpers.h
#pragma once

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>

#define ASSERT_NEAR(actual, expected, tol)                                     \
    do {                                                                       \
        double a_ = (actual), e_ = (expected), t_ = (tol);                    \
        if (std::abs(a_ - e_) > t_) {                                         \
            std::fprintf(stderr, "FAIL: %s:%d: %s = %g, expected %g (tol %g)\n", \
                         __FILE__, __LINE__, #actual, a_, e_, t_);             \
            std::exit(1);                                                      \
        }                                                                      \
    } while (0)

#define ASSERT_TRUE(cond)                                                      \
    do {                                                                       \
        if (!(cond)) {                                                         \
            std::fprintf(stderr, "FAIL: %s:%d: %s\n", __FILE__, __LINE__,      \
                         #cond);                                               \
            std::exit(1);                                                      \
        }                                                                      \
    } while (0)

#define ASSERT_EQ(actual, expected)                                            \
    do {                                                                       \
        auto a_ = (actual);                                                    \
        auto e_ = (expected);                                                  \
        if (a_ != e_) {                                                        \
            std::fprintf(stderr, "FAIL: %s:%d: %s = %d, expected %d\n",        \
                         __FILE__, __LINE__, #actual, (int)a_, (int)e_);        \
            std::exit(1);                                                      \
        }                                                                      \
    } while (0)
