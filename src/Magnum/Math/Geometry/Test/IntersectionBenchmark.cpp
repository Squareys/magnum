/*
    This file is part of Magnum.

    Copyright © 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
              Vladimír Vondruš <mosra@centrum.cz>
    Copyright © 2018 Jonathan Hale <squareys@googlemail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.
*/

#include <Corrade/TestSuite/Tester.h>

#include "Magnum/Math/Geometry/Intersection.h"
#include "Magnum/Math/Angle.h"
#include "Magnum/Math/Matrix4.h"

#include <random>
#include <tuple>
#include <utility>

namespace Magnum { namespace Math { namespace Geometry { namespace Benchmark {

/* @brief Ground truth, slow sphere cone intersection - calculating exact distances,
 *        no optimizations, no precomputations
 * @param sphereCenter Sphere center
 * @param radius Sphere radius
 * @param origin Origin of the cone
 * @param normal Cone normal
 * @param angle Cone opening angle (0 < angle < pi).
 *
 * Returns `true` if the sphere intersects with the cone. */
template<class T> bool sphereConeGT(
        const Math::Vector3<T>& sphereCenter, const T radius,
        const Math::Vector3<T>& origin, const Math::Vector3<T>& normal, const Rad<T> angle) {
    const Math::Vector3<T> diff = sphereCenter - origin;
    const Math::Vector3<T> dir = diff.normalized();
    const Rad<T> halfAngle = angle/T(2);

    /* Compute angle between normal and point */
    const Rad<T> actual = Math::acos(dot(normal, dir));

    /* Distance from cone surface */
    const T distanceFromCone = Math::sin(actual - halfAngle)*diff.length();

    /* Either the sphere center lies in cone, or cone is max radius away from the cone */
    return actual <= halfAngle || distanceFromCone <= radius;
}

template<class T>
Math::Matrix4<T> coneViewFromCone(const Math::Vector3<T>& origin, const Math::Vector3<T>& normal) {
    return Math::Matrix4<T>::lookAt(origin, origin + normal, Math::Vector3<T>::yAxis()).inverted();
}

typedef Math::Vector2<Float> Vector2;
typedef Math::Vector3<Float> Vector3;
typedef Math::Vector4<Float> Vector4;
typedef Math::Matrix4<Float> Matrix4;
typedef Math::Frustum<Float> Frustum;

struct IntersectionBenchmark: Corrade::TestSuite::Tester {
    explicit IntersectionBenchmark();

    void sphereConeNaive();
    void sphereCone();
    void sphereConeView();

    std::vector<std::pair<Vector3, Float>> _spheres;
    std::vector<std::tuple<Vector3, Vector3, Rad<Float>>> _cones;
};

IntersectionBenchmark::IntersectionBenchmark() {
    addBenchmarks({&IntersectionBenchmark::sphereConeNaive,
                   &IntersectionBenchmark::sphereCone,
                   &IntersectionBenchmark::sphereConeView}, 1000);

    /* Generate random data for the benchmarks */
    std::random_device rnd;
    std::mt19937 g(rnd());
    /* Position distribution */
    std::uniform_real_distribution<float> pd(-10.0f, 10.0f);
    /* Radius distribution */
    std::uniform_real_distribution<float> rd(0.0f, 5.0f);
    /* Cone angle distribution */
    std::uniform_real_distribution<float> ad(1.0f, 179.0f);

    _cones.reserve(4);
    for(int i = 0; i < 4; ++i) {
        _cones.emplace_back(Vector3{pd(g), pd(g), pd(g)},
                Vector3{pd(g), pd(g), pd(g)}.normalized(),
                Rad<Float>(Deg<Float>(ad(g))));
    }

    _spheres.reserve(100);
    for(int i = 0; i < 100; ++i) {
        _spheres.emplace_back(Vector3{pd(g), pd(g), pd(g)}, rd(g));
    }
}

void IntersectionBenchmark::sphereConeNaive() {
    volatile bool b = false;
    CORRADE_BENCHMARK(10) {
        for(auto& cone : _cones) {
            for(auto& sphere : _spheres) {
                b = b ^ sphereConeGT<Float>(sphere.first, sphere.second, std::get<0>(cone), std::get<1>(cone), std::get<2>(cone));
            }
        }
    }
}

void IntersectionBenchmark::sphereCone() {
    volatile bool b = false;
    CORRADE_BENCHMARK(10) {
        for(auto& cone : _cones) {
            const Float sinAngle = Math::sin(std::get<2>(cone));
            const Float tanAngle = Math::tan(std::get<2>(cone));
            const Float tanAngleSqPlusOne = tanAngle*tanAngle + 1.0f;
            for(auto& sphere : _spheres) {
                b = b ^ Intersection::sphereCone(sphere.first, sphere.second, std::get<0>(cone), std::get<1>(cone), sinAngle, tanAngleSqPlusOne);
            }
        }
    }
}

void IntersectionBenchmark::sphereConeView() {
    volatile bool b = false;
    CORRADE_BENCHMARK(10) {
        for(auto& cone : _cones) {
            const Matrix4 coneView = coneViewFromCone(std::get<0>(cone), std::get<1>(cone));
            for(auto& sphere : _spheres) {
                b = b ^ Intersection::sphereConeView(sphere.first, sphere.second, coneView, std::get<2>(cone));
            }
        }
    }
}

}}}}

CORRADE_TEST_MAIN(Magnum::Math::Geometry::Benchmark::IntersectionBenchmark)
