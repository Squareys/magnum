/*
    This file is part of Magnum.

    Copyright © 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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
#include <Corrade/Containers/Array.h>

#include <random>
#include <tuple>
#include <utility>

#include "Magnum/Math/Geometry/Intersection.h"
#include "Magnum/Math/Geometry/IntersectionBatch.h"
#include "Magnum/Math/Angle.h"
#include "Magnum/Math/Matrix4.h"

namespace Magnum { namespace Math { namespace Geometry { namespace Test {

template<class T> bool boxFrustumNaive(const Math::Range3D<T>& box, const Math::Frustum<T>& frustum) {
    for(const Math::Vector4<T>& plane: frustum.planes()) {
        bool cornerHit = 0;

        for(UnsignedByte c = 0; c != 8; ++c) {
            const Math::Vector3<T> corner = Math::lerp(box.min(), box.max(), Math::BoolVector<3>{c});

            if(Distance::pointPlaneScaled<T>(corner, plane) >= T(0)) {
                cornerHit = true;
                break;
            }
        }

        /* All corners are outside this plane */
        if(!cornerHit) return false;
    }

    return true;
}

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
        const Math::Vector3<T>& origin, const Math::Vector3<T>& normal, const Math::Rad<T> angle) {
    const Math::Vector3<T> diff = sphereCenter - origin;
    const Math::Vector3<T> dir = diff.normalized();
    const Math::Rad<T> halfAngle = angle/T(2);

    /* Compute angle between normal and point */
    const Math::Rad<T> actual = Math::acos(dot(normal, dir));

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
typedef Math::Range3D<Float> Range3D;
typedef Math::Deg<Float> Deg;
typedef Math::Rad<Float> Rad;

struct IntersectionBenchmark: Corrade::TestSuite::Tester {
    explicit IntersectionBenchmark();

    void boxFrustumNaive();
    void boxFrustum();
    void rangeFrustumBatch();

    void sphereFrustum();
    void sphereConeNaive();
    void sphereCone();
    void sphereConeView();

    const Frustum _frustum{
        {1.0f, 0.0f, 0.0f, 0.0f},
        {-1.0f, 0.0f, 0.0f, 10.0f},
        {0.0f, 1.0f, 0.0f, 0.0f},
        {0.0f, -1.0f, 0.0f, 10.0f},
        {0.0f, 0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, -1.0f, 10.0f}};
    std::tuple<Vector3, Vector3, Rad> _cone;

    std::vector<Range3D> _boxes;
    std::vector<Vector4> _spheres;
};

IntersectionBenchmark::IntersectionBenchmark() {
    addBenchmarks({&IntersectionBenchmark::boxFrustumNaive,
                   &IntersectionBenchmark::boxFrustum,
                   &IntersectionBenchmark::rangeFrustumBatch,

                   &IntersectionBenchmark::sphereFrustum,
                   &IntersectionBenchmark::sphereConeNaive,
                   &IntersectionBenchmark::sphereCone,
                   &IntersectionBenchmark::sphereConeView}, 25);

    /* Generate random data for the benchmarks */
    std::random_device rnd;
    std::mt19937 g(rnd());
    /* Position distribution */
    std::uniform_real_distribution<float> pd(-10.0f, 10.0f);
    /* Radius distribution */
    std::uniform_real_distribution<float> rd(0.0f, 5.0f);
    /* Cone angle distribution */
    std::uniform_real_distribution<float> ad(1.0f, 179.0f);

    _cone = std::make_tuple(Vector3{pd(g), pd(g), pd(g)},
                Vector3{pd(g), pd(g), pd(g)}.normalized(),
                Rad(Deg(ad(g))));

    _boxes.reserve(512);
    for(int i = 0; i < 512; ++i) {
        _boxes.emplace_back(Vector3{pd(g), pd(g), pd(g)},
                            Vector3{pd(g), pd(g), pd(g)});
    }

    _spheres.reserve(512);
    for(int i = 0; i < 512; ++i) {
        _spheres.emplace_back(Vector3{pd(g), pd(g), pd(g)}, rd(g));
    }
}

void IntersectionBenchmark::boxFrustumNaive() {
    volatile bool b = false;
    CORRADE_BENCHMARK(50) for(auto& box: _boxes) {
        b = b ^ Test::boxFrustumNaive<Float>(box, _frustum);
    }
}

void IntersectionBenchmark::boxFrustum() {
    volatile bool b = false;
    CORRADE_BENCHMARK(50) for(auto& box: _boxes) {
        b = b ^ Intersection::boxFrustum(box, _frustum);
    }
}

void IntersectionBenchmark::rangeFrustumBatch() {
    Corrade::Containers::Array<UnsignedLong> results{size_t(Math::ceil(_boxes.size()/64.0f))};
    CORRADE_BENCHMARK(50) {
        Intersection::rangeBatchFrustum(Corrade::Containers::ArrayView<const Range3D>(_boxes.data(), _boxes.size()), _frustum, results);
    }

    volatile long r = 0;
    for(long l : results) {
        r ^= l;
    }
}

void IntersectionBenchmark::sphereFrustum() {
    volatile bool b = false;
    CORRADE_BENCHMARK(50) for(auto& sphere: _spheres) {
        b = b ^ Intersection::sphereFrustum(sphere.xyz(), sphere.w(), _frustum);
    }
}

void IntersectionBenchmark::sphereConeNaive() {
    volatile bool b = false;
    CORRADE_BENCHMARK(50) for(auto& sphere: _spheres) {
        b = b ^ sphereConeGT<Float>(sphere.xyz(), sphere.w(), std::get<0>(_cone), std::get<1>(_cone), std::get<2>(_cone));
    }
}

void IntersectionBenchmark::sphereCone() {
    volatile bool b = false;
    CORRADE_BENCHMARK(50) {
        const Float sinAngle = Math::sin(std::get<2>(_cone));
        const Float tanAngle = Math::tan(std::get<2>(_cone));
        const Float tanAngleSqPlusOne = tanAngle*tanAngle + 1.0f;
        for(auto& sphere: _spheres) {
            b = b ^ Intersection::sphereCone(sphere.xyz(), sphere.w(), std::get<0>(_cone), std::get<1>(_cone), sinAngle, tanAngleSqPlusOne);
        }
    }
}

void IntersectionBenchmark::sphereConeView() {
    volatile bool b = false;
    CORRADE_BENCHMARK(50) {
        const Matrix4 coneView = coneViewFromCone(std::get<0>(_cone), std::get<1>(_cone));
        for(auto& sphere: _spheres) {
            b = b ^ Intersection::sphereConeView(sphere.xyz(), sphere.w(), coneView, std::get<2>(_cone));
        }
    }
}

}}}}

CORRADE_TEST_MAIN(Magnum::Math::Geometry::Test::IntersectionBenchmark)
