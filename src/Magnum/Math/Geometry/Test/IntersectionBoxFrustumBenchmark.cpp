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

#include "Magnum/Math/Geometry/Intersection.h"
#include "Magnum/Math/Angle.h"
#include "Magnum/Math/Matrix4.h"

#include <random>
#include <tuple>
#include <utility>

namespace Magnum { namespace Math { namespace Geometry { namespace Benchmark {

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


typedef Math::Vector3<Float> Vector3;
typedef Math::Frustum<Float> Frustum;
typedef Math::Range3D<Float> Range3D;

struct IntersectionBoxFrustumBenchmark: Corrade::TestSuite::Tester {
    explicit IntersectionBoxFrustumBenchmark();

    void boxFrustumNaive();
    void boxFrustum();

    std::vector<Range3D> _boxes;
    const Frustum _frustum{
        {1.0f, 0.0f, 0.0f, 0.0f},
        {-1.0f, 0.0f, 0.0f, 10.0f},
        {0.0f, 1.0f, 0.0f, 0.0f},
        {0.0f, -1.0f, 0.0f, 10.0f},
        {0.0f, 0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, -1.0f, 10.0f}};
};

IntersectionBoxFrustumBenchmark::IntersectionBoxFrustumBenchmark() {
    addBenchmarks({&IntersectionBoxFrustumBenchmark::boxFrustumNaive,
                   &IntersectionBoxFrustumBenchmark::boxFrustum}, 1000);

    /* Generate random data for the benchmarks */
    std::random_device rnd;
    std::mt19937 g(rnd());
    /* Position distribution */
    std::uniform_real_distribution<float> pd(-10.0f, 10.0f);
    /* Radius distribution */
    std::uniform_real_distribution<float> rd(0.0f, 5.0f);
    /* Cone angle distribution */
    std::uniform_real_distribution<float> ad(1.0f, 179.0f);

    _boxes.reserve(512);
    for(int i = 0; i < 512; ++i) {
        _boxes.emplace_back(Vector3{pd(g), pd(g), pd(g)},
                            Vector3{pd(g), pd(g), pd(g)});
    }
}

void IntersectionBoxFrustumBenchmark::boxFrustumNaive() {
    volatile bool b = false;
    CORRADE_BENCHMARK(10) {
        for(auto& box : _boxes) {
            b = b ^ Benchmark::boxFrustumNaive<Float>(box, _frustum);
        }
    }
}

void IntersectionBoxFrustumBenchmark::boxFrustum() {
    volatile bool b = false;
    CORRADE_BENCHMARK(10) {
        for(auto& box : _boxes) {
            b = b ^ Intersection::boxFrustum(box, _frustum);
        }
    }
}

}}}}

CORRADE_TEST_MAIN(Magnum::Math::Geometry::Benchmark::IntersectionBoxFrustumBenchmark)
