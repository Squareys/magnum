/*
    This file is part of Magnum.

    Copyright © 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
              Vladimír Vondruš <mosra@centrum.cz>
    Copyright © 2016 Jonathan Hale <squareys@googlemail.com>

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

namespace Magnum { namespace Math { namespace Geometry { namespace Test {

struct IntersectionTest: Corrade::TestSuite::Tester {
    explicit IntersectionTest();

    void planeLine();
    void lineLine();

    void pointFrustum();
    void boxFrustum();

    void pointCone();

    void lineCircle();
};

typedef Math::Vector2<Float> Vector2;
typedef Math::Vector3<Float> Vector3;
typedef Math::Vector4<Float> Vector4;
typedef Math::Frustum<Float> Frustum;
typedef Math::Constants<Float> Constants;
typedef Math::Range3D<Float> Range3D;

IntersectionTest::IntersectionTest() {
    addTests({&IntersectionTest::planeLine,
              &IntersectionTest::lineLine,

              &IntersectionTest::pointFrustum,
              &IntersectionTest::boxFrustum,

              &IntersectionTest::pointCone,

              &IntersectionTest::lineCircle});
}

void IntersectionTest::planeLine() {
    const Vector3 planePosition(-1.0f, 1.0f, 0.5f);
    const Vector3 planeNormal(0.0f, 0.0f, 1.0f);

    /* Inside line segment */
    CORRADE_COMPARE(Intersection::planeLine(planePosition, planeNormal,
        {0.0f, 0.0f, -1.0f}, {0.0f, 0.0f, 2.0f}), 0.75f);

    /* Outside line segment */
    CORRADE_COMPARE(Intersection::planeLine(planePosition, planeNormal,
        {0.0f, 0.0f, 1.0f}, {0.0f, 0.0f, 1.0f}), -0.5f);

    /* Line lies on the plane */
    CORRADE_COMPARE(Intersection::planeLine(planePosition, planeNormal,
        {1.0f, 0.5f, 0.5f}, {-1.0f, 0.5f, 0.0f}), Constants::nan());

    /* Line is parallel to the plane */
    CORRADE_COMPARE(Intersection::planeLine(planePosition, planeNormal,
        {1.0f, 0.0f, 1.0f}, {-1.0f, 0.0f, 0.0f}), -Constants::inf());
}

void IntersectionTest::lineLine() {
    const Vector2 p(-1.0f, -1.0f);
    const Vector2 r(1.0, 2.0f);

    /* Inside both line segments */
    CORRADE_COMPARE(Intersection::lineSegmentLineSegment(p, r,
        {0.0f, 0.0f}, {-1.0f, 0.0f}), std::make_pair(0.5f, 0.5f));
    CORRADE_COMPARE(Intersection::lineSegmentLine(p, r,
        {0.0f, 0.0f}, {-1.0f, 0.0f}), 0.5);

    /* Outside both line segments */
    CORRADE_COMPARE(Intersection::lineSegmentLineSegment(p, r,
        {0.0f, -2.0f}, {-1.0f, 0.0f}), std::make_pair(-0.5f, 1.5f));
    CORRADE_COMPARE(Intersection::lineSegmentLine(p, r,
        {0.0f, -2.0f}, {-1.0f, 0.0f}), -0.5f);

    /* Collinear lines */
    const auto tu = Intersection::lineSegmentLineSegment(p, r,
        {0.0f, 1.0f}, {-1.0f, -2.0f});
    CORRADE_COMPARE(tu.first, -Constants::nan());
    CORRADE_COMPARE(tu.second, -Constants::nan());
    CORRADE_COMPARE(Intersection::lineSegmentLine(p, r,
        {0.0f, 1.0f}, {-1.0f, -2.0f}), -Constants::nan());

    /* Parallel lines */
    CORRADE_COMPARE(Intersection::lineSegmentLineSegment(p, r,
        {0.0f, 0.0f}, {1.0f, 2.0f}), std::make_pair(Constants::inf(), Constants::inf()));
    CORRADE_COMPARE(Intersection::lineSegmentLine(p, r,
        {0.0f, 0.0f}, {1.0f, 2.0f}), Constants::inf());
}

void IntersectionTest::pointFrustum() {
    const Frustum frustum{
        {1.0f, 0.0f, 0.0f, 0.0f},
        {-1.0f, 0.0f, 0.0f, 10.0f},
        {0.0f, 1.0f, 0.0f, 0.0f},
        {0.0f, -1.0f, 0.0f, 10.0f},
        {0.0f, 0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, -1.0f, 10.0f}};

    /* Point on edge */
    CORRADE_VERIFY(Intersection::pointFrustum({}, frustum));
    /* Point inside */
    CORRADE_VERIFY(Intersection::pointFrustum({5.0f, 5.0f, 5.0f}, frustum));
    /* Point outside */
    CORRADE_VERIFY(!Intersection::pointFrustum({0.0f, 0.0f, 100.0f}, frustum));
}

void IntersectionTest::boxFrustum() {
    const Frustum frustum{
        {1.0f, 0.0f, 0.0f, 0.0f},
        {-1.0f, 0.0f, 0.0f, 10.0f},
        {0.0f, 1.0f, 0.0f, 0.0f},
        {0.0f, -1.0f, 0.0f, 10.0f},
        {0.0f, 0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, -1.0f, 10.0f}};

    CORRADE_VERIFY(Intersection::boxFrustum({Vector3{1.0f}, Vector3{2.0f}}, frustum));
    /* Bigger than frustum, but still intersects */
    CORRADE_VERIFY(Intersection::boxFrustum(Range3D{Vector3{-100.0f}, Vector3{100.0f}}, frustum));
    /* Outside of frustum */
    CORRADE_VERIFY(!Intersection::boxFrustum(Range3D{Vector3{-10.0f}, Vector3{-5.0f}}, frustum));
}

void IntersectionTest::pointCone() {
    const Vector3 center{0.1f, 0.2f, 0.3f};
    const Vector3 normal{0.0f, 1.0f, 0.0f};
    const Deg<Float> angle{72.0f};

    /* Point on edge */
    CORRADE_VERIFY(Intersection::pointCone(center, center, normal, angle));
    /* Point inside */
    CORRADE_VERIFY(Intersection::pointCone(center + Vector3{0.5f, 1.0f, 0.3f}, center, normal, angle));
    /* Point outside */
    CORRADE_VERIFY(!Intersection::pointCone({3.0f, -10.0f, 100.0f}, center, normal, angle));
    CORRADE_VERIFY(!Intersection::pointCone({0.0f, 0.0f, 0.0f}, center, normal, angle));
}

void IntersectionTest::lineCircle() {
    const Vector2 center{1.0f, 2.0f};
    const Float radius = 2.0f;

    /* One intersection */
    CORRADE_VERIFY(Intersection::twoPointLineCircle(Vector2{3.0f, 2.0f}, Vector2{3.0, 3.0f}, center, radius));
    /* Two intersections */
    CORRADE_VERIFY(Intersection::twoPointLineCircle(Vector2{1.0f, 2.0f}, Vector2{3.0, 4.0f}, center, radius));
    /* No intersection */
    CORRADE_VERIFY(!Intersection::twoPointLineCircle(Vector2{10.0f, 20.0f}, Vector2{30.0, 40.0f}, center, radius));
}

}}}}

CORRADE_TEST_MAIN(Magnum::Math::Geometry::Test::IntersectionTest)
