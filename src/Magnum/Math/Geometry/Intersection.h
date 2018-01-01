#ifndef Magnum_Math_Geometry_Intersection_h
#define Magnum_Math_Geometry_Intersection_h
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

/** @file
 * @brief Namespace @ref Magnum::Math::Geometry::Intersection
 */

#include "Magnum/Math/Frustum.h"
#include "Magnum/Math/Geometry/Distance.h"
#include "Magnum/Math/Range.h"
#include "Magnum/Math/Vector2.h"
#include "Magnum/Math/Vector3.h"

namespace Magnum { namespace Math { namespace Geometry { namespace Intersection {

/**
@brief Intersection of two line segments in 2D
@param p        Starting point of first line segment
@param r        Direction of first line segment
@param q        Starting point of second line segment
@param s        Direction of second line segment

Returns intersection point positions @f$ t @f$, @f$ u @f$ on both lines, NaN if
the lines are collinear or infinity if they are parallel. Intersection point
can be then calculated with @f$ \boldsymbol{p} + t \boldsymbol{r} @f$ or
@f$ \boldsymbol{q} + u \boldsymbol{s} @f$. If @f$ t @f$ is in range
@f$ [ 0 ; 1 ] @f$, the intersection is inside the line segment defined by
@f$ \boldsymbol{p} @f$ and @f$ \boldsymbol{p} + \boldsymbol{r} @f$, if @f$ u @f$
is in range @f$ [ 0 ; 1 ] @f$, the intersection is inside the line segment
defined by @f$ \boldsymbol{q} @f$ and @f$ \boldsymbol{q} + \boldsymbol{s} @f$.

The two lines intersect if @f$ t @f$ and @f$ u @f$ exist such that: @f[
     \boldsymbol p + t \boldsymbol r = \boldsymbol q + u \boldsymbol s
@f]
Crossing both sides with @f$ \boldsymbol{s} @f$, distributing the cross product
and eliminating @f$ \boldsymbol s \times \boldsymbol s = 0 @f$, then solving
for @f$ t @f$ and similarly for @f$ u @f$: @f[
     \begin{array}{rcl}
         (\boldsymbol p + t \boldsymbol r) \times s & = & (\boldsymbol q + u \boldsymbol s) \times s \\
         t (\boldsymbol r \times s) & = & (\boldsymbol q - \boldsymbol p) \times s \\
         t & = & \cfrac{(\boldsymbol q - \boldsymbol p) \times s}{\boldsymbol r \times \boldsymbol s} \\
         u & = & \cfrac{(\boldsymbol q - \boldsymbol p) \times r}{\boldsymbol r \times \boldsymbol s}
     \end{array}
@f]

See also @ref lineSegmentLine() which calculates only @f$ t @f$, useful if you
don't need to test that the intersection lies inside line segment defined by
@f$ \boldsymbol{q} @f$ and @f$ \boldsymbol{q} + \boldsymbol{s} @f$.
    */
template<class T> inline std::pair<T, T> lineSegmentLineSegment(const Vector2<T>& p, const Vector2<T>& r, const Vector2<T>& q, const Vector2<T>& s) {
    const Vector2<T> qp = q - p;
    const T rs = cross(r, s);
    return {cross(qp, s)/rs, cross(qp, r)/rs};
}

/**
@brief Intersection of line segment and line in 2D
@param p        Starting point of first line segment
@param r        Direction of first line segment
@param q        Starting point of second line
@param s        Direction of second line

Returns intersection point position @f$ t @f$ on first line, NaN if the lines
are collinear or infinity if they are parallel. Intersection point can be then
calculated with @f$ \boldsymbol{p} + t \boldsymbol{r} @f$. If returned value is
in range @f$ [ 0 ; 1 ] @f$, the intersection is inside the line segment defined
by @f$ \boldsymbol{p} @f$ and @f$ \boldsymbol{p} + \boldsymbol{r} @f$.

Unlike @ref lineSegmentLineSegment() calculates only @f$ t @f$.
*/
template<class T> inline T lineSegmentLine(const Vector2<T>& p, const Vector2<T>& r, const Vector2<T>& q, const Vector2<T>& s) {
    return cross(q - p, s)/cross(r, s);
}

/**
@brief Intersection of a plane and line
@param planePosition    Plane position
@param planeNormal      Plane normal
@param p                Starting point of the line
@param r                Direction of the line

Returns intersection point position @f$ t @f$ on the line, NaN if the line lies
on the plane or infinity if the intersection doesn't exist. Intersection point
can be then calculated from with @f$ \boldsymbol{p} + t \boldsymbol{r} @f$. If
returned value is in range @f$ [ 0 ; 1 ] @f$, the intersection is inside the
line segment defined by @f$ \boldsymbol{p} @f$ and @f$ \boldsymbol{r} @f$.

First the parameter @f$ f @f$ of parametric equation of the plane is calculated
from plane normal @f$ \boldsymbol{n} @f$ and plane position: @f[
     \begin{pmatrix} n_0 \\ n_1 \\ n_2 \end{pmatrix} \cdot
     \begin{pmatrix} x \\ y \\ z \end{pmatrix} - f = 0
@f]
Using plane normal @f$ \boldsymbol{n} @f$, parameter @f$ f @f$ and line defined
by @f$ \boldsymbol{p} @f$ and @f$ \boldsymbol{r} @f$, value of @f$ t @f$ is
calculated and returned. @f[
     \begin{array}{rcl}
         f & = & \boldsymbol n \cdot (\boldsymbol p + t \boldsymbol r) \\
         \Rightarrow t & = & \cfrac{f - \boldsymbol n \cdot \boldsymbol p}{\boldsymbol n \cdot \boldsymbol r}
     \end{array}
@f]
    */
template<class T> inline T planeLine(const Vector3<T>& planePosition, const Vector3<T>& planeNormal, const Vector3<T>& p, const Vector3<T>& r) {
    const T f = dot(planePosition, planeNormal);
    return (f - dot(planeNormal, p))/dot(planeNormal, r);
}

/**
@brief Intersection of a point and a camera frustum
@param point    Point
@param frustum  Frustum planes with normals pointing outwards

Returns `true` if the point is on or inside the frustum.

Checks for each plane of the frustum whether the point is behind the plane (the
points distance from the plane is negative) using @ref Distance::pointPlaneScaled().
*/
template<class T> bool pointFrustum(const Vector3<T>& point, const Frustum<T>& frustum);

/**
@brief Intersection of an axis-aligned box and a camera frustum
@param box      Axis-aligned box
@param frustum  Frustum planes with normals pointing outwards

Returns `true` if the box intersects with the camera frustum.

Counts for each plane of the frustum how many points of the box lie in front of
the plane (outside of the frustum). If none, the box must lie entirely outside
of the frustum and there is no intersection. Else, the box is considered as
intersecting, even if it is merely corners of the box overlapping with corners
of the frustum, since checking the corners is less efficient.
*/
template<class T> bool boxFrustum(const Range3D<T>& box, const Frustum<T>& frustum);

template<class T> bool pointFrustum(const Vector3<T>& point, const Frustum<T>& frustum) {
    for(const Vector4<T>& plane: frustum.planes()) {
        /* The point is in front of one of the frustum planes (normals point
           outwards) */
        if(Distance::pointPlaneScaled<T>(point, plane) < T(0))
            return false;
    }

    return true;
}

template<class T> bool boxFrustum(const Range3D<T>& box, const Frustum<T>& frustum) {
    for(const Vector4<T>& plane: frustum.planes()) {
        bool cornerHit = 0;

        for(UnsignedByte c = 0; c != 8; ++c) {
            const Vector3<T> corner = Math::lerp(box.min(), box.max(), Math::BoolVector<3>{c});

            if(Distance::pointPlaneScaled<T>(corner, plane) >= T(0)) {
                cornerHit = true;
                break;
            }
        }

        /* All corners are outside this plane */
        if(!cornerHit) return false;
    }

    /** @todo potentially check corners here to avoid false positives */

    return true;
}

/**
@brief Intersection of a point and a cone
@param p        The point
@param origin   Origin of the cone
@param normal   Normal of the cone
@param angle    Apex angle of the cone

Returns `true` if the point is inside the cone.

@see pointCone(Vector3, Vector3, Vector3, T)
@todo: That ref is never going to work...
*/
template<class T> bool pointCone(const Vector3<T>& p, const Vector3<T>& origin, const Vector3<T>& normal, const Deg<T> angle) {
    const T x = T(1)+Math::pow<T>(Math::tan<T>(angle/T(2)), T(2));

    return pointCone(p, origin, normal, x);
}

/**
@brief Faster version of intersection of a point and a cone
@param p        The point
@param origin   Origin of the cone
@param normal   Normal of the cone
@param tanAngleSqaredPlusOne Precomputed portion of the cone intersection equation

Returns `true` if the point is inside the cone.

Uses the result of precomputing @f$x = \tan^2{\theta} + 1@f$.

@todo: I believe the normal is expected to be a unit vector...?
*/
template<class T> bool pointCone(const Vector3<T>& p, const Vector3<T>& origin, const Vector3<T>& normal, const T tanAngleSquaredPlusOne) {
    const Vector3<T> c = p - origin;
    const T lenA = dot(c, normal);

    return c.dot() <= lenA*lenA*tanAngleSquaredPlusOne;
}

/**
@brief Line circle intersection
@param from   First point of the line
@param to     Second point of the line
@param center Center of the circle
@param radius Radius of the circle

Returns `true` if the line intersects the circle.
*/
template<class T> bool twoPointLineCircle(const Vector2<T>& from, const Vector2<T>& to, const Vector2<T>& origin, const T radius) {
    /* Tranform to circle-local coords */
    const Vector2<T> f = from - origin;
    const Vector2<T> t = to - origin;

    return twoPointLineOriginCircle(f, t, radius);
}

template<class T> bool twoPointLineOriginCircle(const Vector2<T>& from, const Vector2<T>& to, const T radius) {
    const Vector2<T> dir = to - from;

    return lineOriginCircle(from, dir, radius);
}

template<class T> bool lineOriginCircle(const Vector2<T>& from, const Vector2<T>& dir, const T radius) {
    return lineOriginCircleRadiusSquared(from, dir, radius*radius);
}

template<class T> bool lineOriginCircleRadiusSquared(const Vector2<T>& from, const Vector2<T>& dir, const T radiusSq) {
    const T a = dir.dot();
    const T b = dot(dir, from);
    const T c = from.dot() - radiusSq;
    const T delta = b*b - T(4)*a*c;

    return delta >= T(0);
}

template<class T> bool lineSphere(const Vector3<T>& origin, const Vector3<T>& dir, const Vector3<T>& center, const T radiusSq) {
    const Vector3<T> diff = origin - center;
    const T x = dot(dir, diff);
    const T = x*x - dir.dot()*(diff.dot()-radiusSq);

    return delta >= T(0);
}

template<class T> bool lineCone(const Vector3<T>& from, const Vector3<T>& dir, const Vector3<T>& origin, const Vector3<T>& normal, const T cosAngleSq) {

    /* Calculate intersection line of the two planes */
    const Vector3<T> o = from; // TODO
    const Vector3<T> p1Normal = normal; // TODO
    const Vector3<T> d = cross(dir, origin); // TODO

    /* Calculate intersection of the resulting line with sphere */
    return false;
}

/**
 * Based on https://www.geometrictools.com/Documentation/IntersectionTriangleCone.pdf
 */
template<class T> bool triangleCone(const Vector3<T>& p0, const Vector3<T>& p1, const Vector3<T>& p2, const Vector3<T>& origin, const Vector3<T>& normal, const T cosAngleSq) {
    bool inFront[4]{false, false, false, false};
    Vector3<T> points[4]{p0, p1, p2, p0}; //constexpr?

    for (int i = 0; i < 3; ++i) {
        const Vector3<T> diff = points[i] - origin;
        const T d = dot(normal, diff);
        inFront[i] = d < T(0);
        if(inFront[i]) {
            if(d*d <= cosAngleSq*diff.dot()) {
                return true;
            }
        } else {
            /* behind the cone */
        }
        ++i;
    }

    if(!inFront[0] && !inFront[1] && !inFront[2]) {
        return false;
    }

    inFront[3] = inFront[1];
    /* If any edge intersects, the triangle intersects, therefore test all of them */
    for(int i = 0; i < 3; ++i) {
        if(!inFront[i] && !inFront[i+1]) {
            /* does not intesect */
        } else {
            /* handle edges fully on the cone side */
            const Vector3<T> dir = points[i+1] - points[i];
            const T d = dot(normal, dir);

            const T c2 = d*d - dir.dot()*cosAngleSq;
            if(c2 < T(0)) {
                const Vector3<T> o = points[i] - origin;
                const T dirDotO = dot(dir, o);
                const T normDotO = dot(normal, o);
                const T c1 = d*normDotO - cosAngleSq*dirDotO;
                if(inFront[i] && inFront[i+1]) {
                    if(T(0) <= c1 && c1 <= c2) {
                        const T c0 = dirDotO*dirDotO - cosAngleSq*o.dot();
                        if(c1*c1 >= c0*c2) {
                            return true;
                        }
                    }
                } else if(inFront[i] && !inFront[i+1]) {
                    if(T(0) <= c1 && c2*normDotO <= c1*dirDotO) {
                        const T c0 = normDotO*normDotO - cosAngleSq*o.dot();
                        if(c1*c1 >= c0*c2) {
                            return true;
                        }
                    }
                } else if(!inFront[i] && inFront[i+1]) {
                    /* <--------> only changed condition below from entire block above */
                    if(c1 <= -c2 && c2*normDotO <= c1*dirDotO) {
                        const T c0 = normDotO*normDotO - cosAngleSq*o.dot();
                        if(c1*c1 >= c0*c2) {
                            return true;
                        }
                    }
                }
            }
        }
    }

    /* plane test */
    const Vector3<T> edge0 = points[1] - points[0];
    const Vector3<T> edge1 = points[2] - points[1];
    const Vector3<T> edge2 = points[0] - points[2];

    const Vector3<T> triangleNormal = cross(edge0, edge1);
    const T dotTriangleConeNormal = dot(triangleNormal, normal);

    const Vector3<T> delta0 = points[0] - origin;
    const T dotTriangleNormalDelta0 = dot(triangleNormal, delta0);

    const Vector3<T> u = dotTriangleNormalDelta0*normal - dotTriangleConeNormal*delta0;
    const Vector3<T> nCrossU = cross(triangleNormal, u);

    if(dotTriangleConeNormal >= T(0)) {
        if(dot(nCrossU, edge0) <= T(0) && dot(nCrossU, edge1) <= T(0)) {
            return dot(nCrossU, edge2) <= dotTriangleConeNormal*triangleNormal.dot();
        }
    } else {
        if(dot(nCrossU, edge0) >= T(0) && dot(nCrossU, edge1) >= T(0)) {
            return dot(nCrossU, edge2) >= dotTriangleConeNormal*triangleNormal.dot();
        }
    }

    return false;
}

}}}}

#endif
