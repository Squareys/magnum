#ifndef Magnum_Math_Geometry_Intersection_h
#define Magnum_Math_Geometry_Intersection_h
/*
    This file is part of Magnum.

    Copyright © 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

Returns @cpp true @ce if the point is on or inside the frustum.

Checks for each plane of the frustum whether the point is behind the plane (the
points distance from the plane is negative) using @ref Distance::pointPlaneScaled().
*/
template<class T> bool pointFrustum(const Vector3<T>& point, const Frustum<T>& frustum);

/**
@brief Intersection of an axis-aligned box and a camera frustum
@param box      Axis-aligned box
@param frustum  Frustum planes with normals pointing outwards

Returns @cpp true @ce if the box intersects with the camera frustum.

Counts for each plane of the frustum how many points of the box lie in front of
the plane (outside of the frustum). If none, the box must lie entirely outside
of the frustum and there is no intersection. Else, the box is considered as
intersecting, even if it is merely corners of the box overlapping with corners
of the frustum, since checking the corners is less efficient.
*/
template<class T> bool boxFrustum(const Range3D<T>& box, const Frustum<T>& frustum);

/**
@brief Intersection of a sphere and a camera frustum
@param center   Sphere center
@param radius   Sphere radius
@param frustum  Frustum planes with normals pointing outwards

Returns @cpp true @ce if the sphere intersects the frustum.

Checks for each plane of the frustum whether the sphere is behind the plane (the
points distance larger than the sphere's radius) using @ref Distance::pointPlaneScaled().
*/
template<class T> bool sphereFrustum(const Vector3<T>& center, T radius, const Frustum<T>& frustum);

/**
@brief Intersection of a point and a cone
@param p        The point
@param origin   Origin of the cone
@param normal   Normal of the cone
@param angle    Apex angle of the cone

Returns @cpp true @ce if the point is inside the cone.

@see pointCone(Vector3, Vector3, Vector3, T) @todo: That ref is never going to work...
*/
template<class T> bool pointCone(const Vector3<T>& p, const Vector3<T>& origin, const Vector3<T>& normal, Deg<T> angle);

/**
@brief Faster version of intersection of a point and a cone
@param p        The point
@param origin   Origin of the cone
@param normal   Normal of the cone
@param tanAngleSqPlusOne Precomputed portion of the cone intersection equation: `Math::tan(angle/2)+1`

Returns @cpp true @ce if the point is inside the cone.

Uses the result of precomputing @f$x = \tan^2{\theta} + 1@f$.
*/
template<class T> bool pointCone(const Vector3<T>& p, const Vector3<T>& origin, const Vector3<T>& normal, T tanAngleSqPlusOne);

/**
@brief Faster version of intersection of a point and a double cone
@param p        The point
@param origin   Origin of the cone
@param normal   Normal of the cone
@param tanAngleSqPlusOne Precomputed portion of the cone intersection equation:
                         @cpp Math::tan(angle/2)+1 @ce

Returns @cpp true @ce if the point is inside the double cone.

Uses the result of precomputing @f$ x = \tan^2{\theta} + 1 @f$.
*/
template<class T> bool pointDoubleCone(const Vector3<T>& p, const Vector3<T>& origin, const Vector3<T>& normal, T tanAngleSqPlusOne);

/**
@brief Intersection of a point and a double cone
@param p        The point
@param origin   Origin of the cone
@param normal   Normal of the cone
@param angle    Apex angle of the cone

Returns @cpp true @ce if the point is inside the double cone.

@see pointDoubleCone(Vector3, Vector3, Vector3, T) @todo: That ref is never going to work...
*/
template<class T> bool pointDoubleCone(const Vector3<T>& p, const Vector3<T>& origin, const Vector3<T>& normal, Deg<T> angle);

/**
@brief Intersection of a sphere and a cone view
@param sphereCenter Center of the sphere
@param radius Radius of the sphere
@param coneView View matrix with translation and rotation of the cone
@param angle Cone opening angle

Returns @cpp true @ce if the sphere intersects the cone.
*/
template<class T> bool sphereConeView(const Vector3<T>& sphereCenter, T radius, const Matrix4<T>& coneView, Rad<T> angle);

template<class T> bool sphereConeView(const Vector3<T>& sphereCenter, T radius, const Matrix4<T>& coneView, T sinAngle, T cosAngle, T tanAngle);

/**
@brief Intersection of a sphere and a cone
@param sphereCenter Sphere center
@param radius Sphere radius
@param origin Origin of the cone
@param normal Cone normal
@param angle Cone opening angle (0 < angle < pi).

Offsets the cone plane by @f$ -r\sin{\theta} \cdot \boldsymbol n @f$ (with  @f$ \theta @f$
the cone's half-angle) which separates two half-spaces:
In front of the plane, in which the sphere cone intersection test is equivalent
to testing the sphere's center against a similarly offset cone (which is equivalent
the cone with surface expanded by @f$ r @f$ in surface normal direction), and
begind the plane, where the test is equivalent to testing whether the origin of
the original cone intersects the sphere.

Returns @cpp true @ce if the sphere intersects with the cone.
*/
template<class T> bool sphereCone(const Vector3<T>& sphereCenter, const T radius, const Vector3<T>& origin, const Vector3<T>& normal, const Rad<T> angle);

template<class T> bool sphereCone(const Vector3<T>& sphereCenter, const T radius, const Vector3<T>& origin, const Vector3<T>& normal, const T sinAngle, const T tanAngleSquaredPlusOne);

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

template<class T> bool sphereFrustum(const Vector3<T>& center, const T radius, const Frustum<T>& frustum) {
    const T radiusSq = radius*radius;

    for(const Vector4<T>& plane: frustum.planes()) {
        /* The sphere is in front of one of the frustum planes (normals point
           outwards) */
        if(Distance::pointPlaneScaled<T>(center, plane) < -radiusSq)
            return false;
    }

    return true;
}


template<class T> bool pointCone(const Vector3<T>& p, const Vector3<T>& origin, const Vector3<T>& normal, const Deg<T> angle) {
    const T x = T(1)+Math::pow<T>(Math::tan<T>(angle/T(2)), T(2));

    return pointCone(p, origin, normal, x);
}

template<class T> bool pointCone(const Vector3<T>& p, const Vector3<T>& origin, const Vector3<T>& normal, const T tanAngleSqPlusOne) {
    const Vector3<T> c = p - origin;
    const T lenA = dot(c, normal);

    return lenA >= 0 && c.dot() <= lenA*lenA*tanAngleSqPlusOne;
}

template<class T> bool pointDoubleCone(const Vector3<T>& p, const Vector3<T>& origin, const Vector3<T>& normal, const Deg<T> angle) {
    const T x = T(1)+Math::pow<T>(Math::tan<T>(angle/T(2)), T(2));

    return pointDoubleCone(p, origin, normal, x);
}

template<class T> bool pointDoubleCone(const Vector3<T>& p, const Vector3<T>& origin, const Vector3<T>& normal, const T tanAngleSqPlusOne) {
    const Vector3<T> c = p - origin;
    const T lenA = dot(c, normal);

    return c.dot() <= lenA*lenA*tanAngleSqPlusOne;
}

template<class T> bool sphereConeView(const Vector3<T>& sphereCenter, const T radius, const Matrix4<T>& coneView, const Rad<T> angle) {
    const Float sinAngle = Math::sin(angle/T(2)); /* precomputable */
    const Float cosAngle = Math::cos(angle/T(2)); /* precomputable */
    const Float tanAngle = Math::tan<T>(angle/T(2));
    const Float tanAngleSquaredPlusOne = T(1)+tanAngle*tanAngle;

    return sphereConeView(sphereCenter, radius, coneView, sinAngle, cosAngle, tanAngle);
}

template<class T> bool sphereConeView(const Vector3<T>& sphereCenter, const T radius, const Matrix4<T>& coneView, const T sinAngle, const T cosAngle, const T tanAngle) {
    /* Axis align cone */
    const Vector3<T> center = coneView.transformPoint(sphereCenter);

    /* Test against plane which determins whether to test against shiften cone or center-sphere */
    if (center.z() > -radius*cosAngle) {
        /* Axis aligned point - cone test shifted so that the cones surface is extended by `radius` */
        const T coneRadius = tanAngle*(center.z() + radius/sinAngle);
        return center.xy().dot() <= coneRadius*coneRadius;
    }

    return false;
}

template<class T> bool sphereCone(
        const Vector3<T>& sCenter, const T radius,
        const Vector3<T>& origin, const Vector3<T>& normal, const Rad<T> angle) {

    const T sinAngle = Math::sin(angle / T(2)); /* precomputable */
    const T tanAngleSquaredPlusOne = T(1)+Math::pow<T>(Math::tan<T>(angle/T(2)), T(2));

    return sphereCone(sCenter, radius, origin, normal, sinAngle, tanAngleSquaredPlusOne);
}

template<class T> bool sphereCone(
        const Vector3<T>& sCenter, const T radius,
        const Vector3<T>& origin, const Vector3<T>& normal, const T sinAngle, const T tanAngleSquaredPlusOne) {

    const Vector3<T> diff = sCenter - origin;

    if (dot(diff, normal) > T(0)) {
        /* point - cone test */
        const Vector3<T> c = sinAngle*diff + (normal*radius);
        const T lenA = dot(c, normal);

        return c.dot() <= lenA*(lenA*tanAngleSquaredPlusOne);
    } else {
        /* Simple sphere plane check */
        return diff.dot() <= radius*radius;
    }
}

}}}}

#endif
