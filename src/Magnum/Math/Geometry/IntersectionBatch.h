#ifndef Magnum_Math_Geometry_IntersectionBatch_h
#define Magnum_Math_Geometry_IntersectionBatch_h
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

/** @file
 * @brief Namespace @ref Magnum::Math::Geometry::Intersection
 */

#include <bitset>

#include "Corrade/Containers/ArrayView.h"

#include "Magnum/Math/Geometry/Intersection.h"

#include "Magnum/Math/Frustum.h"
#include "Magnum/Math/Range.h"
#include "Magnum/Math/Vector2.h"
#include "Magnum/Math/Vector3.h"
#include "Magnum/Math/Matrix4.h"

namespace Magnum { namespace Math { namespace Geometry { namespace Intersection {

/**
@brief Batch range frustum intersection
@param rangeBatch Input ranges
@param frustum Frustum to intersect with
@param outBits Bitset to receive the result

Size of @cpp outBits @ce needs to be at least @cpp Math::ceil(rangeBatch/64.0f) @ce
big to fit the results of intersection.
*/
template<typename T> bool rangeBatchFrustum(Corrade::Containers::ArrayView<const Range3D<T>> rangeBatch, const Frustum<T>& frustum, Corrade::Containers::ArrayView<UnsignedLong> outBits);

template<typename T> bool rangeBatchFrustum(Corrade::Containers::ArrayView<const Range3D<T>> rangeBatch, const Frustum<T>& frustum, Corrade::Containers::ArrayView<UnsignedLong> outBits) {
    /* Plane with abs normal and precomputed 2*w */
    alignas(16) const Vector4<T> absPlane[]{
         {Math::abs(frustum[0].xyz()), T(2)*frustum[0].w()},
         {Math::abs(frustum[1].xyz()), T(2)*frustum[1].w()},
         {Math::abs(frustum[2].xyz()), T(2)*frustum[2].w()},
         {Math::abs(frustum[3].xyz()), T(2)*frustum[3].w()},
         {Math::abs(frustum[4].xyz()), T(2)*frustum[4].w()},
         {Math::abs(frustum[5].xyz()), T(2)*frustum[5].w()}};

    /* Initialize with all bits on */
    for(UnsignedLong& l : outBits) {
        l = ~0L;
    }

    int longIndex = 0;
    UnsignedLong bitMask = 1L;

    for(auto& range: rangeBatch) {
        const Vector3<T> center = range.min() + range.max();
        const Vector3<T> extent = range.max() - range.min();

        for(int i = 0; i < 6; ++i) {
            const Float d = Math::dot(center, frustum.planes()[i].xyz());
            const Float r = Math::dot(extent, absPlane[i].xyz()) + absPlane[i].w();
            if(d + r < 0) {
                /* Not visible */
                outBits[longIndex] ^= bitMask;
                break;
            }
        }

        bitMask <<= 1;
        /* Bit fully pushed out to the left? */
        if(bitMask == 0L) {
            ++longIndex;
            bitMask = 1;
        }
    }

    return true;
}

}}}}

#endif
