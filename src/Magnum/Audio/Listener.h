#ifndef Magnum_Audio_Listener_h
#define Magnum_Audio_Listener_h
/*
    This file is part of Magnum.

    Copyright © 2010, 2011, 2012, 2013, 2014, 2015
              Vladimír Vondruš <mosra@centrum.cz>
    Copyright © 2015 Jonathan Hale <squareys@googlemail.com>

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
 * @brief Class @ref Magnum::Audio::Listener
 */

#include <string>
#include <al.h>

#include <Magnum/Math/Matrix3.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/Quaternion.h>
#include <Magnum/SceneGraph/AbstractFeature.h>
#include <Magnum/SceneGraph/AbstractObject.h>

#include "Magnum/Audio/Renderer.h"

namespace Magnum { namespace Audio {

/**
@brief Listener feature for two dimensional float scenes
 */
class Listener2D: public SceneGraph::AbstractFeature2D {
    public:

        /** @brief Constructor */
        explicit Listener2D(SceneGraph::AbstractObject2D& object):
            SceneGraph::AbstractFeature2D{object}
        {
            setCachedTransformations(SceneGraph::CachedTransformation::Absolute);
        }

        /**
         * @brief setNormalizedRotationAxis
         * @param normalizedAxis Axis to rotate around, expected to be normalized.
         * @return Reference to self (for method chaining)
         */
        Listener2D& setNormalizedRotationAxis(Vector3& normalizedAxis) {
            _rotationAxis = normalizedAxis;
            return *this;
        }

    private:
        void clean(const Matrix3& absoluteTransformationMatrix) override {
            Renderer::setListenerPosition(Vector3::pad(absoluteTransformationMatrix.translation()));

            const Matrix2x2 rotation2D = absoluteTransformationMatrix.rotation();
            Rad rad{Float(atan2(rotation2D[1][2], rotation2D[0][0]))};

            Matrix4 rotation3D = Matrix4::rotation(rad, _rotationAxis);

            constexpr Vector3 fwd{0.0, 0.0, 1.0};
            constexpr Vector3 up{0.0, 1.0, 0.0};
            rotation3D.transformVector(fwd);
            rotation3D.transformVector(up);
            Renderer::setListenerOrientation(fwd, up);

            // TODO: velocity
        }

        Vector3 _rotationAxis;
};

/**
@brief Listener feature for three dimensional float scenes
*/
class Listener3D: public SceneGraph::AbstractFeature3D {
   public:

       /** @brief Constructor */
       explicit Listener3D(SceneGraph::AbstractObject3D& object):
           SceneGraph::AbstractFeature3D{object},
           _fwd{0.0, 0.0, 1.0},
           _up{0.0, 1.0, 0.0}
       {
           setCachedTransformations(SceneGraph::CachedTransformation::Absolute);
       }

       /**
        * @brief Set the vectors to consider as forward and up
        * @param fwd Forward vector
        *
        * Default `{0.0, 1.0}` (fwd) for two dimensional and `{0.0, 0.0, 1.0}` for
        * three dimensional feature.
        * @return Reference to self (for method chaining)
        */
       Listener3D& setForwardAndUp(Vector3& forward, Vector3& up) {
           CORRADE_ASSERT(Math::dot(forward, up) == 0, "Forward and up vectors need to be perpendicular.", *this);
           _fwd = forward;
           _up = up;
           return *this;
       }

   private:
       void clean(const Matrix4& absoluteTransformationMatrix) override {
           Renderer::setListenerPosition(absoluteTransformationMatrix.translation());

           Renderer::setListenerOrientation(absoluteTransformationMatrix.rotation()*_fwd,
                                            absoluteTransformationMatrix.rotation()*_up);

           // TODO: velocity
       }

       Vector3 _fwd;
       Vector3 _up;
};


}}

#endif
