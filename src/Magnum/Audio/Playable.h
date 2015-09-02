#ifndef Magnum_Audio_Playable_h
#define Magnum_Audio_Playable_h
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
 * @brief Class @ref Magnum::Audio::Playable
 */

#include <string>
#include <al.h>

#include <Magnum/SceneGraph/SceneGraph.h>
#include <Magnum/SceneGraph/AbstractFeature.h>

#include "Magnum/Audio/visibility.h"
#include "Magnum/Audio/Source.h"

namespace Magnum { namespace Audio {

/**
@brief Playable feature

Manages the position of a @ref Source for an @ref SceneGraph::Object.

## Usage

Create the feature for an Object
@code
Scene3D scene;
Object3D object{&scene};
Source source;
Playable3D playable{object, source};
@endcode

Position of the source will be updated to the translation of `object` when
@ref SceneGraph::Object::setClean() is called, which is done in
@ref SceneGraph::Camera::draw() for example.

-   @ref Playable2D
-   @ref Playable3D

@see @ref Source, @ref Listener
*/
template<UnsignedInt dimensions> class Playable: public SceneGraph::AbstractFeature<dimensions, Float> {
    public:

        /**
         * @brief Constructor
         */
        explicit Playable(SceneGraph::AbstractObject<dimensions, Float>& object, Source& source):
            SceneGraph::AbstractFeature<dimensions, Float>{object},
            _source{source},
            _fwd{}
        {
            SceneGraph::AbstractFeature<dimensions, Float>::setCachedTransformations(SceneGraph::CachedTransformation::Absolute);
            _fwd[dimensions - 1] = 1;
        }

        /**
         * @brief Source which is managed by this feature
         */
        Source& source() {
            return _source;
        }

        /**
         * @brief Set source to manage
         * @return Reference to self (for method chaining)
         */
        Playable& setSource(Source& source) {
            _source = source;
            return *this;
        }

        /**
         * @brief Set the vector to consider as forward
         * @param fwd Forward vector
         *
         * Default `{0.0, 1.0}` for two dimensional and `{0.0, 0.0, 1.0}` for
         * three dimensional feature.
         * @return Reference to self (for method chaining)
         */
        Playable& setForward(VectorTypeFor<dimensions, Float>& fwd) {
            _fwd = fwd;
            return *this;
        }

    private:
        void clean(const MatrixTypeFor<dimensions, Float>& absoluteTransformationMatrix) override {
            _source.setPosition(Vector3::pad(absoluteTransformationMatrix.translation(), 0));
            _source.setDirection(Vector3::pad(absoluteTransformationMatrix.rotation()*_fwd));

            // TODO: velocity
        }

        Source& _source;
        VectorTypeFor<dimensions, Float> _fwd;
};

/**
 * @brief Playable for two dimensional float scenes
 *
 * @see @ref Playable3D
 */
typedef Playable<2> Playable2D;
/**
 * @brief Playable for three dimensional float scenes
 *
 * @see @ref Playable2D
 */
typedef Playable<3> Playable3D;

#if defined(CORRADE_TARGET_WINDOWS) && !defined(__MINGW32__)
extern template class MAGNUM_AUDIO_EXPORT Playable<2>;
extern template class MAGNUM_AUDIO_EXPORT Playable<3>;
#endif

}}

#endif
