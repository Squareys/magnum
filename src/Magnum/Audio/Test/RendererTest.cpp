/*
    This file is part of Magnum.

    Copyright © 2010, 2011, 2012, 2013, 2014, 2015
              Vladimír Vondruš <mosra@centrum.cz>

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

#include <sstream>
#include <Corrade/TestSuite/Tester.h>

#include "Magnum/Magnum.h"
#include "Magnum/Audio/Context.h"
#include "Magnum/Audio/Renderer.h"

namespace Magnum { namespace Audio { namespace Test {

struct RendererTest: TestSuite::Tester {
    explicit RendererTest();

    void debugError();
    void debugDistanceModel();
    void listenerOrientation();
    void listenerPosition();
    void listenerVelocity();
    void listenerGain();
    void speedOfSound();
    void dopplerFactor();
    void distanceModel();

    Context _context;
};

RendererTest::RendererTest(): _context{} {
    addTests({&RendererTest::debugError,
              &RendererTest::debugDistanceModel,
              &RendererTest::listenerOrientation,
              &RendererTest::listenerPosition,
              &RendererTest::listenerVelocity,
              &RendererTest::listenerGain,
              &RendererTest::speedOfSound,
              &RendererTest::dopplerFactor,
              &RendererTest::distanceModel});
}

void RendererTest::debugError() {
    std::ostringstream out;
    Debug(&out) << Renderer::Error::InvalidOperation;
    CORRADE_COMPARE(out.str(), "Audio::Renderer::Error::InvalidOperation\n");
}

void RendererTest::debugDistanceModel() {
    std::ostringstream out;
    Debug(&out) << Renderer::DistanceModel::Inverse;
    CORRADE_COMPARE(out.str(), "Audio::Renderer::DistanceModel::Inverse\n");
}

void RendererTest::listenerOrientation() {
    constexpr Vector3 up{1, 2, 3}, fwd{3, 2, 1};
    Renderer::setListenerOrientation(fwd, up);
    std::array<Vector3, 2> orientation = Renderer::listenerOrientation();

    CORRADE_COMPARE(orientation[0], fwd);
    CORRADE_COMPARE(orientation[1], up);
}

void RendererTest::listenerPosition() {
    constexpr Vector3 pos{1, 3, 2};
    Renderer::setListenerPosition(pos);

    CORRADE_COMPARE(Renderer::listenerPosition(), pos);
}

void RendererTest::listenerVelocity() {
    constexpr Vector3 vel{1, 3, 2};
    Renderer::setListenerVelocity(vel);

    CORRADE_COMPARE(Renderer::listenerVelocity(), vel);
}

void RendererTest::listenerGain() {
    constexpr Float gain = 0.512;
    Renderer::setListenerGain(gain);

    CORRADE_COMPARE(Renderer::listenerGain(), gain);
}

void RendererTest::speedOfSound() {
    constexpr Float speed = 1.25;
    Renderer::setSpeedOfSound(speed);

    CORRADE_COMPARE(Renderer::speedOfSound(), speed);
}

void RendererTest::dopplerFactor() {
    constexpr Float factor = 0.3335;
    Renderer::setDopplerFactor(factor);

    CORRADE_COMPARE(Renderer::dopplerFactor(), factor);
}

void RendererTest::distanceModel() {
    constexpr Renderer::DistanceModel model = Renderer::DistanceModel::InverseClamped;
    Renderer::setDistanceModel(model);

    CORRADE_COMPARE(Renderer::distanceModel(), model);
}

}}}

CORRADE_TEST_MAIN(Magnum::Audio::Test::RendererTest)
