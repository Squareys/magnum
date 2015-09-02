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

#include <sstream>
#include <Corrade/TestSuite/Tester.h>
#include <Corrade/TestSuite/Comparator.h>

#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>

#include "Magnum/Audio/Buffer.h"
#include "Magnum/Audio/Listener.h"
#include "Magnum/Audio/Context.h"

namespace Magnum { namespace Audio { namespace Test {

typedef SceneGraph::Scene<SceneGraph::MatrixTransformation3D> Scene3D;
typedef SceneGraph::Object<SceneGraph::MatrixTransformation3D> Object3D;

struct ListenerTest: TestSuite::Tester {
    explicit ListenerTest();

    void testFeature();

    Context _context;
};

ListenerTest::ListenerTest(): _context() {
    addTests({&ListenerTest::testFeature});
}

void ListenerTest::testFeature() {
    Scene3D scene;
    Object3D object{&scene};
    Listener3D listener{object};

    constexpr Vector3 offset{1, 2, 3};
    object.translate(offset);
    object.setClean();

    // TODO: Wait until Renderer::listenerPosition() has been implemented.
    CORRADE_COMPARE(Vector3(1, 2, 3), Vector3(1, 2, 3));
}

}}}

CORRADE_TEST_MAIN(Magnum::Audio::Test::ListenerTest)
