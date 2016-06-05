#ifndef Magnum_Vk_ShaderStage_h
#define Magnum_Vk_ShaderStage_h
/*
    This file is part of Magnum.

    Copyright © 2010, 2011, 2012, 2013, 2014, 2015, 2016
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
 * @brief Class @ref Magnum::Vk::ShaderStage
 */

#include "Magnum/Magnum.h"
#include "Magnum/Vk/visibility.h"

#include "Magnum/Vk/Device.h"
#include "Magnum/Vk/RenderPass.h"
#include "Magnum/Vk/Shader.h"
#include "Magnum/Vk/DescriptorSet.h"

#include "Magnum/Math/Vector3.h" // TEMPORARY!!!

#include <Corrade/Containers/Array.h>

#include "vulkan.h"

namespace Magnum { namespace Vk {


enum class ShaderStage: UnsignedInt {
    Vertex = VK_SHADER_STAGE_VERTEX_BIT,
    TesslationControl = VK_SHADER_STAGE_TESSELLATION_CONTROL_BIT,
    TesslationEvaluation = VK_SHADER_STAGE_TESSELLATION_EVALUATION_BIT,
    Geometry = VK_SHADER_STAGE_GEOMETRY_BIT,
    Fragment = VK_SHADER_STAGE_FRAGMENT_BIT,
    Compute = VK_SHADER_STAGE_COMPUTE_BIT,
    AllGraphics = VK_SHADER_STAGE_ALL_GRAPHICS,
    All = VK_SHADER_STAGE_ALL,
};

typedef Containers::EnumSet<ShaderStage> ShaderStageFlags;

CORRADE_ENUMSET_OPERATORS(ShaderStageFlags)

}}

#endif
