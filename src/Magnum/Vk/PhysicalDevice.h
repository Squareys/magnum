#ifndef Magnum_Vk_PhysicalDevice_h
#define Magnum_Vk_PhysicalDevice_h
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
 * @brief Class @ref Magnum::Vk::PhysicalDevice
 */

#include "Magnum/Magnum.h"
#include "Magnum/Vk/visibility.h"

#include <Corrade/Containers/EnumSet.h>

#include "vulkan.h"


namespace Magnum { namespace Vk {

enum class MemoryProperty: UnsignedInt {
    DeviceLocal = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
    HostVisible = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT,
    Coherent = VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
    HostCached = VK_MEMORY_PROPERTY_HOST_CACHED_BIT,
    LazilyAllocated = VK_MEMORY_PROPERTY_LAZILY_ALLOCATED_BIT,
};

typedef Corrade::Containers::EnumSet<MemoryProperty> MemoryProperties;
CORRADE_ENUMSET_OPERATORS(MemoryProperties)

enum class QueueFamily: UnsignedInt {
    Graphics = VK_QUEUE_GRAPHICS_BIT,
    Compute = VK_QUEUE_COMPUTE_BIT,
    Transfer = VK_QUEUE_TRANSFER_BIT,
};

enum class QueueFlag: UnsignedInt {
    SparseBinding = VK_QUEUE_SPARSE_BINDING_BIT,
};


class MAGNUM_VK_EXPORT PhysicalDevice {
    public:

        PhysicalDevice(): _physicalDevice{} {
        }

        PhysicalDevice(const PhysicalDevice& d):
            PhysicalDevice(d._physicalDevice) {
        }

        /** @brief Move constructor */
        PhysicalDevice(PhysicalDevice&& other);

        /**
         * @brief Construct from VkPhysicalDevice
         */
        PhysicalDevice(const VkPhysicalDevice& device):
            _physicalDevice(device)
        {
            _deviceMemoryProperties.memoryHeapCount = 0;
            _deviceMemoryProperties.memoryTypeCount = 0;
        }

        /** @brief Destructor */
        ~PhysicalDevice();

        /** @brief Copying is not allowed */
        PhysicalDevice& operator=(const PhysicalDevice& device){
            _physicalDevice = device._physicalDevice;

            _deviceMemoryProperties.memoryHeapCount = 0;
            _deviceMemoryProperties.memoryTypeCount = 0;

            return *this;
        }

        /** @brief Move assignment is not allowed */
        PhysicalDevice& operator=(PhysicalDevice&&) = delete;

        /**
         * @brief Get the underlying VkPhysicalDevice handle
         */
        operator VkPhysicalDevice() const {
            return _physicalDevice;
        }

        VkFormat getSupportedDepthFormat();

        UnsignedInt getQueueFamilyIndex(QueueFamily family);

        VkPhysicalDeviceProperties getProperties() const {
            VkPhysicalDeviceProperties deviceProperties;
            vkGetPhysicalDeviceProperties(_physicalDevice, &deviceProperties);

            return deviceProperties;
        }

        /**
         * @brief Get the index of the memory type with properties
         * @param supportedTypeBits A bit set with a bit for every memory type index. A
         *                  bit at index i should be set if the memory type at index i
         *                  is supported by the resource for which to get the memory type.
         * @param properties Properties required by the memory type
         * @return Index of the memory type for this device
         */
        UnsignedInt getMemoryType(UnsignedInt supportedTypeBits, MemoryProperties properties);

        VkPhysicalDeviceFeatures getFeatures() const {
            VkPhysicalDeviceFeatures features;
            vkGetPhysicalDeviceFeatures(_physicalDevice, &features);
            return features;
        }

    private:
        VkPhysicalDevice _physicalDevice;
        VkPhysicalDeviceMemoryProperties _deviceMemoryProperties;
};

}}

#endif
