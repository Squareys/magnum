/*
    This file is part of Magnum.

    Copyright © 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

#include "Context.h"

#include <Corrade/Utility/Arguments.h>
#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/String.h>
#include "vulkan.h"

PFN_vkCreateDebugReportCallbackEXT CreateDebugReportCallback;
PFN_vkDestroyDebugReportCallbackEXT DestroyDebugReportCallback;
PFN_vkDebugReportMessageEXT dbgBreakCallback;

VkDebugReportCallbackEXT msgCallback;

namespace Magnum { namespace Vk {

const char *validationLayerNames[] =
{
    // This is a meta layer that enables all of the standard
    // validation layers in the correct order :
    // threading, parameter_validation, device_limits, object_tracker, image, core_validation, swapchain, and unique_objects
    "VK_LAYER_LUNARG_standard_validation"
};

Context* Context::_current = nullptr;

bool Context::hasCurrent() { return _current; }

Context& Context::current() {
    CORRADE_ASSERT(_current, "Context::current(): no current context", *_current);
    return *_current;
}

Context::Context(NoCreateT, Int argc, char** argv, void functionLoader()): _functionLoader{functionLoader}, _version{Version::None} {
    Utility::Arguments args{"magnum"};
    args.parse(argc, argv);
}

Context::Context(Flags flags): _functionLoader{nullptr}, _version{Version::None}, _flags(flags) {
    create();
}

Context::Context(Context&& other): _version{std::move(other._version)},
    _flags{std::move(other._flags)}
{
    if(_current == &other) _current = this;
}

Context::~Context() {
    if(_current == this) _current = nullptr;

    if(_instance == nullptr) {
        Error() << "No VkInstance to destroy";
        return;
    }

    if (_flags >= Flag::EnableValidation) {
        DestroyDebugReportCallback(_instance, msgCallback, nullptr);
    }
    vkDestroyInstance(_instance, nullptr);
}

void Context::create() {
    /* Hard exit if the context cannot be created */
    if(!tryCreate()) std::exit(1);
}

VkBool32 messageCallback(VkDebugReportFlagsEXT flags, VkDebugReportObjectTypeEXT objType,
    uint64_t srcObject,  size_t location, int32_t msgCode, const char* pLayerPrefix, const char* pMsg, void* pUserData) {

    if (flags & VK_DEBUG_REPORT_ERROR_BIT_EXT) {
        Error() << "[" << pLayerPrefix << "] Code" << msgCode << ":" << pMsg;
    } else if (flags & VK_DEBUG_REPORT_WARNING_BIT_EXT) {
        Warning() << "[" << pLayerPrefix << "] Code" << msgCode << ":" << pMsg;
    }

    return false;
}

bool Context::tryCreate() {
    VkApplicationInfo appInfo = {};
    appInfo.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    appInfo.pApplicationName = "Vulkan Example";
    appInfo.pEngineName = "Magnum";
    appInfo.apiVersion = UnsignedInt(Version::Vulkan_1_0);

    std::vector<const char*> enabledExtensions = { VK_KHR_SURFACE_EXTENSION_NAME };

    // Enable surface extensions depending on os
    //enabledExtensions.push_back(VK_KHR_WIN32_SURFACE_EXTENSION_NAME);

    VkInstanceCreateInfo instanceCreateInfo = {};
    instanceCreateInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    instanceCreateInfo.pNext = nullptr;
    instanceCreateInfo.pApplicationInfo = &appInfo;
    if (_flags >= Flag::EnableValidation) {
        enabledExtensions.push_back(VK_EXT_DEBUG_REPORT_EXTENSION_NAME);

        instanceCreateInfo.enabledLayerCount = 1;
        instanceCreateInfo.ppEnabledLayerNames = validationLayerNames;
    }
    if (enabledExtensions.size() > 0) {
        instanceCreateInfo.enabledExtensionCount = enabledExtensions.size();
        instanceCreateInfo.ppEnabledExtensionNames = enabledExtensions.data();
    }

    VkResult ret = vkCreateInstance(&instanceCreateInfo, nullptr, &_instance);
    if(ret != VK_SUCCESS) {
        Error() << "Vulkan Instance creation failed with error" << ret;
        return false;
    }

    _current = this;

    /* setup debugging */
    if (_flags >= Flag::EnableValidation) {
        CreateDebugReportCallback = (PFN_vkCreateDebugReportCallbackEXT)vkGetInstanceProcAddr(_instance, "vkCreateDebugReportCallbackEXT");
        DestroyDebugReportCallback = (PFN_vkDestroyDebugReportCallbackEXT)vkGetInstanceProcAddr(_instance, "vkDestroyDebugReportCallbackEXT");
        dbgBreakCallback = (PFN_vkDebugReportMessageEXT)vkGetInstanceProcAddr(_instance, "vkDebugReportMessageEXT");

        VkDebugReportCallbackCreateInfoEXT dbgCreateInfo = {};
        dbgCreateInfo.sType = VK_STRUCTURE_TYPE_DEBUG_REPORT_CREATE_INFO_EXT;
        dbgCreateInfo.pfnCallback = PFN_vkDebugReportCallbackEXT(messageCallback);
        dbgCreateInfo.flags = VK_DEBUG_REPORT_ERROR_BIT_EXT | VK_DEBUG_REPORT_WARNING_BIT_EXT;

        VkResult err = CreateDebugReportCallback(
            _instance,
            &dbgCreateInfo,
            nullptr,
            &msgCallback);
        assert(!err);
    }
    return true;
}

bool Context::isVersionSupported(Version) const {
    return true;
}

}}
