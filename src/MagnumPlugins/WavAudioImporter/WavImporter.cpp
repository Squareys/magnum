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

#include "WavImporter.h"

#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/Debug.h>
#include <Corrade/Utility/Endianness.h>

#include "MagnumPlugins/WavAudioImporter/WavHeader.h"

namespace Magnum { namespace Audio {

WavImporter::WavImporter() = default;

WavImporter::WavImporter(PluginManager::AbstractManager& manager, std::string plugin): AbstractImporter(manager, std::move(plugin)) {}

auto WavImporter::doFeatures() const -> Features { return Feature::OpenData; }

bool WavImporter::doIsOpened() const { return _data; }

void WavImporter::doOpenData(Containers::ArrayView<const char> data) {
    /* Check file size */
    if(data.size() < sizeof(WavHeaderChunk) + sizeof(WavFormatChunk) + sizeof(RiffChunkHeader)) {
        Error() << "Audio::WavImporter::openData(): the file is too short:" << data.size() << "bytes";
        return;
    }

    /* Get the RIFF/WAV header */
    WavHeaderChunk header(*reinterpret_cast<const WavHeaderChunk*>(data.begin()));

    /* Check RIFF/WAV file signature */
    if(std::strncmp(header.chunk.chunkId, "RIFF", 4) != 0 ||
       std::strncmp(header.format, "WAVE", 4) != 0) {
        Error() << "Audio::WavImporter::openData(): the file signature is invalid";
        return;
    }

    Utility::Endianness::littleEndianInPlace(header.chunk.chunkSize);

    /* Check file size */
    if(header.chunk.chunkSize < 36 || header.chunk.chunkSize + 8 != data.size()) {
        Error() << "Audio::WavImporter::openData(): the file has improper size, expected"
                << header.chunk.chunkSize + 8 << "but got" << data.size();
        return;
    }

    data = data.suffix(sizeof(WavHeaderChunk));

    const RiffChunkHeader* chunkHeader;

    while(data.size() > 0) {
        chunkHeader = reinterpret_cast<const RiffChunkHeader*>(data.prefix(8).data());
        Utility::Endianness::littleEndianInPlace(chunkHeader->chunkSize);

        if(data.size() < sizeof(RiffChunkHeader) + chunkHeader->chunkSize) {
            Error() << "WavImporter::openData(): Unexpected end of file.";
            return;
        }

        if(std::strncmp(chunkHeader->chunkId, "fmt ", 4) == 0) {
            if(chunkHeader->chunkSize + sizeof(RiffChunkHeader) < sizeof(WavFormatChunk)) {
                Error() << "WavImporter::openData(): \"fmt \" chunk was too small.";
            }

            /* Get the WAV format header */
            WavFormatChunk formatChunk(*reinterpret_cast<const WavFormatChunk*>(data.begin()));

            /* Fix endianness on Format chunk */
            Utility::Endianness::littleEndianInPlace(
                formatChunk.chunk.chunkSize, formatChunk.audioFormat, formatChunk.numChannels,
                formatChunk.sampleRate, formatChunk.byteRate, formatChunk.blockAlign,
                formatChunk.bitsPerSample);

            /* Check PCM format */
            if(formatChunk.audioFormat == WAVE_FORMAT_PCM) {
                /* Decide about format */
                if(formatChunk.numChannels == 1 && formatChunk.bitsPerSample == 8)
                    _format = Buffer::Format::Mono8;
                else if(formatChunk.numChannels == 1 && formatChunk.bitsPerSample == 16)
                    _format = Buffer::Format::Mono16;
                else if(formatChunk.numChannels == 2 && formatChunk.bitsPerSample == 8)
                    _format = Buffer::Format::Stereo8;
                else if(formatChunk.numChannels == 2 && formatChunk.bitsPerSample == 16)
                    _format = Buffer::Format::Stereo16;
                else {
                    Error() << "Audio::WavImporter::openData(): unsupported channel count"
                            << formatChunk.numChannels << "with" << formatChunk.bitsPerSample
                            << "bits per sample";
                    return;
                }
            }
            /* Check IEEE Float format */
            else if(formatChunk.audioFormat == WAVE_FORMAT_IEEE_FLOAT) {
                if(formatChunk.numChannels == 1 && formatChunk.bitsPerSample == 32)
                    _format = Buffer::Format::MonoFloat;
                else if(formatChunk.numChannels == 2 && formatChunk.bitsPerSample == 32)
                    _format = Buffer::Format::StereoFloat;
                else if(formatChunk.numChannels == 1 && formatChunk.bitsPerSample == 64)
                    _format = Buffer::Format::MonoDouble;
                else if(formatChunk.numChannels == 2 && formatChunk.bitsPerSample == 64)
                    _format = Buffer::Format::StereoDouble;
                else {
                    Error() << "Audio::WavImporter::openData(): unsupported channel count"
                            << formatChunk.numChannels << "with" << formatChunk.bitsPerSample
                            << "bits per sample";
                    return;
                }
            }
            else if(formatChunk.audioFormat == WAVE_FORMAT_ALAW) {
                if(formatChunk.numChannels == 1)
                    _format = Buffer::Format::MonoALaw;
                else if(formatChunk.numChannels == 2)
                    _format = Buffer::Format::StereoALaw;
                else {
                    Error() << "Audio::WavImporter::openData(): unsupported channel count"
                            << formatChunk.numChannels << "with" << formatChunk.bitsPerSample
                            << "bits per sample";
                    return;
                }
            }
            else if(formatChunk.audioFormat == WAVE_FORMAT_MULAW) {
                if(formatChunk.numChannels == 1)
                    _format = Buffer::Format::MonoMuLaw;
                else if(formatChunk.numChannels == 2)
                    _format = Buffer::Format::StereoMuLaw;
                else {
                    Error() << "Audio::WavImporter::openData(): unsupported channel count"
                            << formatChunk.numChannels << "with" << formatChunk.bitsPerSample
                            << "bits per sample";
                    return;
                }
            }
            else if(formatChunk.audioFormat == WAVE_FORMAT_EXTENSIBLE) {
                Error() << "Audio::WavImporter::openData(): unsupported audio format: extensible not implememented" << formatChunk.audioFormat;
                return;
            }
            else {
                Error() << "Audio::WavImporter::openData(): unsupported audio format" << formatChunk.audioFormat;
                return;
            }

            /* Format sanity checks */
            if(formatChunk.blockAlign != formatChunk.numChannels * formatChunk.bitsPerSample / 8 ||
               formatChunk.byteRate != formatChunk.sampleRate * formatChunk.blockAlign) {
                Error() << "Audio::WavImporter::openData(): the file is corrupted";
                return;
            }

            /* Save frequency */
            _frequency = formatChunk.sampleRate;

        } else if(std::strncmp(chunkHeader->chunkId, "data", 4) == 0) {
            /** @todo Convert the data from little endian too */
            CORRADE_INTERNAL_ASSERT(!Utility::Endianness::isBigEndian());

            /* Copy the data */
            _data = Corrade::Containers::Array<char>{chunkHeader->chunkSize};

            Containers::ArrayView<const char> pureData = data.slice(sizeof(RiffChunkHeader), sizeof(RiffChunkHeader)+chunkHeader->chunkSize);
            std::copy(pureData.begin(), pureData.end(), _data.begin());
        } else if(std::strncmp(chunkHeader->chunkId, "JUNK", 4) == 0) {
            /* JUNK chunk is empty data for padding */
        } else {
            /* Skip unknown chunks */
        }

        /* Offset data to next chunk */
        data = data.suffix(sizeof(RiffChunkHeader) + chunkHeader->chunkSize + (chunkHeader->chunkSize % 2));
    }
    return;
}

void WavImporter::doClose() { _data = nullptr; }

Buffer::Format WavImporter::doFormat() const { return _format; }

UnsignedInt WavImporter::doFrequency() const { return _frequency; }

Containers::Array<char> WavImporter::doData() {
    Containers::Array<char> copy(_data.size());
    std::copy(_data.begin(), _data.end(), copy.begin());
    return copy;
}

}}
