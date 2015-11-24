/*
 * VstPresets.cpp
 * --------------
 * Purpose: VST plugin preset / bank handling
 * Notes  : (currently none)
 * Authors: OpenMPT Devs
 * The OpenMPT source code is released under the BSD license. Read LICENSE for more details.
 */


#include "stdafx.h"

#ifndef NO_VST
#include "../soundlib/Sndfile.h"
#include "Vstplug.h"
#include <vstsdk2.4/pluginterfaces/vst2.x/vstfxstore.h>
#include "VstPresets.h"
#include "../common/FileReader.h"
#include <ostream>


OPENMPT_NAMESPACE_BEGIN


#ifdef NEEDS_PRAGMA_PACK
#pragma pack(push, 1)
#endif

// This part of the header is identical for both presets and banks.
struct PACKED ChunkHeader
{
	VstInt32 chunkMagic;		///< 'CcnK'
	VstInt32 byteSize;			///< size of this chunk, excl. magic + byteSize

	VstInt32 fxMagic;			///< 'FxBk' (regular) or 'FBCh' (opaque chunk)
	VstInt32 version;			///< format version (1 or 2)
	VstInt32 fxID;				///< fx unique ID
	VstInt32 fxVersion;			///< fx version

	// Convert all multi-byte numeric values to current platform's endianness or vice versa.
	void ConvertEndianness()
	{
		SwapBytesBE(chunkMagic);
		SwapBytesBE(byteSize);
		SwapBytesBE(fxMagic);
		SwapBytesBE(version);
		SwapBytesBE(fxID);
		SwapBytesBE(fxVersion);
	}
};

STATIC_ASSERT(sizeof(ChunkHeader) == 24);

#ifdef NEEDS_PRAGMA_PACK
#pragma pack(pop)
#endif

VSTPresets::ErrorCode VSTPresets::LoadFile(FileReader &file, IMixPlugin &plugin)
//------------------------------------------------------------------------------
{
	const bool firstChunk = file.GetPosition() == 0;
	ChunkHeader header;
	if(!file.ReadConvertEndianness(header) || header.chunkMagic != cMagic)
	{
		return invalidFile;
	}
	if(header.fxID != plugin.GetUID())
	{
		return wrongPlugin;
	}
	CVstPlugin *vstPlug = dynamic_cast<CVstPlugin *>(&plugin);

	if(header.fxMagic == fMagic || header.fxMagic == chunkPresetMagic)
	{
		// Program
		PlugParamIndex numParams = file.ReadUint32BE();

		if(vstPlug != nullptr)
		{
			VstPatchChunkInfo info;
			info.version = 1;
			info.pluginUniqueID = header.fxID;
			info.pluginVersion = header.fxVersion;
			info.numElements = numParams;
			MemsetZero(info.future);
			vstPlug->Dispatch(effBeginLoadProgram, 0, 0, &info, 0.0f);
		}
		plugin.BeginSetProgram();

		char prgName[28];
		file.ReadString<mpt::String::maybeNullTerminated>(prgName, 28);
		plugin.SetCurrentProgramName(prgName);

		if(header.fxMagic == fMagic)
		{
			if(plugin.GetNumParameters() != numParams)
			{
				return wrongParameters;
			}
			for(PlugParamIndex p = 0; p < numParams; p++)
			{
				plugin.SetParameter(p, file.ReadFloatBE());
			}
		} else
		{
			uint32 chunkSize = file.ReadUint32BE();
			// Some nasty plugins (e.g. SmartElectronix Ambience) write to our memory block.
			// Directly writing to a memory-mapped file block results in a crash...
			char *chunkData = new (std::nothrow) char[chunkSize];
			if(chunkData)
			{
				file.ReadRaw(chunkData, chunkSize);
				plugin.SetChunk(chunkSize, chunkData, false);
				delete[] chunkData;
			} else
			{
				return outOfMemory;
			}
		}
		plugin.EndSetProgram();
	} else if((header.fxMagic == bankMagic || header.fxMagic == chunkBankMagic) && firstChunk)
	{
		// Bank - only read if it's the first chunk in the file, not if it's a sub chunk.
		uint32 numProgs = file.ReadUint32BE();
		uint32 currentProgram = file.ReadUint32BE();
		file.Skip(124);

		if(vstPlug != nullptr)
		{
			VstPatchChunkInfo info;
			info.version = 1;
			info.pluginUniqueID = header.fxID;
			info.pluginVersion = header.fxVersion;
			info.numElements = numProgs;
			MemsetZero(info.future);
			vstPlug->Dispatch(effBeginLoadBank, 0, 0, &info, 0.0f);
		}

		if(header.fxMagic == bankMagic)
		{
			VstInt32 oldCurrentProgram = plugin.GetCurrentProgram();
			for(uint32 p = 0; p < numProgs; p++)
			{
				plugin.BeginSetProgram(p);
				ErrorCode retVal = LoadFile(file, plugin);
				if(retVal != noError)
				{
					return retVal;
				}
				plugin.EndSetProgram();
			}
			plugin.SetCurrentProgram(oldCurrentProgram);
		} else
		{
			uint32 chunkSize = file.ReadUint32BE();
			// Some nasty plugins (e.g. SmartElectronix Ambience) write to our memory block.
			// Directly writing to a memory-mapped file block results in a crash...
			char *chunkData = new (std::nothrow) char[chunkSize];
			if(chunkData)
			{
				file.ReadRaw(chunkData, chunkSize);
				plugin.SetChunk(chunkSize, chunkData, true);
				delete[] chunkData;
			} else
			{
				return outOfMemory;
			}
		}
		if(header.version >= 2)
		{
			plugin.SetCurrentProgram(currentProgram);
		}
	}

	return noError;
}


bool VSTPresets::SaveFile(std::ostream &f, IMixPlugin &plugin, bool bank)
//-----------------------------------------------------------------------
{
	if(!bank)
	{
		SaveProgram(f, plugin);
	} else
	{
		bool writeChunk = plugin.ProgramsAreChunks();
		ChunkHeader header;
		header.chunkMagic = cMagic;
		header.byteSize = 0; // will be corrected later
		header.version = 2;
		header.fxID = plugin.GetUID();
		header.fxVersion = plugin.GetVersion();

		// Write unfinished header... We need to update the size once we're done writing.
		mpt::IO::WriteConvertEndianness(f, header);

		uint32 numProgs = std::max(plugin.GetNumPrograms(), VstInt32(1)), curProg = plugin.GetCurrentProgram();
		mpt::IO::WriteIntBE(f, numProgs);
		mpt::IO::WriteIntBE(f, curProg);
		uint8 reserved[124];
		MemsetZero(reserved);
		mpt::IO::WriteRaw(f, reserved, sizeof(reserved));

		if(writeChunk)
		{
			char *chunk = nullptr;
			uint32 chunkSize = mpt::saturate_cast<uint32>(plugin.GetChunk(chunk, true));
			if(chunkSize && chunk)
			{
				mpt::IO::WriteIntBE(f, chunkSize);
				mpt::IO::WriteRaw(f, chunk, chunkSize);
			} else
			{
				// The plugin returned no chunk! Gracefully go back and save parameters instead...
				writeChunk = false;
			}
		}
		if(!writeChunk)
		{
			for(uint32 p = 0; p < numProgs; p++)
			{
				plugin.SetCurrentProgram(p);
				SaveProgram(f, plugin);
			}
			plugin.SetCurrentProgram(curProg);
		}

		// Now we know the correct chunk size.
		std::streamoff end = f.tellp();
		header.byteSize = static_cast<VstInt32>(end - 8);
		header.fxMagic = writeChunk ? chunkBankMagic : bankMagic;
		mpt::IO::SeekBegin(f);
		mpt::IO::WriteConvertEndianness(f, header);
	}

	return true;
}


void VSTPresets::SaveProgram(std::ostream &f, IMixPlugin &plugin)
//---------------------------------------------------------------
{
	bool writeChunk = plugin.ProgramsAreChunks();
	ChunkHeader header;
	header.chunkMagic = cMagic;
	header.byteSize = 0; // will be corrected later
	header.version = 1;
	header.fxID = plugin.GetUID();
	header.fxVersion = plugin.GetVersion();

	// Write unfinished header... We need to update the size once we're done writing.
	mpt::IO::Offset start = mpt::IO::TellWrite(f);
	mpt::IO::WriteConvertEndianness(f, header);

	const uint32 numParams = plugin.GetNumParameters();
	mpt::IO::WriteIntBE(f, numParams);

	char name[28];
	mpt::String::Write<mpt::String::maybeNullTerminated>(name, mpt::ToCharset(mpt::CharsetLocale, plugin.GetCurrentProgramName()));
	mpt::IO::WriteRaw(f, name, 28);

	if(writeChunk)
	{
		char *chunk = nullptr;
		uint32 chunkSize = mpt::saturate_cast<uint32>(plugin.GetChunk(chunk, false));
		if(chunkSize && chunk)
		{
			mpt::IO::WriteIntBE(f, chunkSize);
			mpt::IO::WriteRaw(f, chunk, chunkSize);
		} else
		{
			// The plugin returned no chunk! Gracefully go back and save parameters instead...
			writeChunk = false;
		}
	}
	if(!writeChunk)
	{
		for(uint32 p = 0; p < numParams; p++)
		{
			mpt::IO::Write(f, IEEE754binary32BE(plugin.GetParameter(p)));
		}
	}

	// Now we know the correct chunk size.
	mpt::IO::Offset end = mpt::IO::TellWrite(f);
	header.byteSize = static_cast<VstInt32>(end - start - 8);
	header.fxMagic = writeChunk ? chunkPresetMagic : fMagic;
	mpt::IO::SeekAbsolute(f, start);
	mpt::IO::WriteConvertEndianness(f, header);
	mpt::IO::SeekAbsolute(f, end);
}


// Translate error code to string. Returns nullptr if there was no error.
const char *VSTPresets::GetErrorMessage(ErrorCode code)
//-----------------------------------------------------
{
	switch(code)
	{
	case VSTPresets::invalidFile:
		return "This does not appear to be a valid preset file.";
	case VSTPresets::wrongPlugin:
		return "This file appears to be for a different plugin.";
	case VSTPresets::wrongParameters:
		return "The number of parameters in this file is incompatible with the current plugin.";
	case VSTPresets::outOfMemory:
		return "Not enough memory to load preset data.";
	}
	return nullptr;
}

#endif // NO_VST


OPENMPT_NAMESPACE_END
