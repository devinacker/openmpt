/*
* Load_symmod.cpp
* ------------
* Purpose: SymMOD (Symphonie Pro) module loader
* Notes  : Some effects aren't supported.
*          Virtual instruments are partially supported.
*          DSP settings aren't supported (yet?), but instruments which use the DSP will be set to use
*          plugin number 1, so you can easily configure a suitable effect yourself in OpenMPT.
* Authors: Devin Acker
*
* The OpenMPT source code is released under the BSD license. Read LICENSE for more details.
*
* Based in part on Patrick Meng's Java-based Symphonie player and its source.
* Some effect behavior and other things are based on the original Amiga assembly source.
*/

#include "stdafx.h"
#include "Loaders.h"
#include "modsmp_ctrl.h"

#include <vector>
#include <map>

#ifdef LIBOPENMPT_BUILD
// TODO fix
//#define MPT_SYMMOD_USE_REAL_SUBSONGS
#endif

OPENMPT_NAMESPACE_BEGIN

struct SymEvent {
	enum Command
	{
		KeyOn,
		VolSlideUp,
		VolSlideDown,
		PitchSlideUp,
		PitchSlideDown,
		ReplayFrom,
		FromAndPitch,
		SetFromAdd,
		FromAdd,
		SetSpeed,
		AddPitch,
		AddVolume,
		Vibrato,
		Tremolo,
		SampleVib,
		PitchSlideTo,
		Retrig,
		Emphasis,
		AddHalfTone,
		CV,
		CVAdd,

		Filter = 23,
		DSPEcho,
		DSPDelay,
	};
	uint8be command;
	int8be  note;
	enum Volume
	{
		VolCommand = 200,
		StopSample = 254,
		ContSample = 253,
		StartSample = 252, // unused
		KeyOff = 251,
		SpeedDown = 250,
		SpeedUp = 249,
		SetPitch = 248,
		PitchUp = 247,
		PitchDown = 246,
		PitchUp2 = 245,
		PitchDown2 = 244,
		PitchUp3 = 243,
		PitchDown3 = 242
	};
	uint8be param;
	uint8be inst;

	// used to compare DSP events for mapping them to MIDI macro numbers
	bool operator < (const SymEvent& other) const
	{
#define cmp(x) if (x != other.x) return x < other.x
		cmp(command);
		cmp(note);
		cmp(param);
		cmp(inst);
#undef  cmp
		return false;
	}
};

MPT_BINARY_STRUCT(SymEvent, 4);


// Virtual instrument info
// This allows instruments to be created based on a mix of other instruments.
// The sample mixing is done at load time.
struct SymVirtualInst {
	uint32be id; // "ViRT"
	uint8be dummy[10]; 

	uint16be numEvents;
	uint16be maxEvents; // always 20
	uint16be eventSize; // always 4

	SymEvent events[20];
};


// Instrument definition
struct SymInstrument {
	union
	{
		char name[128];
		SymVirtualInst virt;
	};

	enum Type 
	{
		Silent  = -8,
		Kill    = -4,
		Normal  = 0,
		Loop    = 4,
		Sustain = 8
	};
	int8be   type;
	uint8be  loopStart;
	uint8be  loopLen;
	uint8be  numLoops; // for "sustain" instruments
	enum Channels {
		Mono,
		StereoL,
		StereoR,
		LineSrc // virtual mix instrument
	};
	uint8be  channels;
	uint8be  dummy1; // called "automaximize" (normalize?) in amiga source, but unused
	uint8be  volume; // 0-200
	uint8be  dummy2[3]; // info about "parent/child" and sample format
	int8be   finetune; // -128..127 ~= 2 semitones
	int8be   transpose;
	enum SampleFlags 
	{
		PlayReverse   = 1, // reverse sample
		AsQueue       = 2, // "queue" virtual instrument (rendered track)
		MirrorX       = 4, // phase invert sample
		Is16Bit       = 8, // not used, we already know the bit depth of the samples
		NewLoopSystem = 16, // use fine loop start/len values

		MakeNewSample = (PlayReverse | MirrorX)
	};
	uint8be  sampleFlags;
	uint8be  filterType; // maybe same as filter events (0 = none, 1 = lowpass, 2 = highpass, 3 = bandpass)?
	enum InstFlags 
	{
		NoTranspose = 1, // don't apply sequence/position transpose
		NoDSP       = 2, // don't apply DSP effects
		SyncPlay    = 4  // play a stereo instrument pair on consecutive channels
	};
	uint8be  instFlags;
	uint8be  downsample; // downsample factor; affects sample tuning
	uint8be  dummy4[2];  // resonance, "loadflags" (both unused)
	uint8be  info; // ? (only lowest bit used?)
	uint8be  rangeStart; // ?
	uint8be  rangeLen; // ?
	uint8be  dummy5;
	uint16be loopStartFine;
	uint16be loopLenFine;
	uint8be  dummy6[6];
	
	uint8be  filterFlags; // ?
	uint8be  filterPoints; // # of filter envelope points (up to 4, possibly only 1-2 ever actually used)
	struct
	{
		uint8be cutoff;
		uint8be resonance;
	} filter[4];
	
	uint8be  padding[86];

	bool IsVirtual() const
	{
		return virt.id == MAGIC4BE('V', 'i', 'R', 'T')
			&& virt.eventSize == sizeof(SymEvent)
			&& virt.numEvents <= virt.maxEvents
			&& virt.maxEvents <= CountOf(virt.events);
	}

	void ConvertToMPT(ModInstrument &mptInst, ModSample &mptSmp, CSoundFile& module) const
	{
		// name
		if (!IsVirtual())
			mpt::String::Read<mpt::String::maybeNullTerminated>(mptInst.name, name);

		// loop
		// TODO: figure out why loop points not usually being set correctly
		SmpLength nLoopStart = static_cast<SmpLength>(loopStart) << 16;
		SmpLength nLoopLen = (static_cast<SmpLength>(loopLen) << 16);
		if (sampleFlags & NewLoopSystem)
		{
			nLoopStart += loopStartFine;
			nLoopLen += loopLenFine;
		}
		const double loopScale = static_cast<double>(mptSmp.nLength) / (100 << 16);
		nLoopStart = mpt::saturate_cast<SmpLength>(nLoopStart * loopScale);
		nLoopLen = std::min(mptSmp.nLength, mpt::saturate_cast<SmpLength>(nLoopLen * loopScale));
		if (type == Loop)
		{
			mptSmp.uFlags.set(CHN_LOOP);

			mptSmp.nLoopStart = nLoopStart;
			mptSmp.nLoopEnd = nLoopStart + nLoopLen;
		}
		else if (type == Sustain && numLoops > 1)
		{
			// TODO: set up multiple loop
			mptSmp.uFlags.reset(CHN_LOOP);
		}
		else
		{
			mptSmp.uFlags.reset(CHN_LOOP);
		}

		// volume (0-200, default 100)
		// TODO: amplify sample if >100? (Symphonie does)
		if (volume > 0 && volume <= 200)
		{
			mptSmp.nGlobalVol = std::min<uint16>(volume * 64 / 200, 64);
		}
		else
		{
			mptSmp.nGlobalVol = 32u;
		}

		// apply reverse/invert flags
		if (sampleFlags & PlayReverse)
		{
			ctrlSmp::ReverseSample(mptSmp, 0, 0, module);
		}
		if (sampleFlags & MirrorX)
		{
			ctrlSmp::InvertSample(mptSmp, 0, 0, module);
		}

		// tuning info
		// importing a sample of A440 at 44.1kHz and playing it at the closest possible note in Symphonie yields ~427.64Hz
		// which is almost exactly half a semitone lower. if we assume an ideal base sample rate of 44.1kHz, then we should
		// detune it by the same amount in order to correctly match the original tracker
		// TODO: why do some songs sound badly out of tune with the winamp plugin (or maybe all libopenmpt?)
		mptSmp.nC5Speed = 42862u; 
		mptSmp.Transpose(-downsample + (static_cast<double>(transpose) / 12) + (static_cast<double>(finetune) / (128 * 12)));

#ifndef NO_PLUGINS
		// DSP settings
		mptInst.nMixPlug = (instFlags & NoDSP) ? 0 : 1;
#endif

		// filter settings
		if (filterFlags != 0)
		{
			/*
			mptInst.SetCutoff(0x7F, true);
			// just average both resonance values since we can't interpolate them
			mptInst.SetResonance((resonance1 + resonance2) >> 2, true);

			// create filter sweep envelope
			// Symphonie applies this to the actual sample, which means the sweep rate depends
			// on the sample length and the note being played, which we can't really do here
			InstrumentEnvelope &filterEnv = mptInst.GetEnvelope(EnvelopeType::ENV_PITCH);
			filterEnv.dwFlags.set(ENV_ENABLED | ENV_FILTER);

			// TODO: figure out a good way to fake the length
			uint16 sweepLen = 1 + (mptSmp.nLength / 350);
			filterEnv.push_back(EnvelopeNode(0, cutoff1 >> 2));
			filterEnv.push_back(EnvelopeNode(sweepLen, cutoff2 >> 2));
			*/

			mptInst.SetCutoff(filter[0].cutoff >> 1, true);
			mptInst.SetResonance(filter[0].resonance >> 1, true);
		}
		else if (instFlags & NoDSP)
		{
			mptInst.SetCutoff(127, true);
			mptInst.SetResonance(0, true);
		}
	}
};

MPT_BINARY_STRUCT(SymInstrument, 256);

struct SymSequence {
	uint16be start;
	uint16be length;
	uint16be loop;
	int16be  info;
	int16be  transpose;
	
	uint8be  padding[6];
};

MPT_BINARY_STRUCT(SymSequence, 16);

struct SymPosition {
	uint8be  dummy[4];
	uint16be loopNum;
	uint16be loopCount; // unused
	uint16be pattern;
	uint16be start;
	uint16be length;
	uint16be speed;
	int16be  transpose;
	uint16be eventsPerLine; // not used, hopefully

	uint8be  padding[12];

	// used to compare position entries for mapping them to OpenMPT patterns
	bool operator < (const SymPosition& other) const
	{
#define cmp(x) if (x != other.x) return x < other.x
		cmp(pattern);
		cmp(start);
		cmp(length);
		cmp(transpose);
		cmp(speed);
		cmp(loopNum);
#undef  cmp
		return false;
	}
};

MPT_BINARY_STRUCT(SymPosition, 32);


static size_t DecodeChunk(FileReader &file, mpt::byte*& data)
//-----------------------------------------------------------
{
	size_t packedLength, unpackedLength;

	data = 0;
	packedLength = file.ReadUint32BE();
	unpackedLength = 0;

	FileReader chunk = file.ReadChunk(packedLength);

	if (chunk.CanRead(packedLength))
	{
		if (chunk.ReadMagic("PACK\xFF\xFF"))
		{
			// compressed chunk
			unpackedLength = chunk.ReadUint32BE();
			data = new mpt::byte[unpackedLength];

			bool done = false;
			size_t offset = 0;

			while (!done && !chunk.EndOfFile())
			{
				int8 type;
				uint8 len;
				mpt::byte dword[4];

				type = chunk.ReadInt8();

				switch (type)
				{
				case 0:
					// copy raw bytes
					len = chunk.ReadUint8();
					if (offset + len > unpackedLength || !chunk.CanRead(len))
					{
						done = true;
					}
					else
					{
						chunk.ReadRaw(data + offset, len);
						offset += len;
					}
					break;

				case 1:
					// copy a dword multiple times
					len = chunk.ReadUint8();
					if (offset + (len*4) > unpackedLength || !chunk.CanRead(4))
					{
						done = true;
					}
					else
					{
						chunk.ReadArray(dword);
						while (len--)
						{
							data[offset++] = dword[0];
							data[offset++] = dword[1];
							data[offset++] = dword[2];
							data[offset++] = dword[3];
						}
					}
					break;

				case 2:
					// copy a dword twice
					if (offset + 8 > unpackedLength || !chunk.CanRead(4))
					{
						done = true;
					}
					else
					{
						chunk.ReadArray(dword);
						data[offset++] = dword[0];
						data[offset++] = dword[1];
						data[offset++] = dword[2];
						data[offset++] = dword[3];
						data[offset++] = dword[0];
						data[offset++] = dword[1];
						data[offset++] = dword[2];
						data[offset++] = dword[3];
					}
					break;

				case 3:
					len = chunk.ReadUint8();
					if (offset + len > unpackedLength)
					{
						done = true;
					}
					else
					{
						std::memset(data + offset, 0, len);
						offset += len;
					}
					break;

				case -1:
					done = true;
					break;

				default:
					// error
					offset = 0;
					done = true;
					break;
				}
			}

			if (offset < unpackedLength)
			{
				unpackedLength = 0;
			}
		}
		else
		{
			// uncompressed chunk
			chunk.Rewind();

			unpackedLength = packedLength;
			data = new mpt::byte[unpackedLength];
			if (chunk.ReadRaw(data, unpackedLength) != unpackedLength)
			{
				unpackedLength = 0;
			}
		}

		if (!unpackedLength)
		{
			delete[] data;
			data = 0;
		}
	}

	return unpackedLength;
}


template<typename T>
static uint32 DecodeArray(FileReader &file, T*& array)
//----------------------------------------------------
{
	mpt::byte *data;
	size_t size = DecodeChunk(file, data);
	array = reinterpret_cast<T*>(data);
	return (size / sizeof(T));
}


static bool ReadRawSample(ModSample &sample, FileReader &file)
//------------------------------------------------------------------
{
	SampleIO::Bitdepth depth = SampleIO::_8bit;
	SampleIO::Channels channels = SampleIO::mono;

	sample.Initialize();

	if (file.Rewind(), file.ReadMagic("MAESTRO"))
	{
		depth = SampleIO::_16bit;

		file.Seek(12);
		if (file.ReadUint32BE() == 0)
		{
			channels = SampleIO::stereoInterleaved;
			sample.nLength = file.BytesLeft() / 4;
		}
		else
		{
			sample.nLength = file.BytesLeft() / 2;
		}
	}
	else if (file.Rewind(), file.ReadMagic("16BT")) 
	{
		depth = SampleIO::_16bit;

		file.Seek(6);
		sample.nLength = file.BytesLeft() / 2;
	}
	else
	{
		file.Rewind();
		sample.nLength = file.BytesLeft();
	}

	return SampleIO(
		depth, channels,
		SampleIO::bigEndian,
		SampleIO::signedPCM)
		.ReadSample(sample, file) > 0;
}


static bool ReadUnpackedSample(CSoundFile &module, ModSample &sample, SAMPLEINDEX sampleIndex, FileReader &file)
//--------------------------------------------------------------------------------------------------------------
{
	return 
		module.ReadIFFSample(sampleIndex, file) ||
		module.ReadWAVSample(sampleIndex, file) ||
		module.ReadAIFFSample(sampleIndex, file) ||
		ReadRawSample(sample, file);
}


static bool DecodeSample8(CSoundFile &module, ModSample &sample, SAMPLEINDEX sampleIndex, FileReader &file)
//---------------------------------------------------------------------------------------------------------
{
	mpt::byte *data = 0;
	size_t length = DecodeChunk(file, data);

	if (length > 0)
	{
		mpt::byte byte = data[0];

		for (size_t i = 1; i < length; i++)
		{
			byte += data[i];
			data[i] = byte;
		}
	}

	FileReader reader(mpt::const_byte_span(data, length));
	bool ok = ReadUnpackedSample(module, sample, sampleIndex, reader);

	delete[] data;
	return ok;
}


static bool DecodeSample16(CSoundFile &module, ModSample &sample, SAMPLEINDEX sampleIndex, FileReader &file)
//----------------------------------------------------------------------------------------------------------
{
	mpt::byte *data = 0;
	size_t length = DecodeChunk(file, data);
	mpt::byte buf[4096];
	// size of block in 16-bit samples
	const size_t blockSize = sizeof(buf) / 2;

	for (unsigned block = 0; block < length / sizeof(buf); block++)
	{
		mpt::byte byte = 0;
		mpt::byte *blockSrc = data + (block * sizeof(buf));

		// decode LSBs
		for (int i = 0; i < blockSize; i++)
		{
			byte += blockSrc[i];
			buf[i*2+1] = byte;
		}
		// decode MSBs
		for (int i = 0; i < blockSize; i++)
		{
			byte += blockSrc[i + blockSize];
			buf[i * 2] = byte;
		}

		std::memcpy(blockSrc, buf, sizeof(buf));
	}

	FileReader reader(mpt::const_byte_span(data, length));
	bool ok = ReadUnpackedSample(module, sample, sampleIndex, reader);

	delete[] data;
	return ok;
}


static void NormalizeSample(ModSample &sample)
//--------------------------------------------
{
	// pretty much lifted directly from Ctrl_smp.cpp
	SmpLength end = sample.nLength * sample.GetNumChannels();

	if (sample.uFlags[CHN_16BIT])
	{
		int16 *p = sample.pSample16;
		int max = 1;
		for (SmpLength i = 0; i < end; i++)
		{
			if (p[i] > max) max = p[i];
			if (-p[i] > max) max = -p[i];
		}
		if (max < 32767)
		{
			max++;
			for (SmpLength j = 0; j < end; j++)
			{
				int l = (((int)p[j]) << 15) / max;
				p[j] = (int16)l;
			}
		}
	}
	else
	{
		int8 *p = sample.pSample8;
		int max = 1;
		for (SmpLength i = 0; i < end; i++)
		{
			if (p[i] > max) max = p[i];
			if (-p[i] > max) max = -p[i];
		}
		if (max < 127)
		{
			max++;
			for (SmpLength j = 0; j < end; j++)
			{
				int l = (((int)p[j]) << 7) / max;
				p[j] = (int8)l;
			}
		}
	}
}


static void ConvertDSP(const SymEvent *event, char *str)
//----------------------------------------------------------------------
{
	if (event->command == SymEvent::Filter)
	{
		// TODO: this is obviously shit, figure out the real filter params
		// Also check the Java player source to see how they do DSPs

		// 
		// CFILTER_MAXRESO	EQU	185
		// CFILTER_MAXFREQ	EQU	240
		/*
		uint8 cutoff = 0x40 + (event->param >> 2);
		uint8 reso = event->inst >> 1;

		uint8 cutoff = std::min<uint8>(127, event->param * 127 / 240);
		uint8 reso = std::min<uint8>(127, event->inst * 127 / 185);
		*/

		// uint8 cutoff = std::min<uint8>(127, (event->param * 36000 / 28867) >> 1);
		uint8 cutoff = std::min<uint8>(127, event->param * 3 / 2);
		uint8 reso = std::min<uint8>(127, event->inst >> 1);

		if (event->note == 1) // lowpass filter
		{
			// TODO: make sure this is the right way around
			std::sprintf(str, "F0F000%02X F0F001%02X F0F00200", cutoff, reso);
		}
		else if (event->note == 2) // highpass filter
		{
			std::sprintf(str, "F0F000%02X F0F001%02X F0F00210", cutoff, reso);
		}
		else // no filter or unsupported filter type
		{
			std::strcpy(str, "F0F0007F F0F00100");
		}
	}
	// TODO: others
}

bool CSoundFile::ReadSymMOD(FileReader &file, ModLoadingFlags loadFlags)
//----------------------------------------------------------------------
{
	file.Rewind();
	if (!file.ReadMagic("SymM")
		|| file.ReadUint32BE() != 1)
		return false;

	if (loadFlags == onlyVerifyHeader)
		return true;

	InitializeGlobals(MOD_TYPE_SYMMOD);

	// channel panning
	for (CHANNELINDEX nChn = 0; nChn < MAX_BASECHANNELS; nChn++)
	{
		ChnSettings[nChn].nVolume = 64;
		ChnSettings[nChn].dwFlags.reset(CHN_SURROUND);
		ChnSettings[nChn].nPan = (nChn & 1) ? 256 : 0;
	}

	m_SongFlags.set(SONG_LINEARSLIDES);
	// TODO: find a fix for this, this flag gets set again when converting to IT internally so setting it here is useless
	m_playBehaviour.reset(kITPatternLoopTargetReset);

	// TODO: initialize other defaults here

	bool ok = true;

	// Data chunk types
	enum 
	{
		NumChannels     = -1,
		TrackLength     = -2,
		PatternSize     = -3,
		NumInstruments  = -4,
		EventSize       = -5,
		Tempo           = -6,
		ExternalSamples = -7,
		PositionList    = -10,
		SampleFile      = -11,
		EmptySample     = -12,
		PatternEvents   = -13,
		InstrumentList  = -14,
		Sequences       = -15,
		InfoText        = -16,
		SamplePacked    = -17,
		SamplePacked16  = -18,
		InfoType        = -19,
		InfoBinary      = -20,
		InfoString      = -21,
		
		SampleBoost     = 10,
		StereoDetune    = 11,
		StereoPhase     = 12,
	};

	uint32 trackLen = 0; 
	
	SymPosition *positions = 0;
	uint32 numPositions = 0;

	SymSequence *sequences = 0;
	uint32 numSequences = 0;

	SymEvent *patterns = 0;
	uint32 numEvents = 0;

	SymInstrument *instruments = 0;
	uint32 numInstruments = 0;

	// map Symphonie instruments to OpenMPT instruments and samples
	std::vector<SAMPLEINDEX> sampleMap;
	std::vector<INSTRUMENTINDEX> instrumentMap;

	while (ok && !file.EndOfFile())
	{
		int32 chunkType = file.ReadInt32BE();

		switch (chunkType)
		{
		/*
		Simple values
		*/
		case NumChannels:
			// TODO check against max
			m_nChannels = static_cast<CHANNELINDEX>(file.ReadUint32BE());

			if (m_nChannels)
				m_nSamplePreAmp = Clamp(512 / m_nChannels, 16, 128);
			break;

		case TrackLength:
			// TODO check against max
			trackLen = file.ReadUint32BE();
			break;

		case EventSize:
			ok = (file.ReadUint32BE() & 0xFFFF) <= sizeof(SymEvent);
			break;

		case Tempo:
			// not 100% accurate, but...
			m_nDefaultTempo = TEMPO(1.24 * file.ReadUint32BE());
			break;

		/*
		Unused values
		*/
		case NumInstruments: // determined from # of instrument headers instead
		case PatternSize:
		case ExternalSamples:
		case SampleBoost:
		case StereoDetune:
		case StereoPhase:
			file.Skip(4);
			break;

		/*
		Binary chunk types
		*/
		case PositionList:
			if (loadFlags & loadPatternData
				&& positions == 0)
			{
				numPositions = DecodeArray(file, positions);
				ok = numPositions > 0;
			}
			else
			{
				file.Skip(file.ReadUint32BE());
			}
			break;

		case SampleFile:
		case SamplePacked:
		case SamplePacked16:
			if ((loadFlags & loadSampleData) && m_nSamples < MAX_SAMPLES - 1u)
			{
				SAMPLEINDEX sample = GetNextFreeSample();

				if (sample != SAMPLEINDEX_INVALID)
				{
					m_nSamples++;
					sampleMap.push_back(sample);

					if (chunkType == SampleFile)
					{
						FileReader chunk = file.ReadChunk(file.ReadUint32BE());
						ok = ReadUnpackedSample(*this, Samples[sample], sample, chunk);
					}
					else if (chunkType == SamplePacked)
					{
						ok = DecodeSample8(*this, Samples[sample], sample, file);
					}
					else // SamplePacked16
					{
						ok = DecodeSample16(*this, Samples[sample], sample, file);
					}

					if (!ok) break;

					// Symphonie normalizes samples at load time
					NormalizeSample(Samples[sample]);

					// Symphonie represents stereo instruments as two consecutive mono instruments which are
					// automatically played at the same time. if this one uses a stereo sample, split it
					// and map two OpenMPT instruments to the stereo halves to ensure correct playback
					if (Samples[sample].uFlags.test(CHN_STEREO))
					{
						SAMPLEINDEX sampleR = GetNextFreeSample();

						if (sampleR != SAMPLEINDEX_INVALID)
						{ 
							m_nSamples++;
							ReadSampleFromSong(sampleR, *this, sample);
							ctrlSmp::ConvertToMono(Samples[sample], *this, ctrlSmp::onlyLeft);
							ctrlSmp::ConvertToMono(Samples[sampleR], *this, ctrlSmp::onlyRight);
							
							sampleMap.push_back(sampleR);
						}
						else
						{
							// just map the existing stereo sample again
							sampleMap.push_back(sample);
						}
					}
				}
				else
				{
					sampleMap.push_back(0);
				}
			}
			else
			{
				// skip sample
				sampleMap.push_back(0);
				file.Skip(file.ReadUint32BE());
			}
			break;

		case EmptySample:
			sampleMap.push_back(0);
			break;

		case PatternEvents:
			if (loadFlags & loadPatternData
				&& patterns == 0)
			{
				numEvents = DecodeArray(file, patterns);
				ok = numEvents > 0;
			}
			else
			{
				file.Skip(file.ReadUint32BE());
			}
			break;

		case InstrumentList:
			if (loadFlags & loadSampleData
				&& instruments == 0)
			{
				numInstruments = DecodeArray(file, instruments);
				ok = numInstruments > 0;
			}
			else
			{
				file.Skip(file.ReadUint32BE());
			}
			break;

		case Sequences:
			if (loadFlags & loadPatternData
				&& sequences == 0)
			{
				numSequences = DecodeArray(file, sequences);
				ok = numSequences > 0;
			}
			else
			{
				file.Skip(file.ReadUint32BE());
			}
			break;

		case InfoText:
		{
			mpt::byte *text;
			size_t len = DecodeChunk(file, text);

			ok = m_songMessage.Read(text, len, SongMessage::leLF);

			delete[] text;
			break;
		}

		/*
		Unused binary chunks
		*/
		case InfoType:
		case InfoBinary:
		case InfoString:
			file.Skip(file.ReadUint32BE());
			break;

		/*
		Unrecognized chunk/value type
		*/
		default:
			ok = false;
			break;
		}
	}

	// convert instruments
	if (ok && (loadFlags & loadSampleData))
	{
		for (INSTRUMENTINDEX nInst = 0; nInst < numInstruments; nInst++)
		{
			SymInstrument& symInst = instruments[nInst];

			// don't convert empty instruments (all valid instruments have a name)
			if (!symInst.name[0] || symInst.type < 0)
			{
				instrumentMap.push_back(0);
				continue;
			}
			INSTRUMENTINDEX index = GetNextFreeInstrument();
			// create this instrument and assign a sample
			if (index != INSTRUMENTINDEX_INVALID && nInst < sampleMap.size())
			{
				instrumentMap.push_back(index);
				m_nInstruments = std::max<uint16>(m_nInstruments, index);

				SAMPLEINDEX sampleIndex = sampleMap[nInst];
				if (symInst.IsVirtual())
				{
					// "Virtual instruments" allow an arbitrary number of other instruments to be mixed
					// at various pitches/volume levels (to create layered samples, chords, etc.)
					// Theoretically, we could handle most of these by doing some basic mixing here,
					// but for now, let's just use the first base instrument on the list
					SymEvent& event = symInst.virt.events[0];
					INSTRUMENTINDEX nBaseInst = event.inst;

					// we should be able to assume that all the base instruments will come before this one
					if (nBaseInst < instrumentMap.size())
					{
						SymInstrument &baseInst = instruments[nBaseInst];
						SAMPLEINDEX baseSample = GetNextFreeSample();
						if (baseSample != SAMPLEINDEX_INVALID)
						{
							m_nSamples++;
							ReadSampleFromSong(baseSample, *this, sampleMap[nBaseInst]);

							sampleMap[nInst] = sampleIndex = baseSample;
						}
						else
						{
							sampleIndex = sampleMap[nBaseInst];
						}

						// correct pitch
						symInst.downsample += baseInst.downsample;
					}
				}

				ModInstrument *inst = AllocateInstrument(index, sampleIndex);
				ModSample &sample = Samples[sampleIndex];

				symInst.ConvertToMPT(*inst, sample, *this);
			}
			else
			{
				instrumentMap.push_back(0);
			}
		}
	}

	// convert patterns
	if (ok && (loadFlags & loadPatternData) && trackLen > 0)
	{
		// map Symphonie positions to converted patterns
		std::map<SymPosition, PATTERNINDEX> patternMap;
		// map DSP commands to MIDI macro numbers
		std::map<SymEvent, uint8> macroMap;

		const uint32 patternSize = m_nChannels * trackLen;
		PATTERNINDEX numPatterns = mpt::saturate_cast<PATTERNINDEX>(numEvents / patternSize);
		
		Patterns.ResizeArray(numPatterns);
		Order.clear();

#ifdef MPT_SYMMOD_USE_REAL_SUBSONGS
		Order.SetSequence(0);
#endif // MPT_SYMMOD_USE_REAL_SUBSONGS

		for (uint32 nSeq = 0; nSeq < numSequences; nSeq++)
		{
			SymSequence& seq = sequences[nSeq];

			if (seq.info == 1)
				continue;
			if (seq.info == -1)
				break;

			// separate sequences
#ifdef MPT_SYMMOD_USE_REAL_SUBSONGS
			SEQUENCEINDEX seqIndex = Order.AddSequence(false);
			if (seqIndex == SEQUENCEINDEX_INVALID)
				break;
			Order.clear();
			ModSequence &order = Order.GetSequence(static_cast<SEQUENCEINDEX>(seqIndex));
#else
			ModSequence &order = Order;
			if (order.GetLength() > 0)
				order.Append(ModSequence::GetInvalidPatIndex());
#endif // MPT_SYMMOD_USE_REAL_SUBSONGS

			// last note played on a channel
			uint8 lastNote[MAX_CHANNELS] = { 0 };
			// last instrument played on a channel
			int8 lastInst[MAX_CHANNELS] = { 0 };
			// last specified volume of a channel (to avoid excessive Mxx commands)
			uint8 lastVol[MAX_CHANNELS] = { 0 };
			// current volume slide factor of a channel
			double currVolSlide[MAX_CHANNELS] = { 0.0 };
			// cumulative volume slide amount
			double currVolSlideAmt[MAX_CHANNELS] = { 0.0 };
			// current pitch slide factor of a channel
			double currPitchSlide[MAX_CHANNELS] = { 0.0 };
			// cumulative pitch slide amount
			double currPitchSlideAmt[MAX_CHANNELS] = { 0.0 };
			// sample paused or not (affects volume and pitch slides)
			bool chnStopped[MAX_CHANNELS] = { 0 };
			// vibrato and tremolo
			uint8 currVibrato[MAX_CHANNELS] = { 0 };
			uint8 currTremolo[MAX_CHANNELS] = { 0 };

			for (uint16be nPos = seq.start; nPos < seq.start + seq.length && nPos < numPositions; nPos++)
			{
				SymPosition& pos = positions[nPos];
				PATTERNINDEX patternIndex;

				// before checking the map, apply the sequence transpose value
				pos.transpose += seq.transpose;

				// pattern already converted?
				if (patternMap.count(pos))
				{
					patternIndex = patternMap[pos];
				}
				else
				{
					// convert pattern now
					patternIndex = Patterns.InsertAny(pos.length);
					if (patternIndex == PATTERNINDEX_INVALID) 
						break;

					patternMap[pos] = patternIndex;

					// actually convert the pattern
					if (pos.pattern >= numPatterns)
						continue;

					uint8 patternSpeed = static_cast<uint8>(pos.speed);

					SymEvent *srcEvent = patterns + (pos.pattern * patternSize) + (pos.start * m_nChannels);
					for (ROWINDEX nRow = 0; nRow < pos.length; nRow++)
					{
						PatternRow rowBase = Patterns[patternIndex].GetpModCommand(nRow, 0);

						for (CHANNELINDEX chn = 0; chn < m_nChannels; chn++)
						{
							ModCommand &m = rowBase[chn];

							int8 note = srcEvent->note + 24; // TODO: 24 or 25? (see TODO about sample tuning)
							uint8 origInst = srcEvent->inst;
							uint8 mappedInst = 0;
							if (origInst < numInstruments)
							{
								mappedInst = static_cast<uint8>(instrumentMap[origInst]);
								if (!(instruments[origInst].instFlags & SymInstrument::NoTranspose))
									note += static_cast<int8>(pos.transpose);
							}


							switch (srcEvent->command)
							{
							case SymEvent::KeyOn:
								if (srcEvent->param > SymEvent::VolCommand)
								{
									switch (srcEvent->param)
									{
									case SymEvent::StopSample:
										m.volcmd = VOLCMD_VOLUME;
										m.vol = 0;
										chnStopped[chn] = true;
										break;

									case SymEvent::ContSample:
										m.volcmd = VOLCMD_VOLUME;
										m.vol = 64;
										chnStopped[chn] = false;
										break;

									case SymEvent::KeyOff:
										m.volcmd = VOLCMD_OFFSET;
										m.vol = 1;
										break;

									case SymEvent::SpeedDown:
										if (patternSpeed > 1)
										{
											m.command = CMD_SPEED;
											m.param = --patternSpeed;
										}
										break;

									case SymEvent::SpeedUp:
										if (patternSpeed < 0xFF)
										{
											m.command = CMD_SPEED;
											m.param = ++patternSpeed;
										}
										break;

									case SymEvent::SetPitch:
										m.note = lastNote[chn] = note;
										// some modules (like "endless ways") use this even if a prior note hasn't played yet
										// so be sure to set instrument too (but not if inst is 0?)
										if (srcEvent->inst)
											m.instr = lastInst[chn] = mappedInst;
										m.command = CMD_TONEPORTAMENTO;
										m.param = 0xFF;
										currPitchSlide[chn] = 0.0;
										break;

									// fine portamentos with range up to half a semitone
									case SymEvent::PitchUp:
										m.command = CMD_PORTAMENTOUP;
										m.param = 0xF2;
										break;
									case SymEvent::PitchDown:
										m.command = CMD_PORTAMENTODOWN;
										m.param = 0xF2;
										break;
									case SymEvent::PitchUp2:
										m.command = CMD_PORTAMENTOUP;
										m.param = 0xF4;
										break;
									case SymEvent::PitchDown2:
										m.command = CMD_PORTAMENTODOWN;
										m.param = 0xF4;
										break;
									case SymEvent::PitchUp3:
										m.command = CMD_PORTAMENTOUP;
										m.param = 0xF8;
										break;
									case SymEvent::PitchDown3:
										m.command = CMD_PORTAMENTODOWN;
										m.param = 0xF8;
										break;
									}
								}
								else
								{
									if (srcEvent->note >= 0 || srcEvent->param < 100)
									{
										if (srcEvent->note >= 0)
										{
											m.note = lastNote[chn] = note;
											m.instr = lastInst[chn] = mappedInst;
											currPitchSlide[chn] = 0.0;
										}

										if (srcEvent->param > 0)
										{
											uint8 newVol = Util::Round<uint8>(srcEvent->param * 0.64);

											if (lastVol[chn] != newVol || currVolSlide[chn] != 0.0)
											{
												m.command = CMD_CHANNELVOLUME;
												m.param = lastVol[chn] = newVol;
											}
											currVolSlide[chn] = 0.0;
										}
									}
								}

								// key on commands with stereo instruments are played on both channels
								if (srcEvent->note > 0
									&& (instruments[origInst].channels == SymInstrument::StereoL
										|| instruments[origInst].instFlags & SymInstrument::SyncPlay)
									&& (chn < m_nChannels - 1))
								{
									ModCommand &next = rowBase[chn + 1];
									next = m;
									if (instruments[origInst].channels == SymInstrument::StereoL)
										next.instr++;

									lastVol[chn + 1] = lastVol[chn];
									currVolSlide[chn + 1] = currVolSlide[chn];
									currVolSlideAmt[chn + 1] = currVolSlideAmt[chn];
									currPitchSlide[chn + 1] = currPitchSlide[chn];
									currPitchSlideAmt[chn + 1] = currPitchSlideAmt[chn];
								}

								break;

								// volume effects
								// Symphonie has very fine fractional volume slides which are applied at the output sample rate,
								// rather than per tick or per row, so instead let's simulate it based on the pattern speed
								// by keeping track of the volume and using normal volume commands
								// the math here is an approximation which works fine for most songs
							case SymEvent::VolSlideUp:
								currVolSlideAmt[chn] = 0.0;
								currVolSlide[chn] = static_cast<double>(srcEvent->param) * 0.0333;
								break;
							case SymEvent::VolSlideDown:
								currVolSlideAmt[chn] = 0.0;
								currVolSlide[chn] = static_cast<double>(srcEvent->param) * -0.0333;
								break;

							case SymEvent::AddVolume:
								m.command = m.param = 0;
								break;
							case SymEvent::Tremolo:
							{
								// both tremolo speed and depth can go much higher than OpenMPT supports,
								// but modules will probably use pretty sane, supportable values anyway
								// TODO: handle very small nonzero params
								uint8 speed = std::min<uint8>(15, srcEvent->param >> 3);
								uint8 depth = std::min<uint8>(15, srcEvent->inst >> 3);

								currTremolo[chn] = (speed << 4) | depth;
							}
								break;


								// pitch effects
								// Pitch slides have a similar granularity to volume slides, and are approximated
								// the same way here based on a rough comparison against Exx/Fxx slides
							case SymEvent::PitchSlideUp:
								currPitchSlideAmt[chn] = 0.0;
								currPitchSlide[chn] = static_cast<double>(srcEvent->param) * 0.0333;
								break;
							case SymEvent::PitchSlideDown:
								currPitchSlideAmt[chn] = 0.0;
								currPitchSlide[chn] = static_cast<double>(srcEvent->param) * -0.0333;
								break;

							case SymEvent::PitchSlideTo:
								m.command = m.param = 0;
								break;
							case SymEvent::AddPitch:
								/* "The range (-128...127) is about 4 half notes." */
								m.command = m.param = 0;
								break;
							case SymEvent::Vibrato:
							{
								// both vibrato speed and depth can go much higher than OpenMPT supports,
								// but modules will probably use pretty sane, supportable values anyway
								// TODO: handle very small nonzero params
								uint8 speed = std::min<uint8>(15, srcEvent->param >> 3);
								uint8 depth = std::min<uint8>(15, srcEvent->inst);

								currVibrato[chn] = (speed << 4) | depth;
							}
								break;
							case SymEvent::AddHalfTone:
								m.note = (lastNote[chn] += srcEvent->param);
								m.command = CMD_TONEPORTAMENTO;
								m.param = 0xFF;
								break;

								// DSP effects
							case SymEvent::Filter:
#ifndef NO_PLUGINS
							case SymEvent::DSPEcho:
							case SymEvent::DSPDelay:
#endif
								if (macroMap.count(*srcEvent))
								{
									m.command = CMD_MIDI;
									m.param = macroMap[*srcEvent];
								}
								else if (macroMap.size() < 128)
								{
									uint8 param = static_cast<uint8>(macroMap.size());

									ConvertDSP(srcEvent, m_MidiCfg.szMidiZXXExt[param]);

									m.command = CMD_MIDI;
									m.param = macroMap[*srcEvent] = 0x80 | param;
								}
								break;

								// other effects
							case SymEvent::Retrig: // TODO: this takes 2 params apparently?
								m.command = CMD_RETRIG;
								m.param = std::min<uint8>(15, srcEvent->param);
								break;

							case SymEvent::SetSpeed:
								m.command = CMD_SPEED;
								m.param = patternSpeed = srcEvent->param;
								break;

							case SymEvent::CV:
								m.command = m.param = 0;
								break;
							case SymEvent::CVAdd:
								m.command = m.param = 0;
								break;

								// sample effects
							case SymEvent::FromAndPitch:
								lastNote[chn] = note;
								m.instr = lastInst[chn] = mappedInst;
								// fall through
							case SymEvent::ReplayFrom:
								m.note = lastNote[chn];
								// don't always add the command, because often FromAndPitch is used with offset 0
								// to act as a key-on which doesn't cancel volume slides, etc
								if (srcEvent->param > 0)
								{
									m.command = CMD_OFFSETPERCENTAGE;
									m.param = srcEvent->param;
								}
								break;

							}

							srcEvent++;

							// any event which plays a note should re-enable continuous effects
							if (m.note > 0)
								chnStopped[chn] = false;
							else if (chnStopped[chn])
								continue;

							// handle fractional volume slides
							if (currVolSlide[chn] != 0.0)
							{
								currVolSlideAmt[chn] += currVolSlide[chn] * patternSpeed;
								if (!m.command)
								{
									if (patternSpeed > 1 && currVolSlideAmt[chn] >= (patternSpeed - 1))
									{
										uint8 slideAmt = std::min<uint8>(15, static_cast<uint8>(currVolSlideAmt[chn] / (patternSpeed - 1)));
										currVolSlideAmt[chn] -= slideAmt * (patternSpeed - 1);
										// normal slide up
										m.command = CMD_CHANNELVOLSLIDE;
										m.param = slideAmt << 4;
									}
									else if (currVolSlideAmt[chn] >= 1.0)
									{
										uint8 slideAmt = std::min<uint8>(15, static_cast<uint8>(currVolSlideAmt[chn]));
										currVolSlideAmt[chn] -= slideAmt;
										// fine slide up
										m.command = CMD_CHANNELVOLSLIDE;
										m.param = (slideAmt << 4) | 0x0F;
									}
									else if (patternSpeed > 1 && currVolSlideAmt[chn] <= -(patternSpeed - 1))
									{
										uint8 slideAmt = std::min<uint8>(15, static_cast<uint8>(-currVolSlideAmt[chn] / (patternSpeed - 1)));
										currVolSlideAmt[chn] += slideAmt * (patternSpeed - 1);
										// normal slide down
										m.command = CMD_CHANNELVOLSLIDE;
										m.param = slideAmt;
									}
									else if (currVolSlideAmt[chn] <= -1.0)
									{
										uint8 slideAmt = std::min<uint8>(14, static_cast<uint8>(-currVolSlideAmt[chn]));
										currVolSlideAmt[chn] += slideAmt;
										// fine slide down
										m.command = CMD_CHANNELVOLSLIDE;
										m.param = slideAmt | 0xF0;
									}
								}
							}
							// handle fractional pitch slides
							if (currPitchSlide[chn] != 0.0)
							{
								currPitchSlideAmt[chn] += currPitchSlide[chn] * patternSpeed;
								if (!m.command)
								{
									if (patternSpeed > 1 && currPitchSlideAmt[chn] >= (patternSpeed - 1))
									{
										uint8 slideAmt = std::min<uint8>(0xDF, static_cast<uint8>(currPitchSlideAmt[chn] / (patternSpeed - 1)));
										currPitchSlideAmt[chn] -= slideAmt * (patternSpeed - 1);
										// normal slide up
										m.command = CMD_PORTAMENTOUP;
										m.param = slideAmt;
									}
									else if (currPitchSlideAmt[chn] >= 1.0)
									{
										uint8 slideAmt = std::min<uint8>(15, static_cast<uint8>(currPitchSlideAmt[chn]));
										currPitchSlideAmt[chn] -= slideAmt;
										// fine slide up
										m.command = CMD_PORTAMENTOUP;
										m.param = slideAmt | 0xF0;
									}
									else if (patternSpeed > 1 && currPitchSlideAmt[chn] <= -(patternSpeed - 1))
									{
										uint8 slideAmt = std::min<uint8>(0xDF, static_cast<uint8>(-currPitchSlideAmt[chn] / (patternSpeed - 1)));
										currPitchSlideAmt[chn] += slideAmt * (patternSpeed - 1);
										// normal slide down
										m.command = CMD_PORTAMENTODOWN;
										m.param = slideAmt;
									}
									else if (currPitchSlideAmt[chn] <= -1.0)
									{
										uint8 slideAmt = std::min<uint8>(14, static_cast<uint8>(-currPitchSlideAmt[chn]));
										currPitchSlideAmt[chn] += slideAmt;
										// fine slide down
										m.command = CMD_PORTAMENTODOWN;
										m.param = slideAmt | 0xF0;
									}
								}
								// TODO: use volume column if effect column is occupied
								else if (!m.volcmd)
								{
									if (patternSpeed > 1 && currPitchSlideAmt[chn] / 4 >= (patternSpeed - 1))
									{
										uint8 slideAmt = std::min<uint8>(9, static_cast<uint8>(currPitchSlideAmt[chn] / (patternSpeed - 1)) / 4);
										currPitchSlideAmt[chn] -= slideAmt * (patternSpeed - 1) * 4;
										m.volcmd = VOLCMD_PORTAUP;
										m.vol = slideAmt;
									}
									else if (patternSpeed > 1 && currPitchSlideAmt[chn] / 4 <= -(patternSpeed - 1))
									{
										uint8 slideAmt = std::min<uint8>(9, static_cast<uint8>(-currPitchSlideAmt[chn] / (patternSpeed - 1)) / 4);
										currPitchSlideAmt[chn] += slideAmt * (patternSpeed - 1) * 4;
										m.volcmd = VOLCMD_PORTADOWN;
										m.vol = slideAmt;
									}
								}
							}
							// vibrato and tremolo
							if (m.command == 0 && currVibrato[chn] != 0)
							{
								m.command = CMD_VIBRATO;
								m.param = currVibrato[chn];
							}
							if (m.command == 0 && currTremolo[chn] != 0)
							{
								m.command = CMD_TREMOLO;
								m.param = currTremolo[chn];
							}
						}
					}

					// write speed and loop commands
					Patterns[patternIndex].WriteEffect(EffectWriter(CMD_SPEED, static_cast<uint8>(pos.speed)).Row(0));

					if (pos.loopNum > 1)
					{
						Patterns[patternIndex].WriteEffect(EffectWriter(CMD_S3MCMDEX, 0xB0 + std::min(15, pos.loopNum - 1))
							.Row(static_cast<ROWINDEX>(pos.length - 1)));
					}
				}
				/*
				for (int i = 0; i < pos.loopNum; i++)
					*/
					order.Append(patternIndex);
				// undo transpose tweak
				pos.transpose -= seq.transpose;
			}
		}
	}

	delete[] positions;
	delete[] sequences;
	delete[] instruments;
	delete[] patterns;

	return ok;
}

OPENMPT_NAMESPACE_END
