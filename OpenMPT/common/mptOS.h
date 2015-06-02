/*
 * mptOS.h
 * -------
 * Purpose: Operating system version information.
 * Notes  : (currently none)
 * Authors: OpenMPT Devs
 * The OpenMPT source code is released under the BSD license. Read LICENSE for more details.
 */


#pragma once


OPENMPT_NAMESPACE_BEGIN


namespace mpt
{
namespace Windows
{
namespace Version
{


enum Number
{

	// Win9x
	Win98    = 0x040a,
	WinME    = 0x045a,

	// WinNT
	WinNT4   = 0x0400,
	Win2000  = 0x0500,
	WinXP    = 0x0501,
	WinVista = 0x0600,
	Win7     = 0x0601,
	Win8     = 0x0602,
	Win81    = 0x0603,

};


#if defined(MODPLUG_TRACKER)

#if MPT_OS_WINDOWS

// Initializes version information.
// Version information is cached in order to be safely available in audio real-time contexts.
// Checking version information (especially detecting Wine) can be slow.
void Init();

static inline bool IsWindows() { return true; }

bool IsBefore(mpt::Windows::Version::Number version);
bool IsAtLeast(mpt::Windows::Version::Number version);

bool Is9x();
bool IsNT();

bool IsOriginal();
bool IsWine();

#else // !MPT_OS_WINDOWS

static inline void Init() { return; }

static inline bool IsWindows() { return false; }

static inline bool IsBefore(mpt::Windows::Version::Number /* version */ ) { return false; }
static inline bool IsAtLeast(mpt::Windows::Version::Number /* version */ ) { return false; }

static inline bool Is9x() { return false; }
static inline bool IsNT() { return false; }

static inline bool IsOriginal() { return false; }
static inline bool IsWine() { return false; }

#endif // MPT_OS_WINDOWS

#else // !MODPLUG_TRACKER

#if MPT_OS_WINDOWS

static inline bool IsWindows() { return true; }

bool IsBefore(mpt::Windows::Version::Number version);
bool IsAtLeast(mpt::Windows::Version::Number version);

#else // !MPT_OS_WINDOWS

static inline bool IsWindows() { return false; }

static inline bool IsBefore(mpt::Windows::Version::Number /* version */ ) { return false; }
static inline bool IsAtLeast(mpt::Windows::Version::Number /* version */ ) { return false; }

#endif // MPT_OS_WINDOWS

#endif // MODPLUG_TRACKER


} // namespace Version
} // namespace Windows
} // namespace mpt


OPENMPT_NAMESPACE_END
