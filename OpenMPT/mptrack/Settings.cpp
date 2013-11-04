/*
 * Settings.cpp
 * ------------
 * Purpose: Application setting handling framework.
 * Notes  : (currently none)
 * Authors: Joern Heusipp
 *          OpenMPT Devs
 * The OpenMPT source code is released under the BSD license. Read LICENSE for more details.
 */


#include "stdafx.h"

#include "Settings.h"

#include "../common/misc_util.h"
#include "../common/StringFixer.h"
#include "Mptrack.h"

#include <shlwapi.h>

#include <algorithm>
#include <fstream>
#include <iterator>


static const wchar_t EncodeNibble[16] = { L'0', L'1', L'2', L'3', L'4', L'5', L'6', L'7', L'8', L'9', L'A', L'B', L'C', L'D', L'E', L'F' };

static inline bool DecodeByte(uint8 &byte, wchar_t c1, wchar_t c2)
{
	byte = 0;
	if(L'0' <= c1 && c1 <= L'9')
	{
		byte += static_cast<uint8>((c1 - L'0') << 4);
	} else if(L'A' <= c1 && c1 <= L'F')
	{
		byte += static_cast<uint8>((c1 - L'A' + 10) << 4);
	} else
	{
		return false;
	}
	if(L'0' <= c2 && c2 <= L'9')
	{
		byte += static_cast<uint8>(c2 - L'0');
	} else if(L'A' <= c2 && c2 <= L'F')
	{
		byte += static_cast<uint8>(c2 - L'A' + 10);
	} else
	{
		return false;
	}
	return true;
}


std::wstring SettingBinToHex(const std::vector<char> &src)
{
	std::wstring result;
	for(std::size_t i = 0; i < src.size(); ++i)
	{
		uint8 byte = src[i];
		result.push_back(EncodeNibble[(byte&0xf0)>>4]);
		result.push_back(EncodeNibble[byte&0x0f]);
	}
	return result;
}

std::vector<char> SettingHexToBin(const std::wstring &src)
{
	std::vector<char> result;
	for(std::size_t i = 0; i+1 < src.size(); i += 2)
	{
		uint8 byte = 0;
		if(!DecodeByte(byte, src[i*2+0], src[i*2+1]))
		{
			return result;
		}
		result.push_back(byte);
	}
	return result;
}



std::wstring SettingValue::FormatTypeAsString() const
{
	if(GetType() == SettingTypeNone)
	{
		return L"nil";
	}
	std::wstring result;
	switch(GetType())
	{
		case SettingTypeBool:
			result += L"bool";
			break;
		case SettingTypeInt:
			result += L"int";
			break;
		case SettingTypeFloat:
			result += L"float";
			break;
		case SettingTypeString:
			result += L"string";
			break;
		case SettingTypeBinary:
			result += L"binary";
			break;
		case SettingTypeNone:
		default:
			result += L"nil";
			break;
	}
	if(HasTypeTag() && !GetTypeTag().empty())
	{
		result += L":" + mpt::String::Decode(GetTypeTag(), mpt::CharsetUS_ASCII);
	}
	return result;
}


std::wstring SettingValue::FormatValueAsString() const
{
	switch(GetType())
	{
		case SettingTypeBool:
			return StringifyW(valueBool);
			break;
		case SettingTypeInt:
			return StringifyW(valueInt);
			break;
		case SettingTypeFloat:
			return StringifyW(valueFloat);
			break;
		case SettingTypeString:
			return valueString;
			break;
		case SettingTypeBinary:
			return SettingBinToHex(valueBinary);
			break;
		case SettingTypeNone:
		default:
			return std::wstring();
			break;
	}
}


std::wstring SettingValue::FormatAsString() const
{
	return L"(" + FormatTypeAsString() + L")" + FormatValueAsString();
}


void SettingValue::SetFromString(const std::wstring &newVal)
{
	switch(GetType())
	{
		case SettingTypeBool:
			valueBool = ConvertStrTo<bool>(newVal);
			break;
		case SettingTypeInt:
			valueInt = ConvertStrTo<int32>(newVal);
			break;
		case SettingTypeFloat:
			valueFloat = ConvertStrTo<double>(newVal);
			break;
		case SettingTypeString:
			valueString = newVal;
			break;
		case SettingTypeBinary:
			valueBinary = SettingHexToBin(newVal);
			break;
		case SettingTypeNone:
		default:
			break;
	}
}



#if defined(MPT_SETTINGS_CACHE)


SettingValue SettingsContainer::BackendsReadSetting(const SettingPath &path, const SettingValue &def) const
{
	return backend->ReadSetting(path, def);
}

void SettingsContainer::BackendsWriteSetting(const SettingPath &path, const SettingValue &val)
{
	backend->WriteSetting(path, val);
}

void SettingsContainer::BackendsRemoveSetting(const SettingPath &path)
{
	backend->RemoveSetting(path);
}

SettingValue SettingsContainer::ReadSetting(const SettingPath &path, const SettingValue &def, const SettingMetadata &metadata) const
{
	SettingsMap::iterator entry = map.find(path);
	if(entry == map.end())
	{
		entry = map.insert(map.begin(), std::make_pair(path, SettingState(def).assign(BackendsReadSetting(path, def), false)));
		mapMetadata[path] = metadata;
	}
	return entry->second;
}

void SettingsContainer::WriteSetting(const SettingPath &path, const SettingValue &val)
{
	SettingsMap::iterator entry = map.find(path);
	if(entry == map.end())
	{
		map[path] = val;
		entry = map.find(path);
	} else
	{
		entry->second = val;
	}
	NotifyListeners(path);
	if(immediateFlush)
	{
		BackendsWriteSetting(path, val);
		entry->second.Clean();
	}
}

void SettingsContainer::RemoveSetting(const SettingPath &path)
{
	map.erase(path);
	BackendsRemoveSetting(path);
}

void SettingsContainer::NotifyListeners(const SettingPath &path)
{
	const SettingsListenerMap::iterator entry = mapListeners.find(path);
	if(entry != mapListeners.end())
	{
		const std::set<ISettingChanged*>::const_iterator beg = entry->second.begin();
		const std::set<ISettingChanged*>::const_iterator end = entry->second.end();
		for(std::set<ISettingChanged*>::const_iterator it = beg; it != end; ++it)
		{
			(*it)->SettingChanged(path);
		}
	}
}

void SettingsContainer::WriteSettings()
{
	for(SettingsMap::iterator i = map.begin(); i != map.end(); ++i)
	{
		if(i->second.IsDirty())
		{
			BackendsWriteSetting(i->first, i->second);
			i->second.Clean();
		}
	}
}

void SettingsContainer::Flush()
{
	WriteSettings();
}

void SettingsContainer::SetImmediateFlush(bool newImmediateFlush)
{
	if(newImmediateFlush)
	{
		Flush();
	}
	immediateFlush = newImmediateFlush;
}

void SettingsContainer::Register(ISettingChanged *listener, const SettingPath &path)
{
	mapListeners[path].insert(listener);
}

void SettingsContainer::UnRegister(ISettingChanged *listener, const SettingPath &path)
{
	mapListeners[path].erase(listener);
}

SettingsContainer::~SettingsContainer()
{
	WriteSettings();
}

#else // !MPT_SETTINGS_CACHE

SettingValue SettingsContainer::ReadSetting(const SettingPath &path, const SettingValue &def, const SettingMetadata & /*metadata*/ ) const
{
	return backend->ReadSetting(path, def);
}

void SettingsContainer::WriteSetting(const SettingPath &path, const SettingValue &val)
{
	backend->WriteSetting(path, val);
}

void SettingsContainer::RemoveSetting(const SettingPath &path)
{
	backend->RemoveSetting(path);
}

void SettingsContainer::Flush()
{
	return;
}

SettingsContainer::~SettingsContainer()
{
	return;
}

#endif // MPT_SETTINGS_CACHE


SettingsContainer::SettingsContainer(ISettingsBackend *backend)
	: backend(backend)
	, immediateFlush(false)
{
	ASSERT(backend);
}




std::vector<char> IniFileSettingsBackend::ReadSettingRaw(const SettingPath &path, const std::vector<char> &def) const
{
	std::vector<char> result = def;
	::GetPrivateProfileStructW(GetSection(path).c_str(), GetKey(path).c_str(), &result[0], result.size(), filename.c_str());
	return result;
}

std::wstring IniFileSettingsBackend::ReadSettingRaw(const SettingPath &path, const std::wstring &def) const
{
	std::vector<WCHAR> buf(128);
	while(::GetPrivateProfileStringW(GetSection(path).c_str(), GetKey(path).c_str(), def.c_str(), &buf[0], buf.size(), filename.c_str()) == buf.size() - 1)
	{
		buf.resize(buf.size() * 2);
	}
	return &buf[0];
}

double IniFileSettingsBackend::ReadSettingRaw(const SettingPath &path, double def) const
{
	std::vector<WCHAR> buf(128);
	while(::GetPrivateProfileStringW(GetSection(path).c_str(), GetKey(path).c_str(), StringifyW(def).c_str(), &buf[0], buf.size(), filename.c_str()) == buf.size() - 1)
	{
		buf.resize(buf.size() * 2);
	}
	return ConvertStrTo<double>(std::wstring(&buf[0]));
}

int32 IniFileSettingsBackend::ReadSettingRaw(const SettingPath &path, int32 def) const
{
	return (int32)::GetPrivateProfileIntW(GetSection(path).c_str(), GetKey(path).c_str(), (UINT)def, filename.c_str());
}

bool IniFileSettingsBackend::ReadSettingRaw(const SettingPath &path, bool def) const
{
	return ::GetPrivateProfileIntW(GetSection(path).c_str(), GetKey(path).c_str(), def?1:0, filename.c_str()) ? true : false;
}


void IniFileSettingsBackend::WriteSettingRaw(const SettingPath &path, const std::vector<char> &val)
{
	::WritePrivateProfileStructW(GetSection(path).c_str(), GetKey(path).c_str(), (LPVOID)&val[0], val.size(), filename.c_str());
}

void IniFileSettingsBackend::WriteSettingRaw(const SettingPath &path, const std::wstring &val)
{
	::WritePrivateProfileStringW(GetSection(path).c_str(), GetKey(path).c_str(), val.c_str(), filename.c_str());

	if(mpt::String::Decode(mpt::String::Encode(val, mpt::CharsetLocale), mpt::CharsetLocale) != val)
	{
		// Value is not representable in ANSI CP.
		// Now check if the string got stored correctly.
		if(ReadSettingRaw(path, std::wstring()) != val)
		{
			// The ini file is probably ANSI encoded.
			ConvertToUnicode();
			// Re-write non-ansi-representable value.
			::WritePrivateProfileStringW(GetSection(path).c_str(), GetKey(path).c_str(), val.c_str(), filename.c_str());
		}
	}
}

void IniFileSettingsBackend::WriteSettingRaw(const SettingPath &path, double val)
{
	::WritePrivateProfileStringW(GetSection(path).c_str(), GetKey(path).c_str(), StringifyW(val).c_str(), filename.c_str());
}

void IniFileSettingsBackend::WriteSettingRaw(const SettingPath &path, int32 val)
{
	::WritePrivateProfileStringW(GetSection(path).c_str(), GetKey(path).c_str(), StringifyW(val).c_str(), filename.c_str());
}

void IniFileSettingsBackend::WriteSettingRaw(const SettingPath &path, bool val)
{
	::WritePrivateProfileStringW(GetSection(path).c_str(), GetKey(path).c_str(), StringifyW(val?1:0).c_str(), filename.c_str());
}

void IniFileSettingsBackend::RemoveSettingRaw(const SettingPath &path)
{
	::WritePrivateProfileStringW(GetSection(path).c_str(), GetKey(path).c_str(), NULL, filename.c_str());
}


std::wstring IniFileSettingsBackend::GetSection(const SettingPath &path)
{
	return mpt::String::Decode(path.GetSection(), mpt::CharsetLocale);
}
std::wstring IniFileSettingsBackend::GetKey(const SettingPath &path)
{
	return mpt::String::Decode(path.GetKey(), mpt::CharsetLocale);
}



IniFileSettingsBackend::IniFileSettingsBackend(const mpt::PathString &filename)
	: filename(filename)
{
	return;
}

IniFileSettingsBackend::~IniFileSettingsBackend()
{
	return;
}

static std::vector<char> ReadFile(const mpt::PathString &filename)
{
	std::ifstream s(filename.c_str(), std::ios::binary);
	std::vector<char> result;
	while(s)
	{
		char buf[4096];
		s.read(buf, 4096);
		std::streamsize count = s.gcount();
		std::copy(buf, buf + count, std::back_inserter(result));
	}
	return result;
}

static void WriteFileUTF16LE(const mpt::PathString &filename, const std::wstring &str)
{
	STATIC_ASSERT(sizeof(wchar_t) == 2);
	std::ofstream inifile(filename.c_str(), std::ios::binary | std::ios::trunc);
	const uint8 UTF16LE_BOM[] = { 0xff, 0xfe };
	inifile.write(reinterpret_cast<const char*>(UTF16LE_BOM), 2);
	inifile.write(reinterpret_cast<const char*>(str.c_str()), str.length() * sizeof(std::wstring::value_type));
	inifile.flush();
	inifile.close();
}

void IniFileSettingsBackend::ConvertToUnicode(const std::wstring &backupTag)
{
	// Force ini file to be encoded in UTF16.
	// This causes WINAPI ini file functions to keep it in UTF16 encoding
	// and thus support storing unicode strings uncorrupted.
	// This is backwards compatible because even ANSI WINAPI behaves the
	// same way in this case.
	const std::vector<char> data = ReadFile(filename);
	if(!data.empty() && IsTextUnicode(&data[0], data.size(), NULL))
	{
		return;
	}
	const std::wstring backupFilename = backupTag.empty()
		? filename + L".ansi.bak"
		: filename + L".ansi." + backupTag + L".bak"
		;
	MoveFileExW(filename.c_str(), backupFilename.c_str(), MOVEFILE_COPY_ALLOWED|MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH);
	WriteFileUTF16LE(filename, mpt::String::Decode(std::string(&data[0], &data[0] + data.size()), mpt::CharsetLocale));
}

SettingValue IniFileSettingsBackend::ReadSetting(const SettingPath &path, const SettingValue &def) const
{
	switch(def.GetType())
	{
	case SettingTypeBool: return SettingValue(ReadSettingRaw(path, def.as<bool>()), def.GetTypeTag()); break;
	case SettingTypeInt: return SettingValue(ReadSettingRaw(path, def.as<int32>()), def.GetTypeTag()); break;
	case SettingTypeFloat: return SettingValue(ReadSettingRaw(path, def.as<double>()), def.GetTypeTag()); break;
	case SettingTypeString: return SettingValue(ReadSettingRaw(path, def.as<std::wstring>()), def.GetTypeTag()); break;
	case SettingTypeBinary: return SettingValue(ReadSettingRaw(path, def.as<std::vector<char> >()), def.GetTypeTag()); break;
	default: return SettingValue(); break;
	}
}

void IniFileSettingsBackend::WriteSetting(const SettingPath &path, const SettingValue &val)
{
	ASSERT(val.GetType() != SettingTypeNone);
	switch(val.GetType())
	{
	case SettingTypeBool: WriteSettingRaw(path, val.as<bool>()); break;
	case SettingTypeInt: WriteSettingRaw(path, val.as<int32>()); break;
	case SettingTypeFloat: WriteSettingRaw(path, val.as<double>()); break;
	case SettingTypeString: WriteSettingRaw(path, val.as<std::wstring>()); break;
	case SettingTypeBinary: WriteSettingRaw(path, val.as<std::vector<char> >()); break;
	default: break;
	}
}

void IniFileSettingsBackend::RemoveSetting(const SettingPath &path)
{
	RemoveSettingRaw(path);
}





IniFileSettingsContainer::IniFileSettingsContainer(const mpt::PathString &filename)
	: IniFileSettingsBackend(filename)
	, SettingsContainer(this)
{
	return;
}

IniFileSettingsContainer::~IniFileSettingsContainer()
{
	return;
}



DefaultSettingsContainer::DefaultSettingsContainer()
	: IniFileSettingsContainer(mpt::String::Decode(theApp.GetConfigFileName(), mpt::CharsetLocale))
{
	return;
}

DefaultSettingsContainer::~DefaultSettingsContainer()
{
	return;
}


