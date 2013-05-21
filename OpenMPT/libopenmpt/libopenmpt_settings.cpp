/*
 * libopenmpt_settings.cpp
 * -----------------------
 * Purpose: libopenmpt plugin settings
 * Notes  : (currently none)
 * Authors: OpenMPT Devs
 * The OpenMPT source code is released under the BSD license. Read LICENSE for more details.
 */

#include "libopenmpt_internal.h"

#define LIBOPENMPT_BUILD_SETTINGS_DLL
#include "libopenmpt_settings.h"
#include "SettingsForm.h"

namespace openmpt {

registry_settings::registry_settings( const char * subkey_ ) : subkey(subkey_) {
	return;
}

registry_settings::~registry_settings() {
	return;
}

void registry_settings::read_setting( const char * key, int & val ) {
	System::String ^ net_root = "HKEY_CURRENT_USER\\Software\\libopenmpt\\";
	System::String ^ net_path = gcnew System::String( subkey );
	System::String ^ net_key = gcnew System::String( key );
	System::Object ^ retval = Microsoft::Win32::Registry::GetValue( System::String::Concat( net_root, net_path ), net_key, val );
	if ( retval == nullptr ) {
		return;
	} else {
		val = (int)retval;
	}
}

void registry_settings::write_setting( const char * key, int val ) {
	System::String ^ net_root = "HKEY_CURRENT_USER\\Software\\libopenmpt\\";
	System::String ^ net_path = gcnew System::String( subkey );
	System::String ^ net_key = gcnew System::String( key );
	Microsoft::Win32::Registry::SetValue( System::String::Concat( net_root, net_path ), net_key, val );
}

void registry_settings::edit_settings( HWND parent, const char * title ) {
	libopenmpt::SettingsForm ^ form = gcnew libopenmpt::SettingsForm( title, this );
	System::Windows::Forms::IWin32Window ^w = System::Windows::Forms::Control::FromHandle((System::IntPtr)parent);
	form->Show(w);
}

} // namespace openmpt
