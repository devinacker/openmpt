
 project "pugixml"
  uuid "07B89124-7C71-42cc-81AB-62B09BB61F9B"
  language "C++"
  location ( "../../build/" .. mpt_projectpathname .. "/ext" )
  objdir "../../build/obj/pugixml"
  dofile "../../build/premake/premake-defaults-LIB.lua"
  dofile "../../build/premake/premake-defaults.lua"
  dofile "../../build/premake/premake-defaults-winver.lua"
  targetname "openmpt-pugixml"
  includedirs { }
	filter {}
	filter { "action:vs*", "not action:vs2008" }
		characterset "Unicode"
		flags { "Unicode" }
	filter { "action:vs*", "action:vs2008" }
		characterset "MBCS"
	filter {}
  files {
   "../../include/pugixml/src/pugixml.cpp",
  }
  files {
   "../../include/pugixml/src/pugiconfig.hpp",
   "../../include/pugixml/src/pugixml.hpp",
  }
