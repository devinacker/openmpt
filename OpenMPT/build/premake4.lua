
-- premake gets a tiny bit confused if the same project appears in multiple
-- solutions in a single run. premake adds a bogus $projectname path to the
-- intermediate objects directory in that case. work-around using multiple
-- invocations of premake and a custom option to distinguish them.

newoption {
 trigger     = "group",
 value       = "PROJECTS",
 description = "OpenMPT project group",
 allowed = {
  { "libopenmpt-all", "libopenmpt-all" },
  { "libopenmpt_test", "libopenmpt_test" },
  { "libopenmpt", "libopenmpt" },
  { "in_openmpt", "in_openmpt" },
  { "xmp-openmpt", "xmp-openmpt" },
  { "openmpt123", "openmpt123" },
  { "all-externals", "all-externals" }
 }
}

function postprocess_vs2008_main (filename)
	local text
	local infile
	local outfile
	infile = io.open(filename, "r")
	text = infile:read("*all")
	infile:close()
	text = string.gsub(text, "\t\t\t\tEntryPointSymbol=\"mainCRTStartup\"\n", "")
	outfile = io.open(filename, "w")
	outfile:write(text)
	outfile:close()
end

function postprocess_vs2010_mfc (filename)
	local text
	local infile
	local outfile
	infile = io.open(filename, "r")
	text = infile:read("*all")
	infile:close()
	text = string.gsub(text, "<UseOfMfc>Dynamic</UseOfMfc>", "<UseOfMfc>Static</UseOfMfc>")
	outfile = io.open(filename, "w")
	outfile:write(text)
	outfile:close()
end

function postprocess_vs2010_main (filename)
	local text
	local infile
	local outfile
	infile = io.open(filename, "r")
	text = infile:read("*all")
	infile:close()
	text = string.gsub(text, "<EntryPointSymbol>mainCRTStartup</EntryPointSymbol>", "")
	outfile = io.open(filename, "w")
	outfile:write(text)
	outfile:close()
end

newaction {
 trigger     = "postprocess",
 description = "OpenMPT postprocess the project files to mitigate premake problems",
 execute     = function ()
  postprocess_vs2008_main("build/vs2008/libopenmpt_test.vcproj")
  postprocess_vs2008_main("build/vs2008/openmpt123.vcproj")
  postprocess_vs2010_main("build/vs2010/libopenmpt_test.vcxproj")
  postprocess_vs2010_main("build/vs2010/openmpt123.vcxproj")
  postprocess_vs2010_mfc("build/vs2010/in_openmpt.vcxproj")
  postprocess_vs2010_mfc("build/vs2010/xmp-openmpt.vcxproj")
 end
}

if _ACTION == "vs2008" and _PREMAKE_VERSION ~= "4.3" then
 print "Premake 4.3 required"
 os.exit(1)
end

if _ACTION == "vs2010" and _PREMAKE_VERSION ~= "4.3" then
 print "Premake 4.3 required"
 os.exit(1)
end

if _ACTION == "postprocess" and _PREMAKE_VERSION ~= "4.3" then
 print "Premake 4.3 required"
 os.exit(1)
end

if _OPTIONS["group"] == "libopenmpt-all" then

solution "libopenmpt-all"
 location ( "../build/" .. _ACTION )
 configurations { "Debug", "Release" }
 platforms { "x32", "x64" }

 dofile "../build/premake4-win/mpt-libopenmpt_test.premake4.lua"
 dofile "../build/premake4-win/mpt-libopenmpt.premake4.lua"
 dofile "../build/premake4-win/mpt-libopenmpt_examples.premake4.lua"
 dofile "../build/premake4-win/mpt-libopenmptDLL.premake4.lua"
 dofile "../build/premake4-win/mpt-libopenmpt_modplug.premake4.lua"
 dofile "../build/premake4-win/mpt-in_openmpt.premake4.lua"
 dofile "../build/premake4-win/mpt-xmp-openmpt.premake4.lua"
 dofile "../build/premake4-win/mpt-openmpt123.premake4.lua"
 dofile "../build/premake4-win/ext-flac.premake4.lua"
 dofile "../build/premake4-win/ext-miniz.premake4.lua"
 dofile "../build/premake4-win/ext-portaudio.premake4.lua"

end

if _OPTIONS["group"] == "libopenmpt_test" then

solution "libopenmpt_test"
 location ( "../build/" .. _ACTION )
 configurations { "Debug", "Release" }
 platforms { "x32", "x64" }

 dofile "../build/premake4-win/mpt-libopenmpt_test.premake4.lua"
 dofile "../build/premake4-win/ext-miniz.premake4.lua"

end

if _OPTIONS["group"] == "in_openmpt" then

solution "in_openmpt"
 location ( "../build/" .. _ACTION )
 configurations { "Debug", "Release" }
 platforms { "x32" }

 dofile "../build/premake4-win/mpt-libopenmpt.premake4.lua"
 dofile "../build/premake4-win/mpt-in_openmpt.premake4.lua"
 dofile "../build/premake4-win/ext-miniz.premake4.lua"
 dofile "../build/premake4-win/ext-pugixml.premake4.lua"

end

if _OPTIONS["group"] == "xmp-openmpt" then

solution "xmp-openmpt"
 location ( "../build/" .. _ACTION )
 configurations { "Debug", "Release" }
 platforms { "x32" }

 dofile "../build/premake4-win/mpt-libopenmpt.premake4.lua"
 dofile "../build/premake4-win/mpt-xmp-openmpt.premake4.lua"
 dofile "../build/premake4-win/ext-miniz.premake4.lua"
 dofile "../build/premake4-win/ext-pugixml.premake4.lua"

end

-- should stay the last libopenmpt solution in order to overwrite the libopenmpt base project with all possible configurations
if _OPTIONS["group"] == "libopenmpt" then

solution "libopenmpt"
 location ( "../build/" .. _ACTION )
 configurations { "Debug", "Release" }
 platforms { "x32", "x64" }

 dofile "../build/premake4-win/mpt-libopenmpt.premake4.lua"
 dofile "../build/premake4-win/mpt-libopenmpt_examples.premake4.lua"
 dofile "../build/premake4-win/mpt-libopenmptDLL.premake4.lua"
 dofile "../build/premake4-win/mpt-libopenmpt_modplug.premake4.lua"
 dofile "../build/premake4-win/ext-miniz.premake4.lua"
 dofile "../build/premake4-win/ext-miniz-shared.premake4.lua"
 dofile "../build/premake4-win/ext-portaudio.premake4.lua"

end

if _OPTIONS["group"] == "openmpt123" then

solution "openmpt123"
 location ( "../build/" .. _ACTION )
 configurations { "Debug", "Release" }
 platforms { "x32", "x64" }

 dofile "../build/premake4-win/mpt-openmpt123.premake4.lua"
 dofile "../build/premake4-win/mpt-libopenmpt.premake4.lua"
 dofile "../build/premake4-win/ext-flac.premake4.lua"
 dofile "../build/premake4-win/ext-miniz.premake4.lua"
 dofile "../build/premake4-win/ext-ogg.premake4.lua"
 dofile "../build/premake4-win/ext-portaudio.premake4.lua"

end

-- overwrite all external projects once again with the full matrix of possible build config combinations
if _OPTIONS["group"] == "all-externals" then

solution "all-externals"
 configurations { "Debug", "Release", "ReleaseNoLTCG" }
 platforms { "x32", "x64" }

 dofile "../build/premake4-win/ext-flac.premake4.lua"
 dofile "../build/premake4-win/ext-lhasa.premake4.lua"
 dofile "../build/premake4-win/ext-miniz.premake4.lua"
 dofile "../build/premake4-win/ext-miniz-shared.premake4.lua"
 dofile "../build/premake4-win/ext-minizip.premake4.lua"
 dofile "../build/premake4-win/ext-ogg.premake4.lua"
 dofile "../build/premake4-win/ext-portaudio.premake4.lua"
 dofile "../build/premake4-win/ext-portmidi.premake4.lua"
 dofile "../build/premake4-win/ext-pugixml.premake4.lua"
 dofile "../build/premake4-win/ext-r8brain.premake4.lua"
 dofile "../build/premake4-win/ext-smbPitchShift.premake4.lua"
 dofile "../build/premake4-win/ext-soundtouch.premake4.lua"
 dofile "../build/premake4-win/ext-UnRAR.premake4.lua"
 dofile "../build/premake4-win/ext-zlib.premake4.lua"

end
