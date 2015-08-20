/*
 * libopenmpt_example_c_mem.c
 * --------------------------
 * Purpose: libopenmpt C API simple example.
 *          This examples demonstrates how to play a module file from memory in C.
 *          The module file to play can be specified as a command line parameter.
 *          If no parameter is given, a built-in module is loaded from memory instead.
 *          PortAudio is used for audio output.
 * Notes  : (currently none)
 * Authors: OpenMPT Devs
 * The OpenMPT source code is released under the BSD license. Read LICENSE for more details.
 */

#include <memory.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libopenmpt/libopenmpt.h>
#include "chipsim.h"

#include <portaudio.h>

#define BUFFERSIZE 480
#define SAMPLERATE 48000

static int16_t left[BUFFERSIZE];
static int16_t right[BUFFERSIZE];
static int16_t * const buffers[2] = { left, right };

#if (defined(_WIN32) || defined(WIN32)) && (defined(_UNICODE) || defined(UNICODE))
int wmain( int argc, wchar_t * argv[] ) {
#else
int main( int argc, char * argv[] ) {
#endif
	FILE * file = 0;
	size_t size = 0;
	void * data = 0;
	openmpt_module * mod = 0;
	size_t count = 0;
	PaStream * stream = 0;
	PaStreamParameters streamparameters;
	memset( &streamparameters, 0, sizeof( PaStreamParameters ) );
	
	if ( argc > 1 ) {
#if (defined(_WIN32) || defined(WIN32)) && (defined(_UNICODE) || defined(UNICODE))
		file = _wfopen( argv[1], L"rb" );
#else
		file = fopen( argv[1], "rb" );
#endif
		fseek( file, 0, SEEK_END );
		size = ftell( file );
		fseek( file, 0, SEEK_SET );
		data = malloc( size );
		size = fread( data, 1, size, file );
		fclose( file );
	} else {
		/* no command line parameter specified - play an embedded example module. */
		data = chipsim_it;
		size = sizeof(chipsim_it);
	}

	mod = openmpt_module_create_from_memory( data, size, NULL, NULL, NULL );

	if ( argc > 1 ) {
		free( data );
	}

	Pa_Initialize();
	streamparameters.device = Pa_GetDefaultOutputDevice();
	streamparameters.channelCount = 2;
	streamparameters.sampleFormat = paInt16 | paNonInterleaved;
	streamparameters.suggestedLatency = Pa_GetDeviceInfo( streamparameters.device )->defaultHighOutputLatency;
	Pa_OpenStream( &stream, NULL, &streamparameters, SAMPLERATE, paFramesPerBufferUnspecified, 0, NULL, NULL );
	Pa_StartStream( stream );
	while ( 1 ) {
		count = openmpt_module_read_stereo( mod, SAMPLERATE, BUFFERSIZE, left, right );
		if ( count == 0 ) {
			break;
		}
		Pa_WriteStream( stream, buffers, (unsigned long)count );
	}
	Pa_StopStream( stream );
	Pa_CloseStream( stream );
	Pa_Terminate();
	openmpt_module_destroy( mod );
	return 0;
}
