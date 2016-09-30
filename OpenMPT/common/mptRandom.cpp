/*
 * mptRandom.cpp
 * -------------
 * Purpose: PRNG
 * Notes  : (currently none)
 * Authors: OpenMPT Devs
 * The OpenMPT source code is released under the BSD license. Read LICENSE for more details.
 */

#include "stdafx.h"

#include "mptRandom.h"

#include "Endianness.h"
#include "mptCRC.h"

#include <algorithm>

#include <cmath>
#include <ctime>
#include <cstdlib>

#if MPT_OS_WINDOWS
#include <windows.h>
#endif // MPT_OS_WINDOWS


OPENMPT_NAMESPACE_BEGIN


namespace mpt
{


template <typename T>
static T log2(T x)
{
	return std::log(x) / std::log(static_cast<T>(2));
}


static inline int lower_bound_entropy_bits(unsigned int x)
{
	// easy to compile-time evaluate even for stupid compilers
	return
		x >= 0xffffffffu ? 32 :
		x >= 0x7fffffffu ? 31 :
		x >= 0x3fffffffu ? 30 :
		x >= 0x1fffffffu ? 29 :
		x >= 0x0fffffffu ? 28 :
		x >= 0x07ffffffu ? 27 :
		x >= 0x03ffffffu ? 26 :
		x >= 0x01ffffffu ? 25 :
		x >= 0x00ffffffu ? 24 :
		x >= 0x007fffffu ? 23 :
		x >= 0x003fffffu ? 22 :
		x >= 0x001fffffu ? 21 :
		x >= 0x000fffffu ? 20 :
		x >= 0x0007ffffu ? 19 :
		x >= 0x0003ffffu ? 18 :
		x >= 0x0001ffffu ? 17 :
		x >= 0x0000ffffu ? 16 :
		x >= 0x00007fffu ? 15 :
		x >= 0x00003fffu ? 14 :
		x >= 0x00001fffu ? 13 :
		x >= 0x00000fffu ? 12 :
		x >= 0x000007ffu ? 11 :
		x >= 0x000003ffu ? 10 :
		x >= 0x000001ffu ?  9 :
		x >= 0x000000ffu ?  8 :
		x >= 0x0000007fu ?  7 :
		x >= 0x0000003fu ?  6 :
		x >= 0x0000001fu ?  5 :
		x >= 0x0000000fu ?  4 :
		x >= 0x00000007u ?  3 :
		x >= 0x00000003u ?  2 :
		x >= 0x00000001u ?  1 :
		0;
}


template <typename T>
static inline bool is_mask(T x)
{
	STATIC_ASSERT(std::numeric_limits<T>::is_integer);
	typedef typename std::make_unsigned<T>::type unsigned_T;
	unsigned_T ux = static_cast<unsigned_T>(x);
	unsigned_T mask = 0;
	for(std::size_t bits = 0; bits <= (sizeof(unsigned_T) * 8); ++bits)
	{
		mask = (mask << 1) | 1u;
		if(ux == mask)
		{
			return true;
		}
	}
	return false;
}


namespace {
template <typename T> struct default_hash { };
template <> struct default_hash<uint8>  { typedef mpt::checksum::crc16 type; };
template <> struct default_hash<uint16> { typedef mpt::checksum::crc16 type; };
template <> struct default_hash<uint32> { typedef mpt::checksum::crc32c type; };
template <> struct default_hash<uint64> { typedef mpt::checksum::crc64_jones type; };
}

template <typename T>
static T generate_timeseed()
{
	// Note: CRC is actually not that good a choice here, but it is simple and we
	// already have an implementaion available. Better choices for mixing entropy
	// would be a hash function with proper avalanche characteristics or a block
	// or stream cipher with any pre-choosen random key and IV. The only aspect we
	// really need here is whitening of the bits.
	typename mpt::default_hash<T>::type hash;
	
#ifdef MPT_BUILD_FUZZER

	return static_cast<T>(mpt::FUZZER_RNG_SEED);

#else // !MPT_BUILD_FUZZER

	{
		#if MPT_OS_WINDOWS
			FILETIME t;
			MemsetZero(t);
			GetSystemTimeAsFileTime(&t);
		#else // !MPT_OS_WINDOWS
			std::time_t t = std::time(nullptr);
		#endif // MPT_OS_WINDOWS
		mpt::byte bytes[sizeof(t)];
		std::memcpy(bytes, &t, sizeof(t));
		MPT_MAYBE_CONSTANT_IF(mpt::endian_is_little())
		{
			std::reverse(bytes + 0, bytes + CountOf(bytes));
		}
		hash(bytes + 0, bytes + CountOf(bytes));
	}

	{
		std::clock_t c = std::clock();
		mpt::byte bytes[sizeof(c)];
		std::memcpy(bytes, &c, sizeof(c));
		MPT_MAYBE_CONSTANT_IF(mpt::endian_is_little())
		{
			std::reverse(bytes + 0, bytes + CountOf(bytes));
		}
		hash(bytes + 0, bytes + CountOf(bytes));
	}

	return static_cast<T>(hash.result());

#endif // MPT_BUILD_FUZZER

}


#ifdef MODPLUG_TRACKER

namespace rng
{

void crand::reseed(uint32 seed)
{
	std::srand(seed);
}

crand::result_type crand::min()
{
	return 0;
}

crand::result_type crand::max()
{
	return RAND_MAX;
}

int crand::result_bits()
{
	return lower_bound_entropy_bits(RAND_MAX);
}

crand::result_type crand::operator()()
{
	return std::rand();
}

} // namespace rng

#endif // MODPLUG_TRACKER

#if MPT_STD_RANDOM

sane_random_device::sane_random_device()
	: rd_reliable(rd.entropy() > 0.0)
{
	if(!rd_reliable)
	{
		uint64 seed_val = mpt::generate_timeseed<uint64>();
		unsigned int seeds[2];
		seeds[0] = static_cast<uint32>(seed_val >> 32);
		seeds[1] = static_cast<uint32>(seed_val >>  0);
		std::seed_seq seed(seeds + 0, seeds + 2);
		rd_fallback = mpt::make_unique<std::mt19937>(seed);
	}
}

sane_random_device::sane_random_device(const std::string & token)
	: rd(token)
	, rd_reliable(rd.entropy() > 0.0)
{
	if(!rd_reliable)
	{
		uint64 seed_val = mpt::generate_timeseed<uint64>();
		std::vector<unsigned int> seeds;
		seeds.push_back(static_cast<uint32>(seed_val >> 32));
		seeds.push_back(static_cast<uint32>(seed_val >>  0));
		for(std::size_t i = 0; i < token.length(); ++i)
		{
			seeds.push_back(static_cast<unsigned int>(static_cast<unsigned char>(token[i])));
		}
		std::seed_seq seed(seeds.begin(), seeds.end());
		rd_fallback = mpt::make_unique<std::mt19937>(seed);
	}
}

sane_random_device::result_type sane_random_device::min()
{
	return std::numeric_limits<result_type>::min();
}

sane_random_device::result_type sane_random_device::max()
{
	return std::numeric_limits<result_type>::max();
}

int sane_random_device::result_bits()
{
	return sizeof(result_type) * 8;
}

sane_random_device::result_type sane_random_device::operator()()
{
	MPT_LOCK_GUARD<mpt::mutex> l(m);
	result_type result = 0;
	if(rd.min() != 0 || !mpt::is_mask(rd.max()))
	{ // insane std::random_device
		//  This implementation is not exactly uniformly distributed but good enough
		// for OpenMPT.
		double rd_min = static_cast<double>(rd.min());
		double rd_max = static_cast<double>(rd.max());
		double rd_range = rd_max - rd_min;
		double rd_size = rd_range + 1.0;
		double rd_entropy = mpt::log2(rd_size);
		int iterations = static_cast<int>(std::ceil(result_bits() / rd_entropy));
		double tmp = 0.0;
		for(int i = 0; i < iterations; ++i)
		{
			tmp = (tmp * rd_size) + (static_cast<double>(rd()) - rd_min);
		}
		double result_01 = std::floor(tmp / std::pow(rd_size, iterations));
		result = static_cast<result_type>(std::floor(result_01 * (static_cast<double>(max() - min()) + 1.0))) + min();
	} else
	{ // sane std::random_device
		result = 0;
		std::size_t rd_bits = mpt::lower_bound_entropy_bits(rd.max());
		for(std::size_t entropy = 0; entropy < (sizeof(result_type) * 8); entropy += rd_bits)
		{
			if(rd_bits < (sizeof(result_type) * 8))
			{
				result = (result << rd_bits) | static_cast<result_type>(rd());
			} else
			{
				result = result | static_cast<result_type>(rd());
			}
		}
	}
	if(!rd_reliable)
	{ // std::random_device is unreliable
		//  XOR the generated random number with more entropy from the time-seeded
		// PRNG.
		//  Note: This is safe even if the std::random_device itself is implemented
		// as a std::mt19937 PRNG because we are very likely using a different
		// seed.
		result ^= mpt::random<result_type>(*rd_fallback);
	}
	return result;
}

#endif // MPT_STD_RANDOM

prng_random_device_seeder::prng_random_device_seeder()
{
	return;
}

uint8 prng_random_device_seeder::generate_seed8()
{
	return mpt::generate_timeseed<uint8>();
}

uint16 prng_random_device_seeder::generate_seed16()
{
	return mpt::generate_timeseed<uint16>();
}

uint32 prng_random_device_seeder::generate_seed32()
{
	return mpt::generate_timeseed<uint32>();
}

uint64 prng_random_device_seeder::generate_seed64()
{
	return mpt::generate_timeseed<uint64>();
}

#if defined(MODPLUG_TRACKER) && !defined(MPT_BUILD_WINESUPPORT)

static mpt::random_device *g_rd = nullptr;
static mpt::thread_safe_prng<mpt::best_prng> *g_best_prng = nullptr;

void set_global_random_device(mpt::random_device *rd)
{
	g_rd = rd;
}

void set_global_prng(mpt::thread_safe_prng<mpt::best_prng> *prng)
{
	g_best_prng = prng;
}

mpt::random_device & global_random_device()
{
	return *g_rd;
}

mpt::thread_safe_prng<mpt::best_prng> & global_prng()
{
	return *g_best_prng;
}

#else

mpt::random_device & global_random_device()
{
	static mpt::random_device g_rd;
	return g_rd;
}

mpt::thread_safe_prng<mpt::best_prng> & global_prng()
{
	static mpt::thread_safe_prng<mpt::best_prng> g_best_prng(global_random_device());
	return g_best_prng;
}

#endif // MODPLUG_TRACKER && !MPT_BUILD_WINESUPPORT


} // namespace mpt


OPENMPT_NAMESPACE_END
