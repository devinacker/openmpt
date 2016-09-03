/*
 * typedefs.h
 * ----------
 * Purpose: Basic data type definitions and assorted compiler-related helpers.
 * Notes  : (currently none)
 * Authors: OpenMPT Devs
 * The OpenMPT source code is released under the BSD license. Read LICENSE for more details.
 */


#pragma once



OPENMPT_NAMESPACE_BEGIN



// nullptr for ancient gcc
#if MPT_COMPILER_GCC
#if MPT_GCC_BEFORE(4,6,0)
#define nullptr 0
#endif
#endif



//  CountOf macro computes the number of elements in a statically-allocated array.
#if MPT_COMPILER_MSVC
#define MPT_ARRAYCOUNT(x) _countof(x)
#else
#define MPT_ARRAYCOUNT(x) (sizeof((x))/sizeof((x)[0]))
#endif



#if MPT_COMPILER_MSVC && MPT_MSVC_BEFORE(2015,0)
#define MPT_ALIGNOF(type) __alignof( type )
#elif MPT_COMPILER_GCC && MPT_GCC_BEFORE(4,5,0)
#define MPT_ALIGNOF(type) __alignof__( type )
#else
#define MPT_ALIGNOF(type) alignof( type )
#endif



// Advanced inline attributes
#if MPT_COMPILER_MSVC
#define MPT_FORCEINLINE __forceinline
#define MPT_NOINLINE    __declspec(noinline)
#elif MPT_COMPILER_GCC || MPT_COMPILER_CLANG || MPT_COMPILER_MSVCCLANGC2
#define MPT_FORCEINLINE __attribute__((always_inline)) inline
#define MPT_NOINLINE    __attribute__((noinline))
#else
#define MPT_FORCEINLINE inline
#define MPT_NOINLINE
#endif



// Use MPT_RESTRICT to indicate that a pointer is guaranteed to not be aliased.
#if MPT_COMPILER_MSVC || MPT_COMPILER_GCC || MPT_COMPILER_CLANG || MPT_COMPILER_MSVCCLANGC2
#define MPT_RESTRICT __restrict
#else
#define MPT_RESTRICT
#endif



// Some functions might be deprecated although they are still in use.
// Tag them with "MPT_DEPRECATED".
#if MPT_COMPILER_MSVC
#define MPT_DEPRECATED __declspec(deprecated)
#elif MPT_COMPILER_GCC || MPT_COMPILER_CLANG || MPT_COMPILER_MSVCCLANGC2
#define MPT_DEPRECATED __attribute__((deprecated))
#else
#define MPT_DEPRECATED
#endif
#if defined(MODPLUG_TRACKER)
#define MPT_DEPRECATED_TRACKER    MPT_DEPRECATED
#define MPT_DEPRECATED_LIBOPENMPT 
#elif defined(LIBOPENMPT_BUILD)
#define MPT_DEPRECATED_TRACKER    
#define MPT_DEPRECATED_LIBOPENMPT MPT_DEPRECATED
#else
#define MPT_DEPRECATED_TRACKER    MPT_DEPRECATED
#define MPT_DEPRECATED_LIBOPENMPT MPT_DEPRECATED
#endif



// Exception type that is used to catch "operator new" exceptions.
#if defined(_MFC_VER)
typedef CMemoryException * MPTMemoryException;
#else
OPENMPT_NAMESPACE_END
#include <new>
OPENMPT_NAMESPACE_BEGIN
typedef std::bad_alloc & MPTMemoryException;
#endif



OPENMPT_NAMESPACE_END
#include <memory>
#include <utility>
OPENMPT_NAMESPACE_BEGIN

#if MPT_GCC_BEFORE(4,6,0)
#define MPT_SHARED_PTR_NULL(T) std::shared_ptr<T>()
#else
#define MPT_SHARED_PTR_NULL(T) nullptr
#endif



// We provide a simple scoped_ptr which can also be useful in addition to
// std::unique_ptr.
namespace mpt {
template<typename T> class scoped_ptr
{
private:
	// Copying is not supported.
	MPT_DEPRECATED scoped_ptr(const scoped_ptr & other) : m_p(other.m_p) {} // = delete
	// Copying is not supported.
	MPT_DEPRECATED scoped_ptr & operator=(const scoped_ptr & other) { m_p = other.m_p; return *this; } // = delete
public:
	typedef T element_type;
	class initializer_impl
	{
		friend class scoped_ptr<element_type>;
	private:
		element_type * m_p;
	private:
		inline initializer_impl(const initializer_impl & other) : m_p(other.m_p) {}
		inline initializer_impl & operator=(const initializer_impl & other) { m_p = other.m_p; return *this; }
	public:
		explicit initializer_impl(T * p) : m_p(p) {}
	public:
		template <typename Tobj> static inline typename mpt::scoped_ptr<T>::initializer_impl make() { return mpt::scoped_ptr<T>::initializer_impl(new Tobj()); }
		template <typename Tobj, typename T1> static inline typename mpt::scoped_ptr<T>::initializer_impl make(const T1 &x1) { return mpt::scoped_ptr<T>::initializer_impl(new Tobj(x1)); }
		template <typename Tobj, typename T1, typename T2> static inline typename mpt::scoped_ptr<T>::initializer_impl make(const T1 &x1, const T2 &x2) { return mpt::scoped_ptr<T>::initializer_impl(new Tobj(x1, x2)); }
		template <typename Tobj, typename T1, typename T2, typename T3> static inline typename mpt::scoped_ptr<T>::initializer_impl make(const T1 &x1, const T2 &x2, const T3 &x3) { return mpt::scoped_ptr<T>::initializer_impl(new Tobj(x1, x2, x3)); }
		template <typename Tobj, typename T1, typename T2, typename T3, typename T4> static inline typename mpt::scoped_ptr<T>::initializer_impl make(const T1 &x1, const T2 &x2, const T3 &x3, const T4 &x4) { return mpt::scoped_ptr<T>::initializer_impl(new Tobj(x1, x2, x3, x4)); }
	};
	class initializer_impl_result
	{
		friend class scoped_ptr<element_type>;
	private:
		element_type * m_p;
	public:
		inline initializer_impl_result(const initializer_impl_result & other) : m_p(other.m_p) {}
	private:
		MPT_DEPRECATED initializer_impl_result & operator=(const initializer_impl_result & other ) { m_p = other.m_p; return *this; } // = delete
	public:
		inline initializer_impl_result(const initializer_impl & other) : m_p(other.m_p) {}
		inline initializer_impl_result & operator=(const initializer_impl & other) { m_p = other.m_p; return *this; }
	public:
		explicit initializer_impl_result(T * p) : m_p(p) {}
	};
	class initializer
	{
		friend class scoped_ptr<element_type>;
	private:
		element_type * m_p;
	private:
		MPT_DEPRECATED initializer(const initializer & /* other */ ) {} // = delete
		MPT_DEPRECATED initializer & operator=(const initializer & /* other */ ) {return *this;} // = delete
	public:
		inline initializer(const initializer_impl_result & other) : m_p(other.m_p) {}
		inline initializer & operator=(const initializer_impl_result & other) { m_p = other.m_p; return *this; }
	public:
		explicit initializer(T * p) : m_p(p) {}
	};
private:
	element_type * m_p;
public:
	// Creates a scoped_ptr without any owned object.
	scoped_ptr() : m_p(nullptr) {}
	// Creates a scoped_ptr and assumes ownership of p (if p != nulltr).
	explicit scoped_ptr(T * p) : m_p(p) {}
	// Deletes the currently owned object (if any).
	~scoped_ptr() { if(m_p) { delete m_p; m_p = nullptr; } }
	// Deletes the currently owned object (if any), and assumes ownership of the passed object p (if any).
	void reset(T * p = nullptr) { if(m_p) { delete m_p; m_p = nullptr; } m_p = p; }
	// Returns a reference to the owned object. Behaviour is undefined if there is no owned object.
	T & operator*() const { return *m_p; }
	// Returns a pointer to the owned object, or nullptr if there is not any owned object. Ownership is not transferred.
	T * operator->() const { return m_p; }
	// Returns a pointer to the owned object, or nullptr if there is not any owned object. Ownership is not transferred.
	T * get() const { return m_p; }
	// Give away ownership of the owned object.
	T * release() { T * ret = m_p; m_p = nullptr; return ret; }
	// Creates a scoped_ptr and assumes ownership of the passed object init (if any).
	scoped_ptr (const initializer & init) : m_p(init.m_p) {}
	// Deletes the currently owned object (if any), and assumes ownership of the passed object init (if any).
	scoped_ptr & operator=(const initializer & init) { if(m_p) { delete m_p; m_p = nullptr; } m_p = init.m_p; return *this; }
	// Returns true iff the scoped_ptr currently owns an object (i.e. same as (bool)p).
	operator bool() const { return m_p ? true : false; }
};
template <typename T> inline typename mpt::scoped_ptr<T>::initializer_impl_result make_scoped() { typedef typename mpt::scoped_ptr<T>::initializer_impl init; return init::template make<T>(); }
template <typename T, typename T1> inline typename mpt::scoped_ptr<T>::initializer_impl_result make_scoped(const T1 &x1) { typedef typename mpt::scoped_ptr<T>::initializer_impl init; return init::template make<T>(x1); }
template <typename T, typename T1, typename T2> inline typename mpt::scoped_ptr<T>::initializer_impl_result make_scoped(const T1 &x1, const T2 &x2) { typedef typename mpt::scoped_ptr<T>::initializer_impl init; return init::template make<T>(x1, x2); }
template <typename T, typename T1, typename T2, typename T3> inline typename mpt::scoped_ptr<T>::initializer_impl_result make_scoped(const T1 &x1, const T2 &x2, const T3 &x3) { typedef typename mpt::scoped_ptr<T>::initializer_impl init; return init::template make<T>(x1, x2, x3); }
template <typename T, typename T1, typename T2, typename T3, typename T4> inline typename mpt::scoped_ptr<T>::initializer_impl_result make_scoped(const T1 &x1, const T2 &x2, const T3 &x3, const T4 &x4) { typedef typename mpt::scoped_ptr<T>::initializer_impl init; return init::template make<T>(x1, x2, x3, x4); }
} // namespace mpt



#if MPT_COMPILER_MSVC
#define MPT_CONSTANT_IF(x) \
  __pragma(warning(push)) \
  __pragma(warning(disable:4127)) \
  if(x) \
  __pragma(warning(pop)) \
/**/
#define MPT_MAYBE_CONSTANT_IF(x) \
  __pragma(warning(push)) \
  __pragma(warning(disable:4127)) \
  if(x) \
  __pragma(warning(pop)) \
/**/
#endif

#if MPT_COMPILER_GCC
#if MPT_GCC_AT_LEAST(4,6,0)
#define MPT_MAYBE_CONSTANT_IF(x) \
  _Pragma("GCC diagnostic push") \
  _Pragma("GCC diagnostic ignored \"-Wtype-limits\"") \
  if(x) \
  _Pragma("GCC diagnostic pop") \
/**/
#elif MPT_GCC_AT_LEAST(4,5,0)
#define MPT_MAYBE_CONSTANT_IF(x) \
  _Pragma("GCC diagnostic ignored \"-Wtype-limits\"") \
  if(x) \
/**/
#else
// GCC 4.4 does not like _Pragma diagnostic inside functions.
// As GCC 4.4 is one of our major compilers, we do not want a noisy build.
// Thus, just disable this warning globally. (not required for now)
//#pragma GCC diagnostic ignored "-Wtype-limits"
#endif
#endif

#if MPT_COMPILER_CLANG || MPT_COMPILER_MSVCCLANGC2
#define MPT_MAYBE_CONSTANT_IF(x) \
  _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Wunknown-pragmas\"") \
  _Pragma("clang diagnostic ignored \"-Wtype-limits\"") \
  _Pragma("clang diagnostic ignored \"-Wtautological-constant-out-of-range-compare\"") \
  if(x) \
  _Pragma("clang diagnostic pop") \
/**/
#endif

#if !defined(MPT_CONSTANT_IF)
// MPT_CONSTANT_IF disables compiler warnings for conditions that are either always true or always false for some reason (dependent on template arguments for example)
#define MPT_CONSTANT_IF(x) if(x)
#endif

#if !defined(MPT_MAYBE_CONSTANT_IF)
// MPT_MAYBE_CONSTANT_IF disables compiler warnings for conditions that may in some case be either always false or always true (this may turn out to be useful in ASSERTions in some cases).
#define MPT_MAYBE_CONSTANT_IF(x) if(x)
#endif



#if MPT_COMPILER_MSVC
// MSVC warns for the well-known and widespread "do { } while(0)" idiom with warning level 4 ("conditional expression is constant").
// It does not warn with "while(0,0)". However this again causes warnings with other compilers.
// Solve it with a macro.
#define MPT_DO do
#define MPT_WHILE_0 while(0,0)
#endif

#ifndef MPT_DO
#define MPT_DO do
#endif
#ifndef MPT_WHILE_0
#define MPT_WHILE_0 while(0)
#endif



// Static code checkers might need to get the knowledge of our assertions transferred to them.
#define MPT_CHECKER_ASSUME_ASSERTIONS 1
//#define MPT_CHECKER_ASSUME_ASSERTIONS 0

#ifdef MPT_BUILD_ANALYZED

#if MPT_COMPILER_MSVC

#if MPT_CHECKER_ASSUME_ASSERTIONS
#define MPT_CHECKER_ASSUME(x) __analysis_assume(!!(x))
#endif

#elif MPT_COMPILER_CLANG || MPT_COMPILER_MSVCCLANGC2

#if MPT_CHECKER_ASSUME_ASSERTIONS
#ifdef NDEBUG
#error "Builds for static analyzers depend on std::assert being enabled, but the current build has #define NDEBUG. This makes no sense."
#endif
OPENMPT_NAMESPACE_END
#include <cassert>
OPENMPT_NAMESPACE_BEGIN
#define MPT_CHECKER_ASSUME(x) assert(!!(x))
#endif

#endif // MPT_COMPILER

#endif // MPT_BUILD_ANALYZED

#ifndef MPT_CHECKER_ASSUME
#define MPT_CHECKER_ASSUME(x) MPT_DO { } MPT_WHILE_0
#endif



#if defined(_MFC_VER)

#if !defined(ASSERT)
#error "MFC is expected to #define ASSERT"
#endif // !defined(ASERRT)
#define MPT_FRAMEWORK_ASSERT_IS_DEFINED

#if defined(_DEBUG)
 #define MPT_FRAMEWORK_ASSERT_IS_ACTIVE 1
#else // !_DEBUG
 #define MPT_FRAMEWORK_ASSERT_IS_ACTIVE 0
#endif // _DEBUG

// let MFC handle our asserts
#define MPT_ASSERT_USE_FRAMEWORK 1

#else // !_MFC_VER

#if defined(ASSERT)
#define MPT_FRAMEWORK_ASSERT_IS_DEFINED
#if defined(_DEBUG)
 #define MPT_FRAMEWORK_ASSERT_IS_ACTIVE 1
#else // !_DEBUG
 #define MPT_FRAMEWORK_ASSERT_IS_ACTIVE 0
#endif // _DEBUG
#endif // !defined(ASERRT)

// handle assert in our own way without relying on some platform-/framework-specific assert implementation
#define MPT_ASSERT_USE_FRAMEWORK 0

#endif // _MFC_VER


#if defined(MPT_FRAMEWORK_ASSERT_IS_DEFINED) && (MPT_ASSERT_USE_FRAMEWORK == 1)

#define MPT_ASSERT_NOTREACHED()          ASSERT(0)
#define MPT_ASSERT(expr)                 ASSERT((expr))
#define MPT_ASSERT_MSG(expr, msg)        ASSERT((expr) && (msg))
#if (MPT_FRAMEWORK_ASSERT_IS_ACTIVE == 1)
#define MPT_ASSERT_ALWAYS(expr)          ASSERT((expr))
#define MPT_ASSERT_ALWAYS_MSG(expr, msg) ASSERT((expr) && (msg))
#else
#define MPT_ASSERT_ALWAYS(expr)          MPT_DO { if(!(expr)) { AssertHandler(__FILE__, __LINE__, __FUNCTION__, #expr); } MPT_CHECKER_ASSUME(expr); } MPT_WHILE_0
#define MPT_ASSERT_ALWAYS_MSG(expr, msg) MPT_DO { if(!(expr)) { AssertHandler(__FILE__, __LINE__, __FUNCTION__, #expr, msg); } MPT_CHECKER_ASSUME(expr); } MPT_WHILE_0
#ifndef MPT_ASSERT_HANDLER_NEEDED
#define MPT_ASSERT_HANDLER_NEEDED
#endif
#endif

#elif defined(NO_ASSERTS)

#define MPT_ASSERT_NOTREACHED()          MPT_CHECKER_ASSUME(0)
#define MPT_ASSERT(expr)                 MPT_CHECKER_ASSUME(expr)
#define MPT_ASSERT_MSG(expr, msg)        MPT_CHECKER_ASSUME(expr)
#define MPT_ASSERT_ALWAYS(expr)          MPT_DO { if(!(expr)) { AssertHandler(__FILE__, __LINE__, __FUNCTION__, #expr); } MPT_CHECKER_ASSUME(expr); } MPT_WHILE_0
#define MPT_ASSERT_ALWAYS_MSG(expr, msg) MPT_DO { if(!(expr)) { AssertHandler(__FILE__, __LINE__, __FUNCTION__, #expr, msg); } MPT_CHECKER_ASSUME(expr); } MPT_WHILE_0
#ifndef MPT_ASSERT_HANDLER_NEEDED
#define MPT_ASSERT_HANDLER_NEEDED
#endif

#else // !NO_ASSERTS

#define MPT_ASSERT_NOTREACHED()          MPT_DO { MPT_CONSTANT_IF(!(0)) { AssertHandler(__FILE__, __LINE__, __FUNCTION__, "0"); } MPT_CHECKER_ASSUME(0); } MPT_WHILE_0
#define MPT_ASSERT(expr)                 MPT_DO { if(!(expr)) { AssertHandler(__FILE__, __LINE__, __FUNCTION__, #expr); } MPT_CHECKER_ASSUME(expr); } MPT_WHILE_0
#define MPT_ASSERT_MSG(expr, msg)        MPT_DO { if(!(expr)) { AssertHandler(__FILE__, __LINE__, __FUNCTION__, #expr, msg); } MPT_CHECKER_ASSUME(expr); } MPT_WHILE_0
#define MPT_ASSERT_ALWAYS(expr)          MPT_DO { if(!(expr)) { AssertHandler(__FILE__, __LINE__, __FUNCTION__, #expr); } MPT_CHECKER_ASSUME(expr); } MPT_WHILE_0
#define MPT_ASSERT_ALWAYS_MSG(expr, msg) MPT_DO { if(!(expr)) { AssertHandler(__FILE__, __LINE__, __FUNCTION__, #expr, msg); } MPT_CHECKER_ASSUME(expr); } MPT_WHILE_0
#ifndef MPT_ASSERT_HANDLER_NEEDED
#define MPT_ASSERT_HANDLER_NEEDED
#endif

#endif // NO_ASSERTS


#if defined(MPT_ASSERT_HANDLER_NEEDED)
// custom assert handler needed
MPT_NOINLINE void AssertHandler(const char *file, int line, const char *function, const char *expr, const char *msg=nullptr);
#endif // MPT_ASSERT_HANDLER_NEEDED



// Compile time assert.
#define MPT_STATIC_ASSERT(expr) static_assert((expr), "compile time assertion failed: " #expr)



// Macro for marking intentional fall-throughs in switch statements - can be used for static analysis if supported.
#if MPT_COMPILER_MSVC
#define MPT_FALLTHROUGH __fallthrough
#elif MPT_COMPILER_CLANG || MPT_COMPILER_MSVCCLANGC2
#define MPT_FALLTHROUGH [[clang::fallthrough]]
#else
#define MPT_FALLTHROUGH MPT_DO { } MPT_WHILE_0
#endif



OPENMPT_NAMESPACE_END
#include <climits>
#include <cstdint>
#include <stdint.h>
OPENMPT_NAMESPACE_BEGIN

typedef std::int8_t   int8;
typedef std::int16_t  int16;
typedef std::int32_t  int32;
typedef std::int64_t  int64;
typedef std::uint8_t  uint8;
typedef std::uint16_t uint16;
typedef std::uint32_t uint32;
typedef std::uint64_t uint64;

const int8 int8_min     = INT8_MIN;
const int16 int16_min   = INT16_MIN;
const int32 int32_min   = INT32_MIN;
const int64 int64_min   = INT64_MIN;

const int8 int8_max     = INT8_MAX;
const int16 int16_max   = INT16_MAX;
const int32 int32_max   = INT32_MAX;
const int64 int64_max   = INT64_MAX;

const uint8 uint8_max   = UINT8_MAX;
const uint16 uint16_max = UINT16_MAX;
const uint32 uint32_max = UINT32_MAX;
const uint64 uint64_max = UINT64_MAX;


// 24-bit integer wrapper (for 24-bit PCM)
struct int24
{
	uint8 bytes[3];
	int24() { bytes[0] = bytes[1] = bytes[2] = 0; }
	explicit int24(int other)
	{
		#ifdef MPT_PLATFORM_BIG_ENDIAN
			bytes[0] = (static_cast<unsigned int>(other)>>16)&0xff;
			bytes[1] = (static_cast<unsigned int>(other)>> 8)&0xff;
			bytes[2] = (static_cast<unsigned int>(other)>> 0)&0xff;
		#else
			bytes[0] = (static_cast<unsigned int>(other)>> 0)&0xff;
			bytes[1] = (static_cast<unsigned int>(other)>> 8)&0xff;
			bytes[2] = (static_cast<unsigned int>(other)>>16)&0xff;
		#endif
	}
	operator int() const
	{
		#ifdef MPT_PLATFORM_BIG_ENDIAN
			return (static_cast<int8>(bytes[0]) * 65536) + (bytes[1] * 256) + bytes[2];
		#else
			return (static_cast<int8>(bytes[2]) * 65536) + (bytes[1] * 256) + bytes[0];
		#endif
	}
};
MPT_STATIC_ASSERT(sizeof(int24) == 3);
#define int24_min (0-0x00800000)
#define int24_max (0+0x007fffff)


typedef float float32;
MPT_STATIC_ASSERT(sizeof(float32) == 4);

typedef double float64;
MPT_STATIC_ASSERT(sizeof(float64) == 8);



namespace mpt {

MPT_STATIC_ASSERT(CHAR_BIT == 8);

MPT_STATIC_ASSERT(sizeof(char) == 1);

typedef unsigned char byte;
MPT_STATIC_ASSERT(sizeof(mpt::byte) == 1);

} // namespace mpt



#if MPT_COMPILER_MSVC

	#if defined(_M_X64)
		#define MPT_ARCH_BITS 64
		#define MPT_ARCH_BITS_32 0
		#define MPT_ARCH_BITS_64 1
	#elif defined(_M_IX86)
		#define MPT_ARCH_BITS 32
		#define MPT_ARCH_BITS_32 1
		#define MPT_ARCH_BITS_64 0
	#endif

#elif MPT_COMPILER_GCC || MPT_COMPILER_CLANG || MPT_COMPILER_MSVCCLANGC2

	#if defined(__SIZEOF_POINTER__)
		#if (__SIZEOF_POINTER__ == 8)
			#define MPT_ARCH_BITS 64
			#define MPT_ARCH_BITS_32 0
			#define MPT_ARCH_BITS_64 1
		#elif (__SIZEOF_POINTER__ == 4)
			#define MPT_ARCH_BITS 32
			#define MPT_ARCH_BITS_32 1
			#define MPT_ARCH_BITS_64 0
		#endif
	#endif

#endif // MPT_COMPILER

// fallback

#if !defined(MPT_ARCH_BITS)
#include <cstdint>
#include <stdint.h>
MPT_STATIC_ASSERT(sizeof(std::uintptr_t) == sizeof(void*));
#if defined(UINTPTR_MAX)
	#if (UINTPTR_MAX == 0xffffffffffffffffull)
		#define MPT_ARCH_BITS 64
		#define MPT_ARCH_BITS_32 0
		#define MPT_ARCH_BITS_64 1
	#elif (UINTPTR_MAX == 0xffffffffu)
		#define MPT_ARCH_BITS 32
		#define MPT_ARCH_BITS_32 1
		#define MPT_ARCH_BITS_64 0
	#endif
#endif // UINTPTR_MAX
#endif // MPT_ARCH_BITS



#if MPT_COMPILER_GCC || MPT_COMPILER_CLANG || MPT_COMPILER_MSVCCLANGC2
#define MPT_PRINTF_FUNC(formatstringindex,varargsindex) __attribute__((format(printf, formatstringindex, varargsindex)))
#else
#define MPT_PRINTF_FUNC(formatstringindex,varargsindex)
#endif



#if MPT_COMPILER_MSVC && defined(UNREFERENCED_PARAMETER)
#define MPT_UNREFERENCED_PARAMETER(x) UNREFERENCED_PARAMETER(x)
#else
#define MPT_UNREFERENCED_PARAMETER(x) (void)(x)
#endif

#define MPT_UNUSED_VARIABLE(x) MPT_UNREFERENCED_PARAMETER(x)



#if MPT_COMPILER_MSVC
// warning LNK4221: no public symbols found; archive member will be inaccessible
// There is no way to selectively disable linker warnings.
// #pragma warning does not apply and a command line option does not exist.
// Some options:
//  1. Macro which generates a variable with a unique name for each file (which means we have to pass the filename to the macro)
//  2. unnamed namespace containing any symbol (does not work for c++11 compilers because they actually have internal linkage now)
//  3. An unused trivial inline function.
// Option 3 does not actually solve the problem though, which leaves us with option 1.
// In any case, for optimized builds, the linker will just remove the useless symbol.
#define MPT_MSVC_WORKAROUND_LNK4221_CONCAT_DETAIL(x,y) x##y
#define MPT_MSVC_WORKAROUND_LNK4221_CONCAT(x,y) MPT_MSVC_WORKAROUND_LNK4221_CONCAT_DETAIL(x,y)
#define MPT_MSVC_WORKAROUND_LNK4221(x) int MPT_MSVC_WORKAROUND_LNK4221_CONCAT(mpt_msvc_workaround_lnk4221_,x) = 0;
#endif

#ifndef MPT_MSVC_WORKAROUND_LNK4221
#define MPT_MSVC_WORKAROUND_LNK4221(x)
#endif



// legacy
#define CountOf(x) MPT_ARRAYCOUNT(x)
#define forceinline MPT_FORCEINLINE
#define STATIC_ASSERT(x) MPT_STATIC_ASSERT(x)



OPENMPT_NAMESPACE_END
