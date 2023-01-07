#pragma once

#if SEAL_COMPILER == SEAL_COMPILER_MSVC

// Require Visual Studio 2015 or newer and C++14 support
#if (_MSC_VER < 1900) || (_MSVC_LANG < 201402L)
#error "Microsoft Visual Studio 2015 (14.0) or newer required; C++14 support required"
#endif

// Can work with Visual Studio 2015 (C++14) with some limitations
#if _MSVC_LANG == 1900
#undef SEAL_USE_SHARED_MUTEX_FOR_RW_LOCK
#endif

// Deal with Debug mode in Visual Studio
#ifdef _DEBUG
#define SEAL_DEBUG
#endif

// Cannot use std::shared_mutex when compiling with /clr
#ifdef _M_CEE
#undef SEAL_USE_SHARED_MUTEX_FOR_RW_LOCK
#endif

// Try to check presence of additional headers using __has_include
#ifdef __has_include

// Check for MSGSL
#if __has_include(<gsl/gsl>)
#include <gsl/gsl>
#define SEAL_USE_MSGSL
#else
#undef SEAL_USE_MSGSL
#endif //__has_include(<gsl/gsl>)

#endif

// X64
#ifdef _M_X64

// Use compiler intrinsics for better performance
#define SEAL_USE_INTRIN

#ifdef SEAL_USE_INTRIN
#include <intrin.h>

#pragma intrinsic(_addcarry_u64)
#define SEAL_ADD_CARRY_UINT64(operand1, operand2, carry, result) _addcarry_u64(     \
    carry,                                                                          \
    static_cast<unsigned long long>(operand1),                                      \
    static_cast<unsigned long long>(operand2),                                      \
    reinterpret_cast<unsigned long long*>(result))

#pragma intrinsic(_subborrow_u64)
#define SEAL_SUB_BORROW_UINT64(operand1, operand2, borrow, result) _subborrow_u64(  \
    borrow,                                                                         \
    static_cast<unsigned long long>(operand1),                                      \
    static_cast<unsigned long long>(operand2),                                      \
    reinterpret_cast<unsigned long long*>(result))

#pragma intrinsic(_BitScanReverse64)
#define SEAL_MSB_INDEX_UINT64(result, value) _BitScanReverse64(result, value)

#pragma intrinsic(_umul128)
#define SEAL_MULTIPLY_UINT64(operand1, operand2, result128) {                       \
    result128[0] = _umul128(                                                        \
        static_cast<unsigned long long>(operand1),                                  \
        static_cast<unsigned long long>(operand2),                                  \
        reinterpret_cast<unsigned long long*>(result128 + 1));                      \
}
#define SEAL_MULTIPLY_UINT64_HW64(operand1, operand2, hw64) {                       \
    _umul128(                                                                       \
        static_cast<unsigned long long>(operand1),                                  \
        static_cast<unsigned long long>(operand2),                                  \
        reinterpret_cast<unsigned long long*>(hw64));                               \
}

#endif
#else 
#undef SEAL_USE_INTRIN

#endif //_M_X64

#endif