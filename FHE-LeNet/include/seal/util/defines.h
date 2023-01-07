#pragma once

// For security reasons one should never throw when decoding fails due
// to overflow, but in some cases this might help in diagnosing problems.
//#undef SEAL_THROW_ON_DECODER_OVERFLOW

// Multiplication by a plaintext zero should not be allowed, and by
// default SEAL throws an error in this case. For performance reasons
// one might want to undefine this if appropriate checks are performed
// elsewhere.
//#undef SEAL_THROW_ON_MULTIPLY_PLAIN_BY_ZERO

// Expose thread-unsafe memory pools through the MemoryPoolHandle class.
//#undef SEAL_ENABLE_THREAD_UNSAFE_MEMORY_POOL

// Compile for big-endian system (not implemented)
#undef SEAL_BIG_ENDIAN

// Bound on the bit-length of user-defined moduli
#define SEAL_USER_MODULO_BIT_BOUND 60

// Bound on the number of coefficient moduli
#define SEAL_COEFF_MOD_COUNT_BOUND 62

// Bound on polynomial modulus degree
#define SEAL_POLY_MOD_DEGREE_BOUND 65536

// Maximum value for decomposition bit count
#define SEAL_DBC_MAX 60

// Minimum value for decomposition bit count
#define SEAL_DBC_MIN 1

// Debugging help
#define SEAL_ASSERT(condition) { if(!(condition)){ std::cerr << "ASSERT FAILED: "   \
    << #condition << " @ " << __FILE__ << " (" << __LINE__ << ")" << std::endl; } }

// String expansion 
#define _SEAL_STRINGIZE(x) #x
#define SEAL_STRINGIZE(x) _SEAL_STRINGIZE(x)

// Check that double is 64 bits
static_assert(sizeof(double) == 8, "sizeof(double) != 8");

// Use std::shared_mutex (C++17) for reader-writer lock (seal/util/lock.h) 
#define SEAL_USE_SHARED_MUTEX_FOR_RW_LOCK

// Detect compiler
#define SEAL_COMPILER_MSVC 1
#define SEAL_COMPILER_CLANG 2
#define SEAL_COMPILER_GCC 3

#if defined(_MSC_VER)
#define SEAL_COMPILER SEAL_COMPILER_MSVC
#elif defined(__clang__)
#define SEAL_COMPILER SEAL_COMPILER_CLANG
#elif defined(__GNUC__) && !defined(__clang__)
#define SEAL_COMPILER SEAL_COMPILER_GCC
#endif

// MSVC support
#include "seal/util/msvc.h"

// clang support
#include "seal/util/clang.h"

// gcc support
#include "seal/util/gcc.h"

// Use generic functions as (slower) fallback
#ifndef SEAL_ADD_CARRY_UINT64
#define SEAL_ADD_CARRY_UINT64(operand1, operand2, carry, result) add_uint64_generic(operand1, operand2, carry, result)
//#pragma message("SEAL_ADD_CARRY_UINT64 not defined. Using add_uint64_generic (see util/defines.h)")
#endif

#ifndef SEAL_SUB_BORROW_UINT64
#define SEAL_SUB_BORROW_UINT64(operand1, operand2, borrow, result) sub_uint64_generic(operand1, operand2, borrow, result)
//#pragma message("SEAL_SUB_BORROW_UINT64 not defined. Using sub_uint64_generic (see util/defines.h).")
#endif

#ifndef SEAL_MULTIPLY_UINT64
#define SEAL_MULTIPLY_UINT64(operand1, operand2, result128) {                      \
    multiply_uint64_generic(operand1, operand2, result128);                        \
}
//#pragma message("SEAL_MULTIPLY_UINT64 not defined. Using multiply_uint64_generic (see util/defines.h).")
#endif

#ifndef SEAL_MULTIPLY_UINT64_HW64
#define SEAL_MULTIPLY_UINT64_HW64(operand1, operand2, hw64) {                      \
    multiply_uint64_hw64_generic(operand1, operand2, hw64);                        \
}
//#pragma message("SEAL_MULTIPLY_UINT64 not defined. Using multiply_uint64_generic (see util/defines.h).")
#endif

#ifndef SEAL_MSB_INDEX_UINT64
#define SEAL_MSB_INDEX_UINT64(result, value) get_msb_index_generic(result, value)
//#pragma message("SEAL_MSB_INDEX_UINT64 not defined. Using get_msb_index_generic (see util/defines.h).")
#endif
