#ifndef _plu_atomic_h_
#define _plu_atomic_h_

#include <stdint.h>
#include <unistd.h>

/**
 * If the compiler provides atomic primitives we prefer to use
 * them instead of our own atomic assembly.
 */
#if defined(__FUJITSU)
  #undef HAVE_ATOMIC_XLC_32_BUILTINS
#endif
#if defined(HAVE_ATOMIC_XLC_32_BUILTINS)
    #include <stdio.h>
    static inline void plu_mfence( void )
    {
        __sync();
    }
    static inline int plu_atomic_bor_32b( volatile uint32_t* location,
                                            uint32_t value )
    {
        uint32_t old_value = __fetch_and_or(location, value);
        return old_value | value;
    }
    static inline int plu_atomic_band_32b( volatile uint32_t* location,
                                            uint32_t value )
    {
        uint32_t old_value = __fetch_and_and(location, value);
        return old_value & value;
    }
    static inline int plu_atomic_cas_32b( volatile uint32_t* location,
                                            uint32_t old_value,
                                            uint32_t new_value )
    {
        int32_t old = (int32_t)old_value;
        return __compare_and_swap( (volatile int*)location, &old, new_value );
    }
    #if defined(HAVE_ATOMIC_XLC_64_BUILTINS)
    static inline int plu_atomic_cas_64b( volatile uint64_t* location,
                                            uint64_t old_value,
                                            uint64_t new_value )
    {
        int64_t old = (int64_t)old_value;
        return __compare_and_swaplp( (volatile long*)location, &old, new_value );
    }
    #else
    static inline int plu_atomic_cas_64b( volatile uint64_t* location,
                                            uint64_t old_value,
                                            uint64_t new_value )
    {
        fprintf(stderr, "Use of 64b CAS using atomic-xlc without compiler support\n \n");
        return -1;
    }
    #endif
    #define PLU_ATOMIC_HAS_ATOMIC_ADD_32B
    static inline int32_t plu_atomic_add_32b( volatile int32_t *location, int32_t i )
    {
        register int32_t old_val, tmp_val;

        __sync();
        do {
            old_val = __lwarx( (volatile int*)location );
            tmp_val = old_val + i;
        } while( !__stwcx( (volatile int*)location, tmp_val ) );

        return( tmp_val );
    }
    #define PLU_ATOMIC_HAS_ATOMIC_SUB_32B
    static inline int32_t plu_atomic_sub_32b( volatile int32_t *location, int32_t i )
    {
        register int32_t old_val, tmp_val;

        __sync();
        do {
            old_val = __lwarx( (volatile int*)location );
            tmp_val = old_val - i;
        } while( !__stwcx( (volatile int*)location, tmp_val ) );
        return( tmp_val );
    }
#elif defined(MAC_OS_X)
    #ifndef MAC_OS_X
    #error This file should only be included on MAC OS X (Snow Leopard
    #endif
    #include <libkern/OSAtomic.h>
    static inline void plu_mfence( void )
    {
        OSMemoryBarrier();
    }
    static inline int plu_atomic_bor_32b( volatile uint32_t* location,
                                            uint32_t value )
    {
        return OSAtomicOr32( value, location );
    }
    static inline int plu_atomic_band_32b( volatile uint32_t* location,
                                            uint32_t value )
    {
        return OSAtomicAnd32( value, location );
    }
    static inline int plu_atomic_cas_32b( volatile uint32_t* location,
                                            uint32_t old_value,
                                            uint32_t new_value )
    {
        return OSAtomicCompareAndSwap32( old_value, new_value, (volatile int32_t*)location );
    }
    static inline int plu_atomic_cas_64b( volatile uint64_t* location,
                                            uint64_t old_value,
                                            uint64_t new_value )
    {
        return OSAtomicCompareAndSwap64( old_value, new_value, (volatile int64_t*)location );
    }
    #define PLU_ATOMIC_HAS_ATOMIC_INC_32B
    static inline int32_t plu_atomic_inc_32b( volatile int32_t *location )
    {
        return OSAtomicIncrement32( (int32_t*)location );
    }
    #define PLU_ATOMIC_HAS_ATOMIC_DEC_32B
    static inline int32_t plu_atomic_dec_32b( volatile int32_t *location )
    {
        return OSAtomicDecrement32( (int32_t*)location );
    }
    #define PLU_ATOMIC_HAS_ATOMIC_ADD_32B
    static inline int32_t plu_atomic_add_32b( volatile int32_t *location, int32_t i )
    {
        return OSAtomicAdd32( i, location );
    }
    #define PLU_ATOMIC_HAS_ATOMIC_SUB_32B
    static inline int32_t plu_atomic_sub_32b( volatile int32_t *location, int32_t i )
    {
        return OSAtomicAdd32( -i, location );
    }
    #if defined(HAVE_ATOMIC_GCC_128_BUILTINS)
    #define PLU_ATOMIC_HAS_ATOMIC_CAS_128B
    static inline int plu_atomic_cas_128b( volatile __int128_t* location,
                                            __int128_t old_value,
                                            __int128_t new_value )
    {
        return (__sync_bool_compare_and_swap(location, old_value, new_value) ? 1 : 0);
    }
    #endif
#elif defined(PLU_ARCH_PPC)
    #if defined(__bgp__)
        #ifndef __PPC
        #warning This file is only for PowerPC
        #endif  /* __ PPC */

        #ifndef __bgp__
        #warning This file is only for the BG/P
        #endif  /* __bgp__ */

        #ifndef PLU_ATOMIC_BGP_HAS_BEEN_INCLUDED
        #define PLU_ATOMIC_BGP_HAS_BEEN_INCLUDED

        #warning BGP atomic included

        #include <common/bgp_personality.h>
        #include <common/bgp_personality_inlines.h>
        /*#include <bpcore/ppc450_inlines.h>*/
        #include <assert.h>

        static inline void plu_mfence( void )
        {
            _bgp_msync();
        }

        static inline int plu_atomic_bor_32b( volatile uint32_t* location,
                                                uint32_t mask )
        {
            register uint32_t old_val, tmp_val;

            _bgp_msync();
            do {
                old_val = _bgp_LoadReserved( location );
                tmp_val = old_val | mask;
            } while( !_bgp_StoreConditional( location, tmp_val ) );

            return( tmp_val );
        }

        static inline int plu_atomic_band_32b( volatile uint32_t* location,
                                                uint32_t mask )
        {
            register uint32_t old_val, tmp_val;

            _bgp_msync();
            do {
                old_val = _bgp_LoadReserved( location );
                tmp_val = old_val & mask;
            } while( !_bgp_StoreConditional( location, tmp_val ) );

            return( tmp_val );
        }

        static inline int plu_atomic_cas_32b( volatile uint32_t* location,
                                                uint32_t old_value,
                                                uint32_t new_value )
        {
            uint32_t tmp_val;

            do {
                tmp_val = _bgp_LoadReserved( location );
                if( old_value != tmp_val ) {
                    old_value = tmp_val;
                    return( 0 );
                }
            } while( !_bgp_StoreConditional(location, new_value ) );

            return( 1 );
        }

        static inline int plu_atomic_cas_64b( volatile uint64_t* location,
                                                uint64_t old_value,
                                                uint64_t new_value )
        {
            assert(0);  /* Not supported */
        }

        #define PLU_ATOMIC_HAS_ATOMIC_ADD_32B
        static inline int32_t plu_atomic_add_32b( volatile int32_t *location, int32_t i )
        {
            register int32_t old_val, tmp_val;


            _bgp_msync();
            do {
                old_val = _bgp_LoadReserved( location );
                tmp_val = old_val + i;
            } while( !_bgp_StoreConditional( location, tmp_val ) );

            return( tmp_val );
        }

        #define PLU_ATOMIC_HAS_ATOMIC_SUB_32B
        static inline int32_t plu_atomic_sub_32b( volatile int32_t *location, int32_t i )
        {
            register int32_t old_val, tmp_val;

            _bgp_msync();
            do {
                old_val = _bgp_LoadReserved( location );
                tmp_val = old_val - i;
            } while( !_bgp_StoreConditional( location, tmp_val ) );

            return( tmp_val );
        }

        #endif  /* PLU_ATOMIC_BGP_HAS_BEEN_INCLUDED */
    #else
        #ifndef __PPC
        #warning This file is only for PowerPC
        #endif  /* __ PPC */

        static inline void plu_mfence( void )
        {
            __asm__ __volatile__ ("lwsync\n\t":::"memory");
        }

        static inline int plu_atomic_bor_32b( volatile uint32_t* location,
                                                uint32_t mask )
        {
        #if !defined(__IBMC__)
        int32_t old, t;

        __asm__ __volatile__(
                                "1:   lwarx   %0,  0, %3   \n\t"
                                "     or      %1, %0, %2   \n\t"
                                "     stwcx.  %1,  0, %3   \n\t"
                                "     bne-    1b           \n\t"
                                : "=&r" (old), "=&r" (t)
                                : "r" (mask), "r" (location)
                                : "cc", "memory");

        return t;
        #else
        return mask | __fetch_and_or(location, mask);
        #endif  /* !defined(__IBMC__) */
        }

        static inline int plu_atomic_band_32b( volatile uint32_t* location,
                                                uint32_t mask )
        {
        #if !defined(__IBMC__)
        int32_t old, t;

        __asm__ __volatile__(
                                "1:   lwarx   %0,  0, %3   \n\t"
                                "     andc    %1, %0, %2   \n\t"
                                "     stwcx.  %1,  0, %3   \n\t"
                                "     bne-    1b           \n\t"
                                : "=&r" (old), "=&r" (t)
                                : "r" (mask), "r" (location)
                                : "cc", "memory");

        return t;
        #else
        return mask & __fetch_and_and(location, mask);
        #endif  /* !defined(__IBMC__) */
        }

        static inline int plu_atomic_cas_32b( volatile uint32_t* location,
                                                uint32_t old_value,
                                                uint32_t new_value )
        {
        #if !defined(__IBMC__)
        int32_t ret;

        __asm__ __volatile__ (
                                "1: lwarx   %0, 0, %2  \n\t"
                                "   cmpw    0, %0, %3  \n\t"
                                "   bne-    2f         \n\t"
                                "   stwcx.  %4, 0, %2  \n\t"
                                "   bne-    1b         \n\t"
                                "2:"
                                : "=&r" (ret), "=m" (*location)
                                : "r" (location), "r" (old_value), "r" (new_value), "m" (*location)
                                : "cr0", "memory");

        return (ret == old_value);
        #else
        return __compare_and_swap((volatile int*)location, (int*)&old_value, (int)new_value);
        #endif  /* !defined(__IBMC__) */
        }

        static inline int plu_atomic_cas_64b( volatile uint64_t* location,
                                                uint64_t old_value,
                                                uint64_t new_value )
        {
        #if !defined(__IBMC__)
        int64_t ret;

        __asm__ __volatile__ (
                                "1: ldarx   %0, 0, %2  \n\t"
                                "   cmpd    0, %0, %3  \n\t"
                                "   bne-    2f         \n\t"
                                "   stdcx.  %4, 0, %2  \n\t"
                                "   bne-    1b         \n\t"
                                "2:"
                                : "=&r" (ret), "=m" (*location)
                                : "r" (location), "r" (old_value), "r" (new_value), "m" (*location)
                                : "cr0", "memory");

        return (ret == old_value);
        #else
        return __compare_and_swaplp((volatile long*)location, (long*)&old_value, (long)new_value);
        #endif  /* !defined(__IBMC__) */
        }

        #define PLU_ATOMIC_HAS_ATOMIC_INC_32B
        static inline int32_t plu_atomic_inc_32b( volatile int32_t *location )
        {
        #if !defined(__IBMC__)
        int32_t t;

        __asm__ __volatile__(
                                "1:   lwarx   %0, 0, %1    \n\t"
                                "     addic   %0, %0, 1    \n\t"
                                "     stwcx.  %0, 0, %1    \n\t"
                                "     bne-    1b           \n\t"
                                : "=&r" (t)
                                : "r" (location)
                                : "cc", "memory");

        return t;
        #else
        return 1 + __fetch_and_add( (volatile int*)location, 1);
        #endif  /* !defined(__IBMC__) */
        }

        #define PLU_ATOMIC_HAS_ATOMIC_DEC_32B
        static inline int32_t plu_atomic_dec_32b( volatile int32_t *location )
        {
        #if !defined(__IBMC__)
        int32_t t;

        __asm__ __volatile__(
                                "1:   lwarx   %0, 0,%1     \n\t"
                                "     addic   %0,%0,-1     \n\t"
                                "     stwcx.  %0,0,%1      \n\t"
                                "     bne-    1b           \n\t"
                                : "=&r" (t)
                                : "r" (location)
                                : "cc", "memory");

        return t;
        #else
        return __fetch_and_add( (volatile int*)location, -1) - 1;
        #endif  /* !defined(__IBMC__) */
        }

        #define PLU_ATOMIC_HAS_ATOMIC_ADD_32B
        static inline int32_t plu_atomic_add_32b( volatile int32_t *location, int32_t i )
        {
        #if !defined(__IBMC__)
        int32_t t;

        __asm__ __volatile__(
                                "1:   lwarx   %0, 0, %1    \n\t"
                                "     addic   %0, %0,%2    \n\t"
                                "     stwcx.  %0, 0, %1    \n\t"
                                "     bne-    1b           \n\t"
                                : "=&r" (t)
                                : "r" (location), "r" (i)
                                : "cc", "memory");

        return t;
        #else
        return i + __fetch_and_add( (volatile int*)location, i);
        #endif  /* !defined(__IBMC__) */
        }

        #define PLU_ATOMIC_HAS_ATOMIC_SUB_32B
        static inline int32_t plu_atomic_sub_32b( volatile int32_t *location, int32_t i )
        {
        #if !defined(__IBMC__)
        int32_t t;

        __asm__ __volatile__(
                                "1:   lwarx   %0, 0,%1     \n\t"
                                "     subf    %0,%0,%2     \n\t"
                                "     stwcx.  %0,0,%1      \n\t"
                                "     bne-    1b           \n\t"
                                : "=&r" (t)
                                : "r" (location), "r" (i)
                                : "cc", "memory");

        return t;
        #else
        return __fetch_and_add( (volatile int*)location, i) - i;
        #endif  /* !defined(__IBMC__) */
        }
    #endif
#elif defined(HAVE_ATOMIC_GCC_32_BUILTINS)
    #include <stdio.h>
    static inline void plu_mfence( void )
    {
        __sync_synchronize();
    }

    static inline int plu_atomic_bor_32b( volatile uint32_t* location,
                                            uint32_t value )
    {
        uint32_t old_value = __sync_fetch_and_or(location, value);
        return old_value | value;
    }

    static inline int plu_atomic_band_32b( volatile uint32_t* location,
                                            uint32_t value )
    {
        uint32_t old_value = __sync_fetch_and_and(location, value);
        return old_value & value;
    }

    static inline int plu_atomic_cas_32b( volatile uint32_t* location,
                                            uint32_t old_value,
                                            uint32_t new_value )
    {
        return (__sync_bool_compare_and_swap(location, old_value, new_value) ? 1 : 0);
    }

    #if defined(HAVE_ATOMIC_GCC_64_BUILTINS)
    static inline int plu_atomic_cas_64b( volatile uint64_t* location,
                                            uint64_t old_value,
                                            uint64_t new_value )
    {
        return (__sync_bool_compare_and_swap(location, old_value, new_value) ? 1 : 0);
    }
    #else
    static inline int plu_atomic_cas_64b( volatile uint64_t* location,
                                            uint64_t old_value,
                                            uint64_t new_value )
    {
        fprintf(stderr, "Use of 64b CAS using atomic-gcc without __GCC_HAVE_SYNC_COMPARE_AND_SWAP_8 set\n \n");
        (void)location; (void)old_value; (void)new_value;
        return -1;
    }
    #endif

    #if defined(HAVE_ATOMIC_GCC_128_BUILTINS)
    #define PLU_ATOMIC_HAS_ATOMIC_CAS_128B
    static inline int plu_atomic_cas_128b( volatile __uint128_t* location,
                                            __uint128_t old_value,
                                            __uint128_t new_value )
    {
        return (__sync_bool_compare_and_swap(location, old_value, new_value) ? 1 : 0);
    }
    #else
    static inline int plu_atomic_cas_128b( volatile uint64_t* location,
                                            uint64_t old_value,
                                            uint64_t new_value )
    {
        fprintf(stderr, "Use of 128b CAS using atomic-gcc without __GCC_HAVE_SYNC_COMPARE_AND_SWAP_16 set\n \n");
        (void)location; (void)old_value; (void)new_value;
        return -1;
    }
    #endif

    #define PLU_ATOMIC_HAS_ATOMIC_ADD_32B
    static inline int32_t plu_atomic_add_32b(volatile int32_t* location, int32_t i)
    {
        return __sync_add_and_fetch(location, i);
    }

    #define PLU_ATOMIC_HAS_ATOMIC_SUB_32B
    static inline int32_t plu_atomic_sub_32b(volatile int32_t* location, int32_t i)
    {
        return __sync_sub_and_fetch(location, i);
    }
#elif defined(PLU_ARCH_X86)
    static inline void plu_mfence(void)
    {
        __asm__ __volatile__ ("mfence\n\t":::"memory");
    }

    static inline int plu_atomic_cas_32b(volatile uint32_t* location,
                                        uint32_t old_value,
                                        uint32_t new_value)
    {
        unsigned char ret;
        __asm__ __volatile__ (
                            "lock; cmpxchgl %3,%4   \n\t"
                            "sete     %0      \n\t"
                            : "=qm" (ret), "=a" (old_value), "=m" (*location)
                            : "q"(new_value), "m"(*location), "1"(old_value)
                            : "memory", "cc");

        return (int)ret;
    }

    static inline int plu_atomic_bor_32b(volatile uint32_t* location,
                                        uint32_t value)
    {
        uint32_t old_value;

        do {
            old_value = *location;
        } while( !plu_atomic_cas_32b(location, old_value, (old_value|value) ));
        return old_value | value;
    }

    static inline int plu_atomic_band_32b(volatile uint32_t* location,
                                            uint32_t value)
    {
        uint32_t old_value;

        do {
            old_value = *location;
        } while( !plu_atomic_cas_32b(location, old_value, (old_value&value) ));
        return old_value & value;
    }

    #define ll_low(x)	*(((unsigned int *)&(x)) + 0)
    #define ll_high(x)	*(((unsigned int *)&(x)) + 1)

    static inline int plu_atomic_cas_64b(volatile uint64_t* location,
                                        uint64_t old_value,
                                        uint64_t new_value)
    {
    /*
        * Compare EDX:EAX with m64. If equal, set ZF and load ECX:EBX into
        * m64. Else, clear ZF and load m64 into EDX:EAX.
        */
        unsigned char ret;

        __asm__ __volatile__(
                        "push %%ebx            \n\t"
                        "movl %3, %%ebx        \n\t"
                        "lock cmpxchg8b (%4)  \n\t"
                        "sete %0               \n\t"
                        "pop %%ebx             \n\t"
                        : "=qm"(ret),"=a"(ll_low(old_value)), "=d"(ll_high(old_value))
                        : "D"(location), "1"(ll_low(old_value)), "2"(ll_high(old_value)),
                        "r"(ll_low(new_value)), "c"(ll_high(new_value))
                        : "cc", "memory", "ebx");
        return (int) ret;
    }

    #define PLU_ATOMIC_HAS_ATOMIC_INC_32B
    static inline int32_t plu_atomic_inc_32b(volatile int32_t *location)
    {
        __asm__ __volatile__ (
                            "lock; incl %0\n"
                            : "+m" (*(location)));
        return (*location);
    }

    #define PLU_ATOMIC_HAS_ATOMIC_DEC_32B
    static inline int32_t plu_atomic_dec_32b(volatile int32_t *location)
    {
        __asm__ __volatile__ (
                            "lock; decl %0\n"
                            : "+m" (*(location)));
        return (*location);
    }

    #define PLU_ATOMIC_HAS_ATOMIC_ADD_32B
    static inline int32_t plu_atomic_add_32(volatile int32_t* v, int32_t i)
    {
        int ret = i;
    __asm__ __volatile__(
                            "lock; xaddl %1,%0"
                            :"=m" (*v), "+r" (ret)
                            :"m" (*v)
                            :"memory", "cc");
    return (ret+i);
    }

    #define PLU_ATOMIC_HAS_ATOMIC_SUB_32B
    static inline int32_t plu_atomic_sub_32(volatile int32_t* v, int32_t i)
    {
        int ret = -i;
    __asm__ __volatile__(
                            "lock; xaddl %1,%0"
                            :"=m" (*v), "+r" (ret)
                            :"m" (*v)
                            :"memory", "cc");
    return (ret-i);
    }
#elif defined(PLU_ARCH_X86_64)
    static inline void plu_mfence(void)
    {
        __asm__ __volatile__ ("mfence\n\t":::"memory");
    }

    static inline int plu_atomic_cas_32b(volatile uint32_t* location,
                                        uint32_t old_value,
                                        uint32_t new_value)
    {
        unsigned char ret;
        __asm__ __volatile__ (
                            "lock; cmpxchgl %3,%4   \n\t"
                            "      sete     %0      \n\t"
                            : "=qm" (ret), "=a" (old_value), "=m" (*location)
                            : "q"(new_value), "m"(*location), "1"(old_value)
                            : "memory", "cc");

        return (int)ret;
    }

    static inline int plu_atomic_bor_32b(volatile uint32_t* location,
                                        uint32_t value)
    {
        uint32_t old_value;

        do {
            old_value = *location;
        } while( !plu_atomic_cas_32b(location, old_value, (old_value|value) ));
        return old_value | value;
    }

    static inline int plu_atomic_band_32b(volatile uint32_t* location,
                                            uint32_t value)
    {
        uint32_t old_value;

        do {
            old_value = *location;
        } while( !plu_atomic_cas_32b(location, old_value, (old_value&value) ));
        return old_value & value;
    }

    static inline int plu_atomic_cas_64b(volatile uint64_t* location,
                                        uint64_t old_value,
                                        uint64_t new_value)
    {
        unsigned char ret;
        __asm__ __volatile__ (
                            "lock; cmpxchgq %3,%4   \n\t"
                            "      sete     %0      \n\t"
                            : "=qm" (ret), "=a" (old_value), "=m" (*((volatile long*)location))
                            : "q"(new_value), "m"(*((volatile long*)location)), "1"(old_value)
                            : "memory", "cc");

    return (int)ret;
    }

    #define PLU_ATOMIC_HAS_ATOMIC_INC_32B
    static inline int32_t plu_atomic_inc_32b(volatile int32_t *location)
    {
        __asm__ __volatile__ (
                            "lock; incl %0\n"
                            : "+m" (*(location)));
        return (*location);
    }

    #define PLU_ATOMIC_HAS_ATOMIC_DEC_32B
    static inline int32_t plu_atomic_dec_32b(volatile int32_t *location)
    {
        __asm__ __volatile__ (
                            "lock; decl %0\n"
                            : "+m" (*(location)));
        return (*location);
    }

    #define PLU_ATOMIC_HAS_ATOMIC_ADD_32B
    static inline int32_t plu_atomic_add_32(volatile int32_t* v, int32_t i)
    {
        int ret = i;
    __asm__ __volatile__(
                            "lock; xaddl %1,%0"
                            :"=m" (*v), "+r" (ret)
                            :"m" (*v)
                            :"memory", "cc");
    return (ret+i);
    }

    #define PLU_ATOMIC_HAS_ATOMIC_SUB_32B
    static inline int32_t plu_atomic_sub_32(volatile int32_t* v, int32_t i)
    {
        int ret = -i;
    __asm__ __volatile__(
                            "lock; xaddl %1,%0"
                            :"=m" (*v), "+r" (ret)
                            :"m" (*v)
                            :"memory", "cc");
    return (ret-i);
    }
#else
#  error "No safe atomics available"
#endif

#include <assert.h>

static inline int plu_atomic_cas_xxb( volatile void* location,
                                         uint64_t old_value,
                                         uint64_t new_value,
                                         size_t type_size )
{
    switch(type_size){
    case 4:
        return plu_atomic_cas_32b( (volatile uint32_t*)location,
                                     (uint32_t)old_value, (uint32_t)new_value );
    case 8:
        return plu_atomic_cas_64b( (volatile uint64_t*)location,
                                     (uint64_t)old_value, (uint64_t)new_value );
    }
    return 0;
}

static inline uint64_t plu_atomic_bor_xxb( volatile void* location,
                                             uint64_t or_value,
                                             size_t type_size )
{
    assert( 4 == type_size );
    (void)type_size;
    return (uint64_t)plu_atomic_bor_32b( (volatile uint32_t*)location,
                                           (uint32_t)or_value);
}

#define plu_atomic_band(LOCATION, OR_VALUE)  \
    (__typeof__(*(LOCATION)))plu_atomic_band_xxb(LOCATION, OR_VALUE, sizeof(*(LOCATION)) )

#define plu_atomic_bor(LOCATION, OR_VALUE)  \
    (__typeof__(*(LOCATION)))plu_atomic_bor_xxb(LOCATION, OR_VALUE, sizeof(*(LOCATION)) )

#define plu_atomic_cas(LOCATION, OLD_VALUE, NEW_VALUE)               \
    plu_atomic_cas_xxb((volatile void*)(LOCATION),                   \
                         (uint64_t)(OLD_VALUE), (uint64_t)(NEW_VALUE), \
                         sizeof(*(LOCATION)))

#define plu_atomic_set_mask(LOCATION, MASK) plu_atomic_bor((LOCATION), (MASK))
#define plu_atomic_clear_mask(LOCATION, MASK)  plu_atomic_band((LOCATION), ~(MASK))

#ifndef PLU_ATOMIC_HAS_ATOMIC_INC_32B
#define PLU_ATOMIC_HAS_ATOMIC_INC_32B /* We now have it ! */

#ifdef PLU_ATOMIC_HAS_ATOMIC_ADD_32B
#define plu_atomic_inc_32b(l)  plu_atomic_add_32b((int32_t*)l, 1)
#else
static inline int32_t plu_atomic_inc_32b( volatile int32_t *location )
{
    uint32_t l;
    do {
        l = (uint32_t)*location;
    } while( !plu_atomic_cas_32b( location, l, l+1 ) );
    return (int32_t)l+1;
}
#endif  /* PLU_ATOMIC_HAS_ATOMIC_ADD_32B */
#endif  /* PLU_ATOMIC_HAS_ATOMIC_INC_32B */

#ifndef PLU_ATOMIC_HAS_ATOMIC_DEC_32B
#define PLU_ATOMIC_HAS_ATOMIC_DEC_32B /* We now have it ! */

#ifdef PLU_ATOMIC_HAS_ATOMIC_SUB_32B
#define plu_atomic_dec_32b(l)  plu_atomic_sub_32b((int32_t*)l, 1)
#else
static inline int32_t plu_atomic_dec_32b( volatile int32_t *location )
{
    uint32_t l;
    do {
        l = (uint32_t)*location;
    } while( !plu_atomic_cas_32b( location, l, l-1 ) );
    return (int32_t)l-1;
}
#endif  /* PLU_ATOMIC_HAS_ATOMIC_SUB_32B */
#endif  /* PLU_ATOMIC_HAS_ATOMIC_DEC_32B */

#ifndef PLU_ATOMIC_HAS_ATOMIC_ADD_32B
#define PLU_ATOMIC_HAS_ATOMIC_ADD_32B
static inline int32_t plu_atomic_add_32b( volatile int32_t *location, int32_t d )
{
    uint32_t l, n;
    do {
        l = *location;
        n = (uint32_t)((int32_t)l + d);
    } while( !plu_atomic_cas_32b( location, l, n ) );
    return n;
}
#endif /* PLU_ATOMIC_HAS_ATOMIC_ADD_32B */

typedef volatile uint32_t plu_atomic_lock_t;
/**
 * Enumeration of lock states
 */
enum {
    PLU_ATOMIC_UNLOCKED = 0,
    PLU_ATOMIC_LOCKED   = 1
};

static inline void plu_atomic_lock( plu_atomic_lock_t* atomic_lock )
{
    while( !plu_atomic_cas( atomic_lock, 0, 1) )
        /* nothing */;
}

static inline void plu_atomic_unlock( plu_atomic_lock_t* atomic_lock )
{
    plu_mfence();
    *atomic_lock = 0;
}

static inline long plu_atomic_trylock( plu_atomic_lock_t* atomic_lock )
{
    return plu_atomic_cas( atomic_lock, 0, 1 );
}

static inline uint64_t plu_atomic_add_64b( volatile uint64_t *location, uint64_t d )
{
    uint64_t l, n;
    do {
        l = (*location);
        n = (l + d);
    } while( !plu_atomic_cas_64b( location, l, n ) );
    return n;
}
#endif /* _plu_atomic_h_ */
