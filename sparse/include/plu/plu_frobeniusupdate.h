#ifndef _frobeniusupdate_h_
#define _frobeniusupdate_h_

#include <math.h>

/**
 * frobenius_update - Update the couple (scale, sumsq) with one element when
 * computing the Froebnius norm.
 *
 * The frobenius norm is equal to scale * sqrt( sumsq ), this method allows to
 * avoid overflow in the sum square computation.
**/
static inline void
#if defined(PRECISION_d)
frobenius_update( int nb, double *scale, double *sumsq, double *value )
{
    double absval = fabs(*value);
    double ratio;
    if ( absval != 0. ){
        if ( (*scale) < absval ) {
            ratio = (*scale) / absval;
            *sumsq = (double)nb + (*sumsq) * ratio * ratio;
            *scale = absval;
        } else {
            ratio = absval / (*scale);
            *sumsq = (*sumsq) + (double)nb * ratio * ratio;
        }
    }
}
#endif

#endif /* _frobeniusupdate_h_ */
