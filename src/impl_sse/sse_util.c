/* SSE utility routines
 * 
 * Depending on how extensively these are used,
 * might want to roll them in to esl_sse.[h/c]
 */

#include "impl_sse.h"

/* Function:  alt_rightshift_ps()
 * Synopsis:  Shift vector elements to the right.
 * Incept:    SRE, Thu Jul 31 17:13:59 2008 [Janelia]
 * Alt:       DLK, Thu Feb 26 2009
 *
 * Purpose:   Returns a vector containing
 *            <{ b[3] a[0] a[1] a[2] }>:
 *            i.e. shift the values in <a> to the
 *            right, and load the LAST value of
 *            <b> into the first slot.
 */
inline __m128
alt_rightshift_ps(__m128 a, __m128 b)
{
  return _mm_move_ss(_mm_shuffle_ps(a, a, _MM_SHUFFLE(2, 1, 0, 0)), _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 3, 3, 3)));
}

#if 0
/* Function:  sse_leftshift_epi16()
 * Date:      DLK, Fri May 1 2009
 *
 * Purpose:   Returns a vector containing
 *            <{ a[1] a[2] a[3] a[4] a[5] a[6] a[7] b[0] }>:
 *            i.e., shift the values in <a> to the
 *            left, and load in the first value of
 *            <b> into the last slot.
 */
inline __m128i
sse_leftshift_epi16(__m128i a, __m128i b)
{
  register __m128i v = _mm_srli_si128(a, 2); /* now a[1] a[2] a[3] a[4] a[5] a[6] a[7] 0 */
  return _mm_insert_epi16(v,_mm_extract_epi16(b,0),7);
}

/* Function:  sse_rightshift_epi16()
 * Date:      DLK, Fri May 1 2009
 *
 * Purpose:   Returns a vector containing
 *            <{ b[7] a[0] a[1] a[2] a[3] a[4] a[5] a[6] }>:
 *            i.e., shift the values in <a> to the
 *            right, and load in the LAST value of
 *            <b> into the first slot.
 */
inline __m128i
sse_rightshift_epi16(__m128i a, __m128i b)
{
  register __m128i v = _mm_slli_si128(a, 2); /* now 0 a[0] a[1] a[2] a[3] a[4] a[5] a[6] */
  return _mm_insert_epi16(v,_mm_extract_epi16(b,7),0);
}
#endif

/* Function:  sse_setlw_neginfv()
 * Date:      DLK, Tue May 12 2009
 *
 * Purpose:   Returns a vector containing
 *            <{ -32768 a[1] a[2] a[3] a[4] a[5] a[6] a[7] }>
 *
 *            This is designed as a limited-use
 *            replacement for the call:
 *            _mm_insert_epi16(a, -32768, 0)
 *            which suffers from a compiler
 *            bug in gcc 3.4.x
 */
inline __m128i
sse_setlw_neginfv(__m128i a)
{
  __m128i mask = _mm_setr_epi16(0x0000,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff);
  __m128i setv = _mm_setr_epi16(-32768,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000);

  return _mm_or_si128(_mm_and_si128(a,mask),setv);
}

/* Function:  sse_select_si128()
 * Date:      DLK, Tue June 9 2009
 *
 * Purposes:  Bitwise vector select for __m128i
 */
inline __m128i
sse_select_si128(__m128i a, __m128i b, __m128i mask)
{
  b = _mm_and_si128(b, mask);
  a = _mm_andnot_si128(a, mask);
  return _mm_or_si128(a, b);
}

/* Function:  sse_hmax_epi16()
 * Date:      DLK, Tue June 9 2009
 *
 * Purposes:  Returns maximum value of 8 elements in epi16 vector
 */
inline int16_t
sse_hmax_epi16(__m128i a)
{
  a = _mm_max_epi16(a, _mm_srli_si128(a, 8));
  a = _mm_max_epi16(a, _mm_srli_si128(a, 4));
  a = _mm_max_epi16(a, _mm_srli_si128(a, 2));
  return _mm_extract_epi16(a, 0);
}

