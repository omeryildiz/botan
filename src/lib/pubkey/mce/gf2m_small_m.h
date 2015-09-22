/**
 * (C) Copyright Projet SECRET, INRIA, Rocquencourt
 * (C) Bhaskar Biswas and  Nicolas Sendrier
 *
 * (C) 2014 cryptosource GmbH
 * (C) 2014 Falko Strenzke fstrenzke@cryptosource.de
 *
 * Botan is released under the Simplified BSD License (see license.txt)
 *
 */

#ifndef BOTAN_GF2M_SMALL_M_H__
#define BOTAN_GF2M_SMALL_M_H__

#include <vector>
#include <botan/types.h>

namespace Botan {

typedef u16bit gf2m;

class GF2m_Field
   {
   public:
      GF2m_Field(size_t extdeg);

      gf2m gf_mul(gf2m x, gf2m y) const
         {
         return ((x) ? gf_mul_fast(x, y) : 0);
         }

      gf2m gf_square(gf2m x) const
         {
         return ((x) ? m_gf_exp_table[_gf_modq_1(m_gf_log_table[x] << 1)] : 0);
         }

      gf2m square_rr(gf2m x) const
         {
         return _gf_modq_1(x << 1);
         }

      gf2m gf_mul_fast(gf2m x, gf2m y) const
         {
         return ((y) ? m_gf_exp_table[_gf_modq_1(m_gf_log_table[x] + m_gf_log_table[y])] : 0);
         }

      /*
      naming convention of GF(2^m) field operations:
        l logarithmic, unreduced
        r logarithmic, reduced
        n normal, non-zero
        z normal, might be zero
      */

      gf2m gf_mul_lll(gf2m a, gf2m b) const
         {
         return  (a + b);
         }

      gf2m gf_mul_rrr(gf2m a, gf2m b) const
         {
         return (_gf_modq_1(gf_mul_lll(a, b)));
         }

      gf2m gf_mul_nrr(gf2m a, gf2m b) const
         {
         return (gf_exp(gf_mul_rrr(a, b)));
         }

      gf2m gf_mul_rrn(gf2m a, gf2m y) const
         {
         return _gf_modq_1(gf_mul_lll(a, gf_log(y)));
         }

      gf2m gf_mul_rnr(gf2m y, gf2m a) const
         {
         return gf_mul_rrn(a, y);
         }

      gf2m gf_mul_lnn(gf2m x, gf2m y) const
         {
         return (m_gf_log_table[x] + m_gf_log_table[y]);
         }

      gf2m gf_mul_rnn(gf2m x, gf2m y) const
         {
         return _gf_modq_1(gf_mul_lnn(x, y));
         }

      gf2m gf_mul_nrn(gf2m a, gf2m y) const
         {
         return m_gf_exp_table[_gf_modq_1((a) + m_gf_log_table[y])];
         }

      /**
      * zero operand allowed
      */
      gf2m gf_mul_zrz(gf2m a, gf2m y) const
         {
         return ( (y == 0) ? 0 : gf_mul_nrn(a, y) );
         }

      gf2m gf_mul_zzr(gf2m a, gf2m y) const
         {
         return gf_mul_zrz(y, a);
         }

      /**
      * non-zero operand
      */
      gf2m gf_mul_nnr(gf2m y, gf2m a) const
         {
         return gf_mul_nrn( a, y);
         }

      gf2m gf_sqrt(gf2m x) const
         {
         return ((x) ? m_gf_exp_table[_gf_modq_1(m_gf_log_table[x] << (m_gf_extension_degree-1))] : 0);
         }

      gf2m gf_div_rnn(gf2m x, gf2m y) const
         {
         return _gf_modq_1(m_gf_log_table[x] - m_gf_log_table[y]);
         }

      gf2m gf_div_rnr(gf2m x, gf2m b) const
         {
         return _gf_modq_1(m_gf_log_table[x] - b);
         }

      gf2m gf_div_nrr(gf2m a, gf2m b) const
         {
         return m_gf_exp_table[_gf_modq_1(a - b)];
         }

      gf2m gf_div_zzr(gf2m x, gf2m b) const
         {
         return ((x) ? m_gf_exp_table[_gf_modq_1(m_gf_log_table[x] - b)] : 0);
         }

      gf2m gf_inv(gf2m x) const
         {
         return m_gf_exp_table[gf_ord() - m_gf_log_table[x]];
         }
      gf2m gf_inv_rn(gf2m x) const
         {
         return (gf_ord() - m_gf_log_table[x]);
         }

      gf2m gf_square_ln(gf2m x) const
         {
         return m_gf_log_table[x] << 1;
         }

      gf2m gf_square_rr(gf2m a) const
         {
         return a << 1;
         }

      gf2m gf_l_from_n(gf2m x) const
         {
         return m_gf_log_table[x];
         }

      gf2m gf_div(gf2m x, gf2m y) const;

      gf2m gf_exp(gf2m i) const
         {
         return m_gf_exp_table[i]; /* alpha^i */
         }

      gf2m gf_log(gf2m i) const
         {
         return m_gf_log_table[i]; /* return i when x=alpha^i */
         }

      gf2m gf_ord() const
         {
         return m_gf_multiplicative_order;
         }

      gf2m get_extension_degree() const
         {
         return m_gf_extension_degree;
         }

      gf2m get_cardinality() const
         {
         return m_gf_cardinality;
         }

      gf2m gf_pow(gf2m x, int i)  const;

   private:
      gf2m _gf_modq_1(s32bit d) const
         {
         return  (((d) & gf_ord()) + ((d) >> m_gf_extension_degree));
         }

      gf2m m_gf_extension_degree, m_gf_cardinality, m_gf_multiplicative_order;
      const std::vector<gf2m>& m_gf_log_table;
      const std::vector<gf2m>& m_gf_exp_table;
   };

u32bit encode_gf2m(gf2m to_enc, byte* mem);

gf2m decode_gf2m(const byte* mem);

}

#endif
