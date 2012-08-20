
/**
 * scm_mex_core.c
 *
 * This file contains functions for calculation of Spatial channel model for
 * Multiple Input Multiple Output (MIMO) simulations. It also has a mex gateway function
 * for Matlab use.
 *
 * Compilation:
 *
 *    In Matlab type:
 *
 *          mex scm_mex_core.c
 *
 * Instructions for input arguments is found in scm_mec_core.m help file.
 *
 * 
 * @author Jussi Salmi, Helsinki University of Technology, Radio Laboratory, SMARAD Center of Excellence
 * @date 2004/09/29
 * Jari (Apr 17, 2005): Modified power normalization of polarized option.
 *
 * Refs:
 *		[1] 3GPP TR 25.996 V6.1.0 (2003-09)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#define PI 3.14159265358979323846
#define DEFAULT_LUT_POINTS 16384 /* must be a power of two */
#define GENERAL 1
#define POLARIZED 2
#define LOS 3


/**
 * Function complex_multiply calculates the multiplication of two complex numbers.
 * The complex numbers must be in cartesian form of (a + ib)
 * (a + ib)(c + id) = (ac - bd) + i[(a + b)(c + d) - ac - bd] 
 * @param a real part of the first number
 * @param b imaginary part of the first number
 * @param c real part of the second number
 * @param d imaginary part of the second number
 * @param re_ans real part of the answer
 * @param im_ans imaginary part of the answer
 */
void complex_multiply(const double a, const double b, 
   const double c, const double d, double *re_ans, double *im_ans) {
      
   *re_ans = a*c-b*d;
   *im_ans = b*c+a*d;
}

/**
 * Function complex_multiply calculates the multiplication of three complex numbers.
 * @param ar real part of the first number
 * @param ai imaginary part of the first number
 * @param br real part of the second number
 * @param bi imaginary part of the second number
 * @param cr real part of the third number
 * @param ci imaginary part of the third number
 * @param re_ans real part of the answer
 * @param im_ans imaginary part of the answer
 */
void complex_multiply_3(const double ar, const double ai, 
   const double br, const double bi, const double cr, const double ci,
   double *re_ans, double *im_ans) {
      
   *re_ans = ar*br*cr-ar*bi*ci-ai*br*ci-ai*bi*cr;
   *im_ans = ar*br*ci+ar*bi*cr+ai*br*cr-ai*bi*ci;
}

/**
 * Function matrix_multiply calculates open a multiplication of [a][B][c], where
 * a and c are two-element complex vectors and B is a 2-by-2 complex matrix.
 * @param a1_re real part of a1
 * @param a1_im imaginary part of a1
 * @param a2_re real part of a2
 * @param a2_im imaginary part of a2
 * @param c1_re real part of c1
 * @param c1_im imaginary part of c1
 * @param c2_re real part of c2
 * @param c2_im imaginary part of c2
 * @param b11_re real part of b11
 * @param b11_im imaginary part of b11
 * @param b12_re real part of b12
 * @param b12_im imaginary part of b12
 * @param b21_re real part of b21
 * @param b21_im imaginary part of b21
 * @param b22_re real part of b22
 * @param b22_im imaginary part of b22
 * @param ans_real_part pointer to output real part
 * @param ans_imag_part pointer to output imaginary part
 */
void matrix_multiply(double a1_re,double a1_im, double a2_re, double a2_im,
   double c1_re, double c1_im, double c2_re, double c2_im,
   double b11_re, double b11_im, double b12_re, double b12_im,
   double b21_re, double b21_im, double b22_re, double b22_im,
   double *ans_real_part, double *ans_imag_part) {
   
   double real_part, imag_part;
   /* calculating the matrix multiplication open */
   complex_multiply_3(a1_re,a1_im,c1_re,c1_im,b11_re,b11_im,&real_part,&imag_part);
   *ans_real_part = real_part;
   *ans_imag_part = imag_part;
   complex_multiply_3(a2_re,a2_im,c1_re,c1_im,b21_re,b21_im,&real_part,&imag_part);
   *ans_real_part += real_part;
   *ans_imag_part += imag_part;
   complex_multiply_3(a1_re,a1_im,c2_re,c2_im,b12_re,b12_im,&real_part,&imag_part);
   *ans_real_part += real_part;
   *ans_imag_part += imag_part;
   complex_multiply_3(a2_re,a2_im,c2_re,c2_im,b22_re,b22_im,&real_part,&imag_part);
   *ans_real_part += real_part;
   *ans_imag_part += imag_part;
}

/**
 * Function sumloop calculates the repeating summation of subpaths.
 * @param m_i starting index of subpaths
 * @param stop final index +1
 * @param exp_helper table of t-independent terms in complex sinusoidials
 * @param exp_t_coeff table of t-dependent terms
 * @param real_multiplier real part of product of gain terms
 * @param imag_multiplier imaginary part of product of gain terms
 * @param t time instant
 * @param lm order of terms to be summed
 * @param real_sum pointer to real part of the sum
 * @param imag_sum pointer to imaginary part of the sum
 */
void sumloop(int m_i, int stop, const double *exp_helper,const double *exp_t_coeff,
   const double *real_multiplier, const double *imag_multiplier, double t, const double *lm, 
   double *real_sum, double *imag_sum) {

   int lm_i;
   double angle, a, b;

   while (m_i < stop) {
      lm_i = (int)lm[m_i]-1; /* order of terms, lm is indexed in matlab */
      angle = exp_helper[lm_i] + exp_t_coeff[lm_i]* t;
      complex_multiply(real_multiplier[lm_i], imag_multiplier[lm_i],
         cos(angle), sin(angle), &a, &b);
      *real_sum += a; /* real part of the multiplication */
      *imag_sum += b; /* imaginary part of the multiplication */
      m_i++;
   }
}

/**
 * Function sumloop_with_table calculates the repeating summation of subpaths using quantized cosine.
 * @param m_i starting index of subpaths
 * @param stop final index +1
 * @param exp_helper table of t-independent terms in complex sinusoidials
 * @param exp_t_coeff table of t-dependent terms
 * @param real_multiplier real part of product of gain terms
 * @param imag_multiplier imaginary part of product of gain terms
 * @param t time instant
 * @param lm order of terms to be summed
 * @param real_sum pointer to real part of the sum
 * @param imag_sum pointer to imaginary part of the sum
 * @param cos_table table of quantized cosine values
 * @param r2p coefficient for look-up table (num_of_points/(2*PI))
 * @param divider constant used by look-up table
 * @param halfpi (pi/2)
 */
void sumloop_with_table(int m_i, int stop, const double *exp_helper,const double *exp_t_coeff,
   const double *real_multiplier, const double *imag_multiplier, double t, const double *lm, 
   double *real_sum, double *imag_sum, const double *cos_table, const double r2p,
   const long int divider, const double halfpi) {
   
   int lm_i;
   double angle, a, b;

   while (m_i < stop) {
      lm_i = (int)lm[m_i]-1; /* order of terms, lm is indexed in matlab */
      /*printf("lm_i: %i\n", lm_i+1);*/
      angle = exp_helper[lm_i] + exp_t_coeff[lm_i]* t;
      complex_multiply(real_multiplier[lm_i], imag_multiplier[lm_i],
         cos_table[abs(angle*r2p)&divider], cos_table[abs((angle-halfpi)*r2p)&divider],
         &a, &b);
      *real_sum += a; /* real part of the multiplication */
      *imag_sum += b; /* imaginary part of the multiplication */
      m_i++;
   }
}

/**
 * Function scm_sum computes coefficients for 3GPP Spatial channel model [1 5.5.1.
 * The parameter angles must be in radians!
 * @param sin_look_up_points if nonzero, then look-up table for sin/cos is used. Use -1 for default number of points.
 * @param u the number of MS elements
 * @param s the number of BS elements
 * @param n number of multipaths
 * @param l number of midpaths
 * @param m number of subpaths for each n multipath
 * @param k number of individual links
 * @param *re_sq_G_BS real part of square root of BS Gain coefficients (size [k][s][n][m])
 * @param *im_sq_G_BS imaginary part of BS Gain coefficients (size [k][s][n][m])
 * @param *re_sq_G_MS real part of MS Gain coefficients (size [k][u][n][m])
 * @param *im_sq_G_MS imaginary part of MS Gain coefficients (size [k][u][n][m])
 * @param k_CONST wavenumber (2*pi/lambda)
 * @param *d_s distance of BS antenna s from ref. antenna (s=1), (size [s])
 * @param *d_u distance of MS antenna u from ref. antenna (u=1), (size [u])
 * @param *aod Angel of Departure (size [k][n][m])
 * @param *aoa Angel of Arrival (size [k][n][m])
 * @param *phase phase of the mth subpath of the nth path (size [k][n][m])
 * @param *v magnitude of the MS velocity vector (size [k])
 * @param *theta_v angle of the MS velocity vector (size [k])
 * @param *ts vector of time samples (size [k][tn])
 * @param tn number of time samples
 * @param *sq_Pn square root of Pn (size [k][n*l])
 * @param *ln number of subpaths per midpath (size [l])
 * @param *lm subpath indexing order for midpaths, Matlab indexing: 1-m (size [m])
 * @param GainsAreScalar has value 1 if gain is scalar, zero if not
 * @param *re_h pointer to real part of output h, h has to be initialized outside scm function (size [u][s][n][tn][k])
 * @param *im_h pointer to imag part of output h
 * @param *output_SubPathPhase pointer to output phases (size [k][n][m])
 * @return 0 if success, 1 if bad arguments
 */
int scm_sum(const long int sin_look_up_points, const int u, const int s, const int n, const int l, const int m, const int k, 
   const double *re_sq_G_BS, const double *im_sq_G_BS, const double *re_sq_G_MS, const double *im_sq_G_MS,
   const double k_CONST, const double *d_s, const double *d_u, const double *aod, const double *aoa,
   const double *phase, const double *v, const double *theta_v, const double *ts, const int tn,
   const double *sq_Pn, const double *ln, const double *lm, int GainsAreScalar, double *re_h, double *im_h, double *output_SubPathPhase) {

   int i, t_i, n_i, nl_i, l_i, m_i, u_i, s_i, k_i, m_count; /* running index variables */

   /* notation conversion variables from [][] to *(p[0]+var) */
   long int ksnm, kunm, knm, usnltk, kn, knl, kt, usnl, usnltk_i, ksn, nl, kun;

   int msreal, bsreal, both_real; /* boolean variables to check complexity of gains */
   double delta_t;
   double a, b, c, d, real_sum, imag_sum, kds, kdu, kv, gainScalar_re, gainScalar_im;
   double *real_multiplier;
   double *imag_multiplier;
   double *exp_helper; /* part of coefficient for exp(j...) in the sum */
   double *exp_t_coeff;
   double *cos_table;
   double r2p;
   double halfpi;
   double *one_over_sq_ln;
   
   long int num_of_points, divider, pows_of_two;

   real_multiplier = (double*)malloc(m*sizeof(double));
   imag_multiplier = (double*)calloc(m,sizeof(double));
   exp_t_coeff = (double*)malloc(m*sizeof(double));
   exp_helper = (double*)malloc(m*sizeof(double));

   /* check if gains are real valued */
   bsreal = 0;
   msreal = 0;
   both_real = 0;
   if (im_sq_G_BS == NULL)
      bsreal = 1;
   if (im_sq_G_MS == NULL)
      msreal = 1;
   if (bsreal && msreal)
      both_real = 1;

   /* check if gains are scalar */
   if (GainsAreScalar) {
      a = *re_sq_G_BS;
      c = *re_sq_G_MS;
      if (both_real) {
         gainScalar_re = a*c;
         gainScalar_im = 0;
      }
      else {
         if (bsreal)
            b = 0.0;
         else
            b = *im_sq_G_BS;
         if (msreal)
            d = 0.0;
         else
            d = *im_sq_G_MS;
         complex_multiply(a, b, c, d,
            &gainScalar_re, &gainScalar_im);
      }
   }
  
   nl = n*l;
   usnl = u*s*nl; /* for output pointing */
   ksn = k*s*n;  /* for parameter pointing */
   kun = k*u*n;
   
   
   /* calculate coefficient terms */
   one_over_sq_ln = (double*)malloc(l*sizeof(double));
   for(i=0;i<l;i++) {
      one_over_sq_ln[i] = 1/sqrt(ln[i]); 
   }

   /* checks if look-up table is used for sin/cos */
   if (sin_look_up_points) {
      if (sin_look_up_points == -1)
         num_of_points = DEFAULT_LUT_POINTS;
      else {
         pows_of_two = 2;
         while (pows_of_two < sin_look_up_points)
            pows_of_two = pows_of_two*2;
         num_of_points = pows_of_two;
         if (pows_of_two != sin_look_up_points)
            printf("Warning: Number of LU points is not a power-of-2: size changed to %li.\n", num_of_points);
      }

      divider = num_of_points-1;
      r2p = num_of_points/(2*PI);
      halfpi = PI*0.5;
      
      cos_table = (double*)malloc(num_of_points*sizeof(double));

      /* calculate the table */
      for (i=0;i<num_of_points;i++) { 
         cos_table[i] = cos((i+0.5)/r2p);
      }

      /* cycling links */
      for(k_i=0; k_i<k; k_i++) {
         if(tn>1) {
            delta_t = *(ts+k+k_i)-*(ts+k_i);
         }
         kv = k_CONST * v[k_i];
         /* cycling u, u is a given constant*/
         for(u_i=0; u_i < u; u_i++) {
            kdu = k_CONST * *(d_u+u_i);
            /* cyckling s, s is a given constant */
            for(s_i=0;s_i< s; s_i++) {
               kds = k_CONST * *(d_s+s_i);
               /* running through paths, n is a given constant*/
               nl_i = 0;
               for(n_i=0;n_i<n;n_i++) {
                  kn = k*n_i +k_i;
                  /* calculating m-variant t-invariant values */
                  for(m_i=0; m_i<m; m_i++) {                     
                     /* calculation of pointer parameters */
                     knm = k*n*m_i+kn;

                     /* calculation of sqrt(G_BS)*sqrt(G_MS): */
                     if (GainsAreScalar) {
                        real_multiplier[m_i] = gainScalar_re;
                        imag_multiplier[m_i] = gainScalar_im;
                     }
                     /* gains aren't scalar */
                     else {
                        ksnm = ksn*m_i+k*s*n_i+k*s_i+ k_i;
                        kunm = kun*m_i+k*u*n_i+k*u_i+k_i;
                        a = *(re_sq_G_BS + ksnm);
                        c = *(re_sq_G_MS + kunm);
                        if(both_real)
                           real_multiplier[m_i] = a*c;
                        else {
                           if (!bsreal) 
                              b = *(im_sq_G_BS + ksnm);
                           if (!msreal) 
                              d = *(im_sq_G_MS + kunm);
                           complex_multiply(a, b, c, d,
                              &real_multiplier[m_i], &imag_multiplier[m_i]);
                        }
                     }
                     exp_helper[m_i] = kds * cos_table[abs((*(aod+knm)-halfpi)*r2p)&divider]
                        + *(phase+knm) + kdu * cos_table[abs((*(aoa+knm)-halfpi)*r2p)&divider];

                     /* next value should be calculated accurately, because it is the main source of
                     * error when multiplied with large time sample values!
                     */
                     exp_t_coeff[m_i] =  kv * cos(*(aoa+knm)-theta_v[k_i]);	

                     if (tn >1) /* output pphase can be evaluated only if more than one time sample */
                        *(output_SubPathPhase+knm) = exp_t_coeff[m_i]* (*(ts+k*(tn-1)+k_i)+delta_t);
                  }  
                  
                  m_count = 0;
                  for (l_i=0;l_i<l;l_i++ ) {
                     knl = k*nl_i+k_i;
                     usnltk_i = usnl*tn*k_i + u*s*nl_i + u*s_i + u_i;
                     /* cycling the time vector */
                     for (t_i=0; t_i<tn; t_i++) {
                        real_sum = 0;
                        imag_sum = 0;
                        kt = k*t_i + k_i;
                        /* cycling m, now no calculations with m dependable variables is required */
                        m_i = m_count;
                        sumloop_with_table(m_i,m_count+(int)ln[l_i],exp_helper,exp_t_coeff, real_multiplier, imag_multiplier,
                           *(ts+kt), lm, &real_sum, &imag_sum, cos_table, r2p,divider, halfpi);
                        /* setting the output values and multiplying them with proper coefficient*/
                        usnltk = usnltk_i + usnl*t_i;
                        *(re_h+usnltk) = real_sum * *(sq_Pn+knl) * one_over_sq_ln[l_i];
                        *(im_h+usnltk) = imag_sum* *(sq_Pn+knl)* one_over_sq_ln[l_i];
                     }
                     m_count = m_count+(int)ln[l_i];
                     nl_i++;
                  }
               }
            }
         }
      }
      /* freeing memory of temporary tables */
      free(cos_table);
   }

   /* sin_look_up not used */
   else {

      /* cycling links */
      for(k_i=0; k_i<k; k_i++) {
         if(tn>1)
            delta_t = *(ts+k+k_i)-*(ts+k_i);
         kv = k_CONST * v[k_i];
         /* cycling u, u is a given constant*/
         for(u_i=0; u_i < u; u_i++) {
            kdu = k_CONST * *(d_u+u_i);
            /* cyckling s, s is a given constant */
            for(s_i=0;s_i< s; s_i++) {
               kds = k_CONST * *(d_s+s_i);
               /* running through paths, n is a given constant*/
               nl_i = 0;
               for(n_i=0;n_i<n;n_i++) {
                  kn = k*n_i +k_i;
                  /* calculating m-variant t-invariant values */
                  for(m_i=0; m_i<m; m_i++) {                     
                     /* calculation of pointer parameters */
                     knm = k*n*m_i+kn;

                     /* calculation of sqrt(G_BS)*sqrt(G_MS): */
                     if (GainsAreScalar) {
                        real_multiplier[m_i] = gainScalar_re;
                        imag_multiplier[m_i] = gainScalar_im;
                     }
                     /* gains aren't scalar */
                     else {
                        ksnm = ksn*m_i+k*s*n_i+k*s_i+ k_i;
                        kunm = kun*m_i+k*u*n_i+k*u_i+k_i;
                        a = *(re_sq_G_BS + ksnm);
                        c = *(re_sq_G_MS + kunm);
                        if(both_real)
                           real_multiplier[m_i] = a*c;
                        else {
                           if (!bsreal) 
                              b = *(im_sq_G_BS + ksnm);
                           if (!msreal) 
                              d = *(im_sq_G_MS + kunm);
                           complex_multiply(a, b, c, d,
                              &real_multiplier[m_i], &imag_multiplier[m_i]);
                        }
                     }

                     /* calculation of t-invariant part of exp(j...) in [1] 5.4 */
                     exp_helper[m_i] = kds * sin(*(aod+knm)) + *(phase+knm) +
                        kdu * sin(*(aoa+knm));

                     /* the coefficient just in front of t in exp(j...)*/
                     exp_t_coeff[m_i] = kv * cos(*(aoa+knm)-theta_v[k_i]);

                     if (tn >1) /* output pphase can be evaluated only if more than one time sample */
                     *(output_SubPathPhase+knm) = exp_t_coeff[m_i]* (*(ts+k*(tn-1)+k_i)+delta_t);
                  }

                  m_count = 0;
                  for (l_i=0;l_i<l;l_i++ ) {
                     knl = k*nl_i+k_i;
                     usnltk_i = usnl*tn*k_i + u*s*nl_i + u*s_i + u_i;
                     /* cycling the time vector */
                     for (t_i=0; t_i<tn; t_i++) {
                        real_sum = 0;
                        imag_sum = 0;
                        kt = k*t_i + k_i;

                        /* cycling m, now no calculations with m dependable variables is required */
                        m_i = m_count;
                        sumloop(m_i,m_count+(int)ln[l_i],exp_helper,exp_t_coeff, real_multiplier, imag_multiplier,
                           *(ts+kt), lm, &real_sum, &imag_sum);

                        /* setting the output values and multiplying them with proper coefficient*/
                        usnltk = usnltk_i + usnl*t_i;
                        *(re_h+usnltk) = real_sum * *(sq_Pn+knl) * one_over_sq_ln[l_i];
                        *(im_h+usnltk) = imag_sum* *(sq_Pn+knl)* one_over_sq_ln[l_i];
                     }
                     m_count = m_count+(int)ln[l_i];
                     nl_i++;
                  }
               }
            }
         }
      }
   }
   /* freeing memory of temporary tables */
   free(real_multiplier);
   free(imag_multiplier);
   free(exp_t_coeff);
   free(exp_helper);
   free(one_over_sq_ln);
   
   return 0;
}

/**
 * Function scm_pol_sum Creates array for Spatial channel model sums with cross-polarization.
 * The parameter angles must be in radians!
 * @param sin_look_up_points if nonzero, then look-up table for sin/cos is used. Use -1 for default number of points.
 * @param u the number of MS elements
 * @param s the number of BS elements
 * @param n number of multipaths
 * @param l number of midpaths
 * @param m number of subpaths for each n multipath
 * @param k the number of individual links
 * @param *re_X_BS_v real part of BS antenna V-pol component response (size [k][s][n][m])
 * @param *im_X_BS_v imaginary part of BS antenna V-pol component response (size [k][s][n][m])
 * @param *re_X_BS_h real part of BS antenna H-pol component response (size [k][s][n][m])
 * @param *im_X_BS_h imaginary part of BS antenna H-pol component response (size [k][s][n][m])
 * @param *re_X_MS_v real part of MS antenna V-pol component response (size [k][u][n][m])
 * @param *im_X_MS_v imaginary part of MS antenna V-pol component response (size [k][u][n][m])
 * @param *re_X_MS_h real part of MS antenna H-pol component response (size [k][u][n][m])
 * @param *im_X_MS_h imaginary part of MS antenna H-pol component response (size [k][u][n][m])
 * @param k wavenumber (2*pi/lambda)
 * @param *d_s distance of BS antenna s from ref. antenna (s=1), (size [s])
 * @param *d_u distance of MS antenna u from ref. antenna (u=1), (size [u])
 * @param *aod Angel of Departure (size [k][n][m])
 * @param *aoa Angel of Arrival (size [k][n][m])
 * @param *phase_v_v Phase offset of the mth subpath of the nth path between vertical components of BS and MS (size = [k][n][m])
 * @param *phase_v_h Phase offset of the mth subpath of the nth path between vertical components of BS and horizontal components of MS (size = [k][n][m])
 * @param *phase_h_v Phase offset of the mth subpath of the nth path between horizontal components of BS and vertical components of MS (size = [k][n][m])
 * @param *phase_h_h Phase offset of the mth subpath of the nth path between horizontal components of BS and MS (size = [k][n][m])
 * @param *sq_r_n1 square root of power ratio of (v-h)/(v-v) (size [k][n])
 * @param *sq_r_n2 square root of power ratio of (h-v)/(v-v) (size [k][n])
 * @param *v magnitude of the MS velocity vector (size [k])
 * @param *theta_v angle of the MS velocity vector (size [k])
 * @param *ts vectors of time samples (size [k][tn])
 * @param tn number of time samples
 * @param *sq_Pn square root of Pn (size [k][n*l])
 * @param *ln number of subpaths per midpath (size [l])
 * @param *lm subpath indexing order for midpaths, Matlab indexing: 1-m (size [m])
 * @param GainsAreScalar has value 1 if gain is scalar, zero if not
 * @param *re_h pointer to real part output matrice h, must be initialized outside scm function (size [u][s][n][t][k])
 * @param *im_h pointer to imag part output matrice h, must be initialized outside scm function (size [u][s][n][t][k])
 * @param *output_SubPathPhase (size [k][n][m])
 * @return 0 if success, 1 if bad arguments
 */
int scm_pol_sum(const long int sin_look_up_points, const int u, const int s, const int n, const int l, const int m, const int k,
   const double *re_X_BS_v, const double *im_X_BS_v, const double *re_X_BS_h, const double *im_X_BS_h, 
   const double *re_X_MS_v, const double *im_X_MS_v, const double *re_X_MS_h, const double *im_X_MS_h, 
   const double k_CONST, const double *d_s, const double *d_u, const double *aod, const double *aoa,
   const double *phase_v_v, const double *phase_v_h, const double *phase_h_v, const double *phase_h_h,
   const double *r_n1, const double *r_n2,	const double *v, const double *theta_v, const double *ts,
   const int tn, const double *sq_Pn, const double *ln, const double *lm, int GainsAreScalar,
   double *re_h, double *im_h, double *output_SubPathPhase) {

   int i, n_i, nl_i, l_i, m_i, u_i, s_i, t_i, k_i, m_count;	/* running index variables */
   /* notation conversion variables from [][] to *(p[0]+var) */
   long int ksnm, kunm, knm, usnltk, kn, knl, nl, kt, usnl, usnltk_i, ksn, kun;

   int bsreal, msreal, both_real;
   double a1_re, a1_im, a2_re, a2_im, c1_re, c1_im, c2_re, c2_im; /* variables for multiplications*/
   double b11_re, b11_im, b12_re, b12_im, b21_re, b21_im, b22_re, b22_im;
   double a1c1_re, a1c2_re, a2c1_re, a2c2_re;
   double tmp; /* temporary variable to lighten calculations */
   double kds, kdu, kv, real_sum, imag_sum;
   double *real_multiplier;
   double *imag_multiplier;
   double *exp_t_coeff, *exp_helper;
   double delta_t;
   double *cos_table;
   double r2p;
   double halfpi;
   long int num_of_points, divider, pows_of_two;
   double *one_over_sq_ln;
   
   double sq_r_n1, sq_r_n2;     /* (Jari, Apr 17, 2005) */
   
   real_multiplier = (double*)malloc(m*sizeof(double));
   imag_multiplier = (double*)calloc(m,sizeof(double));
   exp_t_coeff = (double*)malloc(m*sizeof(double));
   exp_helper = (double*)malloc(m*sizeof(double));

   /* check if input gains are not complex */
   bsreal = 0;
   msreal = 0;
   both_real = 0;
   if (im_X_BS_v == NULL)
      bsreal = 1;
   if (im_X_MS_v == NULL)
      msreal = 1;
   if(bsreal && msreal)
      both_real = 1;

   if (GainsAreScalar) {
      if (both_real) {
         a1c1_re = *re_X_BS_v * *re_X_MS_v;
         a2c1_re = *re_X_BS_h * *re_X_MS_v;
         a1c2_re = *re_X_BS_v * *re_X_MS_h;
         a2c2_re = *re_X_BS_h * *re_X_MS_h;
      }
   }

   nl = n*l;
   usnl = u*s*nl;
   ksn = k*s*n;
   kun = k*u*n;

   /* calculate coefficient terms */
   one_over_sq_ln = (double*)malloc(l*sizeof(double));
   for(i=0;i<l;i++) {
      one_over_sq_ln[i] = 1/sqrt(ln[i]); 
   }

   /* checks if look-up table is used for sin/cos */
   if (sin_look_up_points) {
      if (sin_look_up_points == -1)
         num_of_points = DEFAULT_LUT_POINTS;
      else {
         pows_of_two = 2;
         while (pows_of_two < sin_look_up_points)
            pows_of_two = pows_of_two*2;
         num_of_points = pows_of_two;
         if (pows_of_two != sin_look_up_points)
            printf("Warning: Number of LU points is not a power-of-2: size changed to %li.\n", num_of_points);
      }

      divider =  num_of_points -1;
      r2p = num_of_points/(2*PI);
      halfpi = PI*0.5;

      cos_table = (double*)malloc(num_of_points*sizeof(double));
      /* calculate the table */
      for (i=0;i<num_of_points;i++) 
         cos_table[i] = cos((i+0.5)/r2p);

      /* cycling links */
      for(k_i=0; k_i<k; k_i++) {
         if(tn>1)
            delta_t = *(ts+k+k_i)-*(ts+k_i);
         kv = k_CONST * v[k_i]; /* value depends only on k */
         /* u is a given constant */
         for (u_i=0; u_i<u; u_i++) {
            kdu = k_CONST * *(d_u+u_i); /* value depends only on u */
            /* s is a given constant */
            for (s_i=0; s_i<s; s_i++) {
               kds = k_CONST * *(d_s+s_i); /* value depends only on s */
               /* running through paths, n is a given constant */
               nl_i=0;
               for (n_i=0; n_i<n; n_i++) {
                  kn = k*n_i +k_i;
                  /* m is a given constant */
                  for (m_i=0; m_i<m; m_i++) {
                     /* calculation of pointer parameters */
                     knm = k*n*m_i+kn;
                     
                     /* Jari (Apr 17, 2005): changed power normalization */
							/* Pekka 28.11.2005, changed to per ray XPRs
                     sq_r_n1 = sqrt(*(r_n1+kn));
                     sq_r_n2 = sqrt(*(r_n2+kn)); */
                     sq_r_n1 = sqrt(*(r_n1+knm));
                     sq_r_n2 = sqrt(*(r_n2+knm));
                     
                     b12_re = sq_r_n1 * cos_table[abs((*(phase_v_h+knm))*r2p)&divider];
                     b12_im = sq_r_n1 * cos_table[abs((*(phase_v_h+knm)-halfpi)*r2p)&divider];
                     b11_re = cos_table[abs((*(phase_v_v+knm))*r2p)&divider];
                     b11_im = cos_table[abs((*(phase_v_v+knm)-halfpi)*r2p)&divider];                        
                     b21_re = sq_r_n2 * cos_table[abs((*(phase_h_v+knm))*r2p)&divider];
                     b21_im = sq_r_n2 * cos_table[abs((*(phase_h_v+knm)-halfpi)*r2p)&divider];
                     b22_re = cos_table[abs((*(phase_h_h+knm))*r2p)&divider];
                     b22_im = cos_table[abs((*(phase_h_h+knm)-halfpi)*r2p)&divider];
                     


                     /* Jari (Apr 17, 2005): Xpd independent power, not used in this version */
                     /*
                     tmp = 1/sqrt(1+*(r_n1+kn));
                     b12_re = tmp * cos_table[abs((*(phase_v_h+knm))*r2p)&divider];
                     b12_im = tmp * cos_table[abs((*(phase_v_h+knm)-halfpi)*r2p)&divider];
                     tmp = sqrt(*(r_n1+kn))*tmp;
                     b11_re =  tmp * cos_table[abs((*(phase_v_v+knm))*r2p)&divider];
                     b11_im = tmp * cos_table[abs((*(phase_v_v+knm)-halfpi)*r2p)&divider];
                     tmp = 1/sqrt(1+*(r_n2+kn));
                     b21_re = tmp * cos_table[abs((*(phase_h_v+knm))*r2p)&divider];
                     b21_im = tmp * cos_table[abs((*(phase_h_v+knm)-halfpi)*r2p)&divider];
                     tmp = sqrt(*(r_n2+kn))*tmp;
                     b22_re = tmp * cos_table[abs((*(phase_h_h+knm))*r2p)&divider];
                     b22_im = tmp * cos_table[abs((*(phase_h_h+knm)-halfpi)*r2p)&divider];
                     */

                     if(both_real) { 
                        if (GainsAreScalar) {
                           real_multiplier[m_i] = a1c1_re *b11_re + a2c1_re * b21_re +
                              a1c2_re*b12_re + a2c2_re*b22_re;
                           imag_multiplier[m_i] = a1c1_re*b11_im + a2c1_re*b21_im +
                              a1c2_re*b12_im + a2c2_re*b22_im;
                        }
                        else { /* aren't scalar, but real*/
                           ksnm = ksn*m_i+k*s*n_i+k*s_i+ k_i;
                           kunm = kun*m_i+k*u*n_i+k*u_i+k_i;
                           a1_re = *(re_X_BS_v+ksnm);
                           a2_re = *(re_X_BS_h+ksnm);
                           c1_re = *(re_X_MS_v+kunm);
                           c2_re = *(re_X_MS_h+kunm);
                           real_multiplier[m_i] = a1_re*c1_re*b11_re + a2_re*c1_re*b21_re + 
                              a1_re*c2_re*b12_re + a2_re*c2_re*b22_re;
                           imag_multiplier[m_i] = a1_re*c1_re*b11_im + a2_re*c1_re*b21_im +
                              a1_re*c2_re*b12_im + a2_re*c2_re*b22_im;
                        }
                     }
                     else { /* aren't real*/
                        ksnm = ksn*m_i+k*s*n_i+k*s_i+ k_i;
                        kunm = kun*m_i+k*u*n_i+k*u_i+k_i;
                        a1_re = *(re_X_BS_v+ksnm);
                        a2_re = *(re_X_BS_h+ksnm);
                        c1_re = *(re_X_MS_v+kunm);
                        c2_re = *(re_X_MS_h+kunm);
                        if (bsreal) {
                           a1_im= 0.0;
                           a2_im = 0.0;
                        }
                        else {
                           a1_im = *(im_X_BS_v+ksnm);
                           a2_im = *(im_X_BS_h+ksnm);
                        }
                        if (msreal) {
                           c1_im = 0.0;
                           c2_im = 0.0;
                        }
                        else {
                           c1_im = *(im_X_MS_v+kunm);
                           c2_im =  *(im_X_MS_h+kunm);
                        }
                        matrix_multiply(a1_re,a1_im,a2_re,a2_im,c1_re,c1_im,c2_re,c2_im,
                           b11_re,b11_im,b12_re,b12_im,b21_re,b21_im,b22_re,b22_im,
                           &real_multiplier[m_i], &imag_multiplier[m_i]);
                     }

                     /* time invariant exponent term */
                     exp_helper[m_i] = kds * cos_table[abs((*(aod+knm)-halfpi)*r2p)&divider]
                        + kdu * cos_table[abs((*(aoa+knm)-halfpi)*r2p)&divider];

                     /* the coefficient just in front of t in exp(j...)
                     * next value should be calculated accurately, because it is the main source of
                     * error when multiplied with large time sample values!
                     */
                     exp_t_coeff[m_i] =  kv * cos(*(aoa+knm)-theta_v[k_i]);	


                     if (tn >1) /* output pphase can be evaluated only if more than one time sample */
                     *(output_SubPathPhase+knm) = exp_t_coeff[m_i]* (*(ts+k*(tn-1)+k_i)+delta_t);
                  }
                  m_count = 0;
                  for (l_i=0;l_i<l;l_i++ ) {
                     knl = k*nl_i+k_i;
                     usnltk_i = usnl*tn*k_i + u*s*nl_i + u*s_i + u_i;
                     /* cycling the time vector */
                     for (t_i=0; t_i<tn; t_i++) {
                        real_sum = 0;
                        imag_sum = 0;
                        kt = k*t_i + k_i;
                        /* cycling m, now no calculations with m dependable variables is required */
                        m_i = m_count;
                        sumloop_with_table(m_i,m_count+(int)ln[l_i],exp_helper,exp_t_coeff, real_multiplier, imag_multiplier,
                           *(ts+kt), lm, &real_sum, &imag_sum, cos_table, r2p,divider, halfpi);
                        /* setting the output values and multiplying them with proper coefficient*/
                        usnltk = usnltk_i + usnl*t_i;
                        *(re_h+usnltk) = real_sum * *(sq_Pn+knl) * one_over_sq_ln[l_i];
                        *(im_h+usnltk) = imag_sum* *(sq_Pn+knl)* one_over_sq_ln[l_i];
                     }
                     m_count = m_count+(int)ln[l_i];
                     nl_i++;
                  }
               }
            }
         }
      }
      free(cos_table);
   }

   /* sin_look_up not used */
   else {

      /* cycling links */
      for(k_i=0; k_i<k; k_i++) {
         if(tn>1) {
            delta_t = *(ts+k+k_i)-*(ts+k_i);
         }
         kv = k_CONST * v[k_i]; /* value depends only on k */
         /* u is a given constant */
         for (u_i=0; u_i<u; u_i++) {
            kdu = k_CONST * *(d_u+u_i); /* value depends only on u */
            /* s is a given constant */
            for (s_i=0; s_i<s; s_i++) {
               kds = k_CONST * *(d_s+s_i); /* value depends only on s */
               /* running through paths, n is a given constant */
               nl_i = 0;
               for (n_i=0; n_i<n; n_i++) {
                  kn = k*n_i +k_i;
                  /* m is a given constant */
                  for (m_i=0; m_i<m; m_i++) {
                     /* calculation of pointer parameters */
                     knm = k*n*m_i+kn;
                     
                     
                     /* Jari (Apr 17, 2005): changed power normalization */
					 		/* Pekka 28.11.2005, changed to per ray XPRs
                     sq_r_n1 = sqrt(*(r_n1+kn));
                     sq_r_n2 = sqrt(*(r_n2+kn)); */
                     sq_r_n1 = sqrt(*(r_n1+knm));
                     sq_r_n2 = sqrt(*(r_n2+knm));

                     
                     b12_re = sq_r_n1 * cos(*(phase_v_h+knm));
                     b12_im = sq_r_n1 * sin(*(phase_v_h+knm));
                     b11_re = cos(*(phase_v_v+knm));
                     b11_im = sin(*(phase_v_v+knm));
                     b21_re = sq_r_n2 * cos(*(phase_h_v+knm));
                     b21_im = sq_r_n2 * sin(*(phase_h_v+knm));
                     b22_re = cos(*(phase_h_h+knm));
                     b22_im = sin(*(phase_h_h+knm));

                     
                     
                     /* Jari (Apr 17, 2005): Xpd independent power normalization, not used in this version */
                     /*
                     tmp = 1/sqrt(1+*(r_n1+kn));
                     b12_re = tmp * cos(*(phase_v_h+knm));
                     b12_im = tmp * sin(*(phase_v_h+knm));
                     tmp = sqrt(*(r_n1+kn))*tmp;
                     b11_re =  tmp * cos(*(phase_v_v+knm));
                     b11_im = tmp * sin(*(phase_v_v+knm));
                     tmp = 1/sqrt(1+*(r_n2+kn));
                     b21_re = tmp * cos(*(phase_h_v+knm));
                     b21_im = tmp * sin(*(phase_h_v+knm));
                     tmp = sqrt(*(r_n2+kn))*tmp;
                     b22_re = tmp * cos(*(phase_h_h+knm));
                     b22_im = tmp * sin(*(phase_h_h+knm));
                     */

                     if(both_real) { 
                        if (GainsAreScalar) {
                           real_multiplier[m_i] = a1c1_re *b11_re + a2c1_re * b21_re +
                              a1c2_re*b12_re + a2c2_re*b22_re;
                           imag_multiplier[m_i] = a1c1_re*b11_im + a2c1_re*b21_im +
                              a1c2_re*b12_im + a2c2_re*b22_im;
                        }
                        else { /* aren't scalar*/
                           ksnm = ksn*m_i+k*s*n_i+k*s_i+ k_i;
                           kunm = kun*m_i+k*u*n_i+k*u_i+k_i;
                           a1_re = *(re_X_BS_v+ksnm);
                           a2_re = *(re_X_BS_h+ksnm);
                           c1_re = *(re_X_MS_v+kunm);
                           c2_re = *(re_X_MS_h+kunm);
                           real_multiplier[m_i] = a1_re*c1_re*b11_re + a2_re*c1_re*b21_re + 
                              a1_re*c2_re*b12_re + a2_re*c2_re*b22_re;
                           imag_multiplier[m_i] = a1_re*c1_re*b11_im + a2_re*c1_re*b21_im +
                              a1_re*c2_re*b12_im + a2_re*c2_re*b22_im;
                        }
                     }
                     else { /* aren't real*/
                        ksnm = ksn*m_i+k*s*n_i+k*s_i+ k_i;
                        kunm = kun*m_i+k*u*n_i+k*u_i+k_i;
                        a1_re = *(re_X_BS_v+ksnm);
                        a2_re = *(re_X_BS_h+ksnm);
                        c1_re = *(re_X_MS_v+kunm);
                        c2_re = *(re_X_MS_h+kunm);
                        if (bsreal) {
                           a1_im= 0.0;
                           a2_im = 0.0;
                        }
                        else {
                           a1_im = *(im_X_BS_v+ksnm);
                           a2_im = *(im_X_BS_h+ksnm);
                        }
                        if (msreal) {
                           c1_im = 0.0;
                           c2_im = 0.0;
                        }
                        else {
                           c1_im = *(im_X_MS_v+kunm);
                           c2_im =  *(im_X_MS_h+kunm);
                        }
                        matrix_multiply(a1_re,a1_im,a2_re,a2_im,c1_re,c1_im,c2_re,c2_im,
                           b11_re,b11_im,b12_re,b12_im,b21_re,b21_im,b22_re,b22_im,
                           &real_multiplier[m_i], &imag_multiplier[m_i]);
                     }

                     /* time invariant exponent term */
                     exp_helper[m_i] = kds * sin(*(aod+knm)) + kdu * sin(*(aoa+knm));

                     /* the coefficient just in front of t in exp(j...)*/
                     exp_t_coeff[m_i] = kv * cos(*(aoa+knm)-theta_v[k_i]);

                     if (tn >1) /* output pphase can be evaluated only if more than one time sample */
                        *(output_SubPathPhase+knm) = exp_t_coeff[m_i]* (*(ts+k*(tn-1)+k_i)+delta_t);
                  }
                  m_count = 0;
                  for (l_i=0;l_i<l;l_i++ ) {
                      knl = k*nl_i+k_i;
                      usnltk_i = usnl*tn*k_i + u*s*nl_i + u*s_i + u_i;
                     /* cycling the time vector */
                     for (t_i=0; t_i<tn; t_i++) {
                        real_sum = 0;
                        imag_sum = 0;
                        kt = k*t_i + k_i;

                        /* cycling m, now no calculations with m dependable variables is required */
                        m_i = m_count;
                        sumloop(m_i,m_count+(int)ln[l_i],exp_helper,exp_t_coeff, real_multiplier, imag_multiplier,
                           *(ts+kt), lm, &real_sum, &imag_sum);

                        /* setting the output values and multiplying them with proper coefficient*/
                        usnltk = usnltk_i + usnl*t_i;
                        *(re_h+usnltk) = real_sum * *(sq_Pn+knl) * one_over_sq_ln[l_i];
                        *(im_h+usnltk) = imag_sum* *(sq_Pn+knl)* one_over_sq_ln[l_i];
                     }
                     m_count = m_count+(int)ln[l_i];
                     nl_i++;
                  }
               }
            }
         }
      }
   }
   /* freeing memory of temporary tables */
   free(real_multiplier);
   free(imag_multiplier);
   free(exp_t_coeff);
   free(exp_helper);
   free(one_over_sq_ln);

   return 0;
}

/**
 * Function scm_los creates a line of sight component.
 * The parameter angles must be in radians!
 * @param u the number of MS elements
 * @param s the number of BS elements
 * @param n number of paths
 * @param k the number of individual links
 * @param *re_G_BS real part of square root of BS Gain coefficients (size [k][s])
 * @param *im_G_BS imaginary part of BS Gain coefficients (size [k][s])
 * @param *re_G_MS real part of MS Gain coefficients (size [k][u])
 * @param *im_G_MS imaginary part of MS Gain coefficients (size [k][u])
 * @param k_CONST wavenumber (2*pi/lambda)
 * @param *d_s distance of BS antenna s from ref. antenna (s=1), (size [s])
 * @param *d_u distance of MS antenna u from ref. antenna (u=1), (size [u])
 * @param *theta_BS Angel of Departure for LOS (size [k])
 * @param *theta_MS Angel of Arrival for LOS (size [k])
 * @param phase_los phase of LOS component (size [k])
 * @param *v magnitude of the MS velocity vector (size [k])
 * @param *theta_v angle of the MS velocity vector (size [k])
 * @param *ts vector of time samples (size [k][tn])
 * @param tn number of time samples
 * @param *k_FACTOR the Ricean K factors (size [k])
 * @return 0 if success, 1 if bad arguments
 * @param *re_h pointer to real part of output los (size[u][s][t])
 * @param *im_h pointer to imag part output matrice h, re_h has to be initialized outside scm function
 * @param *output_los_phase output LOS phase (size [k])

 * @return 0 if success, 1 if bad arguments
 */
int scm_los(const int u, const int s, const int n, const int k, const double *re_G_BS, const double *im_G_BS,
   const double *re_G_MS, const double *im_G_MS,	const double k_CONST,
   const double *d_s, const double *d_u, const double *theta_BS, const double *theta_MS,
   const double *phase_los, const double *v, const double *theta_v, const double *ts, const int tn,
   double *k_FACTOR, double *re_h, double *im_h, double *output_los_phase) {

   int t_i, u_i, s_i, k_i, n_i;	/* running index variables */
   long int usntk_i, usn, usntk, ks, ku, kt; 	/* notation conversion variables from [][] to *(p[0]+var) */
   double delta_t;
   double a, b, angle;
   double real_multiplier, imag_multiplier, re_LOS, im_LOS;
   double exp_helper; /* part of coefficient for exp(j...) in the sum */
   double exp_t_coeff;
   double sq_1_over_K_plus_1, sq_K;


   for(k_i=0;k_i<k;k_i++) {
      /* check for LOS component */
      if (k_FACTOR[k_i] != 0) {
         
         /* the coefficients for terms */
         sq_1_over_K_plus_1 = sqrt(1/(k_FACTOR[k_i]+1));
         sq_K = sqrt(k_FACTOR[k_i]);

         /* the coefficient just in front of t in exp(j...)
         doesn't depend on u or s */
         exp_t_coeff = k_CONST * v[k_i] * cos(theta_MS[k_i]-theta_v[k_i]);

         /* output phase can be evaluated only if more than one time sample */
         if (tn >1) {
            delta_t = *(ts+k*1+k_i)-*(ts+k_i);
            *(output_los_phase+k_i) = exp_t_coeff * (*(ts+k*(tn-1)+k_i)+delta_t);
         }
         usn = u*s*n;

         /* u is a given constant*/
         for(u_i=0; u_i < u; u_i++) {
            /* s is a given constant */
            for(s_i=0;s_i< s; s_i++) {

               for(n_i=0; n_i<n; n_i++) {
                  
                  usntk_i = usn*tn*k_i + u*s*n_i + u*s_i + u_i;
                  /* for all components */
                  for(t_i=0; t_i<tn; t_i++) {
                     usntk = usntk_i + usn*t_i;
                     *(re_h + usntk) = sq_1_over_K_plus_1 * *(re_h + usntk);
                     *(im_h + usntk) = sq_1_over_K_plus_1 * *(im_h + usntk);
                  }

                  /* for LOS components */
                  if (n_i == 0) {

                     ks = k*s_i + k_i;
                     ku = k*u_i + k_i;

                     /* Multiplying gain components (t-invariant) */
                     /* check G_BS and G_MS for being real valued */
                     if (im_G_BS != NULL) {
                        if(im_G_MS != NULL) {
                           complex_multiply(*(re_G_BS+ks), *(im_G_BS+ks),
                              *(re_G_MS+ku), *(im_G_MS+ku), &a, &b);
                        }
                        /* only G_BS is complex */
                        else {
                           complex_multiply(*(re_G_BS+ks), *(im_G_BS+ks),
                              *(re_G_MS+ku), 0.0, &a, &b);
                        }
                     }

                     /* all is real */
                     else {
                        a = *(re_G_BS+ks) * *(re_G_MS+ku);
                        b = 0.0;
                     }

                     /* calculation of t-invariant part of exp(j...) in [1] 5.4 */
                     exp_helper = k_CONST * (*(d_s+s_i) * sin(theta_BS[k_i]) + 
                        *(d_u+u_i) * sin(theta_MS[k_i])) + phase_los[k_i];

                     /* product of previous two */
                     complex_multiply(a, b, cos(exp_helper), sin(exp_helper),
                        &real_multiplier, &imag_multiplier);

                     usntk_i = usn*tn*k_i + u*s*n_i + u*s_i + u_i;
                     for(t_i=0; t_i<tn; t_i++) {
                        usntk = usntk_i + usn*t_i;
                        kt = k*t_i + k_i;
                        angle = exp_t_coeff * *(ts+kt);
                        complex_multiply(real_multiplier, imag_multiplier,
                           cos(angle), sin(angle),
                           &re_LOS, &im_LOS);
                        *(re_h + usntk) += sq_K * sq_1_over_K_plus_1 * re_LOS;
                        *(im_h + usntk) += sq_K * sq_1_over_K_plus_1 * im_LOS;
                     }
                  }
               }
            }
         }
      }
   }
   return 0;
}



/**
 * This is the mex-function.
 *
 * Input arguments are:
 *
 * FOR GENERAL COEFFICIENTS MODE:
 *
 * prhs[0] = mode, either 1 (GENERAL), 2 (POLARIZED) or 3 (LOS) 
 * prhs[1] = [G_BS], complex BS_antenna gains (size [k][s][n][m])
 * prhs[2] = [G_MS], complex MS_antenna gains (size [k][u][n][m])
 * prhs[3] = [aod], angles of departure (size [k][n][m])
 * prhs[4] = [aoa], angles of arrival (size [k][n][m])
 * prhs[5] = [d_s], distances of BS antenna s from ref. antenna (s=1), (size [s])
 * prhs[6] = [d_u], distance of MS antenna u from ref. antenna (u=1), (size [u])
 * prhs[7] = [phase], phase of the mth subpath of the nth path (size [k][n][m])
 * prhs[8] = [ts], time samples (size [k][tn])
 * prhs[9] = k_CONST, wave number
 * prhs[10] = [v], magnitude of the MS velocity vector (size [k])
 * prhs[11] = [theta_v], angle of the MS velocity vector (size [k])
 * prhs[12] = [sq_Pn], square root of Pn (size [k][n*l])
 * prhs[13] = look_up_points, 0 if look-up table not used, -1 if default, otherwise the number of points
 * prhs[14] = u, number of MS antennas
 * prhs[15] = s, number of BS antennas
 * prhs[16] = n, number of multipaths
 * prhs[17] = l, number of midpaths
 * prhs[18] = m, number of subpaths
 * prhs[19] = k, number of links
 * prhs[20] = tn, number of time samples
 * prhs[21] = GainsAreScalar, nonzero if gains are scalar
 * prhs[22] = [lm], indexing vector for subpaths (size [m])
 * prhs[23] = [ln], number of subpaths per midpath (size [l])
 * Output arguments are 
 * plhs[0] = [h], output array of the channel coefficients (size [u][s][n][tn][k])
 * plhs[1] = [output_SubPathPhase] pointer to output phases (size [k][n][m])
 *
 * FOR POLARIZED MODE:
 *
 * prhs[0] = mode, either 1 (GENERAL), 2 (POLARIZED) or 3 (LOS) 
 * prhs[1] = [X_BS_v] BS antenna V-pol component response (size[k][s][n][m])
 * prhs[2] = [X_BS_h] BS antenna H-pol component response (size[k][s][n][m])
 * prhs[3] = [X_MS_v] MS antenna V-pol component response (size[k][u][n][m])
 * prhs[4] = [X_MS_h] MS antenna H-pol component response (size[k][u][n][m])
 * prhs[5] = [aod], angles of departure (size[k][n][m])
 * prhs[6] = [aoa], angles of arrival (size[k][n][m])
 * prhs[7] = [d_s], distances of BS antenna s from ref. antenna (s=1), (size [s])
 * prhs[8] = [d_u], distance of MS antenna u from ref. antenna (u=1), (size [u])
 * prhs[9] = [phase_v_v], Phase offset of the mth subpath of the nth path between vertical components of BS and MS (size = [k][n][m])
 * prhs[10] = [phase_v_h], Phase offset of the mth subpath of the nth path between vertical components of BS and horizontal components of MS (size = [k][n][m])
 * prhs[11] = [phase_h_v], Phase offset of the mth subpath of the nth path between horizontal components of BS and vertical components of MS (size = [k][n][m])
 * prhs[12] = [phase_h_h], Phase offset of the mth subpath of the nth path between horizontal components of BS and MS (size = [k][n][m])
 * prhs[13] = [r_n1] power ratio of (v-h)/(v-v) (size [k][n])
 * prhs[14] = [r_n2] power ratio of (h-v)/(v-v) (size [k][n])
 * prhs[15] = [ts], time sample vector (size [k][tn])
 * prhs[16] = k_CONST, wave number
 * prhs[17] = [v], magnitude of the MS velocity vector (size [k])
 * prhs[18] = [theta_v], angle of the MS velocity vector (size [k])
 * prhs[19] = [sq_Pn], square root of Pn (size [k][n*l])
 * prhs[20] = look_up_points, 0 if look-up table not used, -1 if default, otherwise the number of points
 * prhs[21] = u, number of MS antennas
 * prhs[22] = s, number of BS antennas
 * prhs[23] = n, number of multipaths
 * prhs[24] = l, number of midpaths
 * prhs[25] = m, number of subpaths
 * prhs[26] = k, number of links
 * prhs[27] = tn, number of time samples
 * prhs[28] = GainsAreScalar, nonzero if gains are scalar
 * prhs[29] = [lm], indexing vector for subpaths (size [m])
 * prhs[30] = [ln], number of subpaths per midpath (size [l])
 * Output arguments are 
 * plhs[0] = [h], output array of the channel coefficients (size [u][s][n][tn])
 * plhs[1] = [output_SubPathPhase] pointer to output phases (size [n][m])
 *
 * FOR LOS MODE:
 *
 * prhs[0] = mode, either 1 (GENERAL), 2 (POLARIZED) or 3 (LOS) 
 * prhs[1] = [G_BS], complex BS_antenna gains (size [k][s])
 * prhs[2] = [G_MS], complex MS_antenna gains (size [k][u])
 * prhs[3] = theta_BS, angles of departure (size [k])
 * prhs[4] = theta_MS, angles of arrival (size [k])
 * prhs[5] = [d_s], distances of BS antenna s from ref. antenna (s=1), (size [s])
 * prhs[6] = [d_u], distance of MS antenna u from ref. antenna (u=1), (size [u])
 * prhs[7] = phase_los, phase of LOS component (size [k])
 * prhs[8] = [ts], time sample vector (size [k][tn])
 * prhs[9] = k_CONST, wave number
 * prhs[10] = [v], magnitude of the MS velocity vector (size [k])
 * prhs[11] = [theta_v], angle of the MS velocity vector (size[k])
 * prhs[12] = [h_in], input matrices (size [u][s][n][tn][k])
 * prhs[13] = [input_phase_los], (size [k])
 * prhs[14] = [k_FACTOR], (size [k])
 * prhs[15] = u, number of MS antennas
 * prhs[16] = s, number of BS antennas
 * prhs[17] = n, number of multipaths
 * prhs[18] = k, number of links
 * prhs[19] = tn, number of time samples
 * Output arguments are 
 * plhs[0] = [h], output array of the coefficients with LOS terms (size [u][s][n][tn][k])
 * plhs[1] = [output_phase_los] pointer to output phases (size [k])

 *
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

   /*pointers for input parameters*/

   /* variables for general mode coefficients */
   double *re_G_BS, *im_G_BS, *re_G_MS, *im_G_MS, *phase;

   /* variables for polarized mode coefficients */
   double *re_X_BS_v, *im_X_BS_v, *re_X_BS_h, *im_X_BS_h, *re_X_MS_v, *im_X_MS_v, *re_X_MS_h, *im_X_MS_h,
      *phase_v_v, *phase_v_h, *phase_h_v, *phase_h_h, *r_n1, *r_n2;

   /* variables for los mode */
   double *k_FACTOR, *theta_BS, *theta_MS, *phase_los;

   /* variables for more or less common parameters */
   double k_CONST;
   double *d_s, *d_u, *aod, *aoa, *sq_Pn, *v, *theta_v;
   double *ts;
   double *lm, *ln;  /* subpath ordering */
   int u, s, n, l, m, tn, k; /* parameters to be extracted from inputs' dimensions */
   int h_dims[5]; /* table for output array dimensions */
   int out_Phase_dims[3];
   long int look_up_points;
   int error, mode, GainsAreScalar;

   /*pointers to output*/
   double *re_h, *im_h, *output_SubPathPhase, *output_los_phase;

   /* check which mode is used */
   mode = (int)mxGetScalar(prhs[0]);


   /*Check if general coefficients are wanted */
   if (mode == GENERAL) {

      if(nrhs != 24)
         mexErrMsgTxt("scm_mex_core error: exactly 24 input arguments required for GENERAL option");

      /* ripping values from input */
      re_G_BS = (double*)mxGetPr(prhs[1]);
      im_G_BS = (double*)mxGetPi(prhs[1]);
      re_G_MS = (double*)mxGetPr(prhs[2]);
      im_G_MS = (double*)mxGetPi(prhs[2]);	
      aod = (double*)mxGetPr(prhs[3]);
      aoa = (double*)mxGetPr(prhs[4]);
      d_s = (double*)mxGetPr(prhs[5]);
      d_u = (double*)mxGetPr(prhs[6]);
      phase = (double*)mxGetPr(prhs[7]);
      ts = (double*)mxGetPr(prhs[8]);
      k_CONST = (double)mxGetScalar(prhs[9]);
      v = (double*)mxGetPr(prhs[10]);
      theta_v = (double*)mxGetPr(prhs[11]);
      sq_Pn = (double*)mxGetPr(prhs[12]);
      look_up_points = (long int)mxGetScalar(prhs[13]);
      u = (int)mxGetScalar(prhs[14]);
      s = (int)mxGetScalar(prhs[15]);
      n = (int)mxGetScalar(prhs[16]);
      l = (int)mxGetScalar(prhs[17]);
      m = (int)mxGetScalar(prhs[18]);
      k = (int)mxGetScalar(prhs[19]);
      tn = (int)mxGetScalar(prhs[20]);
      GainsAreScalar = (int)mxGetScalar(prhs[21]);
      lm = (double*)mxGetPr(prhs[22]);
      ln = (double*)mxGetPr(prhs[23]);

      /* now we know the output dimensions */
      h_dims[0] = u;
      h_dims[1] = s;
      h_dims[2] = n*l;
      h_dims[3] = tn;
      h_dims[4] = k;
      out_Phase_dims[0] = k;
      out_Phase_dims[1] = n;
      out_Phase_dims[2] = m;

      /* setting up output */
      plhs[0] = mxCreateNumericArray(5, h_dims, mxDOUBLE_CLASS, mxCOMPLEX);
      plhs[1] = mxCreateNumericArray(3, out_Phase_dims, mxDOUBLE_CLASS, mxREAL);
      re_h = (double*)mxGetPr(plhs[0]); /* Create a C pointer to real values of the output array */
      im_h = (double*)mxGetPi(plhs[0]); /* Create a C pointer to imaginary values of output */
      output_SubPathPhase = (double*)mxGetPr(plhs[1]);

      /* calculating the values */
      error = scm_sum(look_up_points, u, s, n, l, m, k, re_G_BS, im_G_BS, re_G_MS, im_G_MS,
         k_CONST, d_s, d_u, aod, aoa, phase, v, theta_v, ts, tn, sq_Pn, ln, lm, GainsAreScalar,
         re_h, im_h, output_SubPathPhase);

   } /* end general coefficients mode */

   /* check for polarized mode */
   else if (mode == POLARIZED)	{

      if(nrhs != 31)
         mexErrMsgTxt("scmmex error: exactly 31 input arguments required");

      /* ripping values from input */
      re_X_BS_v = (double*)mxGetPr(prhs[1]);
      im_X_BS_v = (double*)mxGetPi(prhs[1]);
      re_X_BS_h = (double*)mxGetPr(prhs[2]);
      im_X_BS_h = (double*)mxGetPi(prhs[2]);	
      re_X_MS_v = (double*)mxGetPr(prhs[3]);
      im_X_MS_v = (double*)mxGetPi(prhs[3]);
      re_X_MS_h = (double*)mxGetPr(prhs[4]);
      im_X_MS_h = (double*)mxGetPi(prhs[4]);
      aod = (double*)mxGetPr(prhs[5]);
      aoa = (double*)mxGetPr(prhs[6]);
      d_s = (double*)mxGetPr(prhs[7]);
      d_u = (double*)mxGetPr(prhs[8]);
      phase_v_v = (double*)mxGetPr(prhs[9]);
      phase_v_h = (double*)mxGetPr(prhs[10]);
      phase_h_v = (double*)mxGetPr(prhs[11]);
      phase_h_h = (double*)mxGetPr(prhs[12]);
      r_n1 = (double*)mxGetPr(prhs[13]);
      r_n2 = (double*)mxGetPr(prhs[14]);
      ts = (double*)mxGetPr(prhs[15]);
      k_CONST = (double)mxGetScalar(prhs[16]);
      v = (double*)mxGetPr(prhs[17]);
      theta_v = (double*)mxGetPr(prhs[18]);
      sq_Pn = (double*)mxGetPr(prhs[19]);
      look_up_points = (long int)mxGetScalar(prhs[20]);
      u = (int)mxGetScalar(prhs[21]);
      s = (int)mxGetScalar(prhs[22]);
      n = (int)mxGetScalar(prhs[23]);
      l = (int)mxGetScalar(prhs[24]);
      m = (int)mxGetScalar(prhs[25]);
      k = (int)mxGetScalar(prhs[26]);
      tn = (int)mxGetScalar(prhs[27]);
      GainsAreScalar = (int)mxGetScalar(prhs[28]);
      lm = (double*)mxGetPr(prhs[29]);
      ln = (double*)mxGetPr(prhs[30]);

      /* now we know the output dimensions */
      h_dims[0] = u;
      h_dims[1] = s;
      h_dims[2] = n*l;
      h_dims[3] = tn;
      h_dims[4] = k;
      out_Phase_dims[0] = k;
      out_Phase_dims[1] = n;
      out_Phase_dims[2] = m;

      /* setting up output */
      plhs[0] = mxCreateNumericArray(5, h_dims, mxDOUBLE_CLASS, mxCOMPLEX);
      plhs[1] = mxCreateNumericArray(3, out_Phase_dims, mxDOUBLE_CLASS, mxREAL);
      re_h = (double*)mxGetPr(plhs[0]); /* Create a C pointer to real values of the output array */
      im_h = (double*)mxGetPi(plhs[0]); /* Create a C pointer to imaginary values of output */
      output_SubPathPhase = (double*)mxGetPr(plhs[1]);

      /* calculating the values */
      error = scm_pol_sum(look_up_points, u, s, n, l, m, k, re_X_BS_v, im_X_BS_v, re_X_BS_h, im_X_BS_h, 
         re_X_MS_v, im_X_MS_v, re_X_MS_h, im_X_MS_h, k_CONST, d_s, d_u, aod, aoa,
         phase_v_v, phase_v_h, phase_h_v, phase_h_h,	r_n1, r_n2,
         v, theta_v, ts, tn, sq_Pn, ln, lm, GainsAreScalar, re_h, im_h, output_SubPathPhase);
   } /* end polarized mode */


   /* check los mode */
   else if (mode == LOS) {

      if(nrhs != 20)
         mexErrMsgTxt("scmlosmex error: exactly 20 input arguments required");

      /* ripping values from input */
      re_G_BS = (double*)mxGetPr(prhs[1]);
      im_G_BS = (double*)mxGetPi(prhs[1]);
      re_G_MS = (double*)mxGetPr(prhs[2]);
      im_G_MS = (double*)mxGetPi(prhs[2]);
      theta_BS = (double*)mxGetPr(prhs[3]);
      theta_MS = (double*)mxGetPr(prhs[4]);
      d_s = (double*)mxGetPr(prhs[5]);
      d_u = (double*)mxGetPr(prhs[6]);
      phase_los = (double*)mxGetPr(prhs[7]);
      ts = (double*)mxGetPr(prhs[8]);
      k_CONST = (double)mxGetScalar(prhs[9]);
      v = (double*)mxGetPr(prhs[10]);
      theta_v = (double*)mxGetPr(prhs[11]);
      re_h = (double*)mxGetPr(prhs[12]);
      im_h = (double*)mxGetPi(prhs[12]);
      output_los_phase = (double*)mxGetPr(prhs[13]);
      k_FACTOR = (double*)mxGetPr(prhs[14]);
      u = (int)mxGetScalar(prhs[15]);
      s = (int)mxGetScalar(prhs[16]);
      n = (int)mxGetScalar(prhs[17]);
      k = (int)mxGetScalar(prhs[18]);
      tn = (int)mxGetScalar(prhs[19]);

      /* setting up output */
      plhs[0] = (mxArray*)prhs[12];
      plhs[1] = (mxArray*)(prhs[13]);

      /* calculating the values */
      error = scm_los(u, s, n, k, re_G_BS, im_G_BS,
         re_G_MS, im_G_MS, k_CONST, d_s, d_u, theta_BS, theta_MS,
         phase_los, v, theta_v, ts, tn, k_FACTOR, re_h, im_h, output_los_phase);

      } /* end los mode */

}
