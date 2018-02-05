/*
Copyright © 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
package math.stats.distribution;

import math.FastGamma;
import math.FastMath;
import math.rng.DefaultRng;
import math.rng.PseudoRandom;
import math.stats.ProbabilityFuncs;

/**
 * StudentT distribution (a.k.a "T distribution").
 * <p>
 * <tt>p(x) = const  *  (1 + x^2/&nu;) ^ -(&nu;+1)/2</tt> where
 * <tt>const = &Gamma;((&nu;+1)/2) / (&radic;(&Pi;*&nu;) * &Gamma;(&nu;/2))</tt>
 * and <tt>&Gamma;(a)</tt> being the Gamma function and <tt>&nu;</tt> being the
 * degrees of freedom.
 * </p>
 * <b>Implementation:</b> This is a port of <A HREF=
 * "http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/Random/RandStudentT.html"
 * >RandStudentT</A> used in <A
 * HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep">CLHEP 1.4.0</A> (C++). CLHEP's
 * implementation, in turn, is based on <tt>tpol.c</tt> from the <A
 * HREF="http://www.cis.tu-graz.ac.at/stat/stadl/random.html">C-RAND /
 * WIN-RAND</A> library. C-RAND's implementation, in turn, is based upon
 * <p>
 * R.W. Bailey (1994): Polar generation of random variates with the
 * t-distribution, Mathematics of Computation 62, 779-781.
 * 
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 */
public class StudentT extends AbstractContinuousDistribution {

    private final int df;
    private final double pdfConst;

    public StudentT(final int df) {
        this(DefaultRng.newPseudoRandom(), df);
    }

    public StudentT(final PseudoRandom prng, final int df) {
        super(prng);
        if (df < 1) {
            throw new IllegalArgumentException("df < 1 : " + df);
        }
        final double tmp = FastGamma.logGamma((df + 1.0) / 2.0)
                - FastGamma.logGamma(df / 2.0);
        this.pdfConst = FastMath.exp(tmp) / Math.sqrt(Math.PI * df);
        this.df = df;
    }

    @Override
    public double pdf(final double x) {
        return pdfConst * FastMath.pow((1.0 + x * x / df), -(df + 1.0) * 0.5);
    }

    @Override
    public double cdf(final double x) {
        return ProbabilityFuncs.studentT(df, x);
    }

    @Override
    public double sample() {
        /*
         * Marsaglia's formulation of the Box/Muller polar method for generating
         * Normal variates is adapted to the Student-t distribution. The two
         * generated variates are not independent and the expected number of
         * uniforms per variate is 2.5464.
         * 
         * Reference:
         * 
         * R.W. Bailey (1994): Polar generation of random variates with the
         * t-distribution, Mathematics of Computation 62, 779-781.
         */
        double u1, u2, q;
        do {
            u1 = 2.0 * prng.nextDouble() - 1.0; // between -1 and 1
            u2 = 2.0 * prng.nextDouble() - 1.0; // between -1 and 1
            q = u1 * u1 + u2 * u2;
        } while (q > 1.0);
        return u1
                * Math.sqrt(df * (FastMath.exp(-2.0 / df * Math.log(q)) - 1.0)
                        / q);
    }

    @Override
    public double mean() {
        if (df == 1) {
            return Double.NaN;
        }
        return 0.0;
    }

    @Override
    public double variance() {
        if (df > 2) {
            return df / ((double) df - 2.0);
        }
        if (df == 2) {
            return Double.POSITIVE_INFINITY;
        }
        return Double.NaN;
    }

    public int getDegreesOfFreedom() {
        return df;
    }
}
