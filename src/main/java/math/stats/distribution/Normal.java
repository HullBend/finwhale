/*
Copyright � 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
package math.stats.distribution;

import math.FastMath;
import math.rng.DefaultRng;
import math.rng.PseudoRandom;
import math.stats.ProbabilityFuncs;
import static math.MathConsts.SQRT_TWO_PI;

/**
 * Normal (a.k.a Gaussian) distribution.
 * 
 * <pre>
 *                 1                       2
 *    pdf(x) = -----------   exp( - (x-mean) / 2v ) 
 *             sqrt(2pi*v)
 * 
 *                            x
 *                            -
 *                 1         | |                 2
 *    cdf(x) = -----------   |    exp( - (t-mean) / 2v ) dt
 *             sqrt(2pi*v) | |
 *                          -
 *                         -inf.
 * </pre>
 * 
 * where <tt>v = variance = standardDeviation^2</tt>.
 * <p/>
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 */
public class Normal extends AbstractContinuousDistribution {

    /** Mean of this distribution */
    private final double mean;

    /** Standard deviation of this distribution */
    private final double stdDev;

    /** Variance of this distribution */
    private final double variance;

    /** 1.0 / (stdDev * sqrt(2 * PI)) */
    private final double factor;

    public Normal() {
        this(0.0, 1.0);
    }

    public Normal(final PseudoRandom prng) {
        this(prng, 0.0, 1.0);
    }

    public Normal(final double mean, final double stdDev) {
        this(DefaultRng.newPseudoRandom(), mean, stdDev);
    }

    public Normal(final PseudoRandom prng, final double mean,
            final double stdDev) {
        super(prng);
        if (stdDev <= 0.0) {
            throw new IllegalArgumentException(
                    "Standard deviation must be positive (" + stdDev + ")");
        }
        this.mean = mean;
        this.stdDev = stdDev;
        this.variance = stdDev * stdDev;
        this.factor = (1.0 / (this.variance * SQRT_TWO_PI));
    }

    @Override
    public double pdf(final double x) {
        double xMinusMu = (x - mean);
        return factor * FastMath.exp(-(xMinusMu * xMinusMu) / (2.0 * variance));
    }

    @Override
    public double cdf(final double x) {
        return ProbabilityFuncs.normal(mean, variance, x);
    }

    @Override
    public double sample() {
        return mean + prng.nextGaussian() * stdDev;
    }

    @Override
    public double mean() {
        return mean;
    }

    @Override
    public double variance() {
        return variance;
    }

    @Override
    public String toString() {
        return getSimpleName(mean, stdDev);
    }
}
