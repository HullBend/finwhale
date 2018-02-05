package math.stats.distribution;

import math.FastMath;
import math.rng.DefaultRng;
import math.rng.PseudoRandom;

public class Exponential extends AbstractContinuousDistribution {

    private static final double BIG = 100.0;

    private final double lambda;

    public Exponential(final double lambda) {
        this(DefaultRng.newPseudoRandom(), lambda);
    }

    public Exponential(final PseudoRandom prng, final double lambda) {
        super(prng);
        if (lambda <= 0.0) {
            throw new IllegalArgumentException("lambda <= 0.0 : " + lambda);
        }
        this.lambda = lambda;
    }

    @Override
    public double pdf(double x) {
        return x < 0.0 ? 0.0 : lambda * FastMath.exp(-lambda * x);
    }

    @Override
    public double cdf(double x) {
        if (x <= 0.0) {
            return 0.0;
        }
        double y = lambda * x;
        if (y >= BIG) {
            return 1.0;
        }
        return -FastMath.expm1(-y);
    }

    @Override
    public double sample() {
        double u;
        do {
            u = prng.nextDouble();
        } while (u == 0.0 || u == 1.0);
        return -Math.log(u) / lambda;
    }

    @Override
    public double mean() {
        return (1.0 / lambda);
    }

    @Override
    public double variance() {
        return (1.0 / (lambda * lambda));
    }
}
