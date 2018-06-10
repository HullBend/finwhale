/*
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * Environment:  Java
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package math.stats.distribution.mle;

import math.FastMath;
import math.GammaFun;
import math.MathConsts;
import math.RootFinder;
import math.function.DoubleUnaryOperator;

/**
 * Provides methods for maximum likelihood estimation of distribution
 * parameters.
 */
public final class MLE {

    private static final double LN_EPS = MathConsts.LN_MIN_NORMAL - MathConsts.LN_2;
    private static final String NO_OBS_MSG = "No observations (x[].length = 0)";

    private static final class GammaMLE implements DoubleUnaryOperator {
        private final int n;
        private final double empiricalMean;
        private final double sumLn;

        GammaMLE(int n, double empiricalMean, double sumLn) {
            this.n = n;
            this.empiricalMean = empiricalMean;
            this.sumLn = sumLn;
        }

        @Override
        public double applyAsDouble(double x) {
            if (x <= 0.0) {
                return 1.0e200;
            }
            return (n * Math.log(empiricalMean / x) + n * GammaFun.digamma(x) - sumLn);
        }
    }

    private static final class WeibullMLE implements DoubleUnaryOperator {
        private final double xi[];
        private final double lnXi[];
        private double sumLnXi = 0.0;

        WeibullMLE(double x[]) {
            xi = x.clone();
            lnXi = new double[x.length];

            for (int i = 0; i < x.length; i++) {
                double lnx;
                if (x[i] > 0.0) {
                    lnx = Math.log(x[i]);
                    lnXi[i] = lnx;
                } else {
                    lnx = LN_EPS;
                    lnXi[i] = lnx;
                }
                sumLnXi += lnx;
            }
        }

        @Override
        public double applyAsDouble(double x) {
            if (x <= 0.0) {
                return 1.0e200;
            }
            double sumXiLnXi = 0.0;
            double sumXi = 0.0;
            for (int i = 0; i < xi.length; i++) {
                double xalpha = FastMath.pow(xi[i], x);
                sumXiLnXi += xalpha * lnXi[i];
                sumXi += xalpha;
            }
            return x * (xi.length * sumXiLnXi - sumLnXi * sumXi) - (xi.length * sumXi);
        }
    }

    /**
     * Estimates the parameters {@code shape} ({@code k}) and {@code scale}
     * (&theta;) of the Gamma distribution from the observations {@code x} using
     * the maximum likelihood method.
     * 
     * @param x
     *            the list of observations to use to evaluate parameters
     * @return returns the parameters {@code k} and &theta;
     */
    public static ParGamma getGammaMLE(double[] x) {
        int n = x.length;
        if (n == 0) {
            throw new IllegalArgumentException(NO_OBS_MSG);
        }

        double sum = 0.0;
        double sumLn = 0.0;

        for (int i = 0; i < n; i++) {
            sum += x[i];
            if (x[i] <= 0.0) {
                sumLn += LN_EPS;
            } else {
                sumLn += Math.log(x[i]);
            }
        }
        double empiricalMean = sum / (double) n;

        sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += (x[i] - empiricalMean) * (x[i] - empiricalMean);
        }

        double alphaMME = (empiricalMean * empiricalMean * (double) n) / sum;
        // left endpoint of initial interval
        double left = alphaMME - 10.0;
        if (left <= 0) {
            left = 1.0e-5;
        }
        // right endpoint of initial interval
        double right = alphaMME + 10.0;

        ParGamma params = new ParGamma();
        params.shape = RootFinder.brentDekker(left, right, new GammaMLE(n, empiricalMean, sumLn), 1e-7);
        params.scale = empiricalMean / params.shape;

        return params;
    }

    /**
     * Estimates the parameters &mu; and &sigma; of the LogNormal distribution
     * from the observations {@code x} using the maximum likelihood method.
     * 
     * @param x
     *            the list of observations to use to evaluate parameters
     * @return returns the parameters &mu; and &sigma;
     */
    public static ParLogNormal getLogNormalMLE(double[] x) {
        int n = x.length;
        if (n == 0) {
            throw new IllegalArgumentException(NO_OBS_MSG);
        }

        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            if (x[i] > 0.0) {
                sum += Math.log(x[i]);
            } else {
                sum += LN_EPS; // log(MIN_NORMAL / 2)
            }
        }

        double mu_hat = sum / n;
        double tmp;
        sum = 0.0;

        for (int i = 0; i < n; i++) {
            if (x[i] > 0.0) {
                tmp = Math.log(x[i]) - mu_hat;
            } else {
                tmp = LN_EPS - mu_hat;
            }
            sum += (tmp * tmp);
        }

        ParLogNormal params = new ParLogNormal();
        params.mu = mu_hat;
        params.sigma = Math.sqrt(sum / n);

        return params;
    }

    /**
     * Estimates the parameters {@code scale} (&lambda;) and {@code shape}
     * ({@code k}) of the Weibull distribution from the observations {@code x}
     * using the maximum likelihood method.
     * 
     * @param x
     *            the list of observations to use to evaluate parameters
     * @return returns the parameters &lambda; and {@code k}
     */
    public static ParWeibull getWeibullMLE(double[] x) {
        int n = x.length;
        if (n == 0) {
            throw new IllegalArgumentException(NO_OBS_MSG);
        }

        double sumLn = 0.0;
        double sumLn2 = 0.0;

        for (int i = 0; i < x.length; i++) {
            double lnxi;
            if (x[i] <= 0.0) {
                lnxi = LN_EPS;
            } else {
                lnxi = Math.log(x[i]);
            }
            sumLn += lnxi;
            sumLn2 += (lnxi * lnxi);
        }

        double alpha0 = Math.sqrt((double) n / ((6.0 / MathConsts.PI_SQUARED) * (sumLn2 - sumLn * sumLn / (double) n)));
        // left endpoint of initial interval
        double left = alpha0 - 20.0;
        if (left <= 0.0) {
            left = 1.0e-5;
        }
        // right endpoint of initial interval
        double right = alpha0 + 20.0;

        double k = RootFinder.brentDekker(left, right, new WeibullMLE(x), 1e-5);

        double sumXalpha = 0.0;
        for (int i = 0; i < x.length; i++) {
            sumXalpha += FastMath.pow(x[i], k);
        }
        double scale = 1.0 / (FastMath.pow((double) n / sumXalpha, 1.0 / k));

        ParWeibull params = new ParWeibull();
        params.shape = k;
        params.scale = scale;
        return params;
    }

    private MLE() {
        throw new AssertionError();
    }
}
