/*
 * Class:        GofFormat
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
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
package math.stats.distribution.fit;

import java.util.Arrays;

import math.stats.distribution.ContinuousDistribution;

/**
 * TODO
 */
public final class GoodnessOfFit {

    /**
     * TODO
     * 
     * @param observations
     * @param distribution
     * @return the {@link UniformTestStatistics.Result} for the given
     *         observations and distribution
     */
    public static UniformTestStatistics.Result computeStatistics(double[] observations,
            ContinuousDistribution distribution) {
        double[] transformed = Transformer.uniform(observations, distribution);
        Arrays.sort(transformed);
        UniformTestStatistics.Result statistics = UniformTestStatistics.compareEmpiricalToUniform(transformed);
        if (observations.length == 1) {
            statistics.MEAN = observations[0];
        }
        return statistics;
    }

    private GoodnessOfFit() {
        throw new AssertionError();
    }
}
