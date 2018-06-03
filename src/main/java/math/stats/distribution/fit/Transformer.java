/*
 * Class:        GofStat
 * Description:  Goodness-of-fit test statistics
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

import math.stats.distribution.ContinuousDistribution;

final class Transformer {

    static double[] uniform(double[] observations, ContinuousDistribution dist) {
        if (dist == null) {
            throw new IllegalArgumentException("dist == null");
        }
        if (observations == null) {
            throw new IllegalArgumentException("observations == null");
        }
        double[] transformed = new double[observations.length];
        for (int i = 0; i < transformed.length; ++i) {
            transformed[i] = dist.cdf(observations[i]);
        }
        return transformed;
    }

    private Transformer() {
        throw new AssertionError();
    }
}
