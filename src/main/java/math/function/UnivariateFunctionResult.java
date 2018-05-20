/*
 * Copyright 2013 SPZ
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
package math.function;

public final class UnivariateFunctionResult {

    public final double point;
    public final double value;
    public final int iterations;
    public final boolean converged;

    public UnivariateFunctionResult(double point, double value, int iterations,
            boolean converged) {
        this.point = point;
        this.value = value;
        this.iterations = iterations;
        this.converged = converged;
    }

    public UnivariateFunctionResult(double point, double value, int iterations) {
        this(point, value, iterations, true);
    }

    public UnivariateFunctionResult(double point, double value) {
        this(point, value, 0);
    }
}
