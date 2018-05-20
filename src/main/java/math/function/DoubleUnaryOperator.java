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

/**
 * Represents an operation on a single {@code double}-valued operand that produces
 * a {@code double}-valued result.
 */
public interface DoubleUnaryOperator {
    /**
     * Applies this operator to the given operand.
     *
     * @param x the operand
     * @return the operator result
     */
    double applyAsDouble(double x);
}
