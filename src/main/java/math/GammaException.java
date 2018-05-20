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
package math;

/**
 * Runtime exception thrown from {@link Gamma} methods.
 */
public final class GammaException extends ArithmeticException {

    private static final long serialVersionUID = 6738428636131168196L;

    public GammaException(final String msg) {
        super(msg);
    }

    public GammaException(final String msg, final int location) {
        super(msg + " : loc:" + location);
    }

    public GammaException(final String msg, final String reason, final int location) {
        super(msg + " : reason = \"" + reason + "\", loc:" + location);
    }

    public GammaException(final String msg, final double value, final String reason, final int location) {
        super(msg + " : [x = " + Double.toString(value) + " ], reason = \"" + reason + "\", loc:" + location);
    }
}
