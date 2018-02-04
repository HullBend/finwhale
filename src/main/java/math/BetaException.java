package math;

/**
 * Runtime exception thrown from {@link Beta} methods.
 */
public final class BetaException extends ArithmeticException {

    private static final long serialVersionUID = -5428195944932663589L;

    public BetaException(final String msg) {
        super(msg);
    }

    public BetaException(final String msg, final int location) {
        super(msg + " : loc:" + location);
    }

    public BetaException(final String msg, final String reason, final int location) {
        super(msg + " : reason = \"" + reason + "\", loc:" + location);
    }

    public BetaException(final String msg, final double value, final String reason, final int location) {
        super(msg + " : [x = " + Double.toString(value) + " ], reason = \"" + reason + "\", loc:" + location);
    }
}
