package math.rng;

public final class Seed {

    // from java.util.Random
    private static long seedUniquifier = 0x1ed8b55fac9decL;

    private static long lastSeed = pseudoRandomSeed();

    // from java.util.Random
    private static long nextSeedUniquifier() {
        // Pierre L'Ecuyer: "Tables of Linear Congruential Generators
        // of Different Sizes and Good Lattice Structure"
        seedUniquifier *= 0x106689d45497fdb5L;
        return seedUniquifier;
    }

    /*
     * Returns a reasonably good (pseudo) random seed
     */
    private static long pseudoRandomSeed() {
        long seed = nextSeedUniquifier() ^ System.nanoTime();

        // apply Austin Appleby's fmix64() hash
        seed ^= seed >>> 33;
        seed *= 0xff51afd7ed558ccdL;
        seed ^= seed >>> 33;
        seed *= 0xc4ceb9fe1a85ec53L;
        seed ^= seed >>> 33;

        return seed;
    }

    /**
     * Returns a reasonably good long random seed.
     * 
     * @return a long random seed.
     */
    public static synchronized long seed() {
        long seed = pseudoRandomSeed();
        while (seed == lastSeed || seed == 0L) {
            seed = pseudoRandomSeed();
        }
        lastSeed = seed;
        return seed;
    }

    private Seed() {
    }
}
