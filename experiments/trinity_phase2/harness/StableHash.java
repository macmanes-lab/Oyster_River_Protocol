import java.util.concurrent.atomic.AtomicInteger;

/**
 * Reproducible hash codes for classes that compare by identity but are used
 * as HashMap/HashSet keys. With the JVM's identity hash, the iteration order
 * of those maps (and so the assembly) shifted whenever unrelated code, such as
 * building a debug message, consumed an identity hash first.
 *
 * The creation counter is scrambled (MurmurHash3 fmix32) so map order stays
 * random-like, as the identity hash order was, instead of following creation
 * order: Butterfly's output is sensitive to that order, and creation order
 * shifted assembly size by several percent in either direction.
 */
public class StableHash {
	private static final AtomicInteger next = new AtomicInteger();

	public static int next() {
		int h = next.getAndIncrement();
		h ^= h >>> 16;
		h *= 0x85ebca6b;
		h ^= h >>> 13;
		h *= 0xc2b2ae35;
		h ^= h >>> 16;
		return h;
	}
}
