import re
classes = ["SimpleEdge", "Path", "PathWithOrig", "PathOverlap", "SimplePathNodeEdge",
           "Path_n_MM_count", "AlignmentStats", "My_DFS"]
open("src/src/StableHash.java", "w").write('''import java.util.concurrent.atomic.AtomicInteger;

/**
 * Hash codes in creation order, for classes that compare by identity but are
 * used as HashMap/HashSet keys. With the JVM's identity hash, the iteration
 * order of those maps (and so the assembly) shifted whenever unrelated code,
 * such as building a debug message, consumed an identity hash first.
 */
public class StableHash {
	private static final AtomicInteger next = new AtomicInteger();

	public static int next() {
		return next.getAndIncrement();
	}
}
''')
for c in classes:
    p = "src/src/%s.java" % c
    s = open(p).read()
    m = re.search(r"public\s+class\s+%s\b[^{]*\{" % c, s)
    assert m, c
    ins = ("\n\tprivate final int _stable_hash = StableHash.next();\n\n"
           "\tpublic int hashCode() {\n\t\treturn _stable_hash;\n\t}\n")
    s = s[:m.end()] + ins + s[m.end():]
    open(p, "w").write(s)
    print("patched", c)
