package mascotdatastreams.benchmark;

import beast.base.parser.XMLParser;
import beast.base.inference.Runnable;
import beast.base.core.BEASTInterface;
import beast.base.util.Randomizer;
import mascotdatastreams.distribution.MascotLogPflag;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Determinism check: parse the same XML once, then call calculateLogP() N times
 * in a row and report whether logP and the doEuler count are stable.
 *
 * Usage: java mascotdatastreams.benchmark.RepeatabilityCheck XML [N=10]
 */
public class RepeatabilityCheck {
    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: RepeatabilityCheck XML [N=10]");
            System.exit(2);
        }
        int n = args.length >= 2 ? Integer.parseInt(args[1]) : 10;

        Randomizer.setSeed(42L);
        XMLParser parser = new XMLParser();
        Runnable runnable = parser.parseFile(new java.io.File(args[0]));
        MascotLogPflag mascot = findMascot(runnable);
        if (mascot == null) throw new IllegalStateException("no MascotLogPflag");

        for (int k = 0; k < n; k++) {
            MascotLogPflag.getAndResetDoEulerCallCount();
            long t0 = System.nanoTime();
            double logP = mascot.calculateLogP();
            long t1 = System.nanoTime();
            long c = MascotLogPflag.getAndResetDoEulerCallCount();
            System.out.printf("call %2d: logP=%.10f  doEuler=%d  ms=%.3f%n",
                    k, logP, c, (t1 - t0) / 1.0e6);
        }
    }

    private static MascotLogPflag findMascot(Object root) {
        Set<Object> visited = new HashSet<>();
        List<Object> stack = new ArrayList<>();
        stack.add(root);
        while (!stack.isEmpty()) {
            Object o = stack.remove(stack.size() - 1);
            if (o == null || !visited.add(o)) continue;
            if (o instanceof MascotLogPflag) return (MascotLogPflag) o;
            if (o instanceof BEASTInterface)
                for (BEASTInterface c : ((BEASTInterface) o).listActiveBEASTObjects()) stack.add(c);
        }
        return null;
    }
}
