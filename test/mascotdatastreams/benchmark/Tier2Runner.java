package mascotdatastreams.benchmark;

import beast.base.parser.XMLParser;
import beast.base.inference.Runnable;
import beast.base.core.BEASTInterface;
import beast.base.util.Randomizer;
import mascotdatastreams.distribution.MascotLogPflag;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Tier-2 timing harness. For each XML:
 *   1. Parse + initialise.
 *   2. Run WARMUP calculateLogP() calls (let JIT settle).
 *   3. Run N timed calls. For each call, time it individually and accumulate
 *      MascotLogPflag.doEulerCallCount.
 *   4. Emit CSV row: xml_path, dataset, mode, parameter, logP,
 *                    n_doEuler_per_call, mean_ms, std_ms, n_calls
 *
 * Configurable via system properties:
 *   -Dtier2.warmup=N     (default 5)
 *   -Dtier2.timed=N      (default 30)
 *
 * Usage:
 *   java mascotdatastreams.benchmark.Tier2Runner OUTPUT_CSV XML [XML ...]
 *
 * The python generator names files like
 *   <dataset>__<mode>__<param>.xml
 * and we extract the three labels from the file name.
 */
public class Tier2Runner {

    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            System.err.println("Usage: Tier2Runner OUTPUT_CSV XML [XML ...]");
            System.exit(2);
        }
        int warmup = Integer.parseInt(System.getProperty("tier2.warmup", "5"));
        int nTimed = Integer.parseInt(System.getProperty("tier2.timed", "30"));
        long seed = Long.parseLong(System.getProperty("tier2.seed", "42"));

        Path outCsv = Paths.get(args[0]);
        Files.createDirectories(outCsv.toAbsolutePath().getParent());

        // Class-loading warmup: parse the first XML once and discard. Lazy
        // class loading on the FIRST parse in a fresh JVM consumes Randomizer
        // values (RandomTree, operator init, etc.), so without this primer
        // the first XML in the sweep gets a different random tree than later
        // XMLs do, even though we Randomizer.setSeed(seed) before every parse.
        // After this primer, all real parses see the same post-class-loading
        // RNG state.
        if (args.length >= 2) {
            System.out.println("[Tier2Runner] class-loading primer parse on " + args[1]);
            Randomizer.setSeed(seed);
            new XMLParser().parseFile(new java.io.File(args[1]));
        }

        try (BufferedWriter w = new BufferedWriter(new FileWriter(outCsv.toFile()))) {
            w.write("xml_path,dataset,baseline,mode,parameter,logP,n_doEuler_per_call,mean_ms,std_ms,n_calls\n");
            for (int i = 1; i < args.length; i++) {
                String xml = args[i];
                System.out.println("[Tier2Runner] " + xml + " (warmup=" + warmup + ", timed=" + nTimed + ", seed=" + seed + ")");
                Result r = runOne(xml, warmup, nTimed, seed);
                w.write(String.format("%s,%s,%s,%s,%s,%.10f,%.2f,%.4f,%.4f,%d%n",
                        xml, r.dataset, r.baseline, r.mode, r.parameter,
                        r.logP, r.nDoEulerPerCall, r.meanMs, r.stdMs, r.nCalls));
                w.flush();
                System.out.printf(
                    "    logP=%.4f  doEuler/call=%.1f  mean=%.3f ms  std=%.3f ms%n",
                    r.logP, r.nDoEulerPerCall, r.meanMs, r.stdMs);
            }
        }
    }

    private static Result runOne(String xmlPath, int warmup, int nTimed, long seed) throws Exception {
        // Reset the RNG seed before EVERY parse so XMLs that use RandomTree (e.g.
        // the SARS XML) get the same tree across all variants in the sweep,
        // making timings/logPs comparable.
        Randomizer.setSeed(seed);
        XMLParser parser = new XMLParser();
        Runnable runnable = parser.parseFile(new java.io.File(xmlPath));
        MascotLogPflag mascot = findMascot(runnable);
        if (mascot == null) {
            throw new IllegalStateException("Could not find a MascotLogPflag inside " + xmlPath);
        }

        // Warm-up
        double logP = 0;
        for (int k = 0; k < warmup; k++) logP = mascot.calculateLogP();

        // Timed
        MascotLogPflag.getAndResetDoEulerCallCount();
        double[] times = new double[nTimed];
        for (int k = 0; k < nTimed; k++) {
            long t0 = System.nanoTime();
            logP = mascot.calculateLogP();
            long t1 = System.nanoTime();
            times[k] = (t1 - t0) / 1.0e6;
        }
        long totalEuler = MascotLogPflag.getAndResetDoEulerCallCount();

        double mean = 0;
        for (double t : times) mean += t;
        mean /= nTimed;
        double var = 0;
        for (double t : times) var += (t - mean) * (t - mean);
        double std = Math.sqrt(var / Math.max(1, nTimed - 1));

        Result r = new Result();
        r.logP = logP;
        r.meanMs = mean;
        r.stdMs = std;
        r.nCalls = nTimed;
        r.nDoEulerPerCall = totalEuler / (double) nTimed;
        r.dataset = parseTagFromName(xmlPath, 0);
        r.baseline = parseTagFromName(xmlPath, 1);
        r.mode = parseTagFromName(xmlPath, 2);
        r.parameter = parseTagFromName(xmlPath, 3);
        return r;
    }

    private static MascotLogPflag findMascot(Object root) {
        Set<Object> visited = new HashSet<>();
        List<Object> stack = new ArrayList<>();
        stack.add(root);
        while (!stack.isEmpty()) {
            Object o = stack.remove(stack.size() - 1);
            if (o == null || !visited.add(o)) continue;
            if (o instanceof MascotLogPflag) return (MascotLogPflag) o;
            if (o instanceof BEASTInterface) {
                BEASTInterface b = (BEASTInterface) o;
                for (BEASTInterface child : b.listActiveBEASTObjects()) stack.add(child);
            }
        }
        return null;
    }

    /**
     * File names: <dataset>__<mode>__<param>.xml
     * which==0 -> dataset, 1 -> mode, 2 -> param.
     */
    private static String parseTagFromName(String path, int which) {
        String name = Paths.get(path).getFileName().toString();
        if (name.endsWith(".xml")) name = name.substring(0, name.length() - 4);
        String[] parts = name.split("__");
        if (which < 0 || which >= parts.length) return "?";
        return parts[which];
    }

    private static class Result {
        double logP;
        double meanMs;
        double stdMs;
        int nCalls;
        double nDoEulerPerCall;
        String dataset;
        String baseline;
        String mode;
        String parameter;
    }
}
