package mascotdatastreams.benchmark;

import beast.base.parser.XMLParser;
import beast.base.inference.Runnable;
import beast.base.core.BEASTInterface;
import beast.base.util.Randomizer;
import mascotdatastreams.distribution.MascotLogPflag;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Tier-1 correctness harness: load each XML, find the MascotLogPflag distribution
 * (by id="Mascot.t:SimDataset" or any class match), call calculateLogP() once,
 * write a CSV row per XML.
 *
 * Usage:
 *   java mascotdatastreams.benchmark.Tier1Runner OUTPUT_CSV XML [XML ...]
 *
 * Output CSV columns:
 *   xml_path,mode,parameter,logP,wall_ms
 *
 * mode/parameter are derived from the XML file name (the python generator
 * encodes them).
 */
public class Tier1Runner {

    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            System.err.println("Usage: Tier1Runner OUTPUT_CSV XML [XML ...]");
            System.exit(2);
        }
        Path outCsv = Paths.get(args[0]);
        Files.createDirectories(outCsv.toAbsolutePath().getParent());

        try (BufferedWriter w = new BufferedWriter(new FileWriter(outCsv.toFile()))) {
            w.write("xml_path,mode,parameter,logP,wall_ms\n");
            for (int i = 1; i < args.length; i++) {
                String xml = args[i];
                System.out.println("[Tier1Runner] " + xml);
                Result r = runOne(xml);
                w.write(String.format("%s,%s,%s,%.10f,%.3f%n",
                        xml, r.mode, r.parameter, r.logP, r.wallMs));
                w.flush();
            }
        }
    }

    private static Result runOne(String xmlPath) throws Exception {
        // Pin the RNG seed so XMLs that use RandomTree initialisation are
        // reproducible across the sweep and across separate JVM invocations.
        Randomizer.setSeed(42L);

        // Silence BEAST's chatty initialisation while still letting fatal errors through.
        PrintStream origOut = System.out;
        PrintStream origErr = System.err;

        XMLParser parser = new XMLParser();
        Runnable runnable = parser.parseFile(new java.io.File(xmlPath));
        // The parser also stores all parsed objects on its registry; walk that.
        MascotLogPflag mascot = findMascot(runnable);
        if (mascot == null) {
            throw new IllegalStateException("Could not find a MascotLogPflag inside " + xmlPath);
        }

        // Warm the JIT a touch so the recorded time isn't first-call cost.
        for (int k = 0; k < 3; k++) mascot.calculateLogP();

        long t0 = System.nanoTime();
        double logP = mascot.calculateLogP();
        long t1 = System.nanoTime();

        Result r = new Result();
        r.logP = logP;
        r.wallMs = (t1 - t0) / 1.0e6;
        r.mode = parseTagFromName(xmlPath, "mode");
        r.parameter = parseTagFromName(xmlPath, "param");
        return r;
    }

    /** Walk a BEASTObject graph (via its outputs/inputs) to find a MascotLogPflag. */
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
                for (BEASTInterface child : b.listActiveBEASTObjects()) {
                    stack.add(child);
                }
            }
        }
        return null;
    }

    /**
     * The python generator names files like 6_2_<mode>__<param>.xml where
     * mode ∈ {old, new}, param is e.g. "grid1x", "grid10x", "max1e-3".
     * We extract these so the CSV is self-describing.
     */
    private static String parseTagFromName(String path, String which) {
        String name = Paths.get(path).getFileName().toString();
        if (name.endsWith(".xml")) name = name.substring(0, name.length() - 4);
        String[] parts = name.split("__");
        // Last two segments: mode, param (in that order)
        if (parts.length < 2) return "?";
        if ("mode".equals(which)) return parts[parts.length - 2];
        if ("param".equals(which)) return parts[parts.length - 1];
        return "?";
    }

    private static class Result {
        double logP;
        double wallMs;
        String mode;
        String parameter;
    }
}
