package mascotdatastreams.benchmark;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.IntervalType;
import beast.base.parser.XMLParser;
import beast.base.util.Randomizer;
import mascot.distribution.StructuredTreeIntervals;
import mascot.dynamics.Dynamics;
import mascotdatastreams.distribution.MascotLogPflag;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Dump the sequence of (length, t_start, t_end, event_at_start) intervals that
 * MascotLogPflag iterates through for a given XML and mode, to a CSV file.
 *
 * "Old" mode = the original calculateLogP loop: outer iterations land on
 * min(nextTreeEvent, nextRateShift). Each row is one ODE step.
 *
 * "New" mode = calculateLogP_maxInterval: each tree-event interval is split
 * into ceil(L / maxInterval) equal sub-intervals.
 *
 * event_at_start labels what kind of event is at t_start of the row:
 *   start      — t_start = 0 (most recent tip side)
 *   sampling   — a sample event happened at t_start (a leaf was added)
 *   coalescent — a coalescent event happened at t_start
 *   interval   — t_start is a grid breakpoint (old) or max-interval breakpoint
 *                (new), i.e. not a tree event
 *   nothing    — IntervalType.NOTHING, rare
 *
 * Usage:
 *   java mascotdatastreams.benchmark.IntervalDumper OUTPUT_CSV XML MODE [maxInterval]
 *   MODE in {old, new}; maxInterval required when MODE=new.
 *
 * Class-loading primer: parses the XML once and discards before the real
 * parse, so the seed-pinned RandomTree is identical across runs (same fix as
 * Tier2Runner).
 */
public class IntervalDumper {
    public static void main(String[] args) throws Exception {
        if (args.length < 3) {
            System.err.println("Usage: IntervalDumper OUTPUT_CSV XML MODE [maxInterval]");
            System.exit(2);
        }
        String outCsv = args[0];
        String xmlPath = args[1];
        String mode = args[2];
        double maxInterval = (args.length > 3) ? Double.parseDouble(args[3]) : Double.POSITIVE_INFINITY;
        long seed = Long.parseLong(System.getProperty("dumper.seed", "42"));

        // Class-loading primer parse (discarded).
        Randomizer.setSeed(seed);
        new XMLParser().parseFile(new File(xmlPath));

        // Real parse.
        Randomizer.setSeed(seed);
        XMLParser parser = new XMLParser();
        beast.base.inference.Runnable runnable = parser.parseFile(new File(xmlPath));
        MascotLogPflag mascot = findMascot(runnable);
        if (mascot == null) throw new IllegalStateException("no MascotLogPflag in " + xmlPath);

        StructuredTreeIntervals ti = mascot.structuredTreeIntervalsInput.get();
        Dynamics dyn = mascot.dynamicsInput.get();
        ti.calculateIntervals();

        // Two-pass dump: collect rows with the t_start label, then fill in
        // each row's t_end label from the *next* row's t_start label.
        List<Row> rows = new ArrayList<>();
        if ("old".equals(mode)) {
            dumpOld(rows, ti, dyn);
        } else if ("new".equals(mode)) {
            if (!Double.isFinite(maxInterval) || maxInterval <= 0) {
                throw new IllegalArgumentException("MODE=new requires a positive maxInterval");
            }
            dumpNew(rows, ti, maxInterval);
        } else {
            throw new IllegalArgumentException("MODE must be 'old' or 'new', got: " + mode);
        }

        try (BufferedWriter w = new BufferedWriter(new FileWriter(outCsv))) {
            w.write("length,t_start,t_end,t_start_type,t_end_type\n");
            for (int i = 0; i < rows.size(); i++) {
                Row r = rows.get(i);
                String tEndType = (i + 1 < rows.size()) ? rows.get(i + 1).startType : "end";
                w.write(String.format("%.12f,%.12f,%.12f,%s,%s%n",
                        r.length, r.tStart, r.tEnd, r.startType, tEndType));
            }
        }
        System.out.println("[IntervalDumper] wrote " + outCsv);
    }

    private static final class Row {
        final double length;
        final double tStart;
        final double tEnd;
        final String startType;
        Row(double length, double tStart, double tEnd, String startType) {
            this.length = length; this.tStart = tStart; this.tEnd = tEnd;
            this.startType = startType;
        }
    }

    /** Mirror calculateLogP's outer loop: min(nextTreeEvent, nextRateShift). */
    private static void dumpOld(List<Row> rows, StructuredTreeIntervals ti, Dynamics dyn) throws Exception {
        int treeInterval = 0, ratesInterval = 0;
        double currentTime = 0.0;
        double nextTreeEvent;
        try { nextTreeEvent = ti.getInterval(treeInterval); }
        catch (Exception e) { return; }
        double nextRateShift = dyn.getInterval(ratesInterval);
        String eventAtStart = "start";

        while (true) {
            double dt = Math.min(nextTreeEvent, nextRateShift);
            double tStart = currentTime;
            double tEnd = currentTime + dt;
            rows.add(new Row(dt, tStart, tEnd, eventAtStart));

            if (nextTreeEvent <= nextRateShift) {
                IntervalType type = ti.getIntervalType(treeInterval);
                eventAtStart = labelOf(type);
                treeInterval++;
                nextRateShift -= nextTreeEvent;
                try { nextTreeEvent = ti.getInterval(treeInterval); }
                catch (Exception e) { break; }
            } else {
                eventAtStart = "interval";
                ratesInterval++;
                nextTreeEvent -= nextRateShift;
                nextRateShift = dyn.getInterval(ratesInterval);
            }
            currentTime = tEnd;
        }
    }

    /** Mirror calculateLogP_maxInterval's outer loop. */
    private static void dumpNew(List<Row> rows, StructuredTreeIntervals ti, double maxInterval) throws Exception {
        int treeInterval = 0;
        double currentTime = 0.0;
        double nextTreeEvent;
        try { nextTreeEvent = ti.getInterval(treeInterval); }
        catch (Exception e) { return; }
        String eventAtStart = "start";

        while (true) {
            double L = nextTreeEvent;

            if (L > 0) {
                int nSub = (Double.isFinite(maxInterval) && L > maxInterval)
                        ? (int) Math.ceil(L / maxInterval) : 1;
                double subL = L / nSub;
                for (int k = 0; k < nSub; k++) {
                    double tStart = currentTime + k * subL;
                    double tEnd = tStart + subL;
                    String thisLabel = (k == 0) ? eventAtStart : "interval";
                    rows.add(new Row(subL, tStart, tEnd, thisLabel));
                }
                currentTime += L;
            } else {
                rows.add(new Row(0.0, currentTime, currentTime, eventAtStart));
            }

            IntervalType type = ti.getIntervalType(treeInterval);
            eventAtStart = labelOf(type);

            treeInterval++;
            try { nextTreeEvent = ti.getInterval(treeInterval); }
            catch (Exception e) { break; }
        }
    }

    private static String labelOf(IntervalType type) {
        if (type == IntervalType.COALESCENT) return "coalescent";
        if (type == IntervalType.SAMPLE) return "sampling";
        return "nothing";
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
                for (BEASTInterface c : ((BEASTInterface) o).listActiveBEASTObjects()) stack.add(c);
            }
        }
        return null;
    }
}
