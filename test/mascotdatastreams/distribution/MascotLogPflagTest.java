package mascotdatastreams.distribution;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

import org.junit.jupiter.api.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import mascot.distribution.StructuredTreeIntervals;
import mascot.dynamics.Constant;

/**
 * Tests that MascotLogPflag returns the actual Mascot log P when compute_likelihood is true
 * and returns 0 (after initialization) when compute_likelihood is false.
 */
public class MascotLogPflagTest {

    @Test
    public void computeLikelihoodTrueReturnsLogP() throws Exception {
        Tree tree = buildMinimalTree();
        StructuredTreeIntervals intervals = buildStructuredTreeIntervals(tree);
        Constant dynamics = buildConstantDynamics();

        MascotLogPflag mascot = new MascotLogPflag();
        mascot.initByName(
                "tree", tree,
                "dynamics", dynamics,
                "structuredTreeIntervals", intervals,
                "compute_likelihood", true,
                "implementation", MascotLogPflag.MascotImplementation.java
        );
        mascot.initAndValidate();

        double logP = mascot.calculateLogP();
        

        assertNotEquals(0.0, logP, 1e-15, "With compute_likelihood=true, calculateLogP() should return the actual Mascot log P, not 0");
    }

    @Test
    public void computeLikelihoodFalseReturnsZero() throws Exception {
        Tree tree = buildMinimalTree();
        StructuredTreeIntervals intervals = buildStructuredTreeIntervals(tree);
        Constant dynamics = buildConstantDynamics();

        MascotLogPflag mascot = new MascotLogPflag();
        mascot.initByName(
                "tree", tree,
                "dynamics", dynamics,
                "structuredTreeIntervals", intervals,
                "compute_likelihood", false,
                "implementation", MascotLogPflag.MascotImplementation.java
        );
        mascot.initAndValidate();

        double logP = mascot.calculateLogP();

        assertEquals(0.0, logP, 0.0, "With compute_likelihood=false, calculateLogP() should return 0");
    }

    @Test
    public void defaultComputeLikelihoodIsTrue() throws Exception {
        Tree tree = buildMinimalTree();
        StructuredTreeIntervals intervals = buildStructuredTreeIntervals(tree);
        Constant dynamics = buildConstantDynamics();

        MascotLogPflag mascot = new MascotLogPflag();
        mascot.initByName(
                "tree", tree,
                "dynamics", dynamics,
                "structuredTreeIntervals", intervals,
                "implementation", MascotLogPflag.MascotImplementation.java
        );
        mascot.initAndValidate();

        double logP = mascot.calculateLogP();

        assertNotEquals(0.0, logP, 1e-15, "Default compute_likelihood is true; log P should not be 0");
    }

    private static Tree buildMinimalTree() throws Exception {
        Sequence s1 = new Sequence();
        Sequence s2 = new Sequence();
        s1.initByName("taxon", "A", "value", "?");
        s2.initByName("taxon", "B", "value", "?");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", s1, "sequence", s2);
        TaxonSet taxa = new TaxonSet();
        taxa.initByName("alignment", alignment);
        TraitSet traitSet = new TraitSet();
        traitSet.initByName("value", "A=0,B=0", "traitname", "type", "taxa", taxa);
        Tree tree = new TreeParser("(A:0.5,B:0.5):0.0;");
        tree.initByName("taxonset", taxa, "trait", traitSet);
        return tree;
    }

    private static StructuredTreeIntervals buildStructuredTreeIntervals(Tree tree) throws Exception {
        StructuredTreeIntervals intervals = new StructuredTreeIntervals();
        intervals.initByName("tree", tree);
        intervals.initAndValidate();
        return intervals;
    }

    /**
     * Build Constant dynamics for a single deme (dimension 1), matching Mascot package setup:
     * Ne plus backwardsMigration (required; dimension 0 when dimension=1).
     * See deps/Mascot/examples/Constant.xml and fxtemplates/MascotConstant.xml.
     */
    private static Constant buildConstantDynamics() throws Exception {
        RealParameter ne = new RealParameter(new Double[] { 1.0 });
        RealParameter backwardsMigration = new RealParameter(new Double[] {}); // dimension 0 for single deme
        Constant constant = new Constant();
        constant.initByName("Ne", ne, "backwardsMigration", backwardsMigration, "dimension", 1);
        constant.initAndValidate();
        return constant;
    }
}
