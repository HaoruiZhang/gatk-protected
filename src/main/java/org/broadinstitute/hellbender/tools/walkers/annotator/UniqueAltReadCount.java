package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;

import java.util.*;
import java.util.stream.Collectors;

/**
 * This annotation counts up the number of ALT reads that have a unique start position and fragment length.
 * We have seen that the majority of false positives in low allele fraction, cell-free DNA have the profile where
 * the evidence for alternate allele comes solely from one, two, or three sets of PCR-duplicates.
 * Reads in such a set have the same read-start and mate-end position (i.e. they come from the same insert),
 * but they have different unique molecular identifiers (UMIs) so UMI-aware Mark Duplicates considers them unique pieces of evidence.
 * Although these reads have different UMIs, we suspect that they really are PCR-duplicates, for two reasons:
 * 1) these sites are false positives and 2) with the library we used, which was not TSCA or NEB,
 * it's highly unlikely that we get multiple fragments with identical start and end positions.
 * The true cause, we suspect, is a barcode swap, which is the process in which two fragments with the same barcode swap UMIs.
 *
 * If evidence for an alternate allele came from a few unique reads, we will filter the variant.
 * Mutect2FilteringEngine::applyDuplicatedAltReadFilter is the accompanying filter.
 */
public class UniqueAltReadCount extends GenotypeAnnotation implements StandardSomaticAnnotation {
    public static final String UNIQUE_ALT_READ_SET_COUNT_KEY = "UNIQ_ALT_READ_COUNT";

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(UNIQUE_ALT_READ_SET_COUNT_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList( new VCFFormatHeaderLine(UNIQUE_ALT_READ_SET_COUNT_KEY, 1, VCFHeaderLineType.Integer,
                "Number of ALT reads with unique start and mate end positions at a variant site"));
    }

    @Override
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        if (g.isHomRef()){
            // skip the normal sample
            return;
        }

        final Allele altAllele = vc.getAlternateAllele(0); // assume single-allelic
        final String tumorSampleName = g.getSampleName();
        Collection<ReadLikelihoods<Allele>.BestAllele> tumorBestAlleles = likelihoods.bestAlleles(tumorSampleName);
        Map<ImmutablePair<Integer, Integer>, Long> duplicateReadMap = tumorBestAlleles.stream()
                .filter(ba -> ba.allele.equals(altAllele) && ba.isInformative())
                .map(ba -> new ImmutablePair<>(ba.read.getStart(), ba.read.getFragmentLength()))
                .collect(Collectors.groupingBy(x -> x, Collectors.counting()));

        gb.attribute(UNIQUE_ALT_READ_SET_COUNT_KEY, duplicateReadMap.size());
    }


}
