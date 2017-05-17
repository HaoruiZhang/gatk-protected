package org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.indexing.INDArrayIndex;

/**
 * A utitlity class for wrapper methods of Nd4j library
 */
public final class Nd4jUtils {

    private Nd4jUtils () {}

    /**
     * Wrapper for the {@link INDArray#get(INDArrayIndex...)}  method to get around the bug described in
     * TODO github/gatk-protected issue #1063
     */
    public static INDArray getNDArrayByIndices(INDArray array, INDArrayIndex indX, INDArrayIndex indY, int sampleNum) {
        if (sampleNum == 1) {
            return array;
        } else {
            return array.get(indX, indY);
        }
    }
}
