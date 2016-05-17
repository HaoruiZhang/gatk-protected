package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.ojalgo.matrix.decomposition.SingularValue;

class OjAlgoSingularValueDecomposer {

    private static final Logger logger = LogManager.getLogger(OjAlgoSingularValueDecomposer.class);

    private OjAlgoSingularValueDecomposer() {
    }

    /**
     * Create a SVD instance using ojAlgo.
     *
     * @param m matrix that is not {@code null}
     * @return SVD instance that is never {@code null}
     */
    public static SVD createSVD(final RealMatrix m) {

        Utils.nonNull(m, "Cannot create SVD on a null matrix.");

        logger.info("Converting Apache RealMatrix to ojAlgo Matrix...");
        //        final PhysicalStore.Factory<Double, PrimitiveDenseStore> tmpFactory = PrimitiveDenseStore.FACTORY;
        //        final PrimitiveDenseStore pds = tmpFactory.rows(m.getData());
        final OjAlgoAdapter pds = OjAlgoAdapter.of(m);

        logger.info("Calculating SVD...");
        final SingularValue<Double> svd = SingularValue.make(pds);
        svd.compute(pds);

        // If this code is too slow, conversion could be implemented by backing a RealMatrix with an ojAlgo Array2D
        logger.info("Converting ojAlgo Matrix to Apache RealMatrix...");
        final RealMatrix pinv = new Array2DRowRealMatrix(svd.getInverse().toRawCopy2D());
        final RealMatrix u = new Array2DRowRealMatrix(svd.getQ1().toRawCopy2D());
        final RealMatrix v = new Array2DRowRealMatrix(svd.getQ2().toRawCopy2D());
        final double[] s = svd.getSingularValues().toRawCopy1D();

        return new SimpleSVD(u, s, v, pinv);
    }
}
