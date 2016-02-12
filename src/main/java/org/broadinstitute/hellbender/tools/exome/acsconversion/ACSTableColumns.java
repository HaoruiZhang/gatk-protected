package org.broadinstitute.hellbender.tools.exome.acsconversion;

import java.util.stream.Stream;

public enum ACSTableColumns {
    CHROMOSOME("Chromosome"), START("Start.bp"), END("End.bp"), NUM_PROBES("n_probes"), LENGTH("length"), NUM_HETS("n_hets"), F("f"), TAU("tau"), SIGMA_TAU("sigma.tau"), MU_MINOR("mu.minor"), SIGMA_MINOR("sigma.minor"), MU_MAJOR("mu.major"), SIGMA_MAJOR("sigma.major"), SEGLABELCNLOH("SegLabelCNLOH");

    private String columnName;  //store the column names

    ACSTableColumns(String columnName) {this.columnName = columnName; }

    @Override
    public String toString() {
        return columnName;
    }

    public static final String[] COLUMN_NAME_ARRAY =
            Stream.of(values()).map(ACSTableColumns::toString).toArray(String[]::new);

}