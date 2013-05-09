package edu.unc.mapseq.pipeline.ncgenes.dx;

import edu.unc.mapseq.pipeline.AbstractPipelineBeanService;

public class NCGenesDXPipelineBeanService extends AbstractPipelineBeanService {

    private String referenceSequence;

    private String icSNPIntervalList;

    private String siteName;

    public NCGenesDXPipelineBeanService() {
        super();
    }

    public String getSiteName() {
        return siteName;
    }

    public void setSiteName(String siteName) {
        this.siteName = siteName;
    }

    public String getReferenceSequence() {
        return referenceSequence;
    }

    public void setReferenceSequence(String referenceSequence) {
        this.referenceSequence = referenceSequence;
    }

    public String getIcSNPIntervalList() {
        return icSNPIntervalList;
    }

    public void setIcSNPIntervalList(String icSNPIntervalList) {
        this.icSNPIntervalList = icSNPIntervalList;
    }

}
