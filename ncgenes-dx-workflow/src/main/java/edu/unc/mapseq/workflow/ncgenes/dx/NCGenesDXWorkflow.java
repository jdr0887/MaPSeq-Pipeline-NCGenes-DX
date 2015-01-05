package edu.unc.mapseq.workflow.ncgenes.dx;

import java.io.File;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.ResourceBundle;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.lang.StringUtils;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.renci.jlrm.condor.CondorJob;
import org.renci.jlrm.condor.CondorJobBuilder;
import org.renci.jlrm.condor.CondorJobEdge;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.commons.ncgenes.dx.RegisterToIRODSRunnable;
import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.Attribute;
import edu.unc.mapseq.dao.model.FileData;
import edu.unc.mapseq.dao.model.MimeType;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.dao.model.Workflow;
import edu.unc.mapseq.dao.model.WorkflowRun;
import edu.unc.mapseq.dao.model.WorkflowRunAttempt;
import edu.unc.mapseq.module.core.ZipCLI;
import edu.unc.mapseq.module.filter.FilterVariantCLI;
import edu.unc.mapseq.module.gatk.GATKApplyRecalibration;
import edu.unc.mapseq.module.gatk.GATKDepthOfCoverageCLI;
import edu.unc.mapseq.module.gatk.GATKDownsamplingType;
import edu.unc.mapseq.module.gatk.GATKPhoneHomeType;
import edu.unc.mapseq.module.gatk.GATKTableRecalibration;
import edu.unc.mapseq.module.picard.PicardSortOrderType;
import edu.unc.mapseq.module.picard.PicardSortSAMCLI;
import edu.unc.mapseq.module.samtools.SAMToolsIndex;
import edu.unc.mapseq.module.samtools.SAMToolsIndexCLI;
import edu.unc.mapseq.module.samtools.SAMToolsViewCLI;
import edu.unc.mapseq.workflow.WorkflowException;
import edu.unc.mapseq.workflow.WorkflowUtil;
import edu.unc.mapseq.workflow.impl.AbstractSampleWorkflow;
import edu.unc.mapseq.workflow.impl.WorkflowJobFactory;

public class NCGenesDXWorkflow extends AbstractSampleWorkflow {

    private final Logger logger = LoggerFactory.getLogger(NCGenesDXWorkflow.class);

    public NCGenesDXWorkflow() {
        super();
    }

    @Override
    public String getName() {
        return NCGenesDXWorkflow.class.getSimpleName().replace("Workflow", "");
    }

    @Override
    public String getVersion() {
        ResourceBundle bundle = ResourceBundle.getBundle("edu/unc/mapseq/workflow/ncgenes/dx/workflow");
        String version = bundle.getString("version");
        return StringUtils.isNotEmpty(version) ? version : "0.0.1-SNAPSHOT";
    }

    @Override
    public Graph<CondorJob, CondorJobEdge> createGraph() throws WorkflowException {
        logger.info("ENTERING createGraph()");

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(
                CondorJobEdge.class);

        int count = 0;
        String version = null;
        String dx = null;
        String summaryCoverageThreshold = "1,2,5,8,10,15,20,30,50";

        Set<Sample> sampleSet = getAggregatedSamples();
        logger.info("sampleSet.size(): {}", sampleSet.size());

        String siteName = getWorkflowBeanService().getAttributes().get("siteName");
        String referenceSequence = getWorkflowBeanService().getAttributes().get("referenceSequence");
        String icSNPIntervalList = getWorkflowBeanService().getAttributes().get("icSNPIntervalList");
        Boolean isIncidental = Boolean.FALSE;
        Workflow ncgenesWorkflow = null;
        try {
            ncgenesWorkflow = getWorkflowBeanService().getMaPSeqDAOBean().getWorkflowDAO().findByName("NCGenes").get(0);
        } catch (MaPSeqDAOException e1) {
            e1.printStackTrace();
        }

        WorkflowRunAttempt attempt = getWorkflowRunAttempt();
        WorkflowRun workflowRun = attempt.getWorkflowRun();

        for (Sample sample : sampleSet) {

            if ("Undetermined".equals(sample.getBarcode())) {
                continue;
            }

            logger.debug(sample.toString());

            File outputDirectory = new File(sample.getOutputDirectory(), getName());
            File tmpDirectory = new File(outputDirectory, "tmp");
            tmpDirectory.mkdirs();

            Set<Attribute> attributeSet = workflowRun.getAttributes();
            if (attributeSet != null && !attributeSet.isEmpty()) {
                Iterator<Attribute> attributeIter = attributeSet.iterator();
                while (attributeIter.hasNext()) {
                    Attribute attribute = attributeIter.next();
                    if ("isIncidental".equals(attribute.getName())) {
                        isIncidental = Boolean.TRUE;
                    }
                    if ("GATKDepthOfCoverage.interval_list.version".equals(attribute.getName())) {
                        version = attribute.getValue();
                    }
                    if ("SAMToolsView.dx.id".equals(attribute.getName())) {
                        dx = attribute.getValue();
                    }
                }
            }
            if (version == null || dx == null) {
                throw new WorkflowException("Both version and DX were null...returning empty dag");
            }

            String format = "/proj/renci/sequence_analysis/annotation/abeast/NCGenes/%1$s/exons_pm_0_v%1$s.interval_list";
            File intervalListByVersionFile = new File(String.format(format, version));
            if (!intervalListByVersionFile.exists()) {
                throw new WorkflowException("Interval list file does not exist: "
                        + intervalListByVersionFile.getAbsolutePath());
            }

            format = "/proj/renci/sequence_analysis/annotation/abeast/NCGenes/%1$s/genes_dxid_%2$s_v_%1$s.interval_list";
            if (isIncidental) {
                format = "/proj/renci/sequence_analysis/annotation/abeast/NCGenes/Incidental/incidental_%2$s_%1$s.interval_list";
            }
            File intervalListByDXAndVersionFile = new File(String.format(format, version, dx));
            if (!intervalListByDXAndVersionFile.exists()) {
                throw new WorkflowException("Interval list file does not exist: "
                        + intervalListByDXAndVersionFile.getAbsolutePath());
            }

            Set<FileData> fileDataSet = sample.getFileDatas();

            File bamFile = null;
            File bamIndexFile = null;

            List<File> potentialBAMFileList = WorkflowUtil.lookupFileByJobAndMimeTypeAndWorkflowId(fileDataSet,
                    getWorkflowBeanService().getMaPSeqDAOBean(), GATKTableRecalibration.class,
                    MimeType.APPLICATION_BAM, ncgenesWorkflow.getId());

            // assume that only one GATKTableRecalibration job exists
            if (!potentialBAMFileList.isEmpty()) {
                bamFile = potentialBAMFileList.get(0);
            }

            File ncgenesDirectory = new File(sample.getOutputDirectory(), "NCGenes");

            for (File file : ncgenesDirectory.listFiles()) {
                if (file.getName().endsWith(".recal.bam")) {
                    bamFile = file;
                    break;
                }
            }

            if (bamFile == null) {
                logger.error("bam file to process was not found");
                throw new WorkflowException("bam file to process was not found");
            }

            List<File> potentialBAMIndexFileList = WorkflowUtil.lookupFileByJobAndMimeTypeAndWorkflowId(fileDataSet,
                    getWorkflowBeanService().getMaPSeqDAOBean(), SAMToolsIndex.class, MimeType.APPLICATION_BAM_INDEX,
                    ncgenesWorkflow.getId());

            if (!potentialBAMIndexFileList.isEmpty()) {
                for (File file : potentialBAMIndexFileList) {
                    if (bamFile != null && bamFile.getName().replace(".bam", ".bai").equals(file.getName())) {
                        bamIndexFile = file;
                        break;
                    }
                }
            }

            if (bamIndexFile == null) {
                logger.error("bam index file was not found");
                throw new WorkflowException("bam index file was not found");
            }

            try {

                // new job
                CondorJobBuilder builder = WorkflowJobFactory.createJob(++count, GATKDepthOfCoverageCLI.class,
                        attempt.getId(), sample.getId()).siteName(siteName);
                String outputPrefix = bamFile.getName().replace(".bam", String.format(".coverage.v%s.gene", version));
                builder.addArgument(GATKDepthOfCoverageCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                        .addArgument(GATKDepthOfCoverageCLI.WORKDIRECTORY, outputDirectory.getAbsolutePath())
                        .addArgument(GATKDepthOfCoverageCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                        .addArgument(GATKDepthOfCoverageCLI.INTERVALMERGING, "OVERLAPPING_ONLY")
                        .addArgument(GATKDepthOfCoverageCLI.REFERENCESEQUENCE, referenceSequence)
                        .addArgument(GATKDepthOfCoverageCLI.VALIDATIONSTRICTNESS, "LENIENT")
                        .addArgument(GATKDepthOfCoverageCLI.OMITDEPTHOUTPUTATEACHBASE)
                        .addArgument(GATKDepthOfCoverageCLI.INPUTFILE, bamFile.getAbsolutePath())
                        .addArgument(GATKDepthOfCoverageCLI.INTERVALS, intervalListByVersionFile.getAbsolutePath())
                        .addArgument(GATKDepthOfCoverageCLI.OUTPUTPREFIX, outputPrefix);
                if (summaryCoverageThreshold.contains(",")) {
                    for (String sct : StringUtils.split(summaryCoverageThreshold, ",")) {
                        builder.addArgument(GATKDepthOfCoverageCLI.SUMMARYCOVERAGETHRESHOLD, sct);
                    }
                }
                CondorJob gatkGeneDepthOfCoverageJob = builder.build();
                logger.info(gatkGeneDepthOfCoverageJob.toString());
                graph.addVertex(gatkGeneDepthOfCoverageJob);

                // new job
                builder = WorkflowJobFactory.createJob(++count, SAMToolsViewCLI.class, attempt.getId(), sample.getId())
                        .siteName(siteName);
                File samtoolsViewOutput = new File(outputDirectory, bamFile.getName().replace(".bam",
                        String.format(".filtered_by_dxid_%s_v%s.bam", dx, version)));
                builder.addArgument(SAMToolsViewCLI.BAMFORMAT)
                        .addArgument(SAMToolsViewCLI.OUTPUT, samtoolsViewOutput.getAbsolutePath())
                        .addArgument(SAMToolsViewCLI.REGIONSFILE, intervalListByDXAndVersionFile.getAbsolutePath())
                        .addArgument(SAMToolsViewCLI.INPUT, bamFile.getAbsolutePath());
                CondorJob samtoolsViewJob = builder.build();
                logger.info(samtoolsViewJob.toString());
                graph.addVertex(samtoolsViewJob);

                // new job
                builder = WorkflowJobFactory
                        .createJob(++count, PicardSortSAMCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
                File picardSortOutput = new File(outputDirectory, samtoolsViewOutput.getName().replace(".bam",
                        ".sorted.bam"));
                builder.addArgument(PicardSortSAMCLI.INPUT, samtoolsViewOutput.getAbsolutePath())
                        .addArgument(PicardSortSAMCLI.OUTPUT, picardSortOutput.getAbsolutePath())
                        .addArgument(PicardSortSAMCLI.SORTORDER,
                                PicardSortOrderType.COORDINATE.toString().toLowerCase());
                CondorJob picardSortSAMJob = builder.build();
                logger.info(picardSortSAMJob.toString());
                graph.addVertex(picardSortSAMJob);
                graph.addEdge(samtoolsViewJob, picardSortSAMJob);

                // new job
                builder = WorkflowJobFactory
                        .createJob(++count, SAMToolsIndexCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
                File picardSortSAMIndexOut = new File(outputDirectory, picardSortOutput.getName().replace(".bam",
                        ".bai"));
                builder.addArgument(SAMToolsIndexCLI.INPUT, picardSortOutput.getAbsolutePath()).addArgument(
                        SAMToolsIndexCLI.OUTPUT, picardSortSAMIndexOut.getAbsolutePath());
                CondorJob picardSortSAMIndexJob = builder.build();
                logger.info(picardSortSAMIndexJob.toString());
                graph.addVertex(picardSortSAMIndexJob);
                graph.addEdge(picardSortSAMJob, picardSortSAMIndexJob);

                // new job
                builder = WorkflowJobFactory.createJob(++count, ZipCLI.class, attempt.getId(), sample.getId())
                        .siteName(siteName);
                File zipOutputFile = new File(outputDirectory, picardSortOutput.getName().replace(".bam", ".zip"));
                builder.addArgument(ZipCLI.ENTRY, picardSortOutput.getAbsolutePath())
                        .addArgument(ZipCLI.ENTRY, picardSortSAMIndexOut.getAbsolutePath())
                        .addArgument(ZipCLI.OUTPUT, zipOutputFile.getAbsolutePath());
                CondorJob zipJob = builder.build();
                logger.info(zipJob.toString());
                graph.addVertex(zipJob);
                graph.addEdge(picardSortSAMIndexJob, zipJob);

                File gatkApplyRecalibrationOut = null;
                List<File> potentialGATKApplyRecalibrationFileList = WorkflowUtil
                        .lookupFileByJobAndMimeTypeAndWorkflowId(fileDataSet, getWorkflowBeanService()
                                .getMaPSeqDAOBean(), GATKApplyRecalibration.class, MimeType.TEXT_VCF, ncgenesWorkflow
                                .getId());
                // assume that only one GATKApplyRecalibrationCLI job exists
                if (potentialGATKApplyRecalibrationFileList.size() > 0) {
                    gatkApplyRecalibrationOut = potentialGATKApplyRecalibrationFileList.get(0);
                }

                if (gatkApplyRecalibrationOut == null) {
                    File ncgenesOutputDirectory = new File(sample.getOutputDirectory(), "NCGenes");
                    List<File> files = Arrays.asList(ncgenesOutputDirectory.listFiles());
                    if (files != null && !files.isEmpty()) {
                        for (File f : files) {
                            if (f.getName().endsWith(".recalibrated.filtered.vcf")) {
                                gatkApplyRecalibrationOut = f;
                                break;
                            }
                        }
                    }
                }

                if (gatkApplyRecalibrationOut == null) {
                    logger.error("gatkApplyRecalibrationOut file to process was not found");
                    throw new WorkflowException("gatkApplyRecalibrationOut file to process was not found");
                }

                // new job
                builder = WorkflowJobFactory
                        .createJob(++count, FilterVariantCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
                File filterVariantOutput = new File(outputDirectory, bamFile.getName().replace(".bam",
                        String.format(".filtered_by_dxid_%s_v%s.vcf", dx, version)));
                builder.addArgument(FilterVariantCLI.INTERVALLIST, intervalListByDXAndVersionFile.getAbsolutePath())
                        .addArgument(FilterVariantCLI.INPUT, gatkApplyRecalibrationOut.getAbsolutePath())
                        .addArgument(FilterVariantCLI.OUTPUT, filterVariantOutput.getAbsolutePath());
                CondorJob filterVariantJob = builder.build();
                logger.info(filterVariantJob.toString());
                graph.addVertex(filterVariantJob);

            } catch (Exception e) {
                throw new WorkflowException(e);
            }

        }

        return graph;
    }

    @Override
    public void postRun() throws WorkflowException {
        logger.info("ENTERING postRun()");

        Set<Sample> sampleSet = getAggregatedSamples();

        String dx = null;
        String version = null;

        WorkflowRunAttempt attempt = getWorkflowRunAttempt();
        WorkflowRun workflowRun = attempt.getWorkflowRun();

        ExecutorService executorService = Executors.newSingleThreadExecutor();

        for (Sample sample : sampleSet) {

            if ("Undetermined".equals(sample.getBarcode())) {
                continue;
            }

            Set<Attribute> attributeSet = workflowRun.getAttributes();
            if (attributeSet != null && !attributeSet.isEmpty()) {
                Iterator<Attribute> attributeIter = attributeSet.iterator();
                while (attributeIter.hasNext()) {
                    Attribute attribute = attributeIter.next();
                    if ("GATKDepthOfCoverage.interval_list.version".equals(attribute.getName())) {
                        version = attribute.getValue();
                    }
                    if ("SAMToolsView.dx.id".equals(attribute.getName())) {
                        dx = attribute.getValue();
                    }
                }
            }

            RegisterToIRODSRunnable runnable = new RegisterToIRODSRunnable();
            runnable.setMaPSeqDAOBean(getWorkflowBeanService().getMaPSeqDAOBean());
            runnable.setMaPSeqConfigurationService(getWorkflowBeanService().getMaPSeqConfigurationService());
            runnable.setSampleId(sample.getId());
            runnable.setDx(dx);
            runnable.setVersion(version);
            executorService.submit(runnable);

        }
        executorService.shutdown();

    }

}
