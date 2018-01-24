package edu.unc.mapseq.workflow.ncgenes.dx;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.collections.CollectionUtils;
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
import edu.unc.mapseq.module.sequencing.filter.FilterVariantCLI;
import edu.unc.mapseq.module.sequencing.gatk.GATKApplyRecalibration;
import edu.unc.mapseq.module.sequencing.gatk.GATKDepthOfCoverageCLI;
import edu.unc.mapseq.module.sequencing.gatk.GATKDownsamplingType;
import edu.unc.mapseq.module.sequencing.gatk.GATKPhoneHomeType;
import edu.unc.mapseq.module.sequencing.gatk.GATKTableRecalibration;
import edu.unc.mapseq.module.sequencing.picard.PicardSortOrderType;
import edu.unc.mapseq.module.sequencing.picard.PicardSortSAMCLI;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsIndex;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsIndexCLI;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsViewCLI;
import edu.unc.mapseq.workflow.WorkflowException;
import edu.unc.mapseq.workflow.core.WorkflowUtil;
import edu.unc.mapseq.workflow.sequencing.AbstractSequencingWorkflow;
import edu.unc.mapseq.workflow.sequencing.SequencingWorkflowJobFactory;
import edu.unc.mapseq.workflow.sequencing.SequencingWorkflowUtil;

public class NCGenesDXWorkflow extends AbstractSequencingWorkflow {

    private static final Logger logger = LoggerFactory.getLogger(NCGenesDXWorkflow.class);

    public NCGenesDXWorkflow() {
        super();
    }

    @Override
    public Graph<CondorJob, CondorJobEdge> createGraph() throws WorkflowException {
        logger.info("ENTERING createGraph()");

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(CondorJobEdge.class);

        int count = 0;
        String version = null;
        String dx = null;
        String summaryCoverageThreshold = "1,2,5,8,10,15,20,30,50";

        Set<Sample> sampleSet = SequencingWorkflowUtil.getAggregatedSamples(getWorkflowBeanService().getMaPSeqDAOBeanService(),
                getWorkflowRunAttempt());
        logger.info("sampleSet.size(): {}", sampleSet.size());

        String siteName = getWorkflowBeanService().getAttributes().get("siteName");
        String referenceSequence = getWorkflowBeanService().getAttributes().get("referenceSequence");
        // String icSNPIntervalList = getWorkflowBeanService().getAttributes().get("icSNPIntervalList");

        Boolean isIncidental = Boolean.FALSE;
        List<Workflow> workflows = new ArrayList<>();

        try {
            List<Workflow> foundWorkflows = getWorkflowBeanService().getMaPSeqDAOBeanService().getWorkflowDAO()
                    .findByName("NCGenesBaseline");
            if (CollectionUtils.isNotEmpty(foundWorkflows)) {
                workflows.addAll(foundWorkflows);
            }
            foundWorkflows = getWorkflowBeanService().getMaPSeqDAOBeanService().getWorkflowDAO().findByName("NCGenesBaselineMem");
            if (CollectionUtils.isNotEmpty(foundWorkflows)) {
                workflows.addAll(foundWorkflows);
            }
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

            File outputDirectory = SequencingWorkflowUtil.createOutputDirectory(sample,
                    getWorkflowRunAttempt().getWorkflowRun().getWorkflow());
            File tmpDirectory = new File(outputDirectory, "tmp");
            tmpDirectory.mkdirs();

            Set<Attribute> attributeSet = workflowRun.getAttributes();
            if (CollectionUtils.isNotEmpty(attributeSet)) {
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

            String dataDirectory = System.getenv("MAPSEQ_DATA_DIRECTORY");
            if (StringUtils.isEmpty(dataDirectory)) {
                dataDirectory = "/projects/mapseq/data";
            }

            String format = "%1$s/resources/annotation/abeast/NCGenes/%2$s/exons_pm_0_v%2$s.interval_list";
            File intervalListByVersionFile = new File(String.format(format, dataDirectory, version));
            if (!intervalListByVersionFile.exists()) {
                throw new WorkflowException("Interval list file does not exist: " + intervalListByVersionFile.getAbsolutePath());
            }

            format = "%1$s/resources/annotation/abeast/NCGenes/%2$s/genes_dxid_%3$s_v_%2$s.interval_list";
            if (isIncidental) {
                format = "%1$s/resources/annotation/abeast/NCGenes/Incidental/incidental_%3$s_%2$s.interval_list";
            }
            File intervalListByDXAndVersionFile = new File(String.format(format, dataDirectory, version, dx));
            if (!intervalListByDXAndVersionFile.exists()) {
                throw new WorkflowException("Interval list file does not exist: " + intervalListByDXAndVersionFile.getAbsolutePath());
            }

            Set<FileData> fileDataSet = sample.getFileDatas();

            File bamFile = null;
            File bamIndexFile = null;
            File gatkApplyRecalibrationOut = null;

            for (Workflow workflow : workflows) {

                bamFile = WorkflowUtil.findFileByJobAndMimeTypeAndWorkflowId(getWorkflowBeanService().getMaPSeqDAOBeanService(),
                        fileDataSet, GATKTableRecalibration.class, MimeType.APPLICATION_BAM, workflow.getId());
                File baselineDirectory = SequencingWorkflowUtil.createOutputDirectory(sample, workflow);

                if (bamFile == null) {
                    List<File> files = Arrays.asList(baselineDirectory.listFiles((a, b) -> {
                        if (b.endsWith(".recal.bam")) {
                            return true;
                        }
                        return false;
                    }));
                    if (CollectionUtils.isNotEmpty(files)) {
                        bamFile = files.get(0);
                    }
                }

                if (bamFile == null) {
                    logger.error("bam file to process was not found");
                    throw new WorkflowException("bam file to process was not found");
                }

                bamIndexFile = WorkflowUtil.findFileByJobAndMimeTypeAndWorkflowId(getWorkflowBeanService().getMaPSeqDAOBeanService(),
                        fileDataSet, SAMToolsIndex.class, MimeType.APPLICATION_BAM_INDEX, workflow.getId());

                if (bamIndexFile == null) {
                    for (File file : baselineDirectory.listFiles()) {
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

                gatkApplyRecalibrationOut = WorkflowUtil.findFileByJobAndMimeTypeAndWorkflowId(
                        getWorkflowBeanService().getMaPSeqDAOBeanService(), fileDataSet, GATKApplyRecalibration.class, MimeType.TEXT_VCF,
                        workflow.getId());

                if (gatkApplyRecalibrationOut == null) {
                    List<File> files = Arrays.asList(baselineDirectory.listFiles((a, b) -> {
                        if (b.endsWith(".recalibrated.filtered.vcf")) {
                            return true;
                        }
                        return false;
                    }));
                    if (CollectionUtils.isNotEmpty(files)) {
                        gatkApplyRecalibrationOut = files.get(0);
                    }
                }

                if (gatkApplyRecalibrationOut == null) {
                    logger.error("gatkApplyRecalibrationOut file to process was not found");
                    throw new WorkflowException("gatkApplyRecalibrationOut file to process was not found");
                }

            }

            try {

                // new job
                CondorJobBuilder builder = SequencingWorkflowJobFactory
                        .createJob(++count, GATKDepthOfCoverageCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
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
                builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsViewCLI.class, attempt.getId(), sample.getId())
                        .siteName(siteName);
                File samtoolsViewOutput = new File(outputDirectory,
                        bamFile.getName().replace(".bam", String.format(".filtered_by_dxid_%s_v%s.bam", dx, version)));
                builder.addArgument(SAMToolsViewCLI.BAMFORMAT).addArgument(SAMToolsViewCLI.OUTPUT, samtoolsViewOutput.getAbsolutePath())
                        .addArgument(SAMToolsViewCLI.REGIONSFILE, intervalListByDXAndVersionFile.getAbsolutePath())
                        .addArgument(SAMToolsViewCLI.INPUT, bamFile.getAbsolutePath());
                CondorJob samtoolsViewJob = builder.build();
                logger.info(samtoolsViewJob.toString());
                graph.addVertex(samtoolsViewJob);

                // new job
                builder = SequencingWorkflowJobFactory.createJob(++count, PicardSortSAMCLI.class, attempt.getId(), sample.getId())
                        .siteName(siteName);
                File picardSortOutput = new File(outputDirectory, samtoolsViewOutput.getName().replace(".bam", ".sorted.bam"));
                builder.addArgument(PicardSortSAMCLI.INPUT, samtoolsViewOutput.getAbsolutePath())
                        .addArgument(PicardSortSAMCLI.OUTPUT, picardSortOutput.getAbsolutePath())
                        .addArgument(PicardSortSAMCLI.SORTORDER, PicardSortOrderType.COORDINATE.toString().toLowerCase());
                CondorJob picardSortSAMJob = builder.build();
                logger.info(picardSortSAMJob.toString());
                graph.addVertex(picardSortSAMJob);
                graph.addEdge(samtoolsViewJob, picardSortSAMJob);

                // new job
                builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, attempt.getId(), sample.getId())
                        .siteName(siteName);
                File picardSortSAMIndexOut = new File(outputDirectory, picardSortOutput.getName().replace(".bam", ".bai"));
                builder.addArgument(SAMToolsIndexCLI.INPUT, picardSortOutput.getAbsolutePath()).addArgument(SAMToolsIndexCLI.OUTPUT,
                        picardSortSAMIndexOut.getAbsolutePath());
                CondorJob picardSortSAMIndexJob = builder.build();
                logger.info(picardSortSAMIndexJob.toString());
                graph.addVertex(picardSortSAMIndexJob);
                graph.addEdge(picardSortSAMJob, picardSortSAMIndexJob);

                // new job
                builder = SequencingWorkflowJobFactory.createJob(++count, ZipCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
                File zipOutputFile = new File(outputDirectory, picardSortOutput.getName().replace(".bam", ".zip"));
                builder.addArgument(ZipCLI.WORKDIR, outputDirectory.getAbsolutePath())
                        .addArgument(ZipCLI.ENTRY, picardSortOutput.getAbsolutePath())
                        .addArgument(ZipCLI.ENTRY, picardSortSAMIndexOut.getAbsolutePath())
                        .addArgument(ZipCLI.OUTPUT, zipOutputFile.getAbsolutePath());
                CondorJob zipJob = builder.build();
                logger.info(zipJob.toString());
                graph.addVertex(zipJob);
                graph.addEdge(picardSortSAMIndexJob, zipJob);

                // new job
                builder = SequencingWorkflowJobFactory.createJob(++count, FilterVariantCLI.class, attempt.getId(), sample.getId())
                        .siteName(siteName);
                File filterVariantOutput = new File(outputDirectory,
                        bamFile.getName().replace(".bam", String.format(".filtered_by_dxid_%s_v%s.vcf", dx, version)));
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

        try {
            ExecutorService es = Executors.newSingleThreadExecutor();

            RegisterToIRODSRunnable runnable = new RegisterToIRODSRunnable(getWorkflowBeanService().getMaPSeqDAOBeanService(),
                    getWorkflowRunAttempt());
            es.submit(runnable);

            es.shutdown();
            es.awaitTermination(1L, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

    }

}
