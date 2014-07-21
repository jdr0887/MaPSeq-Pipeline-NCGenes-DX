package edu.unc.mapseq.workflow.ncgenes.dx;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.ResourceBundle;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.renci.common.exec.BashExecutor;
import org.renci.common.exec.CommandInput;
import org.renci.common.exec.CommandOutput;
import org.renci.common.exec.Executor;
import org.renci.common.exec.ExecutorException;
import org.renci.jlrm.condor.CondorJob;
import org.renci.jlrm.condor.CondorJobBuilder;
import org.renci.jlrm.condor.CondorJobEdge;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.config.RunModeType;
import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.EntityAttribute;
import edu.unc.mapseq.dao.model.FileData;
import edu.unc.mapseq.dao.model.HTSFSample;
import edu.unc.mapseq.dao.model.MimeType;
import edu.unc.mapseq.dao.model.SequencerRun;
import edu.unc.mapseq.dao.model.Workflow;
import edu.unc.mapseq.dao.model.WorkflowRun;
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
import edu.unc.mapseq.workflow.impl.AbstractWorkflow;
import edu.unc.mapseq.workflow.impl.IRODSBean;
import edu.unc.mapseq.workflow.impl.WorkflowJobFactory;

public class NCGenesDXWorkflow extends AbstractWorkflow {

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

        Set<HTSFSample> htsfSampleSet = getAggregateHTSFSampleSet();
        logger.info("htsfSampleSet.size(): {}", htsfSampleSet.size());

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

        WorkflowRun workflowRun = getWorkflowPlan().getWorkflowRun();

        for (HTSFSample htsfSample : htsfSampleSet) {

            if ("Undetermined".equals(htsfSample.getBarcode())) {
                continue;
            }

            SequencerRun sequencerRun = htsfSample.getSequencerRun();
            File outputDirectory = createOutputDirectory(sequencerRun.getName(), htsfSample, getName()
                    .replace("DX", ""), getVersion());

            Set<EntityAttribute> attributeSet = workflowRun.getAttributes();
            if (attributeSet != null && !attributeSet.isEmpty()) {
                Iterator<EntityAttribute> attributeIter = attributeSet.iterator();
                while (attributeIter.hasNext()) {
                    EntityAttribute attribute = attributeIter.next();
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

            Integer laneIndex = htsfSample.getLaneIndex();
            logger.debug("laneIndex = {}", laneIndex);
            Set<FileData> fileDataSet = htsfSample.getFileDatas();

            File bamFile = null;
            File bamIndexFile = null;

            List<File> potentialBAMFileList = WorkflowUtil.lookupFileByJobAndMimeTypeAndWorkflowId(fileDataSet,
                    getWorkflowBeanService().getMaPSeqDAOBean(), GATKTableRecalibration.class,
                    MimeType.APPLICATION_BAM, ncgenesWorkflow.getId());

            // assume that only one GATKTableRecalibration job exists
            if (potentialBAMFileList.size() > 0) {
                bamFile = potentialBAMFileList.get(0);
            }

            if (bamFile == null) {
                logger.error("bam file to process was not found");
                throw new WorkflowException("bam file to process was not found");
            }

            List<File> potentialBAMIndexFileList = WorkflowUtil.lookupFileByJobAndMimeTypeAndWorkflowId(fileDataSet,
                    getWorkflowBeanService().getMaPSeqDAOBean(), SAMToolsIndex.class, MimeType.APPLICATION_BAM_INDEX,
                    ncgenesWorkflow.getId());

            if (potentialBAMIndexFileList.size() > 0) {
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
                        getWorkflowPlan(), htsfSample).siteName(siteName);
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
                builder = WorkflowJobFactory.createJob(++count, SAMToolsViewCLI.class, getWorkflowPlan(), htsfSample)
                        .siteName(siteName);
                File samtoolsViewOutput = new File(outputDirectory, bamFile.getName().replace(".bam", ".filtered.bam"));
                builder.addArgument(SAMToolsViewCLI.BAMFORMAT)
                        .addArgument(SAMToolsViewCLI.OUTPUT, samtoolsViewOutput.getAbsolutePath())
                        .addArgument(SAMToolsViewCLI.REGIONSFILE, intervalListByDXAndVersionFile.getAbsolutePath())
                        .addArgument(SAMToolsViewCLI.INPUT, bamFile.getAbsolutePath());
                CondorJob samtoolsViewJob = builder.build();
                logger.info(samtoolsViewJob.toString());
                graph.addVertex(samtoolsViewJob);

                // new job
                builder = WorkflowJobFactory.createJob(++count, PicardSortSAMCLI.class, getWorkflowPlan(), htsfSample)
                        .siteName(siteName);
                File picardSortOutput = new File(outputDirectory, samtoolsViewOutput.getName().replace(".bam",
                        String.format(".sorted.filtered_by_dxid_%s_v%s.bam", dx, version)));
                builder.addArgument(PicardSortSAMCLI.INPUT, samtoolsViewOutput.getAbsolutePath())
                        .addArgument(PicardSortSAMCLI.OUTPUT, picardSortOutput.getAbsolutePath())
                        .addArgument(PicardSortSAMCLI.SORTORDER,
                                PicardSortOrderType.COORDINATE.toString().toLowerCase());
                CondorJob picardSortSAMJob = builder.build();
                logger.info(picardSortSAMJob.toString());
                graph.addVertex(picardSortSAMJob);
                graph.addEdge(samtoolsViewJob, picardSortSAMJob);

                // new job
                builder = WorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, getWorkflowPlan(), htsfSample)
                        .siteName(siteName);
                File picardSortSAMIndexOut = new File(outputDirectory, picardSortOutput.getName().replace(".bam",
                        ".bai"));
                builder.addArgument(SAMToolsIndexCLI.INPUT, picardSortOutput.getAbsolutePath()).addArgument(
                        SAMToolsIndexCLI.OUTPUT, picardSortSAMIndexOut.getAbsolutePath());
                CondorJob picardSortSAMIndexJob = builder.build();
                logger.info(picardSortSAMIndexJob.toString());
                graph.addVertex(picardSortSAMIndexJob);
                graph.addEdge(picardSortSAMJob, picardSortSAMIndexJob);

                // new job
                builder = WorkflowJobFactory.createJob(++count, ZipCLI.class, getWorkflowPlan(), htsfSample).siteName(
                        siteName);
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
                    logger.error("gatkApplyRecalibrationOut file to process was not found");
                    throw new WorkflowException("gatkApplyRecalibrationOut file to process was not found");
                }

                // new job
                builder = WorkflowJobFactory.createJob(++count, FilterVariantCLI.class, getWorkflowPlan(), htsfSample)
                        .siteName(siteName);
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

        Set<HTSFSample> htsfSampleSet = getAggregateHTSFSampleSet();
        logger.info("htsfSampleSet.size(): {}", htsfSampleSet.size());

        RunModeType runMode = getWorkflowBeanService().getMaPSeqConfigurationService().getRunMode();

        String dx = null;
        String version = null;

        Workflow ncgenesWorkflow = null;
        try {
            ncgenesWorkflow = getWorkflowBeanService().getMaPSeqDAOBean().getWorkflowDAO().findByName("NCGenes").get(0);
        } catch (MaPSeqDAOException e1) {
            e1.printStackTrace();
        }

        WorkflowRun workflowRun = getWorkflowPlan().getWorkflowRun();

        for (HTSFSample htsfSample : htsfSampleSet) {

            if ("Undetermined".equals(htsfSample.getBarcode())) {
                continue;
            }

            SequencerRun sequencerRun = htsfSample.getSequencerRun();
            File outputDirectory = createOutputDirectory(sequencerRun.getName(), htsfSample, getName()
                    .replace("DX", ""), getVersion());
            File tmpDir = new File(outputDirectory, "tmp");
            if (!tmpDir.exists()) {
                tmpDir.mkdirs();
            }

            Set<EntityAttribute> attributeSet = workflowRun.getAttributes();
            if (attributeSet != null && !attributeSet.isEmpty()) {
                Iterator<EntityAttribute> attributeIter = attributeSet.iterator();
                while (attributeIter.hasNext()) {
                    EntityAttribute attribute = attributeIter.next();
                    if ("GATKDepthOfCoverage.interval_list.version".equals(attribute.getName())) {
                        version = attribute.getValue();
                    }
                    if ("SAMToolsView.dx.id".equals(attribute.getName())) {
                        dx = attribute.getValue();
                    }
                }
            }
            if (version == null && dx == null) {
                logger.warn("Both version and DX were null...returning empty dag");
                return;
            }

            Integer laneIndex = htsfSample.getLaneIndex();
            logger.debug("laneIndex = {}", laneIndex);
            Set<FileData> fileDataSet = htsfSample.getFileDatas();

            // assumption: a dash is used as a delimiter between a participantId and the external code
            int idx = htsfSample.getName().lastIndexOf("-");
            String participantId = idx != -1 ? htsfSample.getName().substring(0, idx) : htsfSample.getName();

            File bamFile = null;

            List<File> potentialBAMFileList = WorkflowUtil.lookupFileByJobAndMimeTypeAndWorkflowId(fileDataSet,
                    getWorkflowBeanService().getMaPSeqDAOBean(), GATKTableRecalibration.class,
                    MimeType.APPLICATION_BAM, ncgenesWorkflow.getId());

            // assume that only one GATKTableRecalibration job exists
            if (potentialBAMFileList.size() > 0) {
                bamFile = potentialBAMFileList.get(0);
            }

            if (bamFile == null) {
                logger.error("bam file to process was not found");
                throw new WorkflowException("bam file to process was not found");
            }

            String irodsHome = System.getenv("NCGENES_IRODS_HOME");
            if (StringUtils.isEmpty(irodsHome)) {
                logger.error("irodsHome is not set");
                return;
            }

            String ncgenesIRODSDirectory;

            switch (runMode) {
                case DEV:
                case STAGING:
                    ncgenesIRODSDirectory = String.format("/genomicsDataGridZone/sequence_data/%s/ncgenes/%s/%s",
                            runMode.toString().toLowerCase(), participantId, version);
                    break;
                case PROD:
                default:
                    ncgenesIRODSDirectory = String.format("/genomicsDataGridZone/sequence_data/ncgenes/%s/%s",
                            participantId, version);
                    break;
            }

            CommandOutput commandOutput = null;

            List<CommandInput> commandInputList = new ArrayList<CommandInput>();

            CommandInput commandInput = new CommandInput();
            commandInput.setCommand(String.format("%s/bin/imkdir -p %s", irodsHome, ncgenesIRODSDirectory));
            commandInput.setWorkDir(tmpDir);
            commandInputList.add(commandInput);

            commandInput = new CommandInput();
            commandInput.setCommand(String.format("%s/bin/imeta add -C %s Project NCGENES", irodsHome,
                    ncgenesIRODSDirectory));
            commandInput.setWorkDir(tmpDir);
            commandInputList.add(commandInput);

            commandInput = new CommandInput();
            commandInput.setCommand(String.format("%s/bin/imeta add -C %s ParticipantID %s NCGENES", irodsHome,
                    ncgenesIRODSDirectory, participantId));
            commandInput.setWorkDir(tmpDir);
            commandInputList.add(commandInput);

            List<IRODSBean> files2RegisterToIRODS = new ArrayList<IRODSBean>();

            files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, bamFile.getName().replace(".bam",
                    String.format(".coverage.v%s.gene.sample_cumulative_coverage_counts", version))),
                    "GeneCoverageCount", version, null, runMode));
            files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, bamFile.getName().replace(".bam",
                    String.format(".coverage.v%s.gene.sample_cumulative_coverage_proportions", version))),
                    "GeneCoverageProportions", version, null, runMode));
            files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, bamFile.getName().replace(".bam",
                    String.format(".coverage.v%s.gene.sample_interval_statistics", version))),
                    "GeneIntervalStatistics", version, null, runMode));
            files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, bamFile.getName().replace(".bam",
                    String.format(".coverage.v%s.gene.sample_interval_summary", version))), "GeneIntervalSummary",
                    version, null, runMode));
            files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, bamFile.getName().replace(".bam",
                    String.format(".coverage.v%s.gene.sample_statistics", version))), "GeneSampleStatistics", version,
                    null, runMode));
            files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, bamFile.getName().replace(".bam",
                    String.format(".coverage.v%s.gene.sample_summary", version))), "GeneSampleSummary", version, null,
                    runMode));

            File samtoolsViewOutput = new File(outputDirectory, bamFile.getName().replace(".bam", ".filtered.bam"));
            File picardSortOutput = new File(outputDirectory, samtoolsViewOutput.getName().replace(".bam",
                    String.format(".sorted.filtered_by_dxid_%s_v%s.bam", dx, version)));
            File picardSortSAMIndexOut = new File(outputDirectory, picardSortOutput.getName().replace(".bam", ".bai"));
            files2RegisterToIRODS.add(new IRODSBean(picardSortSAMIndexOut, "FilteredBamIndex", version, dx, runMode));

            File zipOutputFile = new File(outputDirectory, picardSortOutput.getName().replace(".bam", ".zip"));
            files2RegisterToIRODS.add(new IRODSBean(zipOutputFile, "FilteredBamZip", version, dx, runMode));

            File filterVariantOutput = new File(outputDirectory, bamFile.getName().replace(".bam",
                    String.format(".filtered_by_dxid_%s_v%s.vcf", dx, version)));
            files2RegisterToIRODS.add(new IRODSBean(filterVariantOutput, "FilteredVcf", version, dx, runMode));

            for (IRODSBean bean : files2RegisterToIRODS) {

                commandInput = new CommandInput();
                commandInput.setExitImmediately(Boolean.FALSE);

                StringBuilder registerCommandSB = new StringBuilder();
                String registrationCommand = String.format("%s/bin/ireg -f %s %s/%s", irodsHome, bean.getFile()
                        .getAbsolutePath(), ncgenesIRODSDirectory, bean.getFile().getName());
                String deRegistrationCommand = String.format("%s/bin/irm -U %s/%s", irodsHome, ncgenesIRODSDirectory,
                        bean.getFile().getName());
                registerCommandSB.append(registrationCommand).append("\n");
                registerCommandSB.append(String.format("if [ $? != 0 ]; then %s; %s; fi%n", deRegistrationCommand,
                        registrationCommand));
                commandInput.setCommand(registerCommandSB.toString());
                commandInput.setWorkDir(tmpDir);
                commandInputList.add(commandInput);

                commandInput = new CommandInput();
                commandInput.setCommand(String.format("%s/bin/imeta add -d %s/%s ParticipantID %s NCGENES", irodsHome,
                        ncgenesIRODSDirectory, bean.getFile().getName(), participantId));
                commandInput.setWorkDir(tmpDir);
                commandInputList.add(commandInput);

                commandInput = new CommandInput();
                commandInput.setCommand(String.format("%s/bin/imeta add -d %s/%s FileType %s NCGENES", irodsHome,
                        ncgenesIRODSDirectory, bean.getFile().getName(), bean.getType()));
                commandInput.setWorkDir(tmpDir);
                commandInputList.add(commandInput);

                if (StringUtils.isNotEmpty(bean.getDx())) {
                    commandInput = new CommandInput();
                    commandInput.setCommand(String.format("%s/bin/imeta add -d %s/%s DxID %s NCGENES", irodsHome,
                            ncgenesIRODSDirectory, bean.getFile().getName(), bean.getDx()));
                    commandInput.setWorkDir(tmpDir);
                    commandInputList.add(commandInput);
                }

                if (StringUtils.isNotEmpty(bean.getVersion())) {
                    commandInput = new CommandInput();
                    commandInput.setCommand(String.format("%s/bin/imeta add -d %s/%s DxVersion %s NCGENES", irodsHome,
                            ncgenesIRODSDirectory, bean.getFile().getName(), bean.getVersion()));
                    commandInput.setWorkDir(tmpDir);
                    commandInputList.add(commandInput);
                }

                commandInput = new CommandInput();
                commandInput.setCommand(String.format("%s/bin/imeta add -d %s/%s System %s NCGENES", irodsHome,
                        ncgenesIRODSDirectory, bean.getFile().getName(),
                        StringUtils.capitalize(bean.getRunMode().toString().toLowerCase())));
                commandInput.setWorkDir(tmpDir);
                commandInputList.add(commandInput);

            }

            File mapseqrc = new File(System.getProperty("user.home"), ".mapseqrc");
            Executor executor = BashExecutor.getInstance();

            for (CommandInput ci : commandInputList) {
                try {
                    commandOutput = executor.execute(ci, mapseqrc);
                    logger.info("commandOutput.getExitCode(): {}", commandOutput.getExitCode());
                    logger.debug("commandOutput.getStdout(): {}", commandOutput.getStdout());
                } catch (ExecutorException e) {
                    if (commandOutput != null) {
                        logger.warn("commandOutput.getStderr(): {}", commandOutput.getStderr());
                    }
                }
            }

        }

    }

}
