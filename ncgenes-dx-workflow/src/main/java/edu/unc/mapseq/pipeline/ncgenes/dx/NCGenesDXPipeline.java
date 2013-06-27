package edu.unc.mapseq.pipeline.ncgenes.dx;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
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
import edu.unc.mapseq.pipeline.AbstractPipeline;
import edu.unc.mapseq.pipeline.IRODSBean;
import edu.unc.mapseq.pipeline.PipelineException;
import edu.unc.mapseq.pipeline.PipelineJobFactory;
import edu.unc.mapseq.pipeline.PipelineUtil;

public class NCGenesDXPipeline extends AbstractPipeline {

    private final Logger logger = LoggerFactory.getLogger(NCGenesDXPipeline.class);

    public NCGenesDXPipeline() {
        super();
    }

    @Override
    public String getName() {
        return NCGenesDXPipeline.class.getSimpleName().replace("Pipeline", "");
    }

    @Override
    public String getVersion() {
        ResourceBundle bundle = ResourceBundle.getBundle("edu/unc/mapseq/pipeline/ncgenes/dx/pipeline");
        String version = bundle.getString("version");
        return StringUtils.isNotEmpty(version) ? version : "0.0.1-SNAPSHOT";
    }

    @Override
    public Graph<CondorJob, CondorJobEdge> createGraph() throws PipelineException {
        logger.info("ENTERING createGraph()");

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(
                CondorJobEdge.class);

        int count = 0;
        String version = null;
        String dx = null;
        String summaryCoverageThreshold = "1,2,5,8,10,15,20,30,50";

        if (getWorkflowPlan().getSequencerRun() == null && getWorkflowPlan().getHTSFSamples() == null) {
            logger.error("Don't have either sequencerRun and htsfSample");
            throw new PipelineException("Don't have either sequencerRun and htsfSample");
        }

        Set<HTSFSample> htsfSampleSet = new HashSet<HTSFSample>();

        if (getWorkflowPlan().getSequencerRun() != null) {
            logger.info("sequencerRun: {}", getWorkflowPlan().getSequencerRun().toString());
            try {
                htsfSampleSet.addAll(getPipelineBeanService().getMaPSeqDAOBean().getHTSFSampleDAO()
                        .findBySequencerRunId(getWorkflowPlan().getSequencerRun().getId()));
            } catch (MaPSeqDAOException e) {
                e.printStackTrace();
            }
        }

        if (getWorkflowPlan().getHTSFSamples() != null) {
            htsfSampleSet.addAll(getWorkflowPlan().getHTSFSamples());
        }

        logger.info("htsfSampleSet.size(): {}", htsfSampleSet.size());

        String siteName = getPipelineBeanService().getAttributes().get("siteName");
        String referenceSequence = getPipelineBeanService().getAttributes().get("referenceSequence");
        String icSNPIntervalList = getPipelineBeanService().getAttributes().get("icSNPIntervalList");

        for (HTSFSample htsfSample : htsfSampleSet) {

            if ("Undetermined".equals(htsfSample.getBarcode())) {
                continue;
            }

            SequencerRun sequencerRun = htsfSample.getSequencerRun();
            File outputDirectory = createOutputDirectory(sequencerRun.getName(), htsfSample, getName()
                    .replace("DX", ""));

            Set<EntityAttribute> attributeSet = htsfSample.getAttributes();
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

            if (version == null || dx == null) {
                throw new PipelineException("Both version and DX were null...returning empty dag");
            }

            File intervalListByVersionFile = new File(String.format(
                    "/proj/renci/sequence_analysis/annotation/abeast/NCGenes/%1$s/exons_pm_0_v%1$s.interval_list",
                    version));
            if (!intervalListByVersionFile.exists()) {
                throw new PipelineException("Interval list file does not exist: "
                        + intervalListByVersionFile.getAbsolutePath());
            }

            File intervalListByDXAndVersionFile = new File(
                    String.format(
                            "/proj/renci/sequence_analysis/annotation/abeast/NCGenes/%2$s/genes_dxid_%1$s_v_%2$s.interval_list",
                            dx, version));
            if (!intervalListByDXAndVersionFile.exists()) {
                throw new PipelineException("Interval list file does not exist: "
                        + intervalListByDXAndVersionFile.getAbsolutePath());
            }

            Integer laneIndex = htsfSample.getLaneIndex();
            logger.debug("laneIndex = {}", laneIndex);
            Set<FileData> fileDataSet = htsfSample.getFileDatas();

            File bamFile = null;
            File bamIndexFile = null;

            List<File> potentialBAMFileList = PipelineUtil
                    .lookupFileByJobAndMimeType(fileDataSet, getPipelineBeanService().getMaPSeqDAOBean(),
                            GATKTableRecalibration.class, MimeType.APPLICATION_BAM);

            // assume that only one GATKTableRecalibration job exists
            if (potentialBAMFileList.size() > 0) {
                bamFile = potentialBAMFileList.get(0);
            }

            if (bamFile == null) {
                logger.error("bam file to process was not found");
                throw new PipelineException("bam file to process was not found");
            }

            List<File> potentialBAMIndexFileList = PipelineUtil.lookupFileByJobAndMimeType(fileDataSet,
                    getPipelineBeanService().getMaPSeqDAOBean(), SAMToolsIndex.class, MimeType.APPLICATION_BAM_INDEX);

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
                throw new PipelineException("bam index file was not found");
            }

            try {

                // new job
                CondorJob gatkGeneDepthOfCoverageJob = PipelineJobFactory.createJob(++count,
                        GATKDepthOfCoverageCLI.class, getWorkflowPlan(), htsfSample);
                gatkGeneDepthOfCoverageJob.setSiteName(siteName);
                gatkGeneDepthOfCoverageJob.addArgument(GATKDepthOfCoverageCLI.PHONEHOME,
                        GATKPhoneHomeType.NO_ET.toString());
                gatkGeneDepthOfCoverageJob.addArgument(GATKDepthOfCoverageCLI.WORKDIRECTORY,
                        outputDirectory.getAbsolutePath());
                gatkGeneDepthOfCoverageJob.addArgument(GATKDepthOfCoverageCLI.DOWNSAMPLINGTYPE,
                        GATKDownsamplingType.NONE.toString());
                gatkGeneDepthOfCoverageJob.addArgument(GATKDepthOfCoverageCLI.INTERVALMERGING, "OVERLAPPING_ONLY");
                gatkGeneDepthOfCoverageJob.addArgument(GATKDepthOfCoverageCLI.REFERENCESEQUENCE, referenceSequence);
                gatkGeneDepthOfCoverageJob.addArgument(GATKDepthOfCoverageCLI.VALIDATIONSTRICTNESS, "LENIENT");
                gatkGeneDepthOfCoverageJob.addArgument(GATKDepthOfCoverageCLI.OMITDEPTHOUTPUTATEACHBASE);
                if (summaryCoverageThreshold.contains(",")) {
                    for (String sct : StringUtils.split(summaryCoverageThreshold, ",")) {
                        gatkGeneDepthOfCoverageJob.addArgument(GATKDepthOfCoverageCLI.SUMMARYCOVERAGETHRESHOLD, sct);
                    }
                }
                gatkGeneDepthOfCoverageJob.addArgument(GATKDepthOfCoverageCLI.INPUTFILE, bamFile.getAbsolutePath());
                gatkGeneDepthOfCoverageJob.addArgument(GATKDepthOfCoverageCLI.INTERVALS,
                        intervalListByVersionFile.getAbsolutePath());
                String outputPrefix = bamFile.getName().replace(".bam", String.format(".coverage.v%s.gene", version));
                gatkGeneDepthOfCoverageJob.addArgument(GATKDepthOfCoverageCLI.OUTPUTPREFIX, outputPrefix);
                graph.addVertex(gatkGeneDepthOfCoverageJob);

                // new job
                CondorJob samtoolsViewJob = PipelineJobFactory.createJob(++count, SAMToolsViewCLI.class,
                        getWorkflowPlan(), htsfSample);
                samtoolsViewJob.setSiteName(siteName);
                samtoolsViewJob.addArgument(SAMToolsViewCLI.BAMFORMAT);
                File samtoolsViewOutput = new File(outputDirectory, bamFile.getName().replace(".bam", ".filtered.bam"));
                samtoolsViewJob.addArgument(SAMToolsViewCLI.OUTPUT, samtoolsViewOutput.getAbsolutePath());
                samtoolsViewJob.addArgument(SAMToolsViewCLI.REGIONSFILE,
                        intervalListByDXAndVersionFile.getAbsolutePath());
                samtoolsViewJob.addArgument(SAMToolsViewCLI.INPUT, bamFile.getAbsolutePath());
                graph.addVertex(samtoolsViewJob);

                CondorJob picardSortSAMJob = PipelineJobFactory.createJob(++count, PicardSortSAMCLI.class,
                        getWorkflowPlan(), htsfSample);
                picardSortSAMJob.setSiteName(siteName);
                picardSortSAMJob.addArgument(PicardSortSAMCLI.INPUT, samtoolsViewOutput.getAbsolutePath());
                File picardSortOutput = new File(outputDirectory, samtoolsViewOutput.getName().replace(".bam",
                        String.format(".sorted.filtered_by_dxid_%s_v%s.bam", dx, version)));
                picardSortSAMJob.addArgument(PicardSortSAMCLI.OUTPUT, picardSortOutput.getAbsolutePath());
                picardSortSAMJob.addArgument(PicardSortSAMCLI.SORTORDER, PicardSortOrderType.COORDINATE.toString()
                        .toLowerCase());
                graph.addVertex(picardSortSAMJob);
                graph.addEdge(samtoolsViewJob, picardSortSAMJob);

                // new job
                CondorJob picardSortSAMIndexJob = PipelineJobFactory.createJob(++count, SAMToolsIndexCLI.class,
                        getWorkflowPlan(), htsfSample);
                picardSortSAMIndexJob.setSiteName(siteName);
                picardSortSAMIndexJob.addArgument(SAMToolsIndexCLI.INPUT, picardSortOutput.getAbsolutePath());
                File picardSortSAMIndexOut = new File(outputDirectory, picardSortOutput.getName().replace(".bam",
                        ".bai"));
                picardSortSAMIndexJob.addArgument(SAMToolsIndexCLI.OUTPUT, picardSortSAMIndexOut.getAbsolutePath());
                graph.addVertex(picardSortSAMIndexJob);
                graph.addEdge(picardSortSAMJob, picardSortSAMIndexJob);

                // new job
                CondorJob zipJob = PipelineJobFactory.createJob(++count, ZipCLI.class, getWorkflowPlan(), htsfSample);
                zipJob.setSiteName(siteName);
                zipJob.addArgument(ZipCLI.ENTRY, picardSortOutput.getAbsolutePath());
                zipJob.addArgument(ZipCLI.ENTRY, picardSortSAMIndexOut.getAbsolutePath());
                File zipOutputFile = new File(outputDirectory, picardSortOutput.getName().replace(".bam", ".zip"));
                zipJob.addArgument(ZipCLI.OUTPUT, zipOutputFile.getAbsolutePath());
                graph.addVertex(zipJob);
                graph.addEdge(picardSortSAMIndexJob, zipJob);

                File gatkApplyRecalibrationOut = null;
                List<File> potentialGATKApplyRecalibrationFileList = PipelineUtil.lookupFileByJobAndMimeType(
                        fileDataSet, getPipelineBeanService().getMaPSeqDAOBean(), GATKApplyRecalibration.class,
                        MimeType.TEXT_VCF);
                // assume that only one GATKApplyRecalibrationCLI job exists
                if (potentialGATKApplyRecalibrationFileList.size() > 0) {
                    gatkApplyRecalibrationOut = potentialGATKApplyRecalibrationFileList.get(0);
                }

                if (gatkApplyRecalibrationOut == null) {
                    logger.error("gatkApplyRecalibrationOut file to process was not found");
                    throw new PipelineException("gatkApplyRecalibrationOut file to process was not found");
                }

                // new job
                CondorJob filterVariantJob = PipelineJobFactory.createJob(++count, FilterVariantCLI.class,
                        getWorkflowPlan(), htsfSample);
                filterVariantJob.setSiteName(siteName);
                filterVariantJob.addArgument(FilterVariantCLI.INTERVALLIST,
                        intervalListByDXAndVersionFile.getAbsolutePath());
                filterVariantJob.addArgument(FilterVariantCLI.INPUT, gatkApplyRecalibrationOut.getAbsolutePath());
                File filterVariantOutput = new File(outputDirectory, bamFile.getName().replace(".bam",
                        String.format(".filtered_by_dxid_%s_v%s.vcf", dx, version)));
                filterVariantJob.addArgument(FilterVariantCLI.OUTPUT, filterVariantOutput.getAbsolutePath());
                graph.addVertex(filterVariantJob);

            } catch (Exception e) {
                throw new PipelineException(e);
            }

        }

        return graph;
    }

    @Override
    public void postRun() throws PipelineException {
        logger.info("ENTERING postRun()");

        Set<HTSFSample> htsfSampleSet = new HashSet<HTSFSample>();

        if (getWorkflowPlan().getSequencerRun() != null) {
            logger.info("sequencerRun: {}", getWorkflowPlan().getSequencerRun().toString());
            try {
                htsfSampleSet.addAll(getPipelineBeanService().getMaPSeqDAOBean().getHTSFSampleDAO()
                        .findBySequencerRunId(getWorkflowPlan().getSequencerRun().getId()));
            } catch (MaPSeqDAOException e) {
                e.printStackTrace();
            }
        }

        if (getWorkflowPlan().getHTSFSamples() != null) {
            logger.info("htsfSampleSet.size(): {}", htsfSampleSet.size());
            htsfSampleSet.addAll(getWorkflowPlan().getHTSFSamples());
        }

        RunModeType runMode = getPipelineBeanService().getMaPSeqConfigurationService().getRunMode();

        String dx = null;
        String version = null;

        for (HTSFSample htsfSample : htsfSampleSet) {

            if ("Undetermined".equals(htsfSample.getBarcode())) {
                continue;
            }

            SequencerRun sequencerRun = htsfSample.getSequencerRun();
            File outputDirectory = createOutputDirectory(sequencerRun.getName(), htsfSample, getName()
                    .replace("DX", ""));
            File tmpDir = new File(outputDirectory, "tmp");
            tmpDir.mkdirs();

            Set<EntityAttribute> attributeSet = htsfSample.getAttributes();
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

            List<File> potentialBAMFileList = PipelineUtil
                    .lookupFileByJobAndMimeType(fileDataSet, getPipelineBeanService().getMaPSeqDAOBean(),
                            GATKTableRecalibration.class, MimeType.APPLICATION_BAM);

            // assume that only one GATKTableRecalibration job exists
            if (potentialBAMFileList.size() > 0) {
                bamFile = potentialBAMFileList.get(0);
            }

            if (bamFile == null) {
                logger.error("bam file to process was not found");
                throw new PipelineException("bam file to process was not found");
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
