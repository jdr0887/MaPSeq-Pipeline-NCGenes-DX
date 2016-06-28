package edu.unc.mapseq.commons.ncgenes.dx;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.renci.common.exec.BashExecutor;
import org.renci.common.exec.CommandInput;
import org.renci.common.exec.CommandOutput;
import org.renci.common.exec.Executor;
import org.renci.common.exec.ExecutorException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.dao.MaPSeqDAOBeanService;
import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.MimeType;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.module.core.Zip;
import edu.unc.mapseq.module.sequencing.gatk.GATKDepthOfCoverage;
import edu.unc.mapseq.module.sequencing.picard.PicardSortSAM;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsIndex;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsView;
import edu.unc.mapseq.workflow.SystemType;
import edu.unc.mapseq.workflow.sequencing.IRODSBean;
import edu.unc.mapseq.workflow.sequencing.SequencingWorkflowUtil;

public class RegisterToIRODSRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(RegisterToIRODSRunnable.class);

    private MaPSeqDAOBeanService maPSeqDAOBeanService;

    private Long sampleId;

    private String version;

    private String dx;

    private SystemType system;

    public RegisterToIRODSRunnable(MaPSeqDAOBeanService maPSeqDAOBeanService, Long sampleId, String version, String dx, SystemType system) {
        super();
        this.maPSeqDAOBeanService = maPSeqDAOBeanService;
        this.sampleId = sampleId;
        this.version = version;
        this.dx = dx;
        this.system = system;
    }

    @Override
    public void run() {

        Sample sample;
        try {
            sample = maPSeqDAOBeanService.getSampleDAO().findById(sampleId);
        } catch (MaPSeqDAOException e1) {
            e1.printStackTrace();
            return;
        }

        if (sample == null) {
            logger.error("Sample not found");
            return;
        }

        File ncgenesOutputDirectory = new File(sample.getOutputDirectory(), "NCGenesBaseline");
        File tmpDir = new File(sample.getOutputDirectory(), "tmp");
        if (!tmpDir.exists()) {
            tmpDir.mkdirs();
        }

        File outputDirectory = new File(sample.getOutputDirectory(), "NCGenesDX");

        List<File> readPairList = SequencingWorkflowUtil.getReadPairList(sample);

        // assumption: a dash is used as a delimiter between a participantId and
        // the external code
        int idx = sample.getName().lastIndexOf("-");
        String participantId = idx != -1 ? sample.getName().substring(0, idx) : sample.getName();

        // File r1FastqFile = readPairList.get(0);
        // String r1FastqRootName =
        // WorkflowUtil.getRootFastqName(r1FastqFile.getName());

        File r2FastqFile = readPairList.get(1);
        String r2FastqRootName = SequencingWorkflowUtil.getRootFastqName(r2FastqFile.getName());

        String fastqLaneRootName = StringUtils.removeEnd(r2FastqRootName, "_R2");

        String irodsDirectory = String.format("/MedGenZone/sequence_data/ncgenes/%s/%s", participantId, version);

        List<CommandInput> commandInputList = new LinkedList<CommandInput>();

        CommandOutput commandOutput = null;

        CommandInput commandInput = new CommandInput();
        commandInput.setExitImmediately(Boolean.FALSE);
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("$IRODS_HOME/imkdir -p %s%n", irodsDirectory));
        sb.append(String.format("$IRODS_HOME/imeta add -C %s Project NCGENES%n", irodsDirectory));
        sb.append(String.format("$IRODS_HOME/imeta add -C %s ParticipantID %s NCGENES%n", irodsDirectory, participantId));
        commandInput.setCommand(sb.toString());
        commandInput.setWorkDir(tmpDir);
        commandInputList.add(commandInput);

        List<IRODSBean> files2RegisterToIRODS = new ArrayList<IRODSBean>();

        File bwaSAMPairedEndOutFile = new File(ncgenesOutputDirectory, fastqLaneRootName + ".sam");

        File fixRGOutput = new File(ncgenesOutputDirectory, bwaSAMPairedEndOutFile.getName().replace(".sam", ".fixed-rg.bam"));
        File picardMarkDuplicatesOutput = new File(ncgenesOutputDirectory, fixRGOutput.getName().replace(".bam", ".deduped.bam"));
        File indelRealignerOut = new File(ncgenesOutputDirectory, picardMarkDuplicatesOutput.getName().replace(".bam", ".realign.bam"));
        File picardFixMateOutput = new File(ncgenesOutputDirectory, indelRealignerOut.getName().replace(".bam", ".fixmate.bam"));
        File gatkTableRecalibrationOut = new File(ncgenesOutputDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.bam"));

        List<ImmutablePair<String, String>> attributeList = Arrays.asList(new ImmutablePair<String, String>("ParticipantId", participantId),
                new ImmutablePair<String, String>("MaPSeqWorkflowVersion", version),
                new ImmutablePair<String, String>("MaPSeqWorkflowName", "NCGenesDX"),
                new ImmutablePair<String, String>("MaPSeqStudyName", sample.getStudy().getName()),
                new ImmutablePair<String, String>("MaPSeqSampleId", sample.getId().toString()),
                new ImmutablePair<String, String>("MaPSeqSystem", system.getValue()),
                new ImmutablePair<String, String>("MaPSeqFlowcellId", sample.getFlowcell().getId().toString()),
                new ImmutablePair<String, String>("DxID", dx), new ImmutablePair<String, String>("DxVersion", version));

        List<ImmutablePair<String, String>> attributeListWithJob = new ArrayList<>(attributeList);
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKDepthOfCoverage.class.getSimpleName()));
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
        files2RegisterToIRODS
                .add(new IRODSBean(
                        new File(outputDirectory,
                                gatkTableRecalibrationOut.getName().replace(".bam",
                                        String.format(".coverage.v%s.gene.sample_cumulative_coverage_counts", version))),
                        attributeListWithJob));

        attributeListWithJob = new ArrayList<>(attributeList);
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKDepthOfCoverage.class.getSimpleName()));
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
        files2RegisterToIRODS.add(new IRODSBean(
                new File(outputDirectory,
                        gatkTableRecalibrationOut.getName().replace(".bam",
                                String.format(".coverage.v%s.gene.sample_cumulative_coverage_proportions", version))),
                attributeListWithJob));

        attributeListWithJob = new ArrayList<>(attributeList);
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKDepthOfCoverage.class.getSimpleName()));
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
        files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, gatkTableRecalibrationOut.getName().replace(".bam",
                String.format(".coverage.v%s.gene.sample_interval_statistics", version))), attributeListWithJob));

        attributeListWithJob = new ArrayList<>(attributeList);
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKDepthOfCoverage.class.getSimpleName()));
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
        files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory,
                gatkTableRecalibrationOut.getName().replace(".bam", String.format(".coverage.v%s.gene.sample_interval_summary", version))),
                attributeListWithJob));

        attributeListWithJob = new ArrayList<>(attributeList);
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKDepthOfCoverage.class.getSimpleName()));
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
        files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory,
                gatkTableRecalibrationOut.getName().replace(".bam", String.format(".coverage.v%s.gene.sample_statistics", version))),
                attributeListWithJob));

        attributeListWithJob = new ArrayList<>(attributeList);
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKDepthOfCoverage.class.getSimpleName()));
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
        files2RegisterToIRODS.add(new IRODSBean(
                new File(outputDirectory,
                        gatkTableRecalibrationOut.getName().replace(".bam", String.format(".coverage.v%s.gene.sample_summary", version))),
                attributeListWithJob));

        attributeListWithJob = new ArrayList<>(attributeList);
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", SAMToolsView.class.getSimpleName()));
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.APPLICATION_BAM.toString()));
        files2RegisterToIRODS.add(new IRODSBean(
                new File(outputDirectory,
                        gatkTableRecalibrationOut.getName().replace(".bam", String.format(".filtered_by_dxid_%s_v%s.bam", dx, version))),
                attributeListWithJob));

        attributeListWithJob = new ArrayList<>(attributeList);
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", PicardSortSAM.class.getSimpleName()));
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.APPLICATION_BAM.toString()));
        files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory,
                gatkTableRecalibrationOut.getName().replace(".bam", String.format(".filtered_by_dxid_%s_v%s.sorted.bam", dx, version))),
                attributeListWithJob));

        attributeListWithJob = new ArrayList<>(attributeList);
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", SAMToolsIndex.class.getSimpleName()));
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.APPLICATION_BAM_INDEX.toString()));
        files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory,
                gatkTableRecalibrationOut.getName().replace(".bam", String.format(".filtered_by_dxid_%s_v%s.sorted.bai", dx, version))),
                attributeListWithJob));

        attributeListWithJob = new ArrayList<>(attributeList);
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", Zip.class.getSimpleName()));
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.APPLICATION_ZIP.toString()));
        files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory,
                gatkTableRecalibrationOut.getName().replace(".bam", String.format(".filtered_by_dxid_%s_v%s.sorted.zip", dx, version))),
                attributeListWithJob));

        attributeListWithJob = new ArrayList<>(attributeList);
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", Zip.class.getSimpleName()));
        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_VCF.toString()));
        files2RegisterToIRODS.add(new IRODSBean(
                new File(outputDirectory,
                        gatkTableRecalibrationOut.getName().replace(".bam", String.format(".filtered_by_dxid_%s_v%s.vcf", dx, version))),
                attributeListWithJob));

        for (IRODSBean bean : files2RegisterToIRODS) {

            File f = bean.getFile();
            if (!f.exists()) {
                logger.warn("file to register doesn't exist: {}", f.getAbsolutePath());
                continue;
            }

            commandInput = new CommandInput();
            commandInput.setExitImmediately(Boolean.FALSE);

            StringBuilder registerCommandSB = new StringBuilder();
            String registrationCommand = String.format("$IRODS_HOME/ireg -f %s %s/%s", f.getAbsolutePath(), irodsDirectory, f.getName());
            String deRegistrationCommand = String.format("$IRODS_HOME/irm -U %s/%s", irodsDirectory, f.getName());
            registerCommandSB.append(registrationCommand).append("\n");
            registerCommandSB.append(String.format("if [ $? != 0 ]; then %s; %s; fi%n", deRegistrationCommand, registrationCommand));
            commandInput.setCommand(registerCommandSB.toString());
            commandInput.setWorkDir(tmpDir);
            commandInputList.add(commandInput);

            commandInput = new CommandInput();
            commandInput.setExitImmediately(Boolean.FALSE);
            sb = new StringBuilder();
            for (ImmutablePair<String, String> attribute : bean.getAttributes()) {
                sb.append(String.format("$IRODS_HOME/imeta add -d %s/%s %s %s GeneScreen%n", irodsDirectory, bean.getFile().getName(),
                        attribute.getLeft(), attribute.getRight()));
            }
            commandInput.setCommand(sb.toString());
            commandInput.setWorkDir(tmpDir);
            commandInputList.add(commandInput);

        }

        File mapseqrc = new File(System.getProperty("user.home"), ".mapseqrc");
        Executor executor = BashExecutor.getInstance();

        for (CommandInput ci : commandInputList) {
            try {
                commandOutput = executor.execute(ci, mapseqrc);
                if (commandOutput.getExitCode() != 0) {
                    logger.info("commandOutput.getExitCode(): {}", commandOutput.getExitCode());
                    logger.warn("command failed: {}", ci.getCommand());
                }
                logger.debug("commandOutput.getStdout(): {}", commandOutput.getStdout());
            } catch (ExecutorException e) {
                if (commandOutput != null) {
                    logger.warn("commandOutput.getStderr(): {}", commandOutput.getStderr());
                }
            }
        }

        logger.info("FINISHED PROCESSING: {}", sample.toString());

    }

    public MaPSeqDAOBeanService getMaPSeqDAOBeanService() {
        return maPSeqDAOBeanService;
    }

    public void setMaPSeqDAOBeanService(MaPSeqDAOBeanService maPSeqDAOBeanService) {
        this.maPSeqDAOBeanService = maPSeqDAOBeanService;
    }

    public Long getSampleId() {
        return sampleId;
    }

    public void setSampleId(Long sampleId) {
        this.sampleId = sampleId;
    }

    public String getVersion() {
        return version;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getDx() {
        return dx;
    }

    public void setDx(String dx) {
        this.dx = dx;
    }

}
