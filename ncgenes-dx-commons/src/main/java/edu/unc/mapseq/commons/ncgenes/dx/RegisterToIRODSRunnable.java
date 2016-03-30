package edu.unc.mapseq.commons.ncgenes.dx;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.renci.common.exec.BashExecutor;
import org.renci.common.exec.CommandInput;
import org.renci.common.exec.CommandOutput;
import org.renci.common.exec.Executor;
import org.renci.common.exec.ExecutorException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.config.MaPSeqConfigurationService;
import edu.unc.mapseq.config.RunModeType;
import edu.unc.mapseq.dao.MaPSeqDAOBeanService;
import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.Flowcell;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.workflow.impl.IRODSBean;
import edu.unc.mapseq.workflow.impl.SampleWorkflowUtil;

public class RegisterToIRODSRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(RegisterToIRODSRunnable.class);

    private MaPSeqDAOBeanService maPSeqDAOBeanService;

    private MaPSeqConfigurationService maPSeqConfigurationService;

    private Long sampleId;

    private String version;

    private String dx;

    @Override
    public void run() {

        RunModeType runMode = getMaPSeqConfigurationService().getRunMode();

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

        Flowcell flowcell = sample.getFlowcell();

        File ncgenesOutputDirectory = new File(sample.getOutputDirectory(), "NCGenesBaseline");
        File tmpDir = new File(sample.getOutputDirectory(), "tmp");
        if (!tmpDir.exists()) {
            tmpDir.mkdirs();
        }

        File ncgenesDXOutputDirectory = new File(sample.getOutputDirectory(), "NCGenesDX");

        List<File> readPairList = SampleWorkflowUtil.getReadPairList(sample);

        // assumption: a dash is used as a delimiter between a participantId and
        // the external code
        int idx = sample.getName().lastIndexOf("-");
        String participantId = idx != -1 ? sample.getName().substring(0, idx) : sample.getName();

        // File r1FastqFile = readPairList.get(0);
        // String r1FastqRootName =
        // WorkflowUtil.getRootFastqName(r1FastqFile.getName());

        File r2FastqFile = readPairList.get(1);
        String r2FastqRootName = SampleWorkflowUtil.getRootFastqName(r2FastqFile.getName());

        String fastqLaneRootName = StringUtils.removeEnd(r2FastqRootName, "_R2");

        String irodsHome = System.getenv("NCGENES_IRODS_HOME");
        if (StringUtils.isEmpty(irodsHome)) {
            logger.error("irodsHome is not set");
            return;
        }

        String ncgenesIRODSDirectory;

        switch (runMode) {
            case DEV:
            case STAGING:
                ncgenesIRODSDirectory = String.format("/MedGenZone/home/medgenuser/sequence_data/%s/ncgenes/%s/%s",
                        runMode.toString().toLowerCase(), participantId, version);
                break;
            case PROD:
            default:
                ncgenesIRODSDirectory = String.format("/MedGenZone/home/medgenuser/sequence_data/ncgenes/%s/%s", participantId, version);
                break;
        }

        List<CommandInput> commandInputList = new LinkedList<CommandInput>();

        CommandOutput commandOutput = null;

        CommandInput commandInput = new CommandInput();
        commandInput.setExitImmediately(Boolean.FALSE);
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%s/bin/imkdir -p %s%n", irodsHome, ncgenesIRODSDirectory));
        sb.append(String.format("%s/bin/imeta add -C %s Project NCGENES%n", irodsHome, ncgenesIRODSDirectory));
        sb.append(String.format("%s/bin/imeta add -C %s ParticipantID %s NCGENES%n", irodsHome, ncgenesIRODSDirectory, participantId));
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

        files2RegisterToIRODS.add(new IRODSBean(
                new File(ncgenesDXOutputDirectory,
                        gatkTableRecalibrationOut.getName().replace(".bam",
                                String.format(".coverage.v%s.gene.sample_cumulative_coverage_counts", version))),
                "GeneCoverageCount", version, null, runMode));
        files2RegisterToIRODS.add(new IRODSBean(
                new File(ncgenesDXOutputDirectory,
                        gatkTableRecalibrationOut.getName().replace(".bam",
                                String.format(".coverage.v%s.gene.sample_cumulative_coverage_proportions", version))),
                "GeneCoverageProportions", version, null, runMode));
        files2RegisterToIRODS.add(new IRODSBean(
                new File(ncgenesDXOutputDirectory,
                        gatkTableRecalibrationOut.getName().replace(".bam",
                                String.format(".coverage.v%s.gene.sample_interval_statistics", version))),
                "GeneIntervalStatistics", version, null, runMode));
        files2RegisterToIRODS.add(new IRODSBean(
                new File(ncgenesDXOutputDirectory,
                        gatkTableRecalibrationOut.getName().replace(".bam",
                                String.format(".coverage.v%s.gene.sample_interval_summary", version))),
                "GeneIntervalSummary", version, null, runMode));
        files2RegisterToIRODS.add(new IRODSBean(
                new File(ncgenesDXOutputDirectory,
                        gatkTableRecalibrationOut.getName().replace(".bam",
                                String.format(".coverage.v%s.gene.sample_statistics", version))),
                "GeneSampleStatistics", version, null, runMode));
        files2RegisterToIRODS
                .add(new IRODSBean(
                        new File(ncgenesDXOutputDirectory,
                                gatkTableRecalibrationOut.getName().replace(".bam",
                                        String.format(".coverage.v%s.gene.sample_summary", version))),
                        "GeneSampleSummary", version, null, runMode));

        File samtoolsViewOutput = new File(ncgenesDXOutputDirectory,
                gatkTableRecalibrationOut.getName().replace(".bam", String.format(".filtered_by_dxid_%s_v%s.bam", dx, version)));
        File picardSortOutput = new File(ncgenesDXOutputDirectory, samtoolsViewOutput.getName().replace(".bam", ".sorted.bam"));
        File picardSortSAMIndexOut = new File(ncgenesDXOutputDirectory, picardSortOutput.getName().replace(".bam", ".bai"));
        files2RegisterToIRODS.add(new IRODSBean(picardSortSAMIndexOut, "FilteredBamIndex", version, dx, runMode));

        File zipOutputFile = new File(ncgenesDXOutputDirectory, picardSortOutput.getName().replace(".bam", ".zip"));
        files2RegisterToIRODS.add(new IRODSBean(zipOutputFile, "FilteredBamZip", version, dx, runMode));

        File filterVariantOutput = new File(ncgenesDXOutputDirectory,
                gatkTableRecalibrationOut.getName().replace(".bam", String.format(".filtered_by_dxid_%s_v%s.vcf", dx, version)));
        files2RegisterToIRODS.add(new IRODSBean(filterVariantOutput, "FilteredVcf", version, dx, runMode));

        for (IRODSBean bean : files2RegisterToIRODS) {

            File f = bean.getFile();
            if (!f.exists()) {
                logger.warn("file to register doesn't exist: {}", f.getAbsolutePath());
                continue;
            }

            commandInput = new CommandInput();
            commandInput.setExitImmediately(Boolean.FALSE);

            StringBuilder registerCommandSB = new StringBuilder();
            String registrationCommand = String.format("%s/bin/ireg -f %s %s/%s", irodsHome, f.getAbsolutePath(), ncgenesIRODSDirectory,
                    f.getName());
            String deRegistrationCommand = String.format("%s/bin/irm -U %s/%s", irodsHome, ncgenesIRODSDirectory, f.getName());
            registerCommandSB.append(registrationCommand).append("\n");
            registerCommandSB.append(String.format("if [ $? != 0 ]; then %s; %s; fi%n", deRegistrationCommand, registrationCommand));
            commandInput.setCommand(registerCommandSB.toString());
            commandInput.setWorkDir(tmpDir);
            commandInputList.add(commandInput);

            commandInput = new CommandInput();
            commandInput.setExitImmediately(Boolean.FALSE);
            sb = new StringBuilder();
            sb.append(String.format("%s/bin/imeta add -d %s/%s ParticipantID %s NCGENES%n", irodsHome, ncgenesIRODSDirectory, f.getName(),
                    participantId));
            sb.append(String.format("%s/bin/imeta add -d %s/%s FileType %s NCGENES%n", irodsHome, ncgenesIRODSDirectory, f.getName(),
                    bean.getType()));
            sb.append(String.format("%s/bin/imeta add -d %s/%s System %s NCGENES%n", irodsHome, ncgenesIRODSDirectory, f.getName(),
                    StringUtils.capitalize(bean.getRunMode().toString().toLowerCase())));
            commandInput.setCommand(sb.toString());
            commandInput.setWorkDir(tmpDir);
            commandInputList.add(commandInput);

            if (StringUtils.isNotEmpty(bean.getDx())) {
                commandInput = new CommandInput();
                commandInput.setCommand(String.format("%s/bin/imeta add -d %s/%s DxID %s NCGENES", irodsHome, ncgenesIRODSDirectory,
                        f.getName(), bean.getDx()));
                commandInput.setWorkDir(tmpDir);
                commandInputList.add(commandInput);
            }

            if (StringUtils.isNotEmpty(bean.getVersion())) {
                commandInput = new CommandInput();
                commandInput.setCommand(String.format("%s/bin/imeta add -d %s/%s DxVersion %s NCGENES", irodsHome, ncgenesIRODSDirectory,
                        f.getName(), bean.getVersion()));
                commandInput.setWorkDir(tmpDir);
                commandInputList.add(commandInput);
            }

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

    public MaPSeqConfigurationService getMaPSeqConfigurationService() {
        return maPSeqConfigurationService;
    }

    public void setMaPSeqConfigurationService(MaPSeqConfigurationService maPSeqConfigurationService) {
        this.maPSeqConfigurationService = maPSeqConfigurationService;
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
