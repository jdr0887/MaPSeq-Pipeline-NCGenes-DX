package edu.unc.mapseq.commands.ncgenes.dx;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.apache.felix.gogo.commands.Argument;
import org.apache.felix.gogo.commands.Command;
import org.apache.karaf.shell.console.AbstractAction;
import org.renci.common.exec.BashExecutor;
import org.renci.common.exec.CommandInput;
import org.renci.common.exec.CommandOutput;
import org.renci.common.exec.Executor;
import org.renci.common.exec.ExecutorException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.config.MaPSeqConfigurationService;
import edu.unc.mapseq.config.RunModeType;
import edu.unc.mapseq.dao.MaPSeqDAOBean;
import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.HTSFSample;
import edu.unc.mapseq.dao.model.SequencerRun;
import edu.unc.mapseq.workflow.IRODSBean;
import edu.unc.mapseq.workflow.WorkflowUtil;

@Command(scope = "ncgenes-dx", name = "register-to-irods", description = "Register a NCGenesDX sample output to iRODS")
public class RegisterNCGenesDXToIRODSAction extends AbstractAction {

    private final Logger logger = LoggerFactory.getLogger(RegisterNCGenesDXToIRODSAction.class);

    private MaPSeqDAOBean maPSeqDAOBean;

    private MaPSeqConfigurationService maPSeqConfigurationService;

    @Argument(index = 0, name = "htsfSampleId", required = true, multiValued = false)
    private Long htsfSampleId;

    @Argument(index = 1, name = "version", required = true, multiValued = false)
    private String version;

    @Argument(index = 2, name = "dx", required = true, multiValued = false)
    private String dx;

    @Override
    protected Object doExecute() throws Exception {

        RunModeType runMode = getMaPSeqConfigurationService().getRunMode();

        HTSFSample htsfSample;
        try {
            htsfSample = maPSeqDAOBean.getHTSFSampleDAO().findById(htsfSampleId);
        } catch (MaPSeqDAOException e1) {
            e1.printStackTrace();
            return null;
        }

        File mapseqOutputDirectory = new File(System.getenv("MAPSEQ_OUTPUT_DIRECTORY"));
        File baseDir;
        switch (runMode) {
            case DEV:
            case STAGING:
                baseDir = new File(mapseqOutputDirectory, runMode.toString().toLowerCase());
                break;
            case PROD:
            default:
                baseDir = mapseqOutputDirectory;
                break;
        }

        SequencerRun sequencerRun = htsfSample.getSequencerRun();

        File sequencerRunOutputDirectory = new File(baseDir, sequencerRun.getName());
        File workflowDir = new File(sequencerRunOutputDirectory, "NCGenes");
        File outputDirectory = new File(workflowDir, htsfSample.getName());
        File tmpDir = new File(outputDirectory, "tmp");
        if (!tmpDir.exists()) {
            tmpDir.mkdirs();
        }

        List<File> readPairList = WorkflowUtil.getReadPairList(htsfSample.getFileDatas(), sequencerRun.getName(),
                htsfSample.getLaneIndex());

        // assumption: a dash is used as a delimiter between a participantId and
        // the external code
        int idx = htsfSample.getName().lastIndexOf("-");
        String participantId = idx != -1 ? htsfSample.getName().substring(0, idx) : htsfSample.getName();

        // File r1FastqFile = readPairList.get(0);
        // String r1FastqRootName =
        // WorkflowUtil.getRootFastqName(r1FastqFile.getName());

        File r2FastqFile = readPairList.get(1);
        String r2FastqRootName = WorkflowUtil.getRootFastqName(r2FastqFile.getName());

        String fastqLaneRootName = StringUtils.removeEnd(r2FastqRootName, "_R2");

        String irodsHome = System.getenv("NCGENES_IRODS_HOME");
        if (StringUtils.isEmpty(irodsHome)) {
            logger.error("irodsHome is not set");
            return null;
        }

        String ncgenesIRODSDirectory;

        switch (runMode) {
            case DEV:
            case STAGING:
                ncgenesIRODSDirectory = String.format("/genomicsDataGridZone/sequence_data/%s/ncgenes/%s/%s", runMode
                        .toString().toLowerCase(), participantId, version);
                break;
            case PROD:
            default:
                ncgenesIRODSDirectory = String.format("/genomicsDataGridZone/sequence_data/ncgenes/%s/%s",
                        participantId, version);
                break;
        }

        List<CommandInput> commandInputList = new LinkedList<CommandInput>();

        CommandOutput commandOutput = null;

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

        File bwaSAMPairedEndOutFile = new File(outputDirectory, fastqLaneRootName + ".sam");

        File fixRGOutput = new File(outputDirectory, bwaSAMPairedEndOutFile.getName().replace(".sam", ".fixed-rg.bam"));
        File picardMarkDuplicatesOutput = new File(outputDirectory, fixRGOutput.getName().replace(".bam",
                ".deduped.bam"));
        File indelRealignerOut = new File(outputDirectory, picardMarkDuplicatesOutput.getName().replace(".bam",
                ".realign.bam"));
        File picardFixMateOutput = new File(outputDirectory, indelRealignerOut.getName()
                .replace(".bam", ".fixmate.bam"));
        File gatkTableRecalibrationOut = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam",
                ".recal.bam"));

        files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, gatkTableRecalibrationOut.getName().replace(
                ".bam", String.format(".coverage.v%s.gene.sample_cumulative_coverage_counts", version))),
                "GeneCoverageCount", version, null, runMode));
        files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, gatkTableRecalibrationOut.getName().replace(
                ".bam", String.format(".coverage.v%s.gene.sample_cumulative_coverage_proportions", version))),
                "GeneCoverageProportions", version, null, runMode));
        files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, gatkTableRecalibrationOut.getName().replace(
                ".bam", String.format(".coverage.v%s.gene.sample_interval_statistics", version))),
                "GeneIntervalStatistics", version, null, runMode));
        files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, gatkTableRecalibrationOut.getName().replace(
                ".bam", String.format(".coverage.v%s.gene.sample_interval_summary", version))), "GeneIntervalSummary",
                version, null, runMode));
        files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, gatkTableRecalibrationOut.getName().replace(
                ".bam", String.format(".coverage.v%s.gene.sample_statistics", version))), "GeneSampleStatistics",
                version, null, runMode));
        files2RegisterToIRODS.add(new IRODSBean(new File(outputDirectory, gatkTableRecalibrationOut.getName().replace(
                ".bam", String.format(".coverage.v%s.gene.sample_summary", version))), "GeneSampleSummary", version,
                null, runMode));

        File samtoolsViewOutput = new File(outputDirectory, gatkTableRecalibrationOut.getName().replace(".bam",
                ".filtered.bam"));
        File picardSortOutput = new File(outputDirectory, samtoolsViewOutput.getName().replace(".bam",
                String.format(".sorted.filtered_by_dxid_%s_v%s.bam", dx, version)));
        File picardSortSAMIndexOut = new File(outputDirectory, picardSortOutput.getName().replace(".bam", ".bai"));
        files2RegisterToIRODS.add(new IRODSBean(picardSortSAMIndexOut, "FilteredBamIndex", version, dx, runMode));

        File zipOutputFile = new File(outputDirectory, picardSortOutput.getName().replace(".bam", ".zip"));
        files2RegisterToIRODS.add(new IRODSBean(zipOutputFile, "FilteredBamZip", version, dx, runMode));

        File filterVariantOutput = new File(outputDirectory, gatkTableRecalibrationOut.getName().replace(".bam",
                String.format(".filtered_by_dxid_%s_v%s.vcf", dx, version)));
        files2RegisterToIRODS.add(new IRODSBean(filterVariantOutput, "FilteredVcf", version, dx, runMode));

        for (IRODSBean bean : files2RegisterToIRODS) {

            commandInput = new CommandInput();
            commandInput.setExitImmediately(Boolean.FALSE);

            StringBuilder registerCommandSB = new StringBuilder();
            String registrationCommand = String.format("%s/bin/ireg -f %s %s/%s", irodsHome, bean.getFile()
                    .getAbsolutePath(), ncgenesIRODSDirectory, bean.getFile().getName());
            String deRegistrationCommand = String.format("%s/bin/irm -U %s/%s", irodsHome, ncgenesIRODSDirectory, bean
                    .getFile().getName());
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

        System.out.println("FINISHED PROCESSING: " + htsfSample.toString());
        return null;

    }

    public MaPSeqDAOBean getMaPSeqDAOBean() {
        return maPSeqDAOBean;
    }

    public void setMaPSeqDAOBean(MaPSeqDAOBean maPSeqDAOBean) {
        this.maPSeqDAOBean = maPSeqDAOBean;
    }

    public MaPSeqConfigurationService getMaPSeqConfigurationService() {
        return maPSeqConfigurationService;
    }

    public void setMaPSeqConfigurationService(MaPSeqConfigurationService maPSeqConfigurationService) {
        this.maPSeqConfigurationService = maPSeqConfigurationService;
    }

    public Long getHtsfSampleId() {
        return htsfSampleId;
    }

    public void setHtsfSampleId(Long htsfSampleId) {
        this.htsfSampleId = htsfSampleId;
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
