package edu.unc.mapseq.commands.ncgenes.dx;

import java.util.concurrent.Executors;

import org.apache.karaf.shell.commands.Argument;
import org.apache.karaf.shell.commands.Command;
import org.apache.karaf.shell.console.AbstractAction;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.commons.ncgenes.dx.RegisterToIRODSRunnable;
import edu.unc.mapseq.config.MaPSeqConfigurationService;
import edu.unc.mapseq.dao.MaPSeqDAOBean;

@Command(scope = "ncgenes-dx", name = "register-to-irods", description = "Register a NCGenesDX sample output to iRODS")
public class RegisterToIRODSAction extends AbstractAction {

    private final Logger logger = LoggerFactory.getLogger(RegisterToIRODSAction.class);

    private MaPSeqDAOBean maPSeqDAOBean;

    private MaPSeqConfigurationService maPSeqConfigurationService;

    @Argument(index = 0, name = "sampleId", required = true, multiValued = false)
    private Long sampleId;

    @Argument(index = 1, name = "version", required = true, multiValued = false)
    private String version;

    @Argument(index = 2, name = "dx", required = true, multiValued = false)
    private String dx;

    @Override
    protected Object doExecute() throws Exception {
        logger.debug("ENTERING doExecute()");
        RegisterToIRODSRunnable runnable = new RegisterToIRODSRunnable();
        runnable.setMaPSeqDAOBean(maPSeqDAOBean);
        runnable.setMaPSeqConfigurationService(maPSeqConfigurationService);
        runnable.setSampleId(sampleId);
        runnable.setDx(dx);
        runnable.setVersion(version);
        Executors.newSingleThreadExecutor().execute(runnable);
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
