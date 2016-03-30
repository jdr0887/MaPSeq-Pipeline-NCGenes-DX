package edu.unc.mapseq.commands.ncgenes.dx;

import java.util.concurrent.Executors;

import org.apache.karaf.shell.api.action.Action;
import org.apache.karaf.shell.api.action.Argument;
import org.apache.karaf.shell.api.action.Command;
import org.apache.karaf.shell.api.action.lifecycle.Reference;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.commons.ncgenes.dx.RegisterToIRODSRunnable;
import edu.unc.mapseq.config.MaPSeqConfigurationService;
import edu.unc.mapseq.dao.MaPSeqDAOBeanService;

@Command(scope = "ncgenes-dx", name = "register-to-irods", description = "Register a NCGenesDX sample output to iRODS")
public class RegisterToIRODSAction implements Action {

    private final Logger logger = LoggerFactory.getLogger(RegisterToIRODSAction.class);

    @Reference
    private MaPSeqDAOBeanService maPSeqDAOBeanService;

    @Reference
    private MaPSeqConfigurationService maPSeqConfigurationService;

    @Argument(index = 0, name = "sampleId", required = true, multiValued = false)
    private Long sampleId;

    @Argument(index = 1, name = "version", required = true, multiValued = false)
    private String version;

    @Argument(index = 2, name = "dx", required = true, multiValued = false)
    private String dx;

    @Override
    public Object execute() throws Exception {
        logger.debug("ENTERING execute()");
        RegisterToIRODSRunnable runnable = new RegisterToIRODSRunnable();
        runnable.setMaPSeqDAOBeanService(maPSeqDAOBeanService);
        runnable.setMaPSeqConfigurationService(maPSeqConfigurationService);
        runnable.setSampleId(sampleId);
        runnable.setDx(dx);
        runnable.setVersion(version);
        Executors.newSingleThreadExecutor().execute(runnable);
        return null;
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
