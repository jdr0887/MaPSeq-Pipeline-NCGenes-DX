package edu.unc.mapseq.commands.ncgenes.dx;

import java.util.concurrent.Executors;

import org.apache.karaf.shell.commands.Argument;
import org.apache.karaf.shell.commands.Command;
import org.apache.karaf.shell.console.AbstractAction;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.config.MaPSeqConfigurationService;
import edu.unc.mapseq.dao.MaPSeqDAOBean;

@Command(scope = "ncgenes-dx", name = "register-to-irods-migration", description = "re-register an entire flowcell")
public class RegisterToIRODSMigrationAction extends AbstractAction {

    private final Logger logger = LoggerFactory.getLogger(RegisterToIRODSMigrationAction.class);

    private MaPSeqDAOBean maPSeqDAOBean;

    private MaPSeqConfigurationService maPSeqConfigurationService;

    @Argument(index = 0, name = "flowcellId", required = true, multiValued = false)
    private Long flowcellId;

    @Override
    protected Object doExecute() throws Exception {
        logger.debug("ENTERING doExecute()");
        RegisterToIRODSMigrationRunnable runnable = new RegisterToIRODSMigrationRunnable();
        runnable.setMaPSeqDAOBean(maPSeqDAOBean);
        runnable.setMaPSeqConfigurationService(maPSeqConfigurationService);
        runnable.setFlowcellId(flowcellId);
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

    public Long getFlowcellId() {
        return flowcellId;
    }

    public void setFlowcellId(Long flowcellId) {
        this.flowcellId = flowcellId;
    }

}
