package edu.unc.mapseq.commands.ncgenes.dx;

import java.util.concurrent.Executors;

import org.apache.karaf.shell.api.action.Action;
import org.apache.karaf.shell.api.action.Argument;
import org.apache.karaf.shell.api.action.Command;
import org.apache.karaf.shell.api.action.lifecycle.Reference;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.commons.ncgenes.dx.RegisterToIRODSMigrationRunnable;
import edu.unc.mapseq.config.MaPSeqConfigurationService;
import edu.unc.mapseq.dao.MaPSeqDAOBeanService;

@Command(scope = "ncgenes-dx", name = "register-to-irods-migration", description = "re-register an entire flowcell")
public class RegisterToIRODSMigrationAction implements Action {

    private final Logger logger = LoggerFactory.getLogger(RegisterToIRODSMigrationAction.class);

    @Reference
    private MaPSeqDAOBeanService maPSeqDAOBeanService;

    @Reference
    private MaPSeqConfigurationService maPSeqConfigurationService;

    @Argument(index = 0, name = "flowcellId", required = true, multiValued = false)
    private Long flowcellId;

    @Override
    public Object execute() throws Exception {
        logger.debug("ENTERING execute()");
        RegisterToIRODSMigrationRunnable runnable = new RegisterToIRODSMigrationRunnable();
        runnable.setMaPSeqDAOBeanService(maPSeqDAOBeanService);
        runnable.setMaPSeqConfigurationService(maPSeqConfigurationService);
        runnable.setFlowcellId(flowcellId);
        Executors.newSingleThreadExecutor().execute(runnable);
        return null;
    }

    public Long getFlowcellId() {
        return flowcellId;
    }

    public void setFlowcellId(Long flowcellId) {
        this.flowcellId = flowcellId;
    }

}
