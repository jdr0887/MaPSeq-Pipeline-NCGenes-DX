package edu.unc.mapseq.commands.ncgenes.dx;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.karaf.shell.api.action.Action;
import org.apache.karaf.shell.api.action.Command;
import org.apache.karaf.shell.api.action.Option;
import org.apache.karaf.shell.api.action.lifecycle.Reference;
import org.apache.karaf.shell.api.action.lifecycle.Service;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.commons.ncgenes.dx.RegisterToIRODSRunnable;
import edu.unc.mapseq.dao.MaPSeqDAOBeanService;
import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.dao.model.WorkflowRun;

@Command(scope = "ncgenes-dx", name = "register-to-irods", description = "Register a NCGenesDX sample output to iRODS")
@Service
public class RegisterToIRODSAction implements Action {

    private static final Logger logger = LoggerFactory.getLogger(RegisterToIRODSAction.class);

    @Reference
    private MaPSeqDAOBeanService maPSeqDAOBeanService;

    @Option(name = "--sampleId", required = true, multiValued = false)
    private Long sampleId;

    @Option(name = "--workflowRunId", required = true, multiValued = false)
    private Long workflowRunId;

    @Option(name = "--version", required = true, multiValued = false)
    private String version;

    @Option(name = "--dx", required = true, multiValued = false)
    private String dx;

    @Override
    public Object execute() throws Exception {
        logger.debug("ENTERING execute()");
        try {
            ExecutorService es = Executors.newSingleThreadExecutor();
            WorkflowRun workflowRun = maPSeqDAOBeanService.getWorkflowRunDAO().findById(workflowRunId);
            Sample sample = maPSeqDAOBeanService.getSampleDAO().findById(sampleId);
            RegisterToIRODSRunnable runnable = new RegisterToIRODSRunnable(maPSeqDAOBeanService, sample, version, dx, workflowRun);
            es.submit(runnable);
            es.shutdown();
        } catch (MaPSeqDAOException e) {
            logger.error(e.getMessage(), e);
        }
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

    public Long getWorkflowRunId() {
        return workflowRunId;
    }

    public void setWorkflowRunId(Long workflowRunId) {
        this.workflowRunId = workflowRunId;
    }

}
