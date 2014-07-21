package edu.unc.mapseq.executor.ncgenes.dx;

import java.util.Date;
import java.util.List;
import java.util.TimerTask;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.WorkflowDAO;
import edu.unc.mapseq.dao.WorkflowPlanDAO;
import edu.unc.mapseq.dao.WorkflowRunDAO;
import edu.unc.mapseq.dao.model.Workflow;
import edu.unc.mapseq.dao.model.WorkflowPlan;
import edu.unc.mapseq.dao.model.WorkflowRun;
import edu.unc.mapseq.workflow.WorkflowBeanService;
import edu.unc.mapseq.workflow.WorkflowExecutor;
import edu.unc.mapseq.workflow.WorkflowTPE;
import edu.unc.mapseq.workflow.ncgenes.dx.NCGenesDXWorkflow;

public class NCGenesDXWorkflowExecutorTask extends TimerTask {

    private final Logger logger = LoggerFactory.getLogger(NCGenesDXWorkflowExecutorTask.class);

    private final WorkflowTPE threadPoolExecutor = new WorkflowTPE();

    private WorkflowBeanService workflowBeanService;

    public NCGenesDXWorkflowExecutorTask() {
        super();
    }

    @Override
    public void run() {
        logger.info("ENTERING run()");

        threadPoolExecutor.setCorePoolSize(workflowBeanService.getCorePoolSize());
        threadPoolExecutor.setMaximumPoolSize(workflowBeanService.getMaxPoolSize());

        logger.info(String.format("ActiveCount: %d, TaskCount: %d, CompletedTaskCount: %d",
                threadPoolExecutor.getActiveCount(), threadPoolExecutor.getTaskCount(),
                threadPoolExecutor.getCompletedTaskCount()));

        WorkflowDAO workflowDAO = this.workflowBeanService.getMaPSeqDAOBean().getWorkflowDAO();
        WorkflowRunDAO workflowRunDAO = this.workflowBeanService.getMaPSeqDAOBean().getWorkflowRunDAO();
        WorkflowPlanDAO workflowPlanDAO = this.workflowBeanService.getMaPSeqDAOBean().getWorkflowPlanDAO();

        try {
            List<Workflow> workflowList = workflowDAO.findByName("NCGenesDX");
            if (workflowList == null || (workflowList != null && workflowList.isEmpty())) {
                logger.error("No Workflow Found: {}", "NCGenesDX");
                return;
            }
            Workflow workflow = workflowList.get(0);
            List<WorkflowPlan> workflowPlanList = workflowPlanDAO.findEnqueued(workflow.getId());

            if (workflowPlanList != null && workflowPlanList.size() > 0) {

                logger.info("dequeuing {} WorkflowPlans", workflowPlanList.size());
                for (WorkflowPlan workflowPlan : workflowPlanList) {

                    NCGenesDXWorkflow ncGenesDXWorkflow = new NCGenesDXWorkflow();

                    WorkflowRun workflowRun = workflowPlan.getWorkflowRun();
                    workflowRun.setVersion(ncGenesDXWorkflow.getVersion());
                    workflowRun.setDequeuedDate(new Date());
                    workflowRunDAO.save(workflowRun);

                    ncGenesDXWorkflow.setWorkflowBeanService(workflowBeanService);
                    ncGenesDXWorkflow.setWorkflowPlan(workflowPlan);
                    threadPoolExecutor.submit(new WorkflowExecutor(ncGenesDXWorkflow));

                }

            }

        } catch (MaPSeqDAOException e) {
            e.printStackTrace();
        }

    }

    public WorkflowBeanService getWorkflowBeanService() {
        return workflowBeanService;
    }

    public void setWorkflowBeanService(WorkflowBeanService workflowBeanService) {
        this.workflowBeanService = workflowBeanService;
    }

}
