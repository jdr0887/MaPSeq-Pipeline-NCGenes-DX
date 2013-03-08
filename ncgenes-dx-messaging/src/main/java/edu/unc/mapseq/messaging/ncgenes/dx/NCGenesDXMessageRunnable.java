package edu.unc.mapseq.messaging.ncgenes.dx;

import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.Executors;

import org.codehaus.jettison.json.JSONArray;
import org.codehaus.jettison.json.JSONException;
import org.codehaus.jettison.json.JSONObject;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.Account;
import edu.unc.mapseq.dao.model.HTSFSample;
import edu.unc.mapseq.dao.model.SequencerRun;
import edu.unc.mapseq.dao.model.WorkflowPlan;
import edu.unc.mapseq.dao.model.WorkflowRun;
import edu.unc.mapseq.dao.model.WorkflowRunStatusType;
import edu.unc.mapseq.messaging.common.EntityUtil;
import edu.unc.mapseq.pipeline.PipelineExecutor;
import edu.unc.mapseq.pipeline.ncgenes.dx.NCGenesDXPipeline;
import edu.unc.mapseq.pipeline.ncgenes.dx.NCGenesDXPipelineBeanService;

public class NCGenesDXMessageRunnable implements Runnable {

    private final Logger logger = LoggerFactory.getLogger(NCGenesDXMessageListener.class);

    private NCGenesDXPipelineBeanService pipelineBeanService;

    private JSONObject jsonMessage;

    public NCGenesDXMessageRunnable() {
        super();
    }

    @Override
    public void run() {
        logger.debug("ENTERING onMessage(Message)");

        SequencerRun sequencerRun = null;
        Set<HTSFSample> htsfSampleSet = new HashSet<HTSFSample>();
        WorkflowRun workflowRun = null;
        Account account = null;
        NCGenesDXPipeline pipeline = new NCGenesDXPipeline();
        pipeline.setPipelineBeanService(pipelineBeanService);

        try {

            String accountName = jsonMessage.getString("account_name");

            try {
                account = pipelineBeanService.getMaPSeqDAOBean().getAccountDAO().findByName(accountName);
            } catch (MaPSeqDAOException e) {
            }

            if (account == null) {
                logger.error("Must register account first");
                return;
            }

            JSONArray entityArray = jsonMessage.getJSONArray("entities");

            for (int i = 0; i < entityArray.length(); ++i) {

                JSONObject entityJSONObject = entityArray.getJSONObject(i);

                if (entityJSONObject.has("entity_type")) {

                    String entityType = entityJSONObject.getString("entity_type");

                    if ("Sequencer run".equals(entityType) || SequencerRun.class.getSimpleName().equals(entityType)) {
                        sequencerRun = EntityUtil.getSequencerRun(pipelineBeanService.getMaPSeqDAOBean(),
                                entityJSONObject);
                    }

                    if ("HTSF Sample".equals(entityType) || HTSFSample.class.getSimpleName().equals(entityType)) {
                        HTSFSample htsfSample = EntityUtil.getHTSFSample(pipelineBeanService.getMaPSeqDAOBean(),
                                entityJSONObject);
                        htsfSampleSet.add(htsfSample);
                    }

                    if ("Workflow run".equals(entityType) || WorkflowRun.class.getSimpleName().equals(entityType)) {
                        workflowRun = EntityUtil.getWorkflowRun(pipeline, entityJSONObject, account);
                    }

                }

            }
        } catch (JSONException e1) {
            e1.printStackTrace();
            return;
        }

        if (workflowRun == null) {
            logger.warn("Invalid JSON...not running anything");
            return;
        }

        if (sequencerRun == null && htsfSampleSet.size() == 0) {
            logger.warn("Invalid JSON...not running anything");
            workflowRun.setStatus(WorkflowRunStatusType.FAILED);
        }

        try {
            Long workflowRunId = pipelineBeanService.getMaPSeqDAOBean().getWorkflowRunDAO().save(workflowRun);
            workflowRun.setId(workflowRunId);
        } catch (MaPSeqDAOException e) {
            e.printStackTrace();
        }

        try {
            WorkflowPlan workflowPlan = new WorkflowPlan();
            workflowPlan.setWorkflowRun(workflowRun);
            if (htsfSampleSet.size() > 0) {
                workflowPlan.setHTSFSamples(htsfSampleSet);
            }
            if (sequencerRun != null) {
                workflowPlan.setSequencerRun(sequencerRun);
            }
            Long workflowPlanId = pipelineBeanService.getMaPSeqDAOBean().getWorkflowPlanDAO().save(workflowPlan);
            workflowPlan.setId(workflowPlanId);

            pipeline.setWorkflowPlan(workflowPlan);

        } catch (MaPSeqDAOException e) {
            e.printStackTrace();
        }

        if (workflowRun.getStatus() == WorkflowRunStatusType.FAILED) {
            return;
        }

        Executors.newSingleThreadExecutor().execute(new PipelineExecutor(pipeline));

    }

    public NCGenesDXPipelineBeanService getPipelineBeanService() {
        return pipelineBeanService;
    }

    public void setPipelineBeanService(NCGenesDXPipelineBeanService pipelineBeanService) {
        this.pipelineBeanService = pipelineBeanService;
    }

    public JSONObject getJsonMessage() {
        return jsonMessage;
    }

    public void setJsonMessage(JSONObject jsonMessage) {
        this.jsonMessage = jsonMessage;
    }

}
