package edu.unc.mapseq.pipeline.ncgenes.dx;

import java.util.Timer;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class NCGenesDXPipelineExecutorService {

    private final Logger logger = LoggerFactory.getLogger(NCGenesDXPipelineExecutorService.class);

    private final Timer mainTimer = new Timer();

    private NCGenesDXPipelineBeanService pipelineBeanService;

    public void start() throws Exception {
        logger.info("ENTERING stop()");

        long delay = 15 * 1000; // 15 seconds
        long period = 5 * 60 * 1000; // 5 minutes

        NCGenesDXPipelineExecutorTask task = new NCGenesDXPipelineExecutorTask();
        task.setPipelineBeanService(pipelineBeanService);
        mainTimer.scheduleAtFixedRate(task, delay, period);

    }

    public void stop() throws Exception {
        logger.info("ENTERING stop()");
        mainTimer.purge();
        mainTimer.cancel();
    }

    public NCGenesDXPipelineBeanService getPipelineBeanService() {
        return pipelineBeanService;
    }

    public void setPipelineBeanService(NCGenesDXPipelineBeanService pipelineBeanService) {
        this.pipelineBeanService = pipelineBeanService;
    }

}