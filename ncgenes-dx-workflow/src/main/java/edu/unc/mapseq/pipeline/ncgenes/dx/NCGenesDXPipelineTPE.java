package edu.unc.mapseq.pipeline.ncgenes.dx;

import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

public class NCGenesDXPipelineTPE extends ThreadPoolExecutor {

    public NCGenesDXPipelineTPE() {
        super(20, 20, 5L, TimeUnit.MINUTES, new LinkedBlockingQueue<Runnable>());
    }

}
