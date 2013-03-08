package edu.unc.mapseq.messaging.ncgenes.dx;

import java.util.concurrent.Executors;

import javax.jms.JMSException;
import javax.jms.Message;
import javax.jms.MessageListener;
import javax.jms.TextMessage;

import org.apache.commons.lang.StringUtils;
import org.codehaus.jettison.json.JSONException;
import org.codehaus.jettison.json.JSONObject;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.pipeline.ncgenes.dx.NCGenesDXPipelineBeanService;

public class NCGenesDXMessageListener implements MessageListener {

    private final Logger logger = LoggerFactory.getLogger(NCGenesDXMessageListener.class);

    private NCGenesDXPipelineBeanService pipelineBeanService;

    public NCGenesDXMessageListener() {
        super();
    }

    @Override
    public void onMessage(Message message) {
        logger.debug("ENTERING onMessage(Message)");

        String messageValue = null;

        try {
            if (message instanceof TextMessage) {
                logger.debug("received TextMessage");
                TextMessage textMessage = (TextMessage) message;
                messageValue = textMessage.getText();
            }
        } catch (JMSException e2) {
            e2.printStackTrace();
        }

        if (StringUtils.isEmpty(messageValue)) {
            logger.warn("message value is empty");
            return;
        }

        logger.info("messageValue: {}", messageValue);

        JSONObject jsonMessage = null;

        try {
            jsonMessage = new JSONObject(messageValue);
            if (!jsonMessage.has("entities") || !jsonMessage.has("account_name")) {
                logger.error("json lacks entities or account_name");
                return;
            }
        } catch (JSONException e) {
            logger.error("BAD JSON format", e);
            return;
        }

        NCGenesDXMessageRunnable runnable = new NCGenesDXMessageRunnable();
        runnable.setJsonMessage(jsonMessage);
        runnable.setPipelineBeanService(pipelineBeanService);
        Executors.newSingleThreadExecutor().submit(runnable);

    }

    public NCGenesDXPipelineBeanService getPipelineBeanService() {
        return pipelineBeanService;
    }

    public void setPipelineBeanService(NCGenesDXPipelineBeanService pipelineBeanService) {
        this.pipelineBeanService = pipelineBeanService;
    }

}
