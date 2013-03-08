package edu.unc.mapseq.messaging;

import javax.jms.Connection;
import javax.jms.DeliveryMode;
import javax.jms.Destination;
import javax.jms.JMSException;
import javax.jms.MessageProducer;
import javax.jms.Session;

import org.apache.activemq.ActiveMQConnectionFactory;
import org.junit.Test;

public class NCGenesDXMessageTest {

    @Test
    public void testQueue() {
        ActiveMQConnectionFactory connectionFactory = new ActiveMQConnectionFactory(String.format("nio://%s:61616",
                "biodev2.its.unc.edu"));

        Connection connection = null;
        Session session = null;
        try {
            connection = connectionFactory.createConnection();
            session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
            Destination destination = session.createQueue("queue/ncgenes.dx");
            MessageProducer producer = session.createProducer(destination);
            producer.setDeliveryMode(DeliveryMode.NON_PERSISTENT);

            String format = "{\"account_name\":\"%s\",\"entities\":[{\"entity_type\":\"HTSFSample\",\"guid\":\"%d\",\"attributes\":[{\"name\":\"GATKDepthOfCoverage.interval_list.version\",\"value\":\"%d\"},{\"name\":\"SAMToolsView.dx.id\",\"value\":\"%d\"}]},{\"entity_type\":\"WorkflowRun\",\"name\":\"%s\"}]}";
            // producer.send(session.createTextMessage(String.format(format, "rc_renci.svc", 113052, 8, 18,
            // "NCG_00151_V8_Dx18")));

            producer.send(session.createTextMessage(String.format(format, "rc_renci.svc", 86851, 8, 9,
                    "NCG_00049_V8_Dx9")));

        } catch (JMSException e) {
            e.printStackTrace();
        } finally {
            try {
                session.close();
                connection.close();
            } catch (JMSException e) {
                e.printStackTrace();
            }
        }

    }

}
