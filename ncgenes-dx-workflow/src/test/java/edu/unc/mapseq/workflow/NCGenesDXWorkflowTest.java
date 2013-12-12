package edu.unc.mapseq.workflow;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.jgrapht.DirectedGraph;
import org.jgrapht.ext.VertexNameProvider;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.junit.Test;
import org.renci.jlrm.condor.CondorDOTExporter;
import org.renci.jlrm.condor.CondorJob;
import org.renci.jlrm.condor.CondorJobEdge;

import edu.unc.mapseq.module.core.ZipCLI;
import edu.unc.mapseq.module.filter.FilterVariantCLI;
import edu.unc.mapseq.module.gatk.GATKDepthOfCoverageCLI;
import edu.unc.mapseq.module.picard.PicardSortSAMCLI;
import edu.unc.mapseq.module.samtools.SAMToolsIndexCLI;
import edu.unc.mapseq.module.samtools.SAMToolsViewCLI;

public class NCGenesDXWorkflowTest {

    @Test
    public void createDot() {

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(
                CondorJobEdge.class);

        int count = 0;

        // new job
        CondorJob gatkGeneDepthOfCoverageJob = new CondorJob(String.format("%s_%d",
                GATKDepthOfCoverageCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(gatkGeneDepthOfCoverageJob);

        // new job
        CondorJob samtoolsViewJob = new CondorJob(
                String.format("%s_%d", SAMToolsViewCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(samtoolsViewJob);

        // new job
        CondorJob picardSortSAMJob = new CondorJob(String.format("%s_%d", PicardSortSAMCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(picardSortSAMJob);
        graph.addEdge(samtoolsViewJob, picardSortSAMJob);

        // new job
        CondorJob picardSortSAMIndexJob = new CondorJob(String.format("%s_%d", SAMToolsIndexCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(picardSortSAMIndexJob);
        graph.addEdge(picardSortSAMJob, picardSortSAMIndexJob);

        // new job
        CondorJob zipJob = new CondorJob(String.format("%s_%d", ZipCLI.class.getSimpleName(), ++count), null);
        graph.addVertex(zipJob);
        graph.addEdge(picardSortSAMIndexJob, zipJob);

        // new job
        CondorJob filterVariant2Job = new CondorJob(String.format("%s_%d", FilterVariantCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(filterVariant2Job);

        // new job
        CondorJob filterVariant3Job = new CondorJob(String.format("%s_%d", FilterVariantCLI.class.getSimpleName(),
                ++count), null);
        graph.addVertex(filterVariant3Job);

        VertexNameProvider<CondorJob> vnpId = new VertexNameProvider<CondorJob>() {
            @Override
            public String getVertexName(CondorJob job) {
                return job.getName();
            }
        };

        VertexNameProvider<CondorJob> vnpLabel = new VertexNameProvider<CondorJob>() {
            @Override
            public String getVertexName(CondorJob job) {
                return job.getName();
            }
        };

        CondorDOTExporter<CondorJob, CondorJobEdge> dotExporter = new CondorDOTExporter<CondorJob, CondorJobEdge>(
                vnpId, vnpLabel, null, null, null, null);
        File srcSiteResourcesImagesDir = new File("src/site/resources/images");
        if (!srcSiteResourcesImagesDir.exists()) {
            srcSiteResourcesImagesDir.mkdirs();
        }
        File dotFile = new File(srcSiteResourcesImagesDir, "workflow.dag.dot");
        try {
            FileWriter fw = new FileWriter(dotFile);
            dotExporter.export(fw, graph);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}