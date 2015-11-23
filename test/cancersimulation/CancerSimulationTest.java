/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cancersimulation;

import java.util.ArrayList;
import java.util.List;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author malcolmcallis
 */
public class CancerSimulationTest {
    
    public CancerSimulationTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    /**
     * Test of ExponentialKernalFitness, of class CancerSimulation.
     * Tests sample containing no mutations
     */
    @Test
    public void testExponentialKernalFitnessNoMutations() throws Exception {
        System.out.println("main");
        boolean[][] trainData = {{false,false,false,false},{true,false,false,false},{true,true,false,false},{true,true,true,false}};
        List<Integer> mutations = new ArrayList<>();
        double[] geneMutationProbabilities = {.1,.2,.3,.4};
        int[] mutationCount = {0,1,2,3};
        double[][] expectedDivergences = {{0,1,2,3},{1,2,3,4},{2,3,4,5},{3,4,5,6}};
        double[][] standardDeviations = {{1,2,3,4},{2,3,4,5},{3,4,5,6},{4,5,6,7}};
        double[] expectedExponentialKernal = {1,2,3,4};
        int genes = 4;
        int samples = 4;
        double fitness = CancerSimulation.exponentialKernalFitness( trainData, mutations, geneMutationProbabilities, mutationCount, expectedDivergences, standardDeviations, expectedExponentialKernal, genes, samples);
        double expected = .616885;
        //check that the error is less than 1%
        //only works for exponential kernal constant = 3
        assertTrue(Math.abs(fitness - expected) / expected < .01);
    }
    
    
    /**
     * Test of ExponentialKernalFitness, of class CancerSimulation.
     * Tests sample containing no mutations
     */
    
    @Test
    public void testExponentialKernalFitnessOneMutations() throws Exception {
        System.out.println("main");
        boolean[][] trainData = {{false,false,false,false},{true,false,false,false},{true,true,false,false},{true,true,true,false}};
        List<Integer> mutations = new ArrayList<>();
        mutations.add(0);
        double[] geneMutationProbabilities = {.1,.2,.3,.4};
        int[] mutationCount = {0,1,2,3};
        double[][] expectedDivergences = {{0,1,2,3},{1,2,3,4},{2,3,4,5},{3,4,5,6}};
        double[][] standardDeviations = {{1,2,3,4},{2,3,4,5},{3,4,5,6},{4,5,6,7}};
        double[] expectedExponentialKernal = {1,2,3,4};
        int genes = 4;
        int samples = 4;
        double fitness = CancerSimulation.exponentialKernalFitness( trainData, mutations, geneMutationProbabilities, mutationCount, expectedDivergences, standardDeviations, expectedExponentialKernal, genes, samples);
        double expected = .668606667;
        
        //check that the error is less than 1%
        //only works for exponential kernal constant = 3
        assertTrue(Math.abs(fitness - expected) / expected < .01);
    }
    
    
    
}
