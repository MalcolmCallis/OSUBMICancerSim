package cancersimulation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
/**
 * Simulates the accumulation of mutations in cancer cells March, 2014
 *
 * @author Malcolm Callis
 * 
 * GUI calls run function during normal operation
 * 
 */



public final class CancerSimulation{

    /**
     * Default constructor--private to prevent instantiation.
     */
    private CancerSimulation() {
        // no code needed here
    }
    
    //private static final GUI window = new GUI();
    private static boolean stop = false;
    private static boolean inverse = false;
    private static boolean subtractBest = true;
    private static int drugRegime = 0;
    private static boolean regimeChange = false;
    //private static double alpha = .06;
    //private static double epsilon = 1;
    private static double rhoConstant = 1;
    private static int startTime = (int) System.currentTimeMillis();

    private static double nextGaussian(){
        double v1 = 0;
        double v2 = 0;
        double s = 1;
        while (s >=1 ) {
            v1 = 2 * Math.random() - 1;
            v2 = 2 * Math.random() - 1;
            s = v1 * v1 + v2 * v2;
        }
        if (s != 0) {
            s = Math.sqrt(-2 * Math.log(s) / s);
        }
        return v1 * s;
    }
    
    private static long nextPoisson(double meanPoisson) {
            double p = Math.exp(-meanPoisson);
            long n = 0;
            double r = 1.0d;
            double rnd = 1.0d;

            while (n < 1000 * meanPoisson) {
                rnd = Math.random();
                r *= rnd;
                if (r >= p) {
                    n++;
                } else {
                    return n;
                }
            }
            return n;
        }
  
    
     
    private static double[] mScore(boolean sampleData[][], double[] geneMutationProbabilities, boolean[] genome, int mutationCount, int[] sampleDataMutationCounts, int genes, int samples){
      
        double[] mScores = new double[samples];
        
        for (int sample = 0; sample < samples; sample++){
            boolean[] sampleDataGenome = sampleData[sample];
            mScores[sample] = mScore(genome, mutationCount, sampleDataGenome, sampleDataMutationCounts[sample], geneMutationProbabilities, genes);
        }
        
        
        return mScores;
    }
        
    private static double mScore(boolean[] genome1, int mutationCount1, boolean[] genome2, int mutationCount2, double[] geneMutationProbabilities, int genes){
        /*for each gene in each sample, calculate the probability (p for sample1, q for sample2) that it will be mutated
        based on the number of mutations in the sample and the gene lengths
        four possible states:
        sample 1 and 2 are mutated
        sample 1 and 2 are not mutated
        sample 1 is mutated and sample 2 is not mutated
        sample 1 is not mutated and sample 2 is mutated
        for whatever outcome, calculate the probability it would be that outcome based on p and q
        mscore = 1 / (p * q)
        mscore is positive if they both match and negative if they both don't match
        total m score is the sum of all mscores for the sample
        the expected mscore for two random samples is zero
        samples that are more similar will be positive
        samples that are less similar will be negative
       
        
        */
            double totalMScore = 0;
            double totalBestMScore = 0;
            for (int gene = 0; gene < genes; gene++){
                //find the m score for each gene
                double sample1OutcomeProbability = 1 - Math.pow( 1 - geneMutationProbabilities[gene], mutationCount1);
                if (!genome1[gene]){
                    //use 1-p if it is not mutated
                    sample1OutcomeProbability = 1 - sample1OutcomeProbability;
                }
                double sample2OutcomeProbability = 1 - Math.pow(1 - geneMutationProbabilities[gene], mutationCount2);

                if (!genome2[gene]){
                    //use 1-q if it is not mutated
                    sample2OutcomeProbability = 1 - sample2OutcomeProbability;
                }

                double mScore;
                if (inverse){
                    mScore = 1.0 / (sample1OutcomeProbability * sample2OutcomeProbability);
                } else {
                    mScore = Math.log(sample1OutcomeProbability * sample2OutcomeProbability);
                }
                if (genome1[gene] != genome2[gene]){
                    //subtract from the mScore if the mutations do not match up
                    mScore *= -1;
                }
                //the best mscore is the average of the mscore of sample1 to sample1 and the mscore of sample2 to sample2
                //this is because the closest distance is two samples that match up perfectly
                double bestMScore = (Math.log(sample1OutcomeProbability * sample1OutcomeProbability) + Math.log(sample2OutcomeProbability * sample2OutcomeProbability)) / 2;
                totalMScore += mScore;
                totalBestMScore += bestMScore;
            }
            
          
        if (subtractBest){
            return totalMScore - totalBestMScore;
        } else {
            return totalMScore;
        }
    }
    
    private static double[] geneMutationProbabilities(String path, int genes, int otherGenes) throws IOException{
        double[] probabilities = new double[genes];
        //contents of index i is the probability that a mutation will occur in gene i
        //probability of being mutated = genelength / totalgenelength
        //the sum of all probabilities in the array is 1
        List<String> lines = Files.readAllLines(Paths.get(path), StandardCharsets.UTF_8);
        String[] stringLines = new String[lines.size()];
        lines.toArray(stringLines);
        int total = otherGenes;
        for (int gene = 0; gene < genes; gene++){
            probabilities[gene] = Double.parseDouble(stringLines[gene]);
            total += probabilities[gene];
        }
        for (int gene = 0; gene < genes; gene++){
            probabilities[gene] /= total;
        }
        
        return probabilities;
    }
    
    private static double mScoreDivisionRate(double zScores[], int samples, double height, double radius, double nullFitness){
        //calculates a division rate based on a given set of adjusted divergence scores
        double energy = 0;
        double ave = 0;
        double asdf[] = new double[samples];
        for (int sample = 0; sample < samples; sample++){
            ave += zScores[sample];
            //totalScore += Math.log(mScores[sample]);
            double temp = Math.exp( - Math.pow(zScores[sample] / radius, 2) );
            energy += temp;
            asdf[sample] = temp;
            //energy += Math.exp(mScores[sample] / bestMScores[sample]) - 1;
        }
        energy /= samples;
        double fitness = Math.log(energy);
        ave /= samples;
        ave += 0;
        //double ans = .5 + height * (fitness + temp);
        double ans = .5 + height * (Math.log (( energy + .001) / (nullFitness + .001)));
        return ans;
        //return .5 + (energy - samples * Math.exp(-1)) * height;
    }

    private static double nullNormalKernalFitness(boolean trainData[][], double geneMutationProbabilities[], int mutationCount[], int genes, int samples, double radius, boolean inverse, boolean subtractBest){
        
        double nullMScores[] = mScore(trainData, geneMutationProbabilities, new boolean[genes],0, mutationCount, genes, samples);
        
        double energy = 0;
        for (int sample = 0; sample < samples; sample++){
            //totalScore += Math.log(mScores[sample]);
            double temp = Math.exp( - Math.pow(nullMScores[sample] / radius, 2) );
            energy += temp;
            //energy += Math.exp(mScores[sample] / bestMScores[sample]) - 1;
        }
        energy /= samples;
        return energy;
        
    }

    public static double exponentialKernalFitness(boolean trainData[][], boolean[] genome, int mutationCount, double geneMutationProbabilities[], int sampleDataMutationCount[], double expectedDivergences[][], double standardDeviations[][], double expectedExponentialKernal[], int genes, int samples){
        double mScores[] = mScore(trainData, geneMutationProbabilities, genome, mutationCount, sampleDataMutationCount, genes, samples);
        double exponentialKernal;
        if (inverse){
            exponentialKernal = exponentialKernal(mScores, samples);
        } else {
            
            double zScores[] = zScore( mScores, mutationCount, expectedDivergences, standardDeviations, sampleDataMutationCount, samples);
            exponentialKernal = exponentialKernal(zScores, samples);
        }
        
        return exponentialKernal;
    }
  
    public static double observedExpectedRatioFitness(boolean trainData[][], boolean simulatedSampleData[][], boolean[] genome, double[] rho, int mutationCount, double geneMutationProbabilities[], int sampleDataMutationCount[], int genes, int samples, double alpha, double epsilon){
        double mScores[] = mScore(trainData, geneMutationProbabilities, genome, mutationCount, sampleDataMutationCount, genes, samples);
        double expectedMScores[] = mScore(simulatedSampleData, geneMutationProbabilities, genome, mutationCount, sampleDataMutationCount, genes, samples);
        double ret;

        double observed = observed(rho, mScores, mutationCount, sampleDataMutationCount, samples, alpha);
        double expected = observed(rho, expectedMScores, mutationCount, sampleDataMutationCount, samples, alpha);
        ret = (observed + epsilon) / (expected + epsilon);
        
        return ret;
    }
    
    
    public static double similarityTest(boolean testData[][], List<List<Integer>> simulatedClones){
        //average of the smallest distance between each sample int the test data and the simulated clones
        double[] smallestDistances = new double[testData.length];
        for (int sample = 0; sample < testData.length; sample++){
            double distances[] = new double[simulatedClones.size()];
            double pObs = 0;
            for (boolean b : testData[sample]){
                if (b){
                    pObs++;
                }
            }
            pObs /= testData[0].length;
            for (int clone = 0; clone < simulatedClones.size(); clone++){
                double expected = 0;
                double observed = 0;
                for (int gene = 0; gene < testData[0].length; gene++){
                    double pExp = simulatedClones.get(clone).size() / testData[0].length;
                    expected += 1 - (pExp * pObs + (1 - pExp) * (1 - pObs));
                    observed += testData[sample][gene] == simulatedClones.get(clone).contains(gene)? 0 : 1;
                    distances[clone] = Math.min(expected / observed, 2);
                }
            }
            //find the smallest distance
            double smallest = Double.MAX_VALUE;
            for (double x : distances){
                if (x < smallest)
                smallest = x;
            }
            smallestDistances[sample] = smallest;
        }
        //average the distances
        double average = 0;
        for (double x : smallestDistances){
            average += x;
        }
        average /= smallestDistances.length;
        
        return average;
    }
   
    
    
    
    private static double energy(double[] distances, boolean quadratic, double radius, double height){
        //calculates the total energy based on a normal kernal density approximation
        double energy = 0;
        for (int sample = 0; sample < distances.length; sample++){
            //normal kernal distribution = e^(-x^2/radius^2)
            energy += Math.pow(Math.E,- Math.pow(distances[sample]/radius,quadratic? 2 : 1));
        }
        //mutltiply energy by constant factor
        return energy*height;
    }
    
    private static double observed(double[] rho, double divergenceScores[], int mutationCount, int SampleDataMutationCount[], int samples, double alpha){
        double ret = 0;
        for (int sample = 0; sample < samples; sample++){
            ret += rho[mutationCount] / rho[SampleDataMutationCount[sample]] * Math.exp( -alpha * divergenceScores[sample]);
            //ret += Math.exp( -alpha * divergenceScores[sample]);
        }
        return ret;
    }
    
    public static boolean[][] readInputMutationsFile(String path) throws IOException{
        //reads input file

        List<String> lines = Files.readAllLines(Paths.get(path), StandardCharsets.UTF_8);
        

        boolean sampleData[][] = new boolean[lines.size()][]; 
        double total = 0;

        for(int row = 0; row < lines.size(); row++){
            String[] genesString;
            if (path.endsWith(".csv")){
            genesString = lines.get(row).split(",");
            } else{
            genesString = lines.get(row).split("\t");
            }
            boolean[] genesBoolean = new boolean[genesString.length];
            for (int column = 0; column < genesString.length; column++){

                boolean isMutated = genesString[column].compareTo("1") == 0;
                if (isMutated){
                    total++;
                }
                genesBoolean[column] = isMutated;
            }
            sampleData[row] = genesBoolean;
        }
        return sampleData;
    }
    
    private static String[] readInputFile(String path) throws IOException{


        List<String> lines = Files.readAllLines(Paths.get(path), StandardCharsets.UTF_8);
        String[] ret = new String[lines.size()];
        lines.toArray(ret);
        return ret;
        
    }
    
    private static double[][] readInput2dcsv(String path, int genes) throws IOException{
        

        List<String> lines = Files.readAllLines(Paths.get(path), StandardCharsets.UTF_8);
        

        double divergenceScores[][] = new double[genes][genes]; 
        double total = 0;

        for(int row = 0; row < lines.size(); row++){
            String[] genesString;
            genesString = lines.get(row).split(",");
            for (int column = 0; column <= row; column++){
                double temp = Double.parseDouble(genesString[column]);
                
                divergenceScores[row][column] = temp;
                divergenceScores[column][row] = temp;
            }
        }
        return divergenceScores;
    }
    
        private static double[][] readInputSelectiveAdvantage(String path, int genes) throws IOException{
        

        List<String> lines = Files.readAllLines(Paths.get(path), StandardCharsets.UTF_8);
        int columns = lines.get(0).split(",").length;

        double selectiveAdvantages[][] = new double[genes][columns]; 

        for(int row = 0; row < lines.size(); row++){
            String[] elementString;
            elementString = lines.get(row).split(",");
            for (int column = 0; column < columns; column++){
                double temp = Double.parseDouble(elementString[column]);
                
                selectiveAdvantages[row][column] = temp;
            }
        }
        return selectiveAdvantages;
    }
    
    private static double[] readInputNx1Double(String path) throws IOException{
        
        List<String> lines = Files.readAllLines(Paths.get(path), StandardCharsets.UTF_8);
        String[] stringData = new String[lines.size()];
        lines.toArray(stringData);
        double[] doubleData = new double[stringData.length];
        for (int i = 0; i < stringData.length; i++){
            doubleData[i] = Double.parseDouble(stringData[i]);
        }
        return doubleData;
    }
    
    private static int[] readInputTotalMutations(String path) throws IOException{
        List<String> lines = Files.readAllLines(Paths.get(path), StandardCharsets.UTF_8);
        String[] strings = lines.get(0).split(",");
        int[] ret = new int[strings.length];
        for (int i = 0; i < strings.length; i++){
            ret[i] = Integer.parseInt(strings[i]);
        }
        return ret;
    }
    
    private static boolean[][] simulatedSampleData(double[] geneMutationProbabilities,int[] sampleDataMutationCounts, int genes, int samples){
        boolean[][] simulatedSampleData = new boolean[samples][genes];
        
        for (int sample = 0; sample <  samples; sample++){
            int mutationCount = 0;
            while (mutationCount < sampleDataMutationCounts[sample]){
                //add a random mutation
                double rand = Math.random();
                int count = 0;
                while (rand > 0 && count < genes){
                    rand  -= geneMutationProbabilities[count];
                    count++;
                }
                //only add it if it is not mutated already
                if (!simulatedSampleData[sample][count - 1]) {
                    simulatedSampleData[sample][count - 1] = true;
                    mutationCount++;
                } else if (count == genes){
                    //mutation wasn't in genes
                    mutationCount++;
                }
                
            }
        }
        return simulatedSampleData;
    }
    
    private static double[] simulateGeneMutationRate(double geneMutationRateMean, int length){
        double geneMutationRate[] = new double[length];
        //populates the geneMutationRate data based on a normal distribution
       for (int gene = 0; gene < length; gene++){
                 geneMutationRate[gene] = geneMutationRateMean;//mutationRandom.sample();
            }
       return geneMutationRate;
    }
    
    private static void geneInteractions(boolean[][] trainData, boolean[][] simulatedSampleData, String genePath, int genes, int samples, double[] geneMutationProbabilities, int[] sampleDataMutationCount, double nullFitness, double[] rho, double alpha, double epsilon) throws IOException{
        //generates a cell strain for every possible combination of 2 mutations and calculates its division rate, storing it in an excel file
        
        
        List<Double> test = new ArrayList<>();
        File file = new File("Gene Interactions.csv");  
        if ( !file.exists() )
            file.createNewFile();
        String geneNames[] = null;
        try{
            geneNames = readInputFile(genePath);
        } catch (NoSuchFileException err) {
            //window.printError("Gene name file could not be opened.");
            return;
        }
        FileWriter fw = new FileWriter(file);
        Writer writer = new BufferedWriter( fw );
        for (int x = 0; x < genes; x++){
            writer.write("," + geneNames[x]);
        }
        writer.write("\n");
        for (int gene1 = 0; gene1 < genes; gene1++){
            //List<Integer> asdf = new ArrayList<>();
            //asdf.add(gene);
            writer.write(geneNames[gene1] + ",");
            for (int gene2 = 0; gene2 <= gene1; gene2++){
               
                boolean[] genome = new boolean[genes];
                genome[gene1] = true;
                genome[gene2] = true;
                
                int mutationCount;
                
                if (gene1 == gene2){
                        mutationCount = 1;
                } else {
                    mutationCount = 2;
                }
                if (gene1 == 2){
                    int x = 0;
                }
                double fitness = observedExpectedRatioFitness(trainData, simulatedSampleData, genome, rho, mutationCount,  geneMutationProbabilities, sampleDataMutationCount, genes, samples, alpha, epsilon);
        
                double ratio = fitness / nullFitness; 
                
                writer.write(ratio + ",");
                //asdf.remove(1);
            }
            writer.write("\n");
        }
        writer.flush();
        writer.close();
        fw.close();
    }
    
    
    private static void calculateStatistics(Set<CellStrain> strains, int genes,  String genePath, double mutationRate, int generation, long startingCells, long runtime) throws IOException{
            
        //calculate statistics about the final population and display it on the GUI
        if (strains.isEmpty()){
            //window.done(0,"Population extinct",null);
            return;
        }
        
        //calculates some statistics
        double aveRate = 0.0;
        double maxRate = 0.0;
        double minRate = 1.0;
        CellStrain largest = strains.iterator().next(), secondLargest = strains.iterator().next();
        long totalPop = 0;
        List<gene> observedMutations = new ArrayList<>();
        for (int gene = 0; gene < genes; gene++){
            observedMutations.add(new gene(gene));
        }
        
        
        for (CellStrain strain : strains) {
            for (int gene = 0; gene < genes; gene++){
                
                if (strain.genome[gene]){
                    observedMutations.get(gene).mutationFrequency += strain.population;
                }
            }
            totalPop += strain.population;
            double divisionRateStat = strain.divisionRateHistory.get(strain.divisionRateHistory.size() - 1);
            aveRate += divisionRateStat * strain.population;
            if (divisionRateStat > maxRate) {
                maxRate = divisionRateStat;
            } else if (divisionRateStat < minRate) {
                minRate = divisionRateStat;
            }
            if (strain.divisionRateHistory.get(strain.divisionRateHistory.size() - 1) > largest.divisionRateHistory.get(largest.divisionRateHistory.size() - 1)){
                secondLargest = largest;
                largest = strain;
            } else if (strain.divisionRateHistory.get(strain.divisionRateHistory.size() - 1) > secondLargest.divisionRateHistory.get(secondLargest.divisionRateHistory.size() - 1)) {
                secondLargest = strain;
            }
        }
        aveRate /= totalPop;
        
       
        
        //output results
        for (int gene = 0; gene < genes; gene++){
            observedMutations.get(gene).mutationFrequency /= totalPop;
        }
        Collections.sort(observedMutations);
        
        /*
        String[] geneNames;
        try{
            geneNames = readInputFile(genePath);
        } catch (NoSuchFileException err) {
            window.printError("Gene name file could not be opened.");
            return;
        }
        */
        
        StringBuilder sb = new StringBuilder();
        //display mutation frequencies
        boolean exit = false;
        int i = 0;
        List<Integer> significantGenes = new ArrayList<>();
        while (!exit){
            //sb.append("Gene ").append(geneNames[observedMutations.get(i).index]).append(" observed mutation frequency: ").append(observedMutations.get(i).mutationFrequency).append("\n");
            sb.append("Gene ").append(observedMutations.get(i)).append(" observed mutation frequency: ").append(observedMutations.get(i).mutationFrequency).append("\n");
            
            
            i++;
            significantGenes.add(observedMutations.get(i).index);
            if (observedMutations.get(i).mutationFrequency < .1){
                exit = true;
            }        
        }
        sb.append("Runtime in seconds: ").append(runtime / 1000).append("\n");
        sb.append("Gene mutation rate per copy: ").append(mutationRate).append("\n");
        sb.append("Number of strains: ").append(strains.size()).append("\n");
        sb.append("Generations: ").append(generation).append("\n");
        sb.append("Number of starting cells: ").append(startingCells).append("\n");
        sb.append("Number of ending cells: ").append(totalPop).append("\n");
        sb.append("Average division rate: ").append(aveRate).append("\n");
        sb.append("Max division rate: ").append(maxRate).append("\n");
        sb.append("Min division rate: ").append(minRate).append("\n");
        sb.append("Largest Strain:"+ "\n");
        sb.append("Population: ").append(largest.population).append("\n");
        sb.append("Division Rate History: \n");
        
        for(double x : largest.divisionRateHistory){
            sb.append(x).append(" ");
        }
  
        sb.append("\nGene History: \n");
        for(int x = 0; x < genes; x++){
            if (largest.genome[x])
            sb.append(x).append(" ");
        }
        sb.append("\nMutation Times: ");
        for(int x : largest.geneTimes){
            sb.append(x).append(" ");
        }

        
        
        //second largest
        sb.append("Second Largest Strain:\n");
        sb.append("Population: ").append(secondLargest.population).append("\n");
        sb.append("Division Rate History: \n");
        for(double x : secondLargest.divisionRateHistory){
            sb.append(x).append(" ");
        }
        sb.append("\nGene History: \n");
        for(int x = 0; x < genes; x++){
            if (secondLargest.genome[x])
            sb.append(x).append(" ");
        }
        sb.append("\nMutation Times: \n");
        for(int x : secondLargest.geneTimes){
            sb.append(x).append(" ");
        }
        sb.append("--------------------------------------------------------\n");
       
        //window.done(1000,sb.toString(),significantGenes);
        
    }

    private static double calculateMutationFrequencies(boolean[][] trainData, double[] geneMutationFrequency, double[] sampleMutationFrequency, int[] mutationCount, int genes, int samples){
        
        for (int sample = 0; sample < samples; sample++){
            mutationCount[sample] = 0;
        }
        //find the frequency of mutations, used for calculating distance
        double allGeneMutationFrequency = 0.0;
        for (int gene = 0; gene < genes; gene++){
            geneMutationFrequency[gene] = 0.0;
            for (int sample = 0; sample < samples; sample++){
                
                if (trainData[sample][gene] == true){
                    mutationCount[sample]++;
                    allGeneMutationFrequency++;
                    geneMutationFrequency[gene]++;
                    sampleMutationFrequency[sample] += 1.0/genes;
                }
            }
            geneMutationFrequency[gene] /= samples;
        }
        allGeneMutationFrequency /= samples * genes;
        return allGeneMutationFrequency;
    }

    private static double[] zScore(double[] mScores, int strainMutations, double[][] expectedDivergences, double[][] standardDeviations, int[] mutationCount, int samples){
        //turns a mScore into a z score
        //zScore = (mscore - expected mscore) / standard deviation
        double[] zScore = new double[samples];
        for (int sample = 0; sample < samples; sample++){
            zScore[sample] = (mScores[sample] - expectedDivergences[strainMutations][mutationCount[sample]] ) / standardDeviations[strainMutations][mutationCount[sample]];
        }
        return zScore;
    }
    
    private static void generateExpectedDivergence( double[] geneMutationProbabilities, int genes, boolean inverse, boolean subtractBest) throws IOException{
        //generates a table of expected divergences for every combination of samples with any number of mutations
        //expected divergence scores are calculated by generating 10 random samples and calculating the average divergence score
        //also calculates standard deviation
        
        //open up files for writting
        
        File expectedDivergenceFile = new File("expectedDivergence" + genes + ".csv");  
        if ( !expectedDivergenceFile.exists() )
            expectedDivergenceFile.createNewFile();
        File standardDeviationFile = new File("standardDeviation" + genes + ".csv");  
        if ( !standardDeviationFile.exists() )
            standardDeviationFile.createNewFile();
        FileWriter fw = new FileWriter(expectedDivergenceFile);
        Writer expectedWriter = new BufferedWriter( fw );
        FileWriter fw2 = new FileWriter(standardDeviationFile);
        Writer standardDeviationWriter = new BufferedWriter( fw2 );
        
        
        int sampleNumber = 32;
        
        //generate a list of samples
        List<List<boolean[]>> allSamples = new ArrayList<>();
        //first index indicates the number of mutations it will have
        //second index indicates the sample number
        //third index indicates the mutation itself
        for (int gene = 0; gene < genes; gene++){
            List<boolean[]> samples = new ArrayList<>();
            allSamples.add(samples);
            for (int sample = 0; sample < sampleNumber; sample++){
                boolean[] currentSample = new boolean[genes];
                samples.add(currentSample);
                //add genes to the current sample
                int mutationCount = 0;
                while (mutationCount < gene){
                    //add a random mutation
                    double rand = Math.random();
                    int count = 0;
                    while (rand > 0){
                        rand  -= geneMutationProbabilities[count];
                        count++;
                    }
                    //only add it if it is not mutated already
                    if (!currentSample[count - 1]) {
                        currentSample[count - 1] = true;
                        mutationCount++;
                    }
                }
            }
        }
        
        //for each combination of mutations find the average distance
        for (int gene1 = 0; gene1 < genes; gene1++){
            for (int gene2 = 0; gene2 <= gene1; gene2++){
                double average = 0;
                int count = 0;
                double[][] divergences = new double[sampleNumber][sampleNumber];
                for (int sample1 = 0; sample1 < sampleNumber; sample1++){
                    for (int sample2 = 0; sample2 < sampleNumber; sample2++){
                        //don't compare a sample to itself
                        if (gene1 == gene2 && sample1 == sample2){
                            continue;
                        }
                        //get the samples
                        boolean[] genome1 = allSamples.get(gene1).get(sample1);
                        boolean[] genome2 = allSamples.get(gene2).get(sample2);
                        //find the m score for each ecombination of samples
                        
                        double score = mScore (genome1, gene1,genome2, gene2, geneMutationProbabilities, genes);
                        divergences[sample1][sample2] = score;
                        average += score;
                        count++;
                    }
                }
                average /= count;
                expectedWriter.write(average + ",");
                
                //calculate sample standard deviation
                double standardDeviation = 0;
                for (int sample1 = 0; sample1 < sampleNumber; sample1++){
                    for (int sample2 = 0; sample2 < sampleNumber; sample2++){
                        //don't compare a sample to itself
                        if (sample1 == sample2){
                            continue;
                        }
                        standardDeviation += Math.pow(average - divergences[sample1][sample2], 2);
                    }
                }
                
                standardDeviation /= sampleNumber * sampleNumber - 2;
                standardDeviation = Math.sqrt(standardDeviation);
                
                //write to standard deviation file
                standardDeviationWriter.write(standardDeviation + ",");
                
            }
            expectedWriter.write('\n');
            standardDeviationWriter.write('\n');
        }
        
        expectedWriter.flush();
        expectedWriter.close();
        fw.close();
        standardDeviationWriter.flush();
        standardDeviationWriter.close();
        fw2.close();
        
    }
    
    private static void generateExpectedExponentialKernal(boolean sampleData[][], double[] geneMutationProbabilities, int[] numberOfMutations, double[][] expectedDivergences, double[][] standardDeviations, int genes, int samples, boolean inverse, boolean subtractBest) throws IOException{
        //calculates the expected number of kernals that a sample will be in range of
        //given a certain number of random mutations
        
        
        File file = new File("expectedExponentialKernal" + genes + "Const4.csv");  
        if ( !file.exists() )
            file.createNewFile();
        FileWriter fw = new FileWriter(file);
        Writer writer = new BufferedWriter( fw );
        
        
        int sampleNumber = 64;
        
        
        //generate a list of samples
        List<List<boolean[]>> allSamples = new ArrayList<>();
        //first index indicates the number of mutations it will have
        //second index indicates the sample number
        //third index indicates the mutation itself
        for (int gene = 0; gene < genes; gene++){
            List<boolean[]> genomes = new ArrayList<>();
            allSamples.add(genomes);
            for (int sample = 0; sample < sampleNumber; sample++){
                boolean[] currentSample = new boolean[genes];
                genomes.add(currentSample);
                //add genes to the current sample
                int mutationCount = 0;
                while (mutationCount < gene){
                    //add a random mutation
                    double rand = Math.random();
                    int count = 0;
                    while (rand > 0){
                        rand  -= geneMutationProbabilities[count];
                        count++;
                    }
                    //only add it if it is not mutated already
                    if (!currentSample[count - 1]) {
                        currentSample[count - 1] = true;
                        mutationCount++;
                    }
                }
            }
        }
        
        
        for (int gene = 0; gene < genes; gene++){
            double average = 0;
            for (int sample = 0; sample < sampleNumber; sample++){
                boolean[] genome = allSamples.get(gene).get(sample);
                double[] mScores = mScore(sampleData, geneMutationProbabilities,  genome, gene, numberOfMutations, genes, samples);
                //normalize mscores to z scores
                double[] zScore = zScore(mScores, gene, expectedDivergences, standardDeviations,  numberOfMutations, samples);
                
                average += exponentialKernal(zScore, samples);
                
            }
            average /= sampleNumber;
            writer.write(average + "\n");
        }
       
        writer.flush();
        writer.close();
        fw.close();
        
    }
    
    private static double exponentialKernal(double[] zScore, int samples){
        
        // = 1/ (1 + e ^ (c1 + c2 * divergence))
        double exponentialKernal = 0;
        for (int sample = 0; sample < samples; sample++){
            exponentialKernal += 1 / (1 + Math.exp(3 + zScore[sample]));
        }
        return exponentialKernal;
    }
    
    public static double[] optimizeRadius(boolean sampleData[][], double geneMutationProbabilities[], double rho[], int sampleDataMutationCounts[], boolean simulatedSampleData[][],double currentAlpha, double epsilon, double step) throws IOException{
        
        
        
        int genes = sampleData[0].length;
        int samples = sampleData.length;
        final int maxAlphaIndex = 9;
        
        
        boolean[] nullGenome = new boolean[genes];
        
        
        double[] totalFitness = new double[maxAlphaIndex];
        
        
        //for every radius being tested
        double[] nullMScores;
        double[] nullExpectedMScores;
        double[] nullExpected = new double[maxAlphaIndex];
        double[] nullObserved = new double[maxAlphaIndex];
        
        nullMScores = mScore(sampleData, geneMutationProbabilities, nullGenome, 0, sampleDataMutationCounts, genes, samples);
        nullExpectedMScores = mScore(simulatedSampleData, geneMutationProbabilities, nullGenome, 0, sampleDataMutationCounts, genes, samples);

        double alpha;
        
        for (int alphaIndex = 0; alphaIndex < maxAlphaIndex; alphaIndex++){
            alpha = currentAlpha + (alphaIndex - 4) * step;

            

            nullObserved[alphaIndex] = observed(rho, nullMScores, 0, sampleDataMutationCounts, samples, alpha);
            nullExpected[alphaIndex] = observed(rho, nullExpectedMScores, 0, sampleDataMutationCounts, samples, alpha);

            //correct for leave one out
            nullExpected[alphaIndex] *= (samples - 1.0) / samples;

            //for every combination of samples, find the mscore
            for (int sample1 = 0; sample1 < samples; sample1++){
                
                
                //subtract nullMscore since it is a leave one out approach
                double sampleNullFitness = (nullObserved[alphaIndex] - rho[0] / rho[sampleDataMutationCounts[sample1]] * Math.exp(-alpha * nullMScores[sample1])) + epsilon;
        
                sampleNullFitness /= (nullExpected[alphaIndex] + epsilon);
                
                boolean[] genome = sampleData[sample1];
                double fitness;
                double mScores[] = mScore(sampleData, geneMutationProbabilities, genome, sampleDataMutationCounts[sample1], sampleDataMutationCounts, genes, samples);
                double expectedMScores[] = mScore(simulatedSampleData, geneMutationProbabilities, genome, sampleDataMutationCounts[sample1], sampleDataMutationCounts, genes, samples);
                

                double observed = observed(rho, mScores, sampleDataMutationCounts[sample1], sampleDataMutationCounts, samples, alpha);
                double expected = observed(rho, expectedMScores, sampleDataMutationCounts[sample1], sampleDataMutationCounts, samples, alpha);
                
                //correct for leave one out approach
                observed -= 1;
                expected *= (samples - 1.0) / samples;
                
                fitness = (observed + epsilon) / (expected + epsilon);
                
                
              
        
                double score = Math.sqrt(fitness/sampleNullFitness);
                     
                totalFitness[alphaIndex] += score;
                int x = 0;
            }
        }
        
        return totalFitness;
    }
   
    public static double[] optimizeEpsilon(boolean sampleData[][], double geneMutationProbabilities[], double rho[], int sampleDataMutationCounts[], boolean simulatedSampleData[][], double alpha, double currentEpsilon, double step) throws IOException{
        
        
        
        int genes = sampleData[0].length;
        int samples = sampleData.length;
        final int maxEpsilonIndex = 9;
        
        
        boolean[] nullGenome = new boolean[genes];
        
        
        
        double[] totalFitness = new double[maxEpsilonIndex];
        
        
        //for every radius being tested
        double[] nullMScores;
        double[] nullExpectedMScores;
        double[] nullExpected = new double[maxEpsilonIndex];
        double[] nullObserved = new double[maxEpsilonIndex];
        
        nullMScores = mScore(sampleData, geneMutationProbabilities, nullGenome, 0, sampleDataMutationCounts, genes, samples);
        nullExpectedMScores = mScore(simulatedSampleData, geneMutationProbabilities, nullGenome, 0, sampleDataMutationCounts, genes, samples);

        
        
        for (int epsilonIndex = 0; epsilonIndex < maxEpsilonIndex; epsilonIndex++){
            double epsilon = currentEpsilon + (epsilonIndex - 4) * step;

            

            nullObserved[epsilonIndex] = observed(rho, nullMScores, 0, sampleDataMutationCounts, samples, alpha);
            nullExpected[epsilonIndex] = observed(rho, nullExpectedMScores, 0, sampleDataMutationCounts, samples, alpha);

            //correct for leave one out
            nullExpected[epsilonIndex] *= (samples - 1.0) / samples;

            //for every combination of samples, find the mscore
            for (int sample1 = 0; sample1 < samples; sample1++){
                
                
                //subtract nullMscore since it is a leave one out approach
                double sampleNullFitness = (nullObserved[epsilonIndex] - rho[0] / rho[sampleDataMutationCounts[sample1]] * Math.exp(-alpha * nullMScores[sample1])) + epsilon;
        
                sampleNullFitness /= (nullExpected[epsilonIndex] + epsilon);
                
                boolean[] genome = sampleData[sample1];
                double fitness;
                double mScores[] = mScore(sampleData, geneMutationProbabilities, genome, sampleDataMutationCounts[sample1], sampleDataMutationCounts, genes, samples);
                double expectedMScores[] = mScore(simulatedSampleData, geneMutationProbabilities, genome, sampleDataMutationCounts[sample1], sampleDataMutationCounts, genes, samples);
                

                double observed = observed(rho, mScores, sampleDataMutationCounts[sample1], sampleDataMutationCounts, samples, alpha);
                double expected = observed(rho, expectedMScores, sampleDataMutationCounts[sample1], sampleDataMutationCounts, samples, alpha);
                
                //correct for leave one out approach
                observed -= 1;
                expected *= (samples - 1.0) / samples;
                
                fitness = (observed + epsilon) / (expected + epsilon);
                
                
              
        
                double score = Math.sqrt(fitness/sampleNullFitness);
                     
                totalFitness[epsilonIndex] += score;
                int x = 0;
            }
        }
        
        return totalFitness;
    }
    
    public static void sampleDataZScoreCompare(boolean[][] sampleData, double[] geneMutationProbabilities, int[] SampleDataMutationCount, double[][] expectedDivergences, double[][] standardDeviations, int genes, int samples, boolean inverse, boolean subtractBest) throws IOException{
        

        //find and store the zScore between each sample and every other sample
        double[][] zScores = new double[samples][samples];
        for (int sample1 = 0; sample1 < samples; sample1++){
            for (int sample2 = 0; sample2 < sample1; sample2++){
                //get the mScore for each pair of samples
                boolean[] genome1 = sampleData[sample1];
                boolean[] genome2 = sampleData[sample2];
                double mScore = mScore(genome1, SampleDataMutationCount[sample1], genome2, SampleDataMutationCount[sample2], geneMutationProbabilities, genes);
                //convert mScore into zScore
                double zScore = (mScore - expectedDivergences[SampleDataMutationCount[sample1]][SampleDataMutationCount[sample2]]) / standardDeviations[SampleDataMutationCount[sample1]][SampleDataMutationCount[sample2]];
                zScores[sample1][sample2] = zScore;
                zScores[sample2][sample1] = zScore;
            }
        }
        //write the output to a csv file
        
        File zScoresFile = new File("zScores" + genes + ".csv");  
        if ( !zScoresFile.exists() )
            zScoresFile.createNewFile();
        
        FileWriter fw = new FileWriter(zScoresFile);
        Writer writer = new BufferedWriter( fw );
        
        
        for (int sample1 = 0; sample1 < samples; sample1++){
            for (int sample2 = 0; sample2 < samples; sample2++){
                writer.write(zScores[sample1][sample2] + ",");
            }
            writer.write('\n');
        }
        
        writer.flush();
        writer.close();
        fw.close();
    }
    
    private static void writeRho(int mutationCount[], int samples, double constant) throws IOException{
                File rhoFile = new File("rhoConstant" + constant + "Samples" + samples + ".csv");  
        if ( !rhoFile.exists() )
            rhoFile.createNewFile();
        
        FileWriter fw = new FileWriter(rhoFile);
        Writer writer = new BufferedWriter( fw );
        
        double[] rho = new double[9422];
        
        //calculate rho
        for (int mutations = 0; mutations < 9422; mutations++){
            for (int sample = 0; sample < samples; sample++){
                rho[mutations] += 1/(Math.sqrt(2 * Math.PI)) * Math.exp( -1.0 / (2 * constant * constant) * Math.pow(Math.log(mutations + 1) - Math.log(mutationCount[sample] + 1), 2));
            }
            rho[mutations] *= 1.0 / (mutations + 1);
            //write rho
            writer.append(rho[mutations] + "\n");
        }
        
        
        writer.flush();
        writer.close();
        fw.close();
    }
    
    private static void optimizeAlphaBeta(boolean[][] trainData, double[] geneMutationProbabilities, double[] rho, int[] sampleDataMutationCounts, boolean[][] simulatedSampleData) throws IOException{
        double currentAlpha = .85;
        double currentEpsilon = .85;
        double nextAlpha = .85;
        double nextEpsilon = .85;
        int runs = 1;
        for (double step = .2; step > .01; step /= 2){
        

            double[] totalAlphaFitness = new double[9];
            for (int run = 0; run < runs; run++){
                double optimalAlphaFitness[] = optimizeRadius(trainData, geneMutationProbabilities, rho, sampleDataMutationCounts, simulatedSampleData, currentAlpha, currentEpsilon, step);

                //geneMutationProbabilities = geneMutationProbabilities(geneLength,genes,10290030);
                for (int i = 0; i < 9; i++){
                    totalAlphaFitness[i] += optimalAlphaFitness[i];
                }
            }

            //find max alpha
            double max = 0;
            for (int i = 0; i < 9; i++){
                if (totalAlphaFitness[i] > max & currentAlpha + (i - 4) * step > 0){
                    max = totalAlphaFitness[i];
                    nextAlpha = currentAlpha + (i - 4) * step;

                }

            }
            currentAlpha = nextAlpha;
            
            
            double[] totalEpsilonFitness = new double[9];
            for (int run = 0; run < runs; run++){
                double optimalEpsilonFitness[] = optimizeEpsilon( trainData,  geneMutationProbabilities, rho, sampleDataMutationCounts, simulatedSampleData, currentAlpha, currentEpsilon, step);

                //geneMutationProbabilities = geneMutationProbabilities(geneLength,genes,10290030);
                for (int i = 0; i < 9; i++){
                    totalEpsilonFitness[i] += optimalEpsilonFitness[i];
                }
            }
            //find max epsilon
            max = 0;
            for (int i = 0; i < 9; i++){
                if (totalEpsilonFitness[i] > max && currentEpsilon + (i - 4) * step > 0){
                    max = totalEpsilonFitness[i];
                    nextEpsilon = currentEpsilon + (i - 4) * step;

                }

            }
            currentEpsilon = nextEpsilon;
        
        }

        
    }
    
    private static List<List<List<Subdivision>>> setUpSubdivisions(int subdivisionsCount){
        //set subdivisions size for subdivisions
        int maxGridSize = (int) Math.round(Math.pow(subdivisionsCount, 1.0 / 3.0));
        
        List<List<List<Subdivision>>> subdivisions = new ArrayList<>();
        
        //initialize and add all the subdivisions to the subdivisions
        for (int x = 0; x < maxGridSize; x++){
            List<List<Subdivision>> newX = new ArrayList<>();
            subdivisions.add(newX);
            for (int y = 0; y < maxGridSize; y++){
                List<Subdivision> newY = new ArrayList<>();
                subdivisions.get(x).add(newY);
                for (int z = 0; z < maxGridSize; z++){
                    subdivisions.get(x).get(y).add(new Subdivision(x,y,z));
                }
            }
        }
        
        //set up the subdivisions neighbors
        
        for (int x = 0; x < maxGridSize; x++){
            for (int y = 0; y < maxGridSize; y++){
                for (int z = 0; z < maxGridSize; z++){
                    Subdivision subdivision = subdivisions.get(x).get(y).get(z);
                    for (int dimension = 0; dimension < 3; dimension++){
                        for(int offset = -1; offset <=1; offset+=2){
                            int newX = x;
                            int newY = y;
                            int newZ = z;
                            switch (dimension){
                                case 0:
                                    newX += offset;
                                    break;
                                case 1:
                                    newY += offset;
                                    break;
                                case 2:
                                    newZ += offset;
                            }
                            if (newX < 0 || newX == maxGridSize){
                                continue;
                            }
                            if (newY < 0 || newY == maxGridSize){
                                continue;
                            }
                            if (newZ < 0 || newZ == maxGridSize){
                                continue;
                            }
                            
                            subdivision.neighbors.add(subdivisions.get(newX).get(newY).get(newZ));
                        }
                    }
                }
            }
        }
        
        
        return subdivisions;
    }
    
    private static void migrate(Set<CellStrain> strains, List<List<List<Subdivision>>> subdivisions, long maxSubdivisionPopulation, int maxGridSize, double flowRate, int genes){
        
        for (CellStrain strain : strains){
            strain.processed = false;
        }
        
        for (int x = 0; x < maxGridSize; x++){
                
            for (int y = 0; y < maxGridSize; y++){
                    
                for (int z = 0; z < maxGridSize; z++){
                    Subdivision subdivision = subdivisions.get(x).get(y).get(z);
                        double migratePercentage = flowRate / subdivision.neighbors.size();
                        for (CellStrain strain: subdivisions.get(x).get(y).get(z).strains){
                            if (!strain.processed && strain.population > 0){
                                for (Subdivision neighbor: subdivision.neighbors){
                                    long migratePopulation = Math.min(strain.population,(int) Math.round(strain.population * migratePercentage + nextGaussian() * Math.sqrt(strain.population * migratePercentage)));
                                    if (migratePopulation <= 0){
                                        continue;
                                    }
                                    //check if the neighboring subdivision already contains a strain with the same genome
                                    boolean migrate = false;
                                    for (CellStrain neighborStrain: neighbor.strains){
                                        if (strain.genomeID == neighborStrain.genomeID){
                                            strain.population -= migratePopulation;
                                            neighborStrain.population += migratePopulation;
                                            migrate = true;
                                            break;
                                        }
                                    }
                                    if (!migrate){
                                        //make a new strain
                                        strain.population -= migratePopulation;

                                        List<Double> newDivisionRateHistory = new ArrayList<>(strain.divisionRateHistory);
                                         //calculate new divisionRate

                                        //calculate the geneHistory for the new strain
                                        boolean[] newGenome = new boolean[genes];
                                        System.arraycopy(strain.genome, 0, newGenome, 0, genes );


                                        //calculate the geneTimes for the new strain
                                        List<Integer> newGeneTimes = new ArrayList<>(strain.geneTimes);

                                        int newMutationCount = strain.mutationCount;



                                        CellStrain newStrain = new CellStrain(newDivisionRateHistory,migratePopulation,newGenome,newGeneTimes, newMutationCount, neighbor.x, neighbor.y, neighbor.z);
                                        newStrain.genomeID = strain.genomeID;
                                        newStrain.processed = true;
                                        strains.add(newStrain);

                                        neighbor.strains.add(newStrain);
                                    }
                                
                            } 
                        }
                    }
                }
            }
        }
    }
    
    private static double[] generateSimulatedGeneFitness(int genes, double standardDeviation){
        
        double[] geneFitness = new double[genes];
        
        for (int i = 0; i < genes; i++){
            geneFitness[i] = standardDeviation * nextGaussian();
        }
        return geneFitness;
        
    }
    
    private static void writeFitness(int generation,List<List<List<Subdivision>>> subdivisions, int maxGridSize, long subdivisionPop, Writer writer) throws IOException{
       
        
        for (int x = 0; x < maxGridSize; x++){
            for (int y = 0; y < maxGridSize; y++){
                for (int z = 0; z < maxGridSize; z++){
                    Subdivision subdivision = subdivisions.get(x).get(y).get(z);
                    for (CellStrain strain : subdivision.strains){
                        
                        double fitness = strain.divisionRateHistory.get(strain.divisionRateHistory.size() - 1);
                        if (strain.genomeID != 0){
                            double growth =  fitness * subdivision.populationScalingFactor;
                            double pop = strain.population;
                            Random generator = new Random(strain.genomeID + startTime);
                            Random gen = new Random(strain.hashCode());

                            writer.append(generation + "," + (x + gen.nextDouble()) + "," + (y + gen.nextDouble()) + "," + (z + gen.nextDouble()) + "," + pop + "," + generator.nextInt() + "\n");
  
                        }
                    }
                    
                    
                }
            }
        }
        

    }
    
    private static void top3Genes(List<Set<CellStrain>> allStrains, int genes, long pop, int generation) throws IOException{
        
        //find the mutaiton frequencies of the genes
        List<gene> observedMutations = new ArrayList<>();
        for (int gene = 0; gene < genes; gene++){
            observedMutations.add(new gene(gene));
        }  
        
        /*
            for (CellStrain strain : allStrains.get(generation / 10)) {
                for (int gene = 0; gene < genes; gene++){

                    if (strain.genome[gene]){
                        observedMutations.get(gene).mutationFrequency += strain.population;
                    }
                }
            }
        */
                
       //sort them to find the most frequently mutated genes
        for (int gene = 0; gene < genes; gene++){
            observedMutations.get(gene).mutationFrequency /= pop;
        }
        Collections.sort(observedMutations);
        int first = observedMutations.get(0).index;
        int second = observedMutations.get(1).index;
        int third = observedMutations.get(2).index;
        
        
                File file = new File("top3genes.csv");  
        if ( !file.exists() )
            file.createNewFile();
        FileWriter fw = new FileWriter(file);
        Writer writer = new BufferedWriter( fw );
        
        

        for (Set<CellStrain> strains : allStrains){
            double[][][] percentages = new double[2][2][2];
        
            for (CellStrain strain : strains){
                percentages[strain.genome[first] ? 1 : 0][strain.genome[second] ? 1 : 0][strain.genome[third] ? 1 : 0] += ((double) strain.population)/pop;
                
            }
            writer.append(percentages[0][0][0] + "," + percentages[0][0][1]+ "," +percentages[0][1][0]+ "," +percentages[0][1][1]+ "," +percentages[1][0][0]+ "," +percentages[1][0][1]+ "," +percentages[1][1][0]+ "," +percentages[1][1][1] + "\n");
        }
        

        
        writer.flush();
        writer.close();
        fw.close();   
    }
    
    public static void changeDrugRegime(int value){
        drugRegime = value;
        regimeChange = true;
    }
    
    public static double nextNormal(double mean){
        return nextGaussian() * Math.abs(mean - 1.0) / 3.0 + mean;
    }
    
    private static long[][][] count(Set<CellStrain> strains){
        long[][][] counts = new long[2][2][2];
        
        for (CellStrain strain : strains){
            if (strain.genome[0]){
                counts[strain.genome[1]? 1: 0][strain.genome[2]? 1 : 0][strain.genome[3]? 1: 0] += strain.population;
            }
        }
        
        return counts;
    }
    
    public static Tuple run(long startingCells, long runtime, double mutationRate, int subdivisionsCount, double flowRate, int genes, double resistantRate, double sensitiveRate, boolean sequential, boolean smart, int switchEvery, boolean matlab, long cutOff, boolean optimize) throws IOException, InterruptedException{
       
        
        
        //double alpha = .01;
        //double epsilon = 1;
        
        int maxGridSize = (int) Math.round(Math.pow(subdivisionsCount, 1.0 / 3.0));
        int midpoint = maxGridSize / 2;
        List<List<List<Subdivision>>> subdivisions = setUpSubdivisions(subdivisionsCount);
        
        long subdivisionPopulation = startingCells / subdivisionsCount;
        
        long startTime = System.currentTimeMillis();
        //trainData must be final to be used in parallel stream
        //first index is sample number, second index is the gene number, a 1 indicates a mutation

        //double[][] geneFitness = readInputSelectiveAdvantage("src/cancersimulation/selectiveAdvantagesKit.csv", genes);
        
        double[][][][][] fitness = {
            //no drug
            {{{{1.0, 1.0}, {1.0, 1.0}},{{1.0, 1.0}, {1.0, 1.0}}},{{{nextNormal(resistantRate), nextNormal(resistantRate)}, {nextNormal(resistantRate), nextNormal(resistantRate)}},{{nextNormal(resistantRate), nextNormal(resistantRate)}, {nextNormal(resistantRate), nextNormal(resistantRate)}}}},
            //drug 1
            {{{{1.0, 1.0}, {1.0, 1.0}},{{1.0, 1.0}, {1.0, 1.0}}},{{{nextNormal(sensitiveRate), nextNormal(resistantRate)}, {nextNormal(resistantRate), nextNormal(resistantRate)}},{{nextNormal(resistantRate), nextNormal(resistantRate)}, {nextNormal(resistantRate), nextNormal(resistantRate)}}}},
            //drug 2
            {{{{1.0, 1.0}, {1.0, 1.0}},{{1.0, 1.0}, {1.0, 1.0}}},{{{nextNormal(sensitiveRate), nextNormal(resistantRate)}, {nextNormal(resistantRate), nextNormal(resistantRate)}},{{nextNormal(sensitiveRate), nextNormal(resistantRate)}, {nextNormal(resistantRate), nextNormal(resistantRate)}}}},
            //drug 3
            {{{{1.0, 1.0}, {1.0, 1.0}},{{1.0, 1.0}, {1.0, 1.0}}},{{{nextNormal(sensitiveRate), nextNormal(resistantRate)}, {nextNormal(sensitiveRate), nextNormal(resistantRate)}},{{nextNormal(resistantRate), nextNormal(resistantRate)}, {nextNormal(resistantRate), nextNormal(resistantRate)}}}}
        };
        double[] geneMutationProbabilities = {0,.4,.3,.3};
    
        
        drugRegime = 0;
        
        //double[] geneFitness = generateSimulatedGeneFitness(genes, .02);
        
        File file = new File("MATLABData.csv");  
        if ( !file.exists() )
            file.createNewFile();
        FileWriter fw = new FileWriter(file);
        Writer writer = new BufferedWriter( fw );
        boolean exists = true;
        
        File file2 = null;
        Writer writer2 = null;
        FileWriter fw2 = null;
        if(!optimize){
            String name = resistantRate + "," + sensitiveRate + "rates,flowRate=" + flowRate + "trialsSequential.csv";

            if (!sequential){
                name = resistantRate + "," + sensitiveRate + "rates,flowRate=" + flowRate + "trialsAlternating.csv";
            } else if (smart){
                name = resistantRate + "," + sensitiveRate + "rates,flowRate=" + flowRate + "trialsSmart.csv";
            }


            file2 = new File(name);  
            if ( !file2.exists() ){
                file2.createNewFile();
                exists = false;
            }
            fw2 = new FileWriter(file2,true);
            writer2 = new BufferedWriter( fw2 );
            writer2.append("\n");
            if (!exists)

                if (sequential){
                writer2.append("IM regime begins (day), smallest tumor size (g), smallest at (day), tumor size at progression (g), SU regime begins (day), smallest tumor size (g), smallest at (day), tumor size at progression (g),REGO regime begins (day), smallest tumor size (g), smallest at (day), tumor size at progression (g), time of progression (day), 500g tumor at (day), % resistant,divergenceDay\n");

                } else {
                writer2.append("IM regime begins (day), smallest tumor size (g), smallest at (day), tumor size at progression (g), SU/REGO regime begins (day), smallest tumor size (g), smallest at (day), tumor size at progression (g), time of progression (day), 500g tumor at (day), % resistant,divergenceDay\n");

            }
        }
        /*
        
        //for working with kernal density approximation
        
        //tuning constants for 399 genes
        double alpha = .042;
        double epsilon = .63;
        
        final boolean[][] trainDataf = trainData;
        
        
        //the number of samples and genes respectively
        int samples = trainDataf.length;
        int genes = trainDataf[0].length;
        
        //read in input files
        double rho[] = readInputNx1Double("src/cancersimulation/rhoConstant" + rhoConstant + "Samples" + samples + ".csv");
        int[] sampleDataMutationCounts = readInputTotalMutations("src/cancersimulation/totalmutations");
        final double geneMutationProbabilities[] = geneMutationProbabilities(geneLength,genes,10290030);
        boolean simulatedSampleData[][] = simulatedSampleData(geneMutationProbabilities, sampleDataMutationCounts, genes, samples);
        
        boolean[] genome = new boolean[genes];
        double nullFitness = observedExpectedRatioFitness(trainData, simulatedSampleData, genome, rho, 0, geneMutationProbabilities, sampleDataMutationCounts, genes, samples, alpha, epsilon);
       
        */
        
        
        int genomeID = 1;
        int generation = 0;
        long prevPopulation = startingCells;
        long startingMutated = Math.min(startingCells / subdivisionsCount, 10000);
        long mutated = startingMutated;
        
        List<Set<CellStrain>> allStrains = new ArrayList<>();
        Set<CellStrain> strains = new HashSet<>();

        for (int x = 0; x < maxGridSize; x++){
            for (int y = 0; y < maxGridSize; y++){
                for (int z = 0; z < maxGridSize; z++){
                    
                    boolean[] genome = new boolean[genes];
                    CellStrain startingStrain = null;
                    List<Double> divisionRate = new ArrayList<>();
                    
                    
                    divisionRate.add(1.0);

                    
                        List<Integer> times = new ArrayList<>();
                        times.add(0);
                    startingStrain = new CellStrain(divisionRate,subdivisionPopulation, genome, times, 0, x, y, z);
                    startingStrain.genomeID = 0;
                    
                    
                    if (x == midpoint && y == midpoint && z == midpoint){
                    
                        boolean[] genome2 = new boolean[genes];
                        genome2[0] = true;
                        CellStrain startingStrain2 = null;

                        List<Double> divisionRate2 = new ArrayList<>();


                        divisionRate2.add(fitness[0][1][0][0][0]);

                        List<Integer> times2 = new ArrayList<>();
                        times2.add(0);
                        startingStrain2 = new CellStrain(divisionRate2,startingMutated, genome2, times2, 1, x, y, z);
                        startingStrain2.genomeID = -1;

                        startingStrain.population -= startingMutated;
                        strains.add(startingStrain2);
                        subdivisions.get(x).get(y).get(z).add(startingStrain2);
                    }
                        
                    //add the original unmutated strain to the set of strains
                    strains.add(startingStrain);
                    subdivisions.get(x).get(y).get(z).add(startingStrain);
                }
            }
        }
        //geneInteractions( trainData, simulatedSampleData, genePath, genes, samples, geneMutationProbabilities, sampleDataMutationCounts, nullFitness, rho, alpha, epsilon);
        
        //window.setStopEnabled(true);
        //keep looping through generations until the simulation has timed out
        double divergenceDay = 0.0;
        double difference = 9999.0;
        long lowestPopulation = mutated;
        int lowestGeneration = 0;
        int printed = 0;
        int halvingGen = 0;
        long halvingPop = 0;
        double predictedPop = (double) startingCells;
        double doublingTime = 0.0;
        while (System.currentTimeMillis() - startTime < runtime && !stop && mutated > 0 && mutated < startingCells * .75) {
        //writer2.append("IM regime begins (week), smallest tumor size (g), smallest at (week), tumor size at progression (g), SU regime begins (week), smallest tumor size (g), smallest at (week), tumor size at progression (g), REG regime begins (week), smallest tumor size (g), smallest at (week), tumor size at progression (g), 100 g tumor at (week)\n");
        
            
            if (!sequential){
                if (drugRegime == 2 && generation % (2 * switchEvery) == 0){
                    drugRegime = 3;
                    regimeChange = true;
                } else if (drugRegime == 3 && generation % (2 * switchEvery) == switchEvery){
                    drugRegime = 2;
                    regimeChange = true;
                }

            }
            
            if (mutated < lowestPopulation){
                lowestPopulation = mutated;
                lowestGeneration = generation;
            }
            if (mutated > cutOff && drugRegime == 0){
                       
                //reached 100g tumor
                if(!optimize)
                writer2.append(Double.toString(generation) + ",");
                lowestPopulation = mutated;
                lowestGeneration = generation;
                halvingGen = generation;
                halvingPop = mutated;
                drugRegime = 1;
                regimeChange = true;
                doublingTime = Math.log(2)*halvingGen/(Math.log(((double)halvingPop) / startingMutated));//1.267 * (generation - lowestGeneration);
                    
                
            } else if (mutated > lowestPopulation * 1.728 && mutated > 100000000 && drugRegime > 0 && printed < 2){
                if (drugRegime == 1){
                    
                    divergenceDay = generation;
                    //calculate halving time
                    //double halvingTime = Math.log(.5) * (lowestGeneration - halvingGen) / Math.log( ((double)lowestPopulation) / halvingPop);
                    double progressionTime = generation - halvingGen;
                    System.out.println("Progression time: " + progressionTime);
                    //calculate doubling time
                    
                    System.out.println("Doubling time: " + doublingTime);
                    
                    
                    //double doublingTime = halvingGen / (Math.log(halvingPop / startingMutated) / Math.log(2));
                    //System.out.println("Doubling time: " + doublingTime);
                    if (optimize){
                        return new Tuple(doublingTime,progressionTime);
                    }
                }
                
                if(!optimize){
                writer2.append(Double.toString(lowestPopulation / 100000000.0) + ",");
                    
                writer2.append(Double.toString(lowestGeneration) + ",");
                    
                writer2.append(Double.toString(mutated / 100000000.0) + ",");
                
                
                writer2.append(Double.toString(generation) + ",");
                
                }
                lowestPopulation = mutated;
                lowestGeneration = generation;
                
                
                if (sequential){
                    
                       if (smart){
                           switch (drugRegime){
                               case 1:
                                   
                                long[][][] counts = count(strains);
                                drugRegime = counts[1][1][0] > counts[1][0][1]? 2: 3;
                                
                                regimeChange = true;
                                   break;
                               case 2:
                                   drugRegime = 3;
                                   printed++;
                                regimeChange = true;
                                   break;
                               case 3:
                                   drugRegime = 2;
                                   printed++;
                                regimeChange = true;
                                   break;
                           }
                           
                            
                       } else {
                           if (drugRegime == 3){
                                printed = 2;
                            } else {
                            drugRegime++;
                            
                        regimeChange = true;
                       }
                    }
                } else {
                    if (drugRegime == 1){
                        
                        drugRegime = 2;
                        regimeChange = true;
                    } else {
                        
                        printed = 2;
                        
                    }
                }
                
                    
                 
            } else if (mutated > 5 * cutOff && printed > 1 && !optimize){
                       
                //reached 500g tumor
                writer2.append(Double.toString(generation) + ",");
                writer2.append(Double.toString(((double)(strains.stream().mapToLong(x -> (x.genome[3])? x.population : 0).sum())) /  mutated) + ",");
                writer2.append(Double.toString(generation - divergenceDay) + ",");
                difference = generation - divergenceDay;
                stop = true;
                
            } 
            
            if (regimeChange){
                for(CellStrain strain : strains){
                    double divisionRate = fitness[drugRegime][strain.genome[0]? 1 :0][strain.genome[1]? 1 :0][strain.genome[2]? 1 :0][strain.genome[3]? 1 :0];

                    
                    strain.divisionRateHistory.add(divisionRate);
                    
                }
                
                regimeChange = false;
            }
            
            
            if (matlab && generation % 10 == 0){
                Set<CellStrain> tempSet = new HashSet<>();
                for (CellStrain strain: strains){
                    boolean[] newGenome = new boolean[genes];
                    System.arraycopy(strain.genome, 0, newGenome, 0, genes );
                            
                    CellStrain newStrain = new CellStrain(new ArrayList<Double>(strain.divisionRateHistory),strain.population,newGenome,new ArrayList<Integer>(strain.geneTimes),strain.mutationCount,strain.x,strain.y,strain.z);
                    tempSet.add(newStrain);
                }
                
                    allStrains.add(tempSet);
                writeFitness(generation, subdivisions,  maxGridSize, subdivisionPopulation,writer);
            }
            
            
            generation++;
            final int generationf = generation;
            final long prevPopulationf = prevPopulation;
            final double predictedPopf = predictedPop;
            //go through each strain
            strains.stream().parallel().forEach((CellStrain strain) -> {
                //calculate the new population for the strain
                double cellDivisionRate = strain.divisionRateHistory.get(strain.divisionRateHistory.size() - 1);
                long newPop = 0;
                long oldPop = strain.population;
                double sd;
                double r;
                double scalingFactor = subdivisions.get(strain.x).get(strain.y).get(strain.z).populationScalingFactor;
                //use more accurate poisson distribution for small populations and the faster normal distributions for large populations
                if (strain.population * cellDivisionRate>0 && strain.population <= 20 ){
                    //small pop
                    if (generationf % 5 == strain.geneTimes.get(strain.geneTimes.size() - 1) % 5){
                        //only calculate every 5 days
                    newPop =  2 * nextPoisson(strain.population * Math.pow(scalingFactor * cellDivisionRate, 5.0) / 2.0);
                    } else {
                        newPop = strain.population;
                    }
                    
                } else if (strain.population * cellDivisionRate > 0) {
                    sd = Math.sqrt( scalingFactor * strain.population*cellDivisionRate);
                    r = nextGaussian();
                    newPop = Math.round(scalingFactor * strain.population * cellDivisionRate +  r*sd );
                    
                    if (newPop == 0){
                        int x = 0;
                    }
                }
                strain.population = newPop;
                
                if (newPop > 0 && strain.genomeID != 0 && (oldPop > 20 || generationf % 5 == strain.geneTimes.get(strain.geneTimes.size() - 1) % 5)){
                    //calculate any mutations
                    double odds;
                    if (oldPop <= 20){
                        odds = Math.pow(1.0-mutationRate,newPop / 2.0);
                    } else {
                        odds = Math.pow(1.0-mutationRate,(newPop - oldPop * 4.0/5.0) / 2.0);
                    }
                        if (strain.mutationCount < genes && Math.random() > odds){
                             
                            
                            //a unique mutation has occured                            
                            newPop--;
                            List<Double> newDivisionRateHistory = new ArrayList<>(strain.divisionRateHistory);
                             //calculate new divisionRate

                            //calculate the geneHistory for the new strain
                            boolean[] newGenome = new boolean[genes];
                            System.arraycopy(strain.genome, 0, newGenome, 0, genes );
                            
                            
                            int count = 0;
                            do {
                            //add a random mutation
                            double rand = Math.random();
                            
                            count = 0;
                            while (rand > 0 && count < genes){
                                rand  -= geneMutationProbabilities[count];
                                count++;
                            }
                            
                            } while(newGenome[count - 1]);
                            
                            newGenome[count - 1] = true;
                            
                            
                            
                            //calculate the geneTimes for the new strain
                            List<Integer> newGeneTimes = new ArrayList<>(strain.geneTimes);
                            newGeneTimes.add(generationf);
                            int newMutationCount = strain.mutationCount + 1;
                            
                            /*
                            double fitness = observedExpectedRatioFitness(trainData, simulatedSampleData, genome, rho, newMutationCount,  geneMutationProbabilities, sampleDataMutationCounts, genes, samples, alpha, epsilon);
        
                            double newDivisionRate = 1 + .01 * Math.log(fitness/nullFitness);
                            */
                            
                            
                            double newDivisionRate = fitness[drugRegime][newGenome[0]? 1 :0][newGenome[1]? 1 :0][newGenome[2]? 1 :0][newGenome[3]? 1 :0];
                            
                            newDivisionRateHistory.add(newDivisionRate);
                            
                            CellStrain newStrain = new CellStrain(newDivisionRateHistory,1,newGenome,newGeneTimes, newMutationCount, strain.x, strain.y, strain.z);

                            //add the new strain as the child of the parent strain
                            strain.children.add(newStrain);
                        }    
                    

                }
            });
            
        Set<CellStrain> newStrains = new LinkedHashSet<>();

        prevPopulation = 0;
        //add all the strains into a new set
        for(CellStrain strain: strains){
            for (CellStrain child: strain.children){
                //update the genomeID
                child.genomeID = genomeID;
                genomeID++;
                newStrains.add(child);
                subdivisions.get(child.x).get(child.y).get(child.z).strains.add(child);
            }
            if (strain.population > 0){
                strain.children.clear();
                newStrains.add(strain);
            } else {
                subdivisions.get(strain.x).get(strain.y).get(strain.z).strains.remove(strain);
            }
        }
        
        //reset strains
        strains = newStrains;
        
        //migrate cells between subdivisions
        migrate( strains, subdivisions, subdivisionPopulation, maxGridSize, flowRate, genes);
        
        
        //limit the population growth to the amount specified by user
        //limitGrowth(strains, newStrains, maxDivisionRate, prevPopulationf, startingCells);
        
        //reset all population scaling factors
        for (int x = 0; x < maxGridSize; x++){
            for (int y = 0; y < maxGridSize; y++){
                for (int z = 0; z < maxGridSize; z++){
                  subdivisions.get(x).get(y).get(z).populationScalingFactor = 0;
                }
            }
        }
        predictedPop = 0.0;
        //add up the predicted population for the next generation for each subdivision
        for (CellStrain strain: strains){
            double temp = strain.population * strain.divisionRateHistory.get(strain.divisionRateHistory.size() -1);
            predictedPop += temp;
            subdivisions.get(strain.x).get(strain.y).get(strain.z).populationScalingFactor += temp;
        }
        
        //adjust scaling factor
        for (int x = 0; x < maxGridSize; x++){
            for (int y = 0; y < maxGridSize; y++){
                for (int z = 0; z < maxGridSize; z++){
                  subdivisions.get(x).get(y).get(z).populationScalingFactor = Math.min(1.0,subdivisionPopulation / subdivisions.get(x).get(y).get(z).populationScalingFactor);
                }
            }
        }
        
        
        prevPopulation = strains.stream().mapToLong(x -> x.population).sum();
        double maxDivRate = strains.stream().mapToDouble(x -> x.divisionRateHistory.get(x.divisionRateHistory.size() - 1)).max().getAsDouble();
        
        mutated = strains.stream().mapToLong(x -> (x.genomeID != 0)? x.population : 0).sum();
        long first = strains.stream().mapToLong(x -> (x.genome[0] && x.mutationCount == 1)? x.population : 0).sum();
        long second = strains.stream().mapToLong(x -> (x.genome[1])? x.population : 0).sum();
        long third = strains.stream().mapToLong(x -> (x.genome[2])? x.population : 0).sum();
        long fourth = strains.stream().mapToLong(x -> (x.genome[3])? x.population : 0).sum();
    
        int y = 0;
        int brafSubdivisions = strains.stream().mapToInt(x -> (x.genomeID == -1)? 1 : 0).sum();

        if (generation % 20 == 0){
            int x = 0;
        }
        //update the output to the gui
        //window.update(Integer.toString(generation),NumberFormat.getIntegerInstance().format(prevPopulation),Integer.toString(strains.size()),Long.toString((startTime + runtime - System.currentTimeMillis()) / 1000),(int)(1000.0 * (System.currentTimeMillis() - startTime) / runtime),maxDivRate,first, second, third, fourth, mutated);
    }
         
        //top3Genes(allStrains,genes,startingCells,generation);
        //calculates some statistics
        //calculateStatistics( strains, genes,"", mutationRate,generation, startingCells,runtime);
        
        //window.setStopEnabled(false);
        //window.setRunEnabled(true);
        stop = false;
 
        writer.flush();
        writer.close();
        fw.close();  
        if(!optimize){
            writer2.flush();
            writer2.close();
            fw2.close();
            
        return new Tuple(difference,0.0);
        } else {
            return new Tuple(doublingTime,9999);
        }
        /*
        Thread thread = new Thread(){
        @Override
        public void run(){
            try {
                CancerSimulation.run(startingCells, runtime, mutationRate, subdivisionsCount, flowRate,  genes);
            } catch (IOException ex) {
                Logger.getLogger(CancerSimulation.class.getName()).log(Level.SEVERE, null, ex);
            } catch (InterruptedException ex) {
                Logger.getLogger(CancerSimulation.class.getName()).log(Level.SEVERE, null, ex);
            }
       
        
        }
        };
        thread.start();
        */
        
    }

    public static void stop(){
        stop = true;
    }
    
    /**
     * Main method.
     *
     * @param args
     *            the command line arguments; unused here
     * @throws java.io.FileNotFoundException
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {
        long startingCells = 100000000000l;
        long runtime = 600000;
        double mutationRate = 0.000000075;
        int subdivisionsCount = 1000;
        double flowRate = .1;
        int genes = 4;
        /*for 1000 subdivisions*/
        double sensitiveRate = .97;
        double resistantRate = 1.0015;
        
        /*
        double sensitiveRate = .998975;
        double resistantRate = 1.0185;
        */int trials = 100;
        boolean sequential = true;
        boolean smart = false;
        int switchEvery = 3;
        boolean matlab = false;
        long cutOff = 10000000000L;
        boolean optimize = true;
        boolean flip = false;
        if (args.length == 15){
            startingCells = Long.parseLong(args[0]);
            runtime = Long.parseLong(args[1]);
            mutationRate = Double.parseDouble(args[2]);
            subdivisionsCount = Integer.parseInt(args[3]);
            flowRate = Double.parseDouble(args[4]);
           
            
            resistantRate = Double.parseDouble(args[5]);
            sensitiveRate = Double.parseDouble(args[6]);
            trials = Integer.parseInt(args[7]);
            sequential = Boolean.parseBoolean(args[8]);
            smart = Boolean.parseBoolean(args[9]);
            switchEvery = Integer.parseInt(args[10]);
            matlab = Boolean.parseBoolean(args[11]);
            cutOff = Long.parseLong(args[12]);
            optimize = Boolean.parseBoolean(args[13]);
            flip = Boolean.parseBoolean(args[14]);
        } 
        
        System.out.println("starting cells: " + startingCells + " runtime: " + runtime + " mutation rate: " + mutationRate + " subdivisions count" + subdivisionsCount + " flowrate: " + flowRate + " resistant rate: " + resistantRate + " sensitive rate: " + sensitiveRate + " trials: " + trials + " sequential: " + sequential + " matlab: " + matlab + " cutoff: " + cutOff);
        if (optimize){
            int resistanceHigh = 0;
            int resistanceLow = 0;
            int sensitiveHigh = 0;
            int sensitiveLow = 0;
            double resistanceTotal = 0.0;
            double resistanceMean = 378.0;
            boolean stop = false;
            int inARow = 0;
            double median = 0.0;
            double[] sensitives = new double[100];
            double[] resistives = new double[100];
            double[] progressions = new double[100];
                List<Double> divergences = new ArrayList<>();
                
                List<Double> sorted = new ArrayList<>();
            for (int i = 0; i < 100 && !stop; i++){
            
                Tuple t = run(startingCells, runtime, mutationRate, subdivisionsCount, flowRate, genes,resistantRate, sensitiveRate, sequential, smart, switchEvery, matlab, cutOff, optimize);
                resistanceTotal += t.doublingTime;
                resistanceMean = resistanceMean * .8 + t.doublingTime * .2;//resistanceTotal / (i + 1);
                divergences.add(t.halvingTime);
                sensitives[i] = sensitiveRate;
                progressions[i] = t.halvingTime;
                resistives[i] = resistantRate;
                if(divergences.size() > 5){
                    //only use 5 most recent to calculate median
                    divergences.remove(0);
                }
                sorted.clear();
                sorted.addAll(divergences);
                Collections.sort(sorted);
                median = sorted.get(sorted.size() / 2);
                double targetMean = 279.0;
                double targetMedian = 720.0;
                if (Math.abs(resistanceMean - targetMean)/targetMean > 0.1){// 
                    if (resistanceMean > targetMean){
                        resistantRate = 1.0 + (resistantRate - 1.0) * (1.0 + .5/ (1.0+i/5.0));
                        resistanceLow++;
                    } else {
                        resistantRate = 1.0 + (resistantRate - 1.0) / (1.0 + .5/ (1.0+i/5.0));
                        resistanceHigh++;
                    }
                    
                }
                if (Math.abs(median - targetMedian) / targetMedian > .1){
                    if (median < targetMedian ^ flip) {
                        sensitiveRate = Math.max(.75, 1.0 + (sensitiveRate - 1.0) * (1.0 + .5/ (1.0+i/5.0)));
                        sensitiveHigh++;
                    } else {
                        sensitiveRate = Math.max(.75, 1.0 + (sensitiveRate - 1.0) / (1.0 + .5/ (1.0+i/5.0)));
                        sensitiveLow++;
                    }
                }
                
                
                if (i > 30 && Math.abs(resistanceMean - targetMean)/targetMean < 0.1  && Math.abs(median - targetMedian) / targetMedian < 0.1){
                    inARow++;
                    if (inARow == 5){
                        stop = true;
                    }
                } else {
                    inARow = Math.max(0, inARow - 1);
                }
                
                            
            }
                //write the results to a file
                File file = new File("flowRate=" + flowRate + "flip="+flip+".csv");  
                if ( !file.exists() )
                    file.createNewFile();
                FileWriter fw = new FileWriter(file);
                Writer writer = new BufferedWriter( fw );

                writer.append(resistantRate + "," + sensitiveRate+"\n");
                writer.append("Average Doubling time, " + resistanceMean + ", Median time to progression, "+ median);
                writer.flush();
                writer.close();
                fw.close();
        
            } else {
            
                //read in rates
            
                
                List<String> lines = Files.readAllLines(Paths.get("flowRate=" + flowRate + ".csv"), StandardCharsets.UTF_8);
                String[] rates = lines.get(0).split(",");
                resistantRate = Double.parseDouble(rates[0]);
                sensitiveRate = Double.parseDouble(rates[1]);
            
                List<Double> divergences = new ArrayList<>();
                for (int i = 0; i < trials; i++){

                long startTime = System.currentTimeMillis();
                Tuple t;
                do {
                t = run(startingCells, runtime, mutationRate, subdivisionsCount, flowRate, genes,resistantRate, sensitiveRate, sequential, smart, switchEvery, matlab, cutOff, optimize);
                divergences.add(t.doublingTime);
                } while (t.doublingTime == 0.0);

                System.out.println(sequential + " in " + (System.currentTimeMillis() - startTime) + " milliseconds.");
                }
                Collections.sort(divergences);
                double median = divergences.get(divergences.size() / 2);


                String name = resistantRate + "," + sensitiveRate + "rates,flowRate=" + flowRate + "trialsSequential.csv";

                if (!sequential){
                    name = resistantRate + "," + sensitiveRate + "rates,flowRate=" + flowRate + "trialsAlternating.csv";
                } else if (smart){
                    name = resistantRate + "," + sensitiveRate + "rates,flowRate=" + flowRate + "trialsSmart.csv";
                }
                // File (or directory) with old name
                File file = new File(name);

                // File (or directory) with new name
                File file2 = new File("median=" + median + "," + name);

                if (file2.exists())
                   throw new java.io.IOException("file exists");

                // Rename file (or directory)
                boolean success = file.renameTo(file2);




        //window.setVisible(true);
        //window.showBar();
        }
    }


}
