package cancersimulation;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


public class CellStrain {

    //each successive index is the next cellStrain
    //so the divisionRate of the original cellStrain is in index 0
    //and the divisionRate of the current cellStrain is in index length-1
    List<Double> divisionRateHistory = new ArrayList<>();
    int mutationCount;
    long population;
    boolean[] genome;
    int genomeID;
    boolean processed = false;
    
    //position coordinates
    int x;
    int y;
    int z;
    
    //the generation that each cellStrain split off from the parent strain
    List<Integer> geneTimes = new ArrayList<>();
    //contains any children the strain might have
    Set<CellStrain> children = new HashSet<>();

    //default constructor
    public CellStrain(List divisionRateHistory, long population, boolean[] genome, List geneTimes, int mutationCount, int x, int y, int z) {
        this.divisionRateHistory = divisionRateHistory;
        this.population = population;
        this.genome = genome;
        this.geneTimes = geneTimes;
        this.mutationCount = mutationCount;
        this.x = x;
        this.y = y;
        this.z = z;
    }
    
}
