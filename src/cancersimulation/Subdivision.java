/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cancersimulation;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 *
 * @author malcolmcallis
 */


public class Subdivision {
    
    List<Subdivision> neighbors = new ArrayList<>();
    List<CellStrain> strains = new ArrayList<>();
    double populationScalingFactor = 1;
    
    //position in compartment grid
    int x;
    int y;
    int z;
    
    public Subdivision(int x, int y, int z){
        neighbors.clear();
        strains.clear();
        
        this.x = x;
        this.y = y;
            this.z = z;
        
    }
    
    public void add(CellStrain strain){
        strains.add(strain);
    }
    
}
