/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cancersimulation;

import java.util.Comparator;
/**
 *
 * @author Malcolm
 */
public class gene implements Comparator<gene>, Comparable<gene>{
    int index;
    double mutationFrequency = 0.0;
    
    public gene(int index){
        this.index = index;
    }
   
    public int compareTo(gene g){
        int ret = 0;
        if (this.mutationFrequency > g.mutationFrequency){
            ret =  -1;
        }
        if (this.mutationFrequency < g.mutationFrequency){
            ret =  1;
        }
        
        return ret;
    }
    
    
    public int compare(gene g1, gene g2){
        int ret = 0;
        if (g1.mutationFrequency > g2.mutationFrequency){
            ret =  -1;
        }
        if (g1.mutationFrequency < g2.mutationFrequency){
            ret =  1;
        }
        
        return ret;
    }
    
    
}
