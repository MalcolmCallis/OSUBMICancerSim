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
public class sample implements Comparator<sample>, Comparable<sample>{
    int index;
    double distance = 0.0;
    
    public sample(int index, double distance){
        this.index = index;
        this.distance = distance;
    }
   
    public int compareTo(sample g){
        int ret = 0;
        if (this.distance > g.distance){
            ret =  -1;
        }
        if (this.distance < g.distance){
            ret =  1;
        }
        
        return ret;
    }
    
    
    public int compare(sample g1, sample g2){
        int ret = 0;
        if (g1.distance > g2.distance){
            ret =  -1;
        }
        if (g1.distance < g2.distance){
            ret =  1;
        }
        
        return ret;
    }
    
    
}
