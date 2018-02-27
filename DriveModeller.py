#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      I am
#
# Created:     27/02/2018
# Copyright:   (c) I am 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import numpy
import scipy
from matplotlib import pyplot



def ProbOddCrossovers(d):
    return (1-numpy.exp(-2*d/100.))/2


def ConvertPositionToGeneticDistance(pos, chrom):
    chrom_coef={'2L':(-.01, .2, 2.59,-1.56),
            '2R':(-.007, .35, -1.43, 56.91),
            '3L':(-.006, .092, 2.94, -2.9),
            '3R': (.004,.24, -1.63,50.26),
            'X':(-.01,.3,1.15, -1.87)}

    chrom_bounds={'2L':numpy.array([.53, 18.87]),
            '2R':numpy.array([1.87, 20.86]),
            '3L':numpy.array([.75, 19.02]),
            '3R': numpy.array([2.58,27.44]),
            'X':numpy.array(([1.22, 21.21]))}

    genetic_distance= Polynomial(pos, *chrom_coef[chrom])

    #Set meiotic recombination rate to zero in heterochromatin
    min_dist, max_dist =Polynomial(chrom_bounds[chrom], *chrom_coef[chrom])
    genetic_distance[pos<chrom_bounds[chrom][0]]=min_dist
    genetic_distance[genetic_distance<0]=0
    genetic_distance[pos>chrom_bounds[chrom][1]]=max_dist

    return genetic_distance



def Polynomial(x,a,b,c,d):
    return a*x**3+b*x**2+c*x+d

def RecombinationDecay(distance, drive):
    prob_odd=ProbOddCrossovers(distance)
    return 1- (.5+.5*((drive*(1-prob_odd))+(1-drive)*prob_odd))

def main():
    pass

if __name__ == '__main__':
    main()


