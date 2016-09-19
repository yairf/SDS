#!/home/yairf/Develop/python/loopy/bin/python

import os, sys, re, math
import argparse
import numpy as np

import simuOpt
simuOpt.setOptions(quiet=True)

import simuPOP as sim
import math
from simuPOP.utils import *


#-----------------#
# Parse arguments #
#-----------------#

parser = argparse.ArgumentParser(description="Backward simulation of selection using the simuPOP package. The input includes the demographic model, the current derived allele frequency, the selection coefficient, and the backward time window at which selection was active.")


parser.add_argument("-g", "--number_of_generations", 
                    type=int,
                    default=10000,
                    help="The number of generations simulated \
                          (default = %(default)s).")

parser.add_argument("-te", "--selection_end_time", 
                    type=int,
                    default=0,
                    help="The generation (backward in time) at which selection ended \
                          (default = %(default)s, i.e., current time).")

parser.add_argument("-ts", "--selection_start_time", 
                    type=int,
                    default=100,
                    help="The generation (backward in time) at which selection started \
                          (default = %(default)s).")

parser.add_argument("-ps", "--selection_sigma_paramter",
                    type=float,
                    default=0.0,
                    help="The selection \"sigma\" paramter \
                          (default = %(default)s).")

parser.add_argument("-ph", "--selection_dominance_paramter",
                    type=float,
                    default=0.5,
                    help="The selection h paramter \
                          (default = %(default)s).")

parser.add_argument("-f", "--current_snp_frequency",
                    type=float,
                    default=0.6,
                    help="The frequency of the snp at current time \
                          (default = %(default)s).")

parser.add_argument("-sNe", "--scale_Ne", 
                    type=float,
                    default=1.0,
                    help="A scale factor for Ne, the effective population_size, \
                          (default = %(default)s).")


parser.add_argument("-m", "--population_size_model",
                    choices=["Tennessen_CEU","Gravel_AFR","Gravel_CEU","Gravel_CHB","Const_1E4"],
                    default="Tennessen_CEU",
                    help="Define the population size history model. \
                          (default = %(default)s).")



args = parser.parse_args()

#-------------------------------
#
#-------------------------------

OUTPUT = sys.stdout

num_of_genereations = args.number_of_generations-1


#

def genBackwards(gen_forward):

    return int(abs(gen_forward-(num_of_genereations+1)))

#

def Nt(gen):


    if gen==0:
        return 1

    tmp_baseNt = baseNt(gen)
    return int(args.scale_Ne * tmp_baseNt)


def baseNt(gen):

    if args.population_size_model == "Tennessen_CEU":
        return baseNt_Tennessen_CEU(gen)
    elif args.population_size_model == "Gravel_AFR":
        return baseNt_Gravel_AFR(gen)
    elif args.population_size_model == "Gravel_CEU":
        return baseNt_Gravel_CEU(gen)
    elif args.population_size_model == "Gravel_CHB":
        return baseNt_Gravel_CHB(gen)
    elif args.population_size_model == "Const_1E4":
        return baseNt_Const_1E4(gen)


def baseNt_Tennessen_CEU(gen):

    my_gen_backwards = genBackwards(gen)
    if my_gen_backwards > 5920:
        return 7310
    if my_gen_backwards > 2040:
        return 14474
    if my_gen_backwards > 920:
        return 1861
    if my_gen_backwards == 920:
        return 1032
    if my_gen_backwards >= 205:
        return int(1032*(1.0*9300/1032)**(1.0*(920-my_gen_backwards)/(920-205)))
    if my_gen_backwards >= 1:
        return int(9300*(1.0*512000/9300)**(1.0*(205-my_gen_backwards)/(205-1)))


def baseNt_Gravel_AFR(gen):

    my_gen_backwards = genBackwards(gen)
    if my_gen_backwards > 5920:
        return 7310
    if my_gen_backwards >= 1:
        return 14474


def baseNt_Gravel_CEU(gen):

    my_gen_backwards = genBackwards(gen)
    if my_gen_backwards > 5920:
        return 7310
    if my_gen_backwards > 2040:
        return 14474
    if my_gen_backwards > 920:
        return 1861
    if my_gen_backwards == 920:
        return 1032
    if my_gen_backwards >= 1:
        return int(1032*(1.0038)**(1.0*(920-my_gen_backwards)))


def baseNt_Gravel_CHB(gen):

    my_gen_backwards = genBackwards(gen)
    if my_gen_backwards > 5920:
        return 7310
    if my_gen_backwards > 2040:
        return 14474
    if my_gen_backwards > 920:
        return 1861
    if my_gen_backwards == 920:
        return 550
    if my_gen_backwards >= 1:
        return int(550*(1.0048)**(1.0*(920-my_gen_backwards)))


def baseNt_Const_1E4(gen):

    return 10000



#

def Ft(gen, subPop):

    my_gen_backwards = genBackwards(gen)

    if my_gen_backwards < args.selection_end_time:
        return (1.0, 1.0, 1.0)
    if my_gen_backwards > args.selection_start_time:
        return (1.0, 1.0, 1.0)
    else:
        #older version had a bug... fixed: Sep 4, 2015
        #return (1.0, 1.0+(2.0*args.selection_sigma_paramter*args.selection_dominance_paramter), 1.0+(2.0*args.selection_sigma_paramter))
        return (1.0-args.selection_sigma_paramter, 1.0-(args.selection_sigma_paramter*args.selection_dominance_paramter), 1.0)



traj = simulateBackwardTrajectory(N=Nt, fitness=Ft, endGen=num_of_genereations, endFreq=args.current_snp_frequency)
#traj = simulateBackwardTrajectory(N=Nt, fitness=[1,1,1], endGen=num_of_genereations, endFreq=0.6)

#for i in range(0,1+num_of_genereations):
for i in range(num_of_genereations,0,-1):
    my_gen = str(genBackwards(i))
    my_popsize = str(Nt(i))
    my_fitness = Ft(i,0)
    my_fitness0 = str(my_fitness[0])
    my_fitness1 = str(my_fitness[1])
    my_fitness2 = str(my_fitness[2])
    my_freq = str(traj.freq(i,0)[0])
    OUTPUT.write( my_gen + "\t" + my_popsize + "\t" + my_fitness0 + "\t" + my_fitness1 + "\t" + my_fitness2 + "\t" + my_freq  +"\n")
