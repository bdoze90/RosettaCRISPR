"""Algorithm for processing the data of truncation into a valuable formula for predicting activity."""

import numpy
import matplotlib
import random, operator
import math

class EvoAlg:

    def __init__(self,filename):
        self.generations = 100
        self.mutation_rate = 0.01
        self.carrying_capacity = 50
        self.gene_pool = [(0.4,5,10,51,8),
                          (5,-0.25,57,21,2),
                          (51,1,24,246,9),
                          (21,12,0.051,24.01,8),
                          (1.1,0.369,2.1,-9.2,1),
                          (514,0.045,-8,3.1,44),
                          (21,4.45,-2.4,74.1,4)]# This is an example gene pool
        self.operators = "+-*/"

        # This is the current population of individuals.  The first value is the Pearson coefficient
        # and the second value is the genome.
        self.population = list()

        # This is the values that the algorithm is trying to replicate.
        self.ideal_values = [4,2,8,1.2,10]

        #self.import_data(filename)

    # This function is used when there is a large dataset
    def import_data(self,file):
        f = open(file)
        for line in f:
            line[:-1].split("\t")


    # This function initiates the population
    def spontaneous_generation(self):
        for i in range(self.carrying_capacity):
            self.combine(self.gene_pool)

    # This function is the main function in which the promising functions are mated together.  Those with higher
    # regression values get a higher chance of pieces of their genome getting passed along.
    def mating_season(self):
        # Get the top mates by sorting the population:
        self.population = sorted(self.population, key=operator.itemgetter[1])
        pop_cutoff = int(len(self.population)*0.1)
        for individual in self.population[:pop_cutoff]:
            partner = self.population[random.randint(0,pop_cutoff)]
            self.reproduction(individual,partner)

    def reproduction(self,org1, org2):




    # This function uses random "epigenetic" markers and operators to create an organism
    def combine(self, genes):
        genome = str()
        for i in range(len(genes)):
            a = str(round(random.uniform(0.01,4), 3))
            b = str(round(random.uniform(0.01,4), 3))
            epigene = "(" + a + "*" + "self.gene_pool[" + str(i) + "]" + "**" + b + ")"
            operator = self.operators[random.randint(0,3)]
            genome += epigene + operator
        print(genome)
        print(eval(genome[:-1]))


    def splice(self,genome):
        splice_loc = genome.find(self.operators[random.randint(0,3)])  #picks a random operator and uses it as splice
        first_half = genome[:splice_loc]
        second_half = genome[splice_loc:]
        if bool(random.getrandbits(1)):
            return first_half
        else:
            return second_half

E = EvoAlg("/Users/brianmendoza/Desktop/EvoLearningData.txt")
E.spontaneous_generation()





