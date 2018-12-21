"""Algorithm for processing the data of truncation into a valuable formula for predicting activity."""

import numpy, re
import matplotlib
import random, operator
import math

class EvoAlg:

    def __init__(self,filename):
        self.generations = 100
        self.mutation_rate = 0.01  # value must not be bigger than 0.5
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
        if bool(random.getrandbits(1)):
            org1 = self.splice(org1)
        if bool(random.getrandbits(1)):
            org2 = self.splice(org2)
        self.combine([org1,org2])



    # This function uses random "epigenetic" markers and operators to create an organism
    def combine(self, genes):
        genome = str()
        org_size = random.randint(1,len(genes))  # the size of the organism is generated randomly
        for i in range(org_size):
            if org_size == 2:
                gid = i
            else:
                gid = random.randint(0,len(genes))
            a = str(round(random.uniform(0.01,4), 3))
            b = str(round(random.uniform(0.01,4), 3))
            epigene = "(" + a + "*" + "self.gene_pool[" + str(gid) + "]" + "**" + b + ")"
            myoperator = self.operators[random.randint(0,3)]
            genome += epigene + myoperator
        self.mutate(genome)
        print(genome)
        #print(eval(genome[:-1]))


    def mutate(self,genome):
        mutation_key = random.randint(0,1/self.mutation_rate)
        # get the locations of the numbers:
        numbers = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",genome)
        number_locs = list()
        for number in numbers:
            number_locs.append(genome.find(number))
        # Go through every character in the genome and figure out if it needs to be mutated
        for i in range(len(genome)):
            # skip the parentheses:
            if genome[i] == "(" or genome[i] == ")":
                continue
            mutator = random.randint(0,1/self.mutation_rate)
            if mutator == mutation_key:
                print("Mutation created.")
                # mutate an operator
                if self.operators.find(genome[i]) != -1:
                    if genome[i+1] != "*" or genome[i-1] != "*":  # make sure its not an exponent
                        genome = genome[i:].replace(genome[i],self.operators[random.randint(0,3)])
                # mutate a number
                elif i in number_locs:
                    num_mut = random.uniform(-2.0,2.0)
                    # how to see if you landed on the number in general not necessarily the first digit?






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





