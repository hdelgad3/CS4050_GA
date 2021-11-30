"""
For this assignment there is no automated testing. You will instead submit
your *.py file in Canvas. I will download and test your program from Canvas.
Name: Hector Delgado
"""

import time
import random
from queue import PriorityQueue


def adjMatFromFile(filename):
    """ Create an adj/weight matrix from a file with verts, neighbors, and weights. """
    f = open(filename, "r")
    n_verts = int(f.readline())
    print(f" n_verts = {n_verts}")
    adjmat = [[None] * n_verts for i in range(n_verts)]
    for i in range(n_verts):
        adjmat[i][i] = 0
    for line in f:
        int_list = [int(i) for i in line.split()]
        vert = int_list.pop(0)
        assert len(int_list) % 2 == 0
        n_neighbors = len(int_list) // 2
        neighbors = [int_list[n] for n in range(0, len(int_list), 2)]
        distances = [int_list[d] for d in range(1, len(int_list), 2)]
        for i in range(n_neighbors):
            adjmat[vert][neighbors[i]] = distances[i]
    f.close()
    return adjmat


def TSPwGenAlgo(g, max_num_generations=500, population_size=100,
        mutation_rate=0.01, explore_rate=0.6):
    """ A genetic algorithm to attempt to find an optimal solution to TSP  """

    def perform_mutation(arr):
        first = int(random.random() * len(arr))
        second = int(random.random() * len(arr))
        arr[first], arr[second] = arr[second], arr[first]
        return arr

    def get_path_sum(arr)-> int:
        fir = 0
        sec = 1
        total = 0
        while sec < len(arr):
            a = arr[fir]
            b = arr[sec]
            cost = g[a][b]
            total += cost
            fir += 1
            sec += 1
        return total

    def crossover(parent1, parent2):
        # print("Parent1: ",parent1)
        # print("Parent2: ",parent2)
        print("parent1 cost: ", get_path_sum(parent1))
        print("parent2 cost: ", get_path_sum(parent2))
        ch1_start = []
        ch2_start = []
        child1 = []
        child2 = []

        first_random = random.random()
        second_random = random.random()
        gene1 = int( first_random * len(parent1) * 0.5)
        gene2 = int(second_random * len(parent1))
        start = min(gene1, gene2)
        end = max(gene1, gene2)
        for i in range(start, end):
            ch1_start.append(parent1[i])
        ch1_end = [chrome for chrome in parent2 if chrome not in ch1_start]
        child1 = ch1_start + ch1_end

        gene3 = int(first_random * len(parent2) * 0.5)
        gene4 = int(second_random * len(parent2))
        start = min(gene3, gene4)
        end = max(gene3, gene4)
        for i in range(start, end):
            ch2_start.append(parent2[i])
        ch2_end = [chrome for chrome in parent1 if chrome not in ch2_start]
        child2 = ch2_start + ch2_end
        print("child1 cost: ", get_path_sum(child1))
        print("child2 cost: ", get_path_sum(child2))


        return child1, child2


    solution_cycle_distance = None # the distance of the final solution cycle/path
    solution_cycle_path = [] # the sequence of vertices representing final sol path to be returned
    shortest_path_each_generation = [] # store shortest path found in each generation

    # create individual members of the population
    for i in g:
        print(i)

    population = []
    alphabet = [i for i in range(len(g))]
    for r in range(population_size):
        random.shuffle(alphabet)
        population.append(list(alphabet))
    # print(population, "population")

    solutions = []
    # initialize individuals to an initial 'solution'
    for ind in population:
        copy = ind.copy()
        copy.append(copy[0])
        solutions.append(copy)
    print(solutions)



    # loop for x number of generations (with possibly other early-stopping criteria)
    for x in range(max_num_generations):
        # calculate fitness of each individual in the population
        fitness_sums = []
        sorted_ind = []

        for ind in solutions:
            ind_cost = get_path_sum(ind)
            fitness_sums.append(ind_cost)
            sorted_ind.append((ind_cost,ind))
        print(fitness_sums)
        sorted_ind.sort(key = lambda x: x[0])
        print(sorted_ind)


        # (and append distance of the 'fittest' to shortest_path_each_generation)
        shortest_path_each_generation.append(sorted_ind[0])

        # select the individuals to be used to spawn the generation, then create
        # individuals of the new generation (using some form of crossover)
        new_generation = []
        sorted_population = []
        for i in sorted_ind:
            sorted_population.append(i[1])

        couples = sorted_population[:int(len(sorted_population) * (1-explore_rate))]

        while couples:
            parent1 = couples.pop(0)
            parent2 = couples.pop(0)
            if parent2 == None:
                parent2 = parent1.reverse()
            child1, child2 = crossover(parent1, parent2)
            new_generation.append(child1)
            new_generation.append(child2)
        print(new_generation)

        # allow for mutations (this should not happen too often)



        # ...
    # calculate and verify final solution, and update solution_cycle_distance,
    # solution_path, etc.

    # ...

    return {
            'solution': solution_cycle_path,
            'solution_distance': solution_cycle_distance,
            'evolution': shortest_path_each_generation
           }


def TSPwDynProg(g):
    """ (10pts extra credit) A dynamic programming approach to solve TSP """
    solution_cycle_distance = None # the distance of the final solution cycle/path
    solution_cycle_path = [] # the sequence of vertices representing final sol path to be returned

    #...

    return {
            'solution': solution_cycle_path,
            'solution_distance': solution_cycle_distance,
           }


def TSPwBandB(g):
    """ (10pts extra credit) A branch and bound approach to solve TSP """
    solution_cycle_distance = None # the distance of the final solution cycle/path
    solution_cycle_path = [] # the sequence of vertices representing final sol path to be returned

    #...

    return {
            'solution': solution_cycle_path,
            'solution_distance': solution_cycle_distance,
           }


def assign05_main():
    """ Load the graph (change the filename when you're ready to test larger ones) """
    g = adjMatFromFile("complete_graph_n08.txt")

    # Run genetic algorithm to find best solution possible
    start_time = time.time()
    res_ga = TSPwGenAlgo(g)
    elapsed_time_ga = time.time() - start_time
    print(f"GenAlgo runtime: {elapsed_time_ga:.2f}")
    print(f"  sol dist: {res_ga['solution_distance']}")
    print(f"  sol path: {res_ga['solution']}")

    # (Try to) run Dynamic Programming algorithm only when n_verts <= 10
    if len(g) <= 10:
        start_time = time.time()
        res_dyn_prog = TSPwDynProg(g)
        elapsed_time = time.time() - start_time
        if len(res_dyn_prog['solution']) == len(g) + 1:
            print(f"Dyn Prog runtime: {elapsed_time:.2f}")
            print(f"  sol dist: {res_dyn_prog['solution_distance']}")
            print(f"  sol path: {res_dyn_prog['solution']}")

    # (Try to) run Branch and Bound only when n_verts <= 10
    if len(g) <= 10:
        start_time = time.time()
        res_bnb = TSPwBandB(g)
        elapsed_time = time.time() - start_time
        if len(res_bnb['solution']) == len(g) + 1:
            print(f"Branch & Bound runtime: {elapsed_time:.2f}")
            print(f"  sol dist: {res_bnb['solution_distance']}")
            print(f"  sol path: {res_bnb['solution']}")


# Check if the program is being run directly (i.e. not being imported)
if __name__ == "__main__":
    assign05_main()

