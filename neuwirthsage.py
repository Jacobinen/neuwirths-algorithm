import itertools
import neuwirth as nw
import sage.all # Needed to avoid circular dependencies
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.combinat.permutation import Permutation
#from sage.rings.integer import Integer

# Import SAGE tools here without circular conflict

def stream_permutations_mod_cyclic(lists):
        # Create a generator for the Cartesian product of adjusted permutations
        for perm in itertools.product(*([[lst[0]] + list(p) for p in itertools.permutations(lst[1:])] for lst in lists)):
            yield perm

def tcycles(P):
    #returns the number of orbits/transitivity cycles of a permutation
    return len(P.cycle_type())

def count_orbits(perm1, perm2):
    # Create the permutation group from the two permutations
    G = PermutationGroup([perm1, perm2])
    
    # Create a set to keep track of the orbits
    orbits = G.orbits()

    # Return the number of orbits
    return len(orbits)

def count_orbits_3(perm1, perm2, perm3):
    # Create the permutation group from the three permutations
    G = PermutationGroup([perm1, perm2, perm3])
    
    # Create a set to keep track of the orbits
    orbits = G.orbits()

    # Return the number of orbits
    return len(orbits)

def L_One_Case(A,B,C):
    #first we check if |A,C|=1, i.e. that L=1   
    assert count_orbits(A,C)==1

    # Compute the permutations AC and BC, respectively. Note: C is applied last.
    AC=A.right_action_product(C); BC=B.right_action_product(C)


    # Check the two conditions in the L=1 case of Neuwirth's algorithm
    if count_orbits(AC,BC)==1 and tcycles(A)-tcycles(C)+2==tcycles(AC):
        return True

    # If it is not true, then we return false
    return False

def construct_F(A,C):
    G = PermutationGroup([A, C])

    # We start writing down the only non-trivial cycle in F
    F_build='('

    for orbit in G.orbits():
        # Pick out the largest integer in each cycle of (A,C) - this will guarantee the correct size of F
        F_build+=str(orbit[-1])+','

    # Delete the extra comma and close the bracket
    F_build = F_build[:-1]+')'

    return Permutation(F_build)

def L_More_Than_One_Case(A,B,C):
    # Compute L
    L=count_orbits(A,C)

    # Compute the permutations AC and BC, respectively. Note: C is applied last.
    AC=A.right_action_product(C); BC=B.right_action_product(C)

    # Build the additional permutation F
    F=construct_F(A,C)

    # Check the two conditions in the L>1 case of Neuwirth's algorithm
    if count_orbits_3(AC,BC,F)==1 and tcycles(A)-tcycles(C)+2*L==tcycles(AC):
        return True

    # If it is not true, then we return false
    return False    

def bar_reflection(reordering):
    # bar_reflection might be a bad name for this function...
    
    # Takes a list of integers as input, i.e. value table
    # Check the length to know how much we need to shift
    n=len(reordering)
    
    # Return the input concatenated with a shifted inverse of itself
    # This is a value table
    return reordering + [n+e for e in list(Permutation(reordering))] # NEED TO REMOVE .inverse()

def first_half(lst):
    # Make sure that the list has an even number of entries
    assert len(lst)%2==0
    
    # Return the first half of the list
    return lst[:len(lst)//2]

def exists_2_spine(A,B,C):
    # We make a list of the unbarred cycles of C (a partition of the unbarred elements)
    PC=[list(p) for p in first_half(C.to_cycles())]

    # Loop through all permutations that stay in their partition blocks - now we mod out cyclicity
    for config in stream_permutations_mod_cyclic(PC):
        # In value table form - removing "inner brackets" of the tuple and making it into a list
        v_tab=[e for block in config for e in block]

        # Extend to a permutation that is compatible on the unbarred and barred elements
        re_ordering=bar_reflection(v_tab)

        # Use the permutation class in SAGE
        p=Permutation(re_ordering)

        # Check if this particular datum yields a thickening to a 3-manifold
        # We are conjugating A with p to get another admissible A
        
        try: 
            if L_One_Case(p.right_action_product(A).right_action_product(p.inverse()),B,C):
                # Print the result
                print('THE FOLLOWING CONFIGURATION YIELDS A SPINE:\n')

                # Print which reordering we used and show it
                print('Reordering permutation: ', p)
                p.show()

                # Print the involution A we used
                print('Involution A: ', p.right_action_product(A).right_action_product(p.inverse()))
                p.right_action_product(A).right_action_product(p.inverse()).show()

                return True

        except:
            if L_More_Than_One_Case(p.right_action_product(A).right_action_product(p.inverse()),B,C):
                # Print the result
                print('THE FOLLOWING CONFIGURATION YIELDS A SPINE:\n')

                # Print which reordering we used and show it
                print('Reordering permutation: ', p)
                p.show()

                # Print the involution A we used
                print('Involution A: ', p.right_action_product(A).right_action_product(p.inverse()))
                p.right_action_product(A).right_action_product(p.inverse()).show()

                return True

    # If no configuration works then there is no thickening to a 3-manifold
    return False

def check_presentation(gen, rel):
    """We take a balanced presentation as input. Ex: gen=[['a','A'],['b','B']] where A:=a^-1,
    rel=['aaBBB','ABAbab']. Need the inputs to be arrays of strings."""

    prev = nw.prevalences(gen, rel)

    A = Permutation(nw.admissible(gen,rel))
    B = Permutation(nw.B(sum(prev)))
    C = Permutation(nw.C(prev,len(prev)))

    return exists_2_spine(A,B,C)


def main():
    #A=Permutation('(12,2)(13,17)(6,18)(7,19)(8,1)(3,20)(9,15)(4,10)(21,5)(16,11)(22,14)')
    #B=Permutation('(1,12)(2,13)(3,14)(4,15)(5,16)(6,17)(7,18)(8,19)(9,20)(10,21)(11,22)')
    #C=Permutation('(1,2,3,4,5)(16,15,14,13,12)(6,7,8,9,10,11)(22,21,20,19,18,17)')

    print(check_presentation([['a','A'],['b','B']],['aaBBB','ABAbab']))
    print(check_presentation([['a','A'],['b','B']],['aBBBa', 'BAbabA']))
    print(check_presentation([['a','A'],['b','B']],['aBBBa', 'BAbaBBa']))
    print(check_presentation([['a','A'],['b','B']],['aBBBa', 'AbaBBaB']))
    print(check_presentation([['a','A'],['b','B']],['aBBaBBaB', 'AbaBBaB']))
    print(check_presentation([['a','A'],['b','B']],['aBBaBBaB', 'BabAbbA']))
    print(check_presentation([['a','A'],['b','B']],['aBBaBBaB', 'BaBaB']))
    print(check_presentation([['a','A'],['b','B']],['BBaBBaBa', 'AbAbb']))
    print(check_presentation([['a','A'],['b','B']],['BBa', 'AbAbb']))
    print(check_presentation([['a','A'],['b','B']],['BBa', 'Ab']))
    print(check_presentation([['a','A'],['b','B']],['B', 'Ab']))
    print(check_presentation([['a','A'],['b','B']],['B', 'A']))




if __name__=='__main__':
    main()



