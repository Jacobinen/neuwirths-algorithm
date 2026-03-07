#Want to build function that generates an admissible permutation A from a given (balanced) group presentation


#first, we start by defining a buch of convenient functions for the exact data format we will be using

import itertools
from sage.groups.perm_gps.permgroup import PermutationGroup, SymmetricGroup
from sage.groups.perm_gps.permgroup_element import Permutation
#from sage.rings.integer import Integer

def read_j(elem):
    return int(elem.split(',')[1][:-1])

def read_i(elem):
    return int(elem.split(',')[0].split('(')[1])

def bar_unbar(txt):
    if txt[1]=='_':
        return txt[0]+txt[2:]
    else:
        return txt[0]+'_'+txt[1:]

def i_index(symb, gens):
    """determines which index the generator has"""
    for i, elem in enumerate(gens):
        if symb in elem:
            return i+1

def first_j(A,i,bar=True):

    j=1

    if bar:
        while True:
            candidate='p_('+str(i)+','+str(j)+')'
            if candidate not in A and bar_unbar(candidate) not in A:
                return candidate
            j+=1

    j=1

    while True:
        candidate='p('+str(i)+','+str(j)+')'
        if candidate not in A and bar_unbar(candidate) not in A:
            return candidate
        j+=1

def translate_to_integer(txt, prevalences):

    barred=False
    if txt[1]=='_':
        barred=True

    number=sum(prevalences[:read_i(txt)-1])+read_j(txt)

    if barred:
        return number+sum(prevalences)
    return number


def to_e_notation(txt):
    if txt[1]=='_':
        return '\\bar{e}_{'+str(read_i(txt))+'}^{'+str(read_j(txt))+'}'
    else:
        return 'e_{'+str(read_i(txt))+'}^{'+str(read_j(txt))+'}'



def prevalences(gen, rel):
    prevalences=[]

    for a in gen:
        count=0
        for r in rel:
            for e in r:
                if e==a[0] or e==a[1]:
                    count+=1
        prevalences.append(count)

    return prevalences

def to_p_notation(num, prevalences):
    # Count how many unbarred elements of E there are
    n=sum(prevalences)
    
    if num>n:
        bar=True
        # Shift num by n
        num=num-n
    else:
        bar=False

    for k in range(len(prevalences)):
        if num<= sum(prevalences[:k+1]):
            i=k+1
            j=num-sum(prevalences[:k])

            if bar:
                return 'p_('+str(i)+','+str(j)+')'
            return 'p('+str(i)+','+str(j)+')'



def admissible(gen,rel):
    """We take a balanced presentation as input. Ex: gen=[['a','A'],['b','B']] where A:=a^-1,
    rel=['aaBBB','ABAbab']. Need the inputs to be arrays of strings."""

    #n=sum([len(r) for r in rel]) #what do we use this for?

    augmented_rel=[word+word[0] for word in rel] # preparation for cyclic looping

    capitals=[a[1] for a in gen] #need to be able to check if we have an inverse
    #non_capitals=[a[0] for a in gen] #actually a bit superfluous

    prevalences=[] #keep track of how many slots we need to have for the "orbits of C"

    for a in gen:
        count=0
        for r in rel:
            for e in r:
                if e==a[0] or e==a[1]:
                    count+=1
        prevalences.append(count)

    #print('Relators:', rel)

    #print('Prevalence of generators in relators (total):', prevalences)

    g=len(prevalences) #half the number of nodes in our Gamma graph... maybe not useful here

    #offset=sum(prevalences) #how much we need to shift the barred stuff

    E=[]

    for i in range(1,g+1):
        for j in range(1,prevalences[i-1]):
            E.append('p('+str(i)+','+str(j)+')') #first index corresponds to the generator (handle), second to the order of the sheets around it
            E.append('p_('+str(i)+','+str(j)+')')

    #print(E,'\n') #these are the vertices that we are dealing with -> e_i^j better nomenclature?



    A=[] #start boulding the involution (just list the elements and then it's implicit that they are paired up) WARNING!

    for r in augmented_rel:
        #here is where we start reading stuff off

        for e, (previous, current) in enumerate(zip(r, r[1:])):
            #print(previous, current)

            if previous in capitals:
                bar=False
            else:
                bar=True

            i=i_index(previous,gen)

            if e==0:
                final=bar_unbar(first_j(A,i,bar))
                A.append(first_j(A,i,bar))
            else:
                A.append(bar_unbar(A[-1]))



            if current in capitals:
                bar=True
            else:
                bar=False

            i=i_index(current,gen)

            if e==len(r)-2:
                A.append(final)
            else:
                A.append(first_j(A,i,bar))


    #print(A) #now we essentiall have the permutation. Todo: transform it into cycle notation with integers so that Mathematica can handle it.
        

    Output="A=Permutation('"

    for idx, entry in enumerate(A):
        if idx%2==0:
            if idx!=0:
                Output+=''
            Output+='('+str(translate_to_integer(entry, prevalences))+','
        else:
            Output+=str(translate_to_integer(entry, prevalences))+')'
    


    print(Output+"')")

    
    Output=''

    for idx, entry in enumerate(A):
        if idx%2==0:
            Output+='('+to_e_notation(entry)
        else:
            Output+=to_e_notation(entry)+')'

    print(Output)
    
    """

    




    Pseudocode:

    loop through pairs

    aa->e_bar e
    aA->e_bar e_bar
    Aa->e e
    AA->e e_bar

    take "first" available pair in E that fits

    augment to list of transpositions in A

    remember which one we need to start the next transposition with

    when we are at the end of each relator, remember to "tie back" to the start

    """


def B(g):
    """This writes out a mathematica compatible cycle notation for the outer arc involutions 
    where we have 2g symbols in E"""

    Output=''

    for i in range(g):
        if i==0:
            Output+='('+str(i+1)+','+str(i+g+1)+')'
        else:
            Output+='('+str(i+1)+','+str(i+g+1)+')'


    return Output

def C(prevalences, nr_gens):
    """The barred and unbarred stuff has to be compatible"""

    offset=sum(prevalences)

    Output=''

    for i in range(nr_gens):
        start=sum(prevalences[:i])

        #if i!=0:
        #    Output+=','

        for j in range(prevalences[i]):
            if j==0:
                Output+='('+str(start+j+1)
            else:
                Output+=','+str(start+j+1)

        Output+=')'

        for j in range(prevalences[i]):
            if j==0:
                Output+='('+str(start+offset+prevalences[i]-j)
            else:
                Output+=','+str(start+offset+prevalences[i]-j)
        
        Output+=')'


    return Output


    #for e, handle in enumerate(prevalences):
        #first the unbarred stuff

    #    Output 



        #then the barred stuff

    # number=sum(prevalences[:read_i(txt)-1])+read_j(txt)




def stream_permutations(lists):
    # Create a generator for the Cartesian product of permutations for each list
    for perm in itertools.product(*(itertools.permutations(lst) for lst in lists)):
        yield perm

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
    AC=A.right_action_product(C)
    BC=B.right_action_product(C)


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
    AC=A.right_action_product(C)
    BC=B.right_action_product(C)

    # Build the additional permutation F
    F=construct_F(A,C)

    # Check the two conditions in the L>1 case of Neuwirth's algorithm
    if count_orbits_3(AC,BC,F)==1 and tcycles(A)-tcycles(C)+2*L==tcycles(AC):
        return True

    # If it is not true, then we return false
    return False   

def AllPermutationsOfPartition(partition, items=[]):
    # WE DO NOT WANT TO USE THIS - SAVING IN A LIST REQUIRES TOO MUCH RAM SOMETIMES
    """
    Generate all permutations of a given partition [x_1,x_2,...,x_g].
    
    Args:
    partition: A list of integers representing the sizes of parts in the partition.
    
    Returns:
    A list of all distinct permutations of the partition.
    """
    #from sage.groups.perm import Permutation, PermutationGroup
    #from sage.all import symmetric_group

    #X=[list(range(sum(partition[:idx])+1, sum(partition[:idx])+n_i+1)) for idx,n_i in enumerate(partition)]

    if partition==[]:
        return [list(e) for e in items]

    else:
        collection=[]
        x_top=partition.pop()
        S_n=SymmetricGroup(len(x_top))
        for perm in S_n:
            if items==[]:
                collection.append(tuple(perm(x_top)))
            else:
                for item in items:
                    collection.append(tuple(perm(x_top))+item)
        return AllPermutationsOfPartition(partition, collection)
    
def bar_reflection(reordering):
    # bar_reflection might be a bad name for this function...
    
    # Takes a list of integers as input, i.e. value table
    # Check the length to know how much we need to shift
    n=len(reordering)
    
    # Return the input concatenated with a shifted inverse of itself
    # This is a value table
    return reordering + [n+e for e in list(Permutation(reordering))]

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
                print('Reordering permutation: ', p);p.show()

                # Print the involution A we used
                print('Involution A: ', p.right_action_product(A).right_action_product(p.inverse()))
                p.right_action_product(A).right_action_product(p.inverse()).show()

                return True

        except:
            if L_More_Than_One_Case(p.right_action_product(A).right_action_product(p.inverse()),B,C):
                # Print the result
                print('THE FOLLOWING CONFIGURATION YIELDS A SPINE:\n')

                # Print which reordering we used and show it
                print('Reordering permutation: ', p);p.show()

                # Print the involution A we used
                print('Involution A: ', p.right_action_product(A).right_action_product(p.inverse()))
                p.right_action_product(A).right_action_product(p.inverse()).show()

                return True

    # If no configuration works then there is no thickening to a 3-manifold
    return False

def check(A,B,C):
    # This may require too much RAM
    
    PartCycles=[list(p) for p in first_half(C.to_cycles())] #this might be problematic if we have cycles of len=1 -> (1,),...
    Parts=AllPermutationsOfPartition(PartCycles)
    reordering=[bar_reflection(half) for half in Parts] #admissible re-orderings

    for perm in reordering:
        p=Permutation(perm)

        if L_One_Case(p.right_action_product(A).right_action_product(p.inverse()),B,C):
            print('THE FOLLOWING CONFIGURATION YIELDS A SPINE:\n')
            print('Reordering permutation: ', p);p.show()
            print('Involution A: ', p.right_action_product(A).right_action_product(p.inverse()))
            p.right_action_product(A).right_action_product(p.inverse()).show()

            return True
            
    return False

def copy_paste_to_nb(gen, rel):
    info="#("+rel[0]+","+rel[1]+")"
    print(info)
    admissible(gen, rel)
    prev=prevalences(gen,rel)
    print("B=Permutation('"+B(sum(prev))+"')")
    print("C=Permutation('"+C(prev,len(prev))+"')")




def main():
    #(abAA,BBB) - dummy test
    A=Permutation('(8,4)(11,9)(2,10)(3,1)(5,13)(6,14)(7,12)')
    B=Permutation('(1,8)(2,9)(3,10)(4,11)(5,12)(6,13)(7,14)')
    C=Permutation('(1,2,3)(10,9,8)(4,5,6,7)(14,13,12,11)')

    exists_2_spine(A,B,C)


    for p in stream_permutations_mod_cyclic([[1,2,3,4],[5,6,7]]):
        print(p)



if __name__=='__main__':
    main()
    """
    copy_paste_to_nb([['a','A'],['b','B']],['BBa', 'Ab'])
    copy_paste_to_nb([['a','A'],['b','B']],['abAA', 'BBB'])


    copy_paste_to_nb([['a','A'],['b','B']],['baBAA', 'AbbaBBB']) #0
    copy_paste_to_nb([['a','A'],['b','B']],['baBAA', 'bbaBBBA'])
    copy_paste_to_nb([['a','A'],['b','B']],['baBAA', 'abbbABB'])
    copy_paste_to_nb([['a','A'],['b','B']],['baBAA', 'BabbbAB'])
    copy_paste_to_nb([['a','A'],['b','B']],['baBAA', 'BabbAA'])
    copy_paste_to_nb([['a','A'],['b','B']],['bABaa', 'BAbbaa']) #5
    copy_paste_to_nb([['a','A'],['b','B']],['bABaa', 'AbbaaB'])
    copy_paste_to_nb([['a','A'],['b','B']],['bABaa', 'AbbaBaa'])
    copy_paste_to_nb([['a','A'],['b','B']],['bABaa', 'abbaB'])
    copy_paste_to_nb([['a','A'],['b','B']],['bABaa', 'abaa'])
    copy_paste_to_nb([['a','A'],['b','B']],['bABaa', 'AABA']) #10
    copy_paste_to_nb([['a','A'],['b','B']],['bABaa', 'BA'])
    copy_paste_to_nb([['a','A'],['b','B']],['bABaa', 'ab'])
    copy_paste_to_nb([['a','A'],['b','B']],['AAbaB', 'ab'])
    copy_paste_to_nb([['a','A'],['b','B']],['AAbaB', 'ba'])
    copy_paste_to_nb([['a','A'],['b','B']],['AAbaa', 'ba']) #15
    copy_paste_to_nb([['a','A'],['b','B']],['b', 'ba'])
    copy_paste_to_nb([['a','A'],['b','B']],['b', 'AB'])
    copy_paste_to_nb([['a','A'],['b','B']],['b', 'A'])
    copy_paste_to_nb([['a','A'],['b','B']],['b', 'a'])

    print('---NEW SEQUENCE---')

    copy_paste_to_nb([['a','A'],['b','B']],['AbbaBBB', 'aaBAB']) #0
    copy_paste_to_nb([['a','A'],['b','B']],['bbaBBBA', 'aaBAB'])
    copy_paste_to_nb([['a','A'],['b','B']],['bbaBBBaBAB', 'aaBAB'])
    copy_paste_to_nb([['a','A'],['b','B']],['bbaBBBaBAB', 'BaaBA'])
    copy_paste_to_nb([['a','A'],['b','B']],['aBBBaBAb', 'BaaBA'])
    copy_paste_to_nb([['a','A'],['b','B']],['aBBBaBaBA', 'BaaBA']) #5
    copy_paste_to_nb([['a','A'],['b','B']],['aBBBaBaBA', 'aBABa'])
    """

    










