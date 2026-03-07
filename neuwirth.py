#Want to build function that generates an admissible permutation A from a given (balanced) group presentation
#first, we start by defining a buch of convenient functions for the exact data format we will be using

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
        

    Output=""

    for idx, entry in enumerate(A):
        if idx%2==0:
            if idx!=0:
                Output+=''
            Output+='('+str(translate_to_integer(entry, prevalences))+','
        else:
            Output+=str(translate_to_integer(entry, prevalences))+')'
    

    if __name__=='__main__':
        print("A=Permutation('"+Output+"')")
    return Output



    
    Output=''

    for idx, entry in enumerate(A):
        if idx%2==0:
            Output+='('+to_e_notation(entry)
        else:
            Output+=to_e_notation(entry)+')'

    #print(Output)
    
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











 

    


def copy_paste_to_nb(gen, rel):
    info="#("+rel[0]+","+rel[1]+")"
    print(info)
    admissible(gen, rel)
    prev=prevalences(gen,rel)
    print("B=Permutation('"+B(sum(prev))+"')")
    print("C=Permutation('"+C(prev,len(prev))+"')")




def main():
    # Replace with the group presentation that you want to check
    copy_paste_to_nb([['a','A'],['b','B']],['aaBBB', 'ABAbab'])





if __name__=='__main__':
    main()


    """
    # This is the trivialisation sequence that we checked in the paper
    copy_paste_to_nb([['a','A'],['b','B']],['aaBBB', 'ABAbab']) #0
    copy_paste_to_nb([['a','A'],['b','B']],['aBBBa', 'BAbabA'])
    copy_paste_to_nb([['a','A'],['b','B']],['aBBBa', 'BAbaBBa'])
    copy_paste_to_nb([['a','A'],['b','B']],['aBBBa', 'AbaBBaB'])
    copy_paste_to_nb([['a','A'],['b','B']],['aBBaBBaB', 'AbaBBaB'])
    copy_paste_to_nb([['a','A'],['b','B']],['aBBaBBaB', 'BabAbbA'])
    copy_paste_to_nb([['a','A'],['b','B']],['aBBaBBaB', 'BaBaB'])
    copy_paste_to_nb([['a','A'],['b','B']],['BBaBBaBa', 'AbAbb'])
    copy_paste_to_nb([['a','A'],['b','B']],['BBa', 'AbAbb'])
    copy_paste_to_nb([['a','A'],['b','B']],['BBa', 'Ab'])
    copy_paste_to_nb([['a','A'],['b','B']],['B', 'Ab'])
    copy_paste_to_nb([['a','A'],['b','B']],['B', 'A'])
    """

    










