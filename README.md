This repository contains an code that is able to check if the presentation complex of a balanced group presentation is the spine of a closed orientable connected 3-manifold using Neuwirth's algorithm.

HOW TO USE: For a balanced group presentation $\langle a_1,...,a_g \mid r_1,...,r_g\rangle$ that you want to check, start by adopting the naming convention that the i-th generator is the i-th letter of the English alphabet. If $g=4$, input the list [['a','A'],['b','B'],['c','C'],['d','D']] into the first argument of the function copy_paste_to_nb function (and similarly for other values of $g$). For the second argument of the copy_paste_to_nb function you input a list of strings that consists of your relations written as a sequence of your specified letters, where the capitalised version is the inverse of the small letter. 

Example: for the balanced group persentation $\langle a, b \mid b^{-2}a, a^{-1}b\rangle$, you paste copy_paste_to_nb([['a','A'],['b','B']],['BBa', 'Ab']) into the main() function of neuwirth.py script and run the code. 

Then, open neuwirthsalgorithm.ipynb notebook using a SAGE kernel. Paste the result of the copy_paste_to_nb function into a cell beneath the first 23 cells. Finally, write exists_2_spine(A,B,C) in the subsequent cell and execute the entire notebook. 


