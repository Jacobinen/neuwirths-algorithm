This repository contains an code that is able to check if the presentation complex of a balanced group presentation is the spine of a closed orientable connected 3-manifold using Neuwirth's algorithm.

HOW TO USE: For a balanced group presentation $\langle a_1,...,a_g \mid r_1,...,r_g\rangle$ that you want to check, start by adopting a naming convention such that the i-th generator is the i-th letter of the English alphabet. If $g=4$, input the list [['a','A'],['b','B'],['c','C'],['d','D']] into the first argument of the function copy_paste_to_nb function (and similarly for other values of $g$). For the second argument of the copy_paste_to_nb function you input a list of strings that consists of your relations written as a sequence of your specified letters, where the capitalised version is the inverse of the small letter. 

and so on depending on how many you have. 

Example: for the balanced group persentation $\langle a, b \mid b^{-2}a, a^{-1}b\rangle$, you paste copy_paste_to_nb([['a','A'],['b','B']],['BBa', 'Ab']) into the main() function of neuwirth.py and run the code. 

Then, open neuwirthsalgorithm.ipynb using a SAGE kernel and paste the result into a code cell after the first X cells.


