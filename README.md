This repository contains an code that is able to check if the presentation complex of a balanced group presentation is the spine of a closed orientable connected 3-manifold using Neuwirth's algorithm.

REQUIREMENTS:

Python 3 only (no third-party packages). The permutation and orbit routines that used to come from SageMath are provided by permutation.py, which reproduces the behaviour of the SageMath/GAP functions the algorithm relies on.

HOW TO USE: 

Option 1 (running the script): For a balanced group presentation $\langle a_1,...,a_g \mid r_1,...,r_g\rangle$ that you want to check, start by adopting the naming convention that the i-th generator is the i-th letter of the English alphabet. If $g=4$, input the list [['a','A'],['b','B'],['c','C'],['d','D']] into the first argument of the check_presentation function (and similarly for other values of $g$). For the second argument of the check_presentation function you input a list of strings that consists of your relations written as a sequence of your specified letters, where the capitalised version is the inverse of the small letter. 

Example: for the balanced group persentation $\langle a, b \mid b^{-2}a, a^{-1}b\rangle$, you place check_presentation([['a','A'],['b','B']],['BBa', 'Ab']) in the into the main() function of the neuwirthcheck.py script and run `python neuwirthcheck.py`. 

Option 2 (using a Jupyter notebook): Alternatively, you can start by running copy_paste_to_nb() with the same input as check_presentation() would receive in neuwirth.py and copy pasting the result into the neuwirthsalgorithm.ipynb notebook in the indicated code cell and execute the entire notebook using a Python 3 kernel started in this directory (so that permutation.py can be imported).
