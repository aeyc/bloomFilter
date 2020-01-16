Ayca Begum Tascioglu
21600907

To run the code:

1) Please go to the directory of where source code(bloomFilter.py) located.
 
2) Add query and reference files that you want to use, in the same file with the source code.

3) After, in the terminal go to the directory where bloomFilter.py and desired input files( i.e reference.fasta and query.fasta) are located.

4) type 'make'.
	>>make
	python bloomFilter.py program.dat

5) after you see the program compiled, write the parameters as given in the assignment. 
Example:
	>>bloomFilter --ref reference.fasta --query query.fasta --kmer 3 --bloomsize 100

6)The program will output the results as:
	(2) Number of distinct 3-mers indexed in reference(cannonicals are included).
	(4) Number of distinct 3-mers scanned in query(cannonicals are NOT included)
	(1) Number of distinct 3-mers from query found in the reference.
where 2 =  Number of distinct 3-mers indexed in reference(cannonicals are included)
      4 =  Number of distinct 3-mers scanned in query(cannonicals are NOT included)
      1 =  Number of distinct 3-mers from query found in the reference.


Example Run:

(base) Ayca-MacBook-Pro:~ Ayca$ cd Desktop/Tascioglu_Ayca_hw6
(base) Ayca-MacBook-Pro:Tascioglu_Ayca_hw6 Ayca$ ls
Makefile	README.txt	bloomFilter.py	query.fasta	reference.fasta
(base) Ayca-MacBook-Pro:Tascioglu_Ayca_hw6 Ayca$ make
python bloomFilter.py program.dat
bloomFilter --ref reference.fasta --query query.fasta --kmer 3 --bloomsize 100
(2) Number of distinct 3-mers indexed in reference(cannonicals are included).
(4) Number of distinct 3-mers scanned in query(cannonicals are NOT included)
(1) Number of distinct 3-mers from query found in the reference.