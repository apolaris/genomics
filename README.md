# CBMF4761 genomics project

Files:
	main.py is main program.
	multiple.py is the CLUSTAL module.
	multiple2.py is the TCOFFEE module.
	multiple3.py is the quick TCOFFEE module.
	output_example1.txt is the output of 20 real_long_data by TCOFFEE.
	output_example2.txt is the output of 50 real_short_data by CLUSTAL.
	real_long_data.txt is the long real sequences in nomal format.
	real_long_fasta.txt is the long real sequences in FASTA format.
	real_short_data.txt is the long real sequences in nomal format.
	real_short_fasta.txt is the long real sequences in FASTA format.
	(for FASTA, the <input_type> is 0, for nomal format, it's 1.)
	simulate_data.txt is simulated sequences in nomal format.
	makedata.c is c program to simulate some sequences by simulating phylogeny tree and produce sequences.
	



The alignment program is implemented by python2 and it only imports time module. The makedata.c compiles in windows by gcc 4.7.

The main program is main.py.

usage of main.py: python main.py <read_file_name> <write_file_name> <opengap_penalty> <extendgap_penalty> <match_score> <mismatch _score> <seq_num> <option> <refine_time> <input_type>

	<read_file_name>: The file from which to read sequence.

	<write_file_name>: The file to which the program output the alignment.

	<opengap_penalty>: The penalty for opening a gap, should be positive.

	<extendgap_penalty>: The penalty for extending a gap, should be positive.
	
	<match_score>: The score when the base is matched.

	<mismatch_score>: The score when the base is not matched.

	<seq_num>: The number of sequences you want to align. If you want to align all the sequences, just use 0 as the seq_num.

	<option>: Determine which algorithm you want to use. If option is 1, it will run CLUSTAL. If option is 2, it will run TCOFFEE. If option is 3, it will run Quick TCOFFEE.

	<refine_time>: Times you want to refine. 0 is not refine. Besides, in the Quick TCOFFEE algorithm, we didn't implement the refinement.

	<input_type>: The type of the read file. This file support two types read file. First type is .FASTA. Second type is the file that just stores sequences and stores one sequence in one line. If input_type is 0, the read file is .FASTA. If input_type is 1, the read file is second type.

	Example:  python main.py short_real_data.txt output.txt 3 2 -1 2 0 1 10 1

	This command means that the input file is real_short_data.txt, the output file is output.txt, the open gap is 3, extend gap is 2, match score is -1, mismatch score is 2, align all the sequences in the file, use CLUSTAL algorithm, refine 10 times, the input file is second type.

	After running the program, it will output the SP_score to the screen. The performance is better when the SP_score is smaller. And we suggest that the open gap is 4, extend gap is 3, match score is 0, mismatch score is 3.





usage of evalue.py: python evalue.py <filename>

	<filename>: The sequences you want to evalue. But it only accept the output file by main.py. 

	The score and penalty in this program: match score is 0, mismatch score is 3, open gap penalty is 4, extend gap penalty is 3.





usage of makedata.c: Just compile it and run it. Then input the length and the number of individuals. The number of sequences it output is random. And the sequences is made of 20 characters.

	
