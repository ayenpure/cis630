CIS630, Spring 2017, Term Project I, Due date: May 3, 2017

Name : Abhishek Yenpure
Student ID : 951495864


What programming language did you use to write your code?

	C++


Does your program compile on ix (CIS department server)?

	Yes

How should we compile your program on ix? (please provide a makefile)

	execute "make" to compile the program
	execute "make run" to run the program


Does your program run on ix?

	Yes


Does your program calculate the credit values accurately?

	Yes

Does your program have a limit for the number of nodes it can handle in the input graph?

	Yes


If yes, what is the limit on graph size that your program can handle?

	Then number of nodes are bound by the maximum number allowed by the int datatype.
	Also if you try to run the program over around 1000 rounds, the program fails because it cannot allocate enough memory.


How long does it take for your program with two partitions to read the Flickr input files, perform 5 round and write
the output of each round in the output files on ix-dev?

	~ 10-12 seconds


Does your program run correctly with four partitions on ix-dev?

	Yes


Does you program end gracefully after completing the specified number of rounds?

	Yes
