# CourseworkSIMD
Bewith180997
Coursework involving Multithreading and SIMD. We were provided with a basic program and were tasked to make it more efficient by implementing multithreading and Single Instruction Multiple Data (SIMD). Of a class average of 55, this scored 70. The general feedback was that multithreading was properly implemented, and the workload between threads was balanced. More marks could have been acheived by improved by better basic and advanced SIMD optimisations.
The following is a copy of the analysis I was also tasked to provide when submitting this coursework to demonstrate running times and understanding.


Results:
Computer/processor details:

>3 year old HP Pavilion Notebook
Intel Core i5-6200U CPU @ 2.30GHz
	2 Cores
	4 Logical Processors
	L1 cache: 128KB
	L2 cache: 512KB
	L3 cache: 3.0MB
Memory: 8GB @ 2133MHz


#define NB_THREADS 2

CPU 1 Thread :	2.062 s

CPU MT		 :	1.020 s
	Core 1 : 1.015 s
	Core 2 : 1.015 s
	
SIMD 1 Thread:	0.472 s

SIMD MT		 :	0.277 s
	Core 1 : 0.235 s
	Core 2 : 0.235 s



#define NB_THREADS 4

CPU 1 Thread :	1.501 s

CPU MT		 :	0.559 s
	Core 1 : 0.553 s
	Core 2 : 0.553 s
	Core 3 : 0.555 s
	Core 4 : 0.555 s
	
SIMD 1 Thread:	0.394 s

SIMD MT		 :	0.231 s
	Core 1 : 0.198 s
	Core 2 : 0.198 s
	Core 3 : 0.198 s
	Core 4 : 0.198 s


Conclusion:
There would be some fluctuation, but using the simple method sits somewhere around 4x faster for SIMD compared to without.
Single core SIMD is only slightly faster than four MT cores running the simple calculation.
This is likely from the single core SIMD running more calculations than necessary. Particularly with 
the method MulC8(). Whereas, the simple MT would run only the necessary number of repetitions for each complex struct.
There is a section where it would be more beneficial to using the instristic's horizontal_maximum method, which has been
specified at the section.

Roughly 7x speedup from CPU 1 Thread to SIMD MT. 6x When using 4 cores for SIMD MT


