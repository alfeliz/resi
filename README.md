# resi

Octave script to read radial and electrical ALEX short data and transfrom them into voltage, resistivty and other values.

It is a executable script with three parameters given: 
- first parameter, complete name of the wire parameters file. It should have the shot name follow by a "-" and whatever.
- second parameter, compelte name of the radial data.
- third parametr, number of scope channels  to read.

Produced output are the files:
- "SHOTNAME_res.dat" Datasheet with resistivity boundaries as a function of time and temperature.
- "SHOTNAME_voltages.dat" Datasheet with voltages recorded, see first line of the file for explanation of them.
- "SHOTNAME_combo.dat" Datasheet with voltage and resistance data.
- "SHOTNAME_param.txt" Text file with parameters extracted of the calculation (L and R of the anode, time duration of the dark pause, etc.)

To execute it in a linux system, just call him with the parameters in order and with spaces between them.

###Needed Octave functions:
It needs these octave functions, alse present in this repository:
+ **trans.m** 	To scale variables.
+ **deri.m**		To derivate vectors with a lot of noise.
+	**peakdet.m**	To detect peaks with octave under 4.0 version.(http://www.billauer.co.il/peakdet.html)
+ **display_rounded_matrix**	To nicely print numerical data.


