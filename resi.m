#! /usr/bin/octave -q
#La línea de arriba señala que es un script que se le pasará a este programa. Recuerda que el archivo tiene que ser ejecutable. 
#La opción -q es para evitar el mensajillo de saludo al principio.
## Copyright (C) 2016 Gonzalo Rodríguez Prieto <gonzalo.rprieto@uclm.es>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.
###########################################################################################
#
#  OCTAVE SCRIPT TO OBTAIN THE RESISTIVITY DATA FROM ALEX VOLTAGE AND RADIAL DATA
#    Made by Gonzalo Rodríguez Prieto
#              Version 1.0
#
###########################################################################################
#	IT USES THE FUNCTIONS (they should be in the same folder than this program):
#	trans.m 	To scale variables.
#	deri.m		To derivate vectors with a lot of noise.
#	peakdet.m	To detect peaks with octave under 4.0 version.(http://www.billauer.co.il/peakdet.html)
#	display_rounded_matrix	To nicely print numerical data.
###########################################################################################

tic; #Control total computing time.

#Empty memory from user created variables:
clear;

###########################################################################################
# General script parameters nad unit conversions:
###########################################################################################
pi = 3.141592; %Pi value
mu = 4*pi*1e-7; %Vacuum permeability in S.I. units

mm = 1e-3;
us = 1e-6;

###########################################################################################
# Transforming initial data pass through command line:
###########################################################################################
arg_list = argv(); %Arguments passed to the script

printf("Argumentos: \n");
for i = 1:nargin %nargin es el número de argumentos, incluidos el comando de entrada.
  printf (" %s \n", arg_list{i});
endfor
printf ("\n");

#String with the parameters file name:
%datafile = "ALEX345-data.txt"; %Testing purposes.
datafile = arg_list{1};

#file with radial data:
%radfile = "example-radial.dat"; %Testing purposes
radfile = arg_list{2};

#Number of channels recorded with electrical signals:
num_chan = 3; #Testing purposes
%num_chan = str2num(arg_list{3}); %It records the arguments as strings/ch


###########################################################################################
# Reading the parameters file with wire paramters and magnitudes:
# It must be written with parametr name up, and the next line the value and a blank line after.
###########################################################################################
[dfile, msg] = fopen(datafile, "r");
if (dfile == -1) 
   error ("resi script: Unable to open data file name: %s, %s",datafile, msg); 
endif;

lines = fskipl (dfile, Inf); %Counting the lines in the data file
frewind (dfile); %Repositiong the reading counter to the file beginning

r_rows =  "%s \n"; %1 row to read, but three things, including the TAB in between columns.
if feof(dfile)==0 #Read the file until it finds the EndOfFile character.
   data = textscan(dfile, r_rows, "collectoutput", 1);
endif;

fclose(dfile);

i = 1;
while ( i<= prod(size(data{1})) )
 [dat,cool] = str2num(cell2mat(data{1}(i))); %Convert whatever in numbers always
 
	 if cool==0 %We are with a data label and convert data in the next line.
		str = cell2mat(data{1}(i));
		switch (str) %"If" sentence cannot handle strings well.
		 case "dew" %Wire diameter to transform in radius
			radius = str2num(cell2mat(data{1}(i+1)));
			scale_rad = trans(cell2mat(data{1}(i+2)));
			radius = radius * scale_rad* 0.5; %To transform in radius the diameter
			i = i + 2;
		 case "lew" %Wire length
			lwire = str2num(cell2mat(data{1}(i+1)));
			scale_len = trans(cell2mat(data{1}(i+2)));
			lwire = lwire * scale_len;
			i = i + 2;
		 case "distwire" %Wire distance for inductance
			dwire = str2num(cell2mat(data{1}(i+1)));
			scale_len = trans(cell2mat(data{1}(i+2)));
			dwire = dwire * scale_len;
			i = i + 2;
		 case "Cv" %Specific heat in J * mol / K
			Cv = str2num(cell2mat(data{1}(i+1)));
			scale_len = trans(cell2mat(data{1}(i+2)));
			Cv = Cv * scale_len;
			i = i + 2;
		  case "density" %density
			den = str2num(cell2mat(data{1}(i+1)));
			scale_len = trans(cell2mat(data{1}(i+2)));
			den = den * scale_len;
			i = i + 2;	
		  case "mol.weight" %molecular weight
			weight = str2num(cell2mat(data{1}(i+1)));
			scale_len = trans(cell2mat(data{1}(i+2)));
			weight = weight * scale_len;
			i = i + 2;
		  case "tboil" %boiling temperature
			Tboil = str2num(cell2mat(data{1}(i+1)));
			scale_len = trans(cell2mat(data{1}(i+2)));
			Tboil = Tboil * scale_len;
			i = i + 2;	
		  case "tion" %ionization temperature
			Tion = str2num(cell2mat(data{1}(i+1)));
			scale_len = trans(cell2mat(data{1}(i+2)));
			Tion = Tion * scale_len;
			i = i + 2;															
		endswitch;
    
	 endif; 
i = i + 1; 
endwhile;

###########################################################################################
# Reading the data files with electrical channels information:
###########################################################################################
 chan_rows =  "%f %f"; %Rows definition on the CSV channel datafiles
for i=1:num_chan
 cha_name = ["ALL_CH" num2str(i) ".csv"]; %Making the name of each channel present
 cha_file = fopen(cha_name, "r"); %Opening
 data = textscan(cha_file, chan_rows, "endofline", "\n"); %Arra. 1, time in seconds, 2, signal in volts
 fclose(cha_file);	 
	 if i==1 %Ch1 Rogowsky (ALTER THIS IF SEQUENCE WHEN THE CHANNEL CONFIGURATION CHANGES!!!!!)
	  rog = 9.12*1e9*data{2}; %I-dot values (A/s)
	 elseif i==2 %Ch2 2-Resistive divider 
	  res2 = data{2}*1358.91; %Volts
	 elseif i==3 %Ch3 3-resisitive divider
	  res3 = data{2}*2389.13; %Volts
	 endif
endfor;
 
tvec = data{1}; %Time vector. All signals synchronized.
 
###########################################################################################
# Reading the data files with radial data information:
#It assumes that radial data are time(µs) rad(mm) format. Change otherwise.
###########################################################################################
chan_rows =  "%f %f"; %Rows definition on the CSV channel datafiles
[t_rad, rad_data] = textread(radfile, chan_rows, "endofline", "\n", "headerlines",1); %
t_rad = t_rad * us; %In seconds now.
rad_data = rad_data * mm; %In meters now.




###########################################################################################
# Calculating parameters:
###########################################################################################

#Current:
curr = cumtrapz(tvec,rog);

#Measured voltage:
vmeas = res2 - res3 ;

#Finding L2 and R2, the lumped parameters for the anode:
[min_vol, pos] = min(vmeas); %Position of the first minimum peak. To perform the approximation to find L2 and R2


X(:,1) = zeros(length(rog(pos:end)),1); %Theoretically, it does forces the pass fro zero, but it works better...
X(:,2) = rog(pos:end);
X(:,3) = curr(pos:end);

[anode,anode_error] = regress(vmeas(pos:end),X); %Multiple linear regreession.

anode_error = anode_error(:,2)-anode_error(:,1); %L2 and R2 errors for the confidence intervals.

L2 = anode(2); %Henrios

R2 = anode(3); %Ohms

vlast = L2*rog + R2*curr; %Voltage across the last part of the circuit

vew = vmeas - vlast; %Voltage on the wire

#Calculating border inductance:
lbor = ((mu*lwire)/(2*pi)) * (1./(cosh(dwire/rad_data)));
lbor = lbor';
#Derivative of the inductanece:
derlbor = deri(lbor, abs(t_rad(1)-t_rad(2)))';


%~ derlbor =  [0'; %Adding extra terms to account for the ones that we lose on the derivative

idx = ( tvec >= t_rad(1) ) & ( tvec <= t_rad(end) ); %Index with the time of the appearance of the radial data.

#Calculating inductive voltage, adjusting by length of rog(idx):
if length(rog(idx)) > length(lbor)
 lbor = [ lbor; zeros(length(rog(idx))-length(lbor),1)];
elseif length(rog(idx)) < length(lbor)
 lbor = lbor(1:length(rog(idx)));
endif;

if length(rog(idx)) > length(derlbor)
 derlbor = [ derlbor; zeros(length(rog(idx))-length(derlbor),1)];
elseif length(rog(idx)) < length(derlbor)
 derlbor = derlbor(1:length(rog(idx)));
endif;

van = lbor.*rog(idx) + derlbor.*curr(idx);

vanode = zeros(length(vmeas),1); %Vector of same length that other voltage measurements

vanode(idx) = van; %Putting the values of inductive voltage substraction in its correct place

vres = vew - vanode; %Real resistive voltage across the wire.

#Calculating resitive resistance:
Rres = vres./curr; %Ohms and cool!!!!

#Putting the radial expansion in an adequate size and synchronization vector:
rad = zeros(length(vmeas),1);

if length(rad(idx)) > length(rad_data)
 rad(idx) = [ rad_data ; zeros( length(rad(idx)) - length(rad_data) , 1) ];
elseif length(rad(idx)) < length(rad_data)
 rad(idx) = rad_data(1:length(rad(idx)));
elseif length(rad(idx)) == length(rad_data)
 rad(idx) = rad_data;
endif;


##
#Gas resistivity lower limit:
##
resgasdown =  (Rres .* pi .* rad.^2)./lwire; %Ohm*m


#Gas resistivity lower limit initial parameters calculations:
moles = (den * lwire * pi * radius^2) / weight; %Mass in moles

#time interval findings:
#Peak detection with peakdet(http://www.billauer.co.il/peakdet.html)
dvres = abs(deri(vres, abs(tvec(1)-tvec(2)))); %First, I derive the resistive voltage and transform it in a positive value.
peaks = peakdet(dvres,max(dvres)*0.25,tvec(1:length(dvres))); %Now, I find the peaks that are larger than 1/4 of maximum value.
[peak_max,pos] = max(peaks(:,2)); %Position of the maximum peak (the first timing value is here)

t1 = peaks(pos,1); %time of first position in seconds. (boiling moment)

t2 = peaks(end,1); %Last peak is the most away and signals the end of voltage...(ionization moment)

delta_t = t2-t1; 

if delta_t < 0
 error("script resi: time interval for dark pause negative, something wrong.");
endif;

#Linear relation of time with temperature:
alpha = (Tboil-Tion)/(t1-t2);
bet = Tboil - ( alpha * t1);

idx2 = ( tvec >= t1 ) & ( tvec <= t2 ); %Index with the time of the appearance of the radial data.

T =  alpha.*tvec(idx2) + bet; %Sintetic Temperature vector(Kelvin). (VALID ONLY BETWEEN t1 and t2)

temp = zeros(length(vmeas),1);
temp(idx2) = T;

delta_T = Tion-Tboil;

##
#Gas resistivity upper limit:
##
resgasup = (vres.^2 ./ (moles * Cv) ) .* ( (pi .* rad.^2)./lwire ) .* (delta_t/delta_T); %Ohms*m

###########################################################################################
# Removing data garbage from data files:
###########################################################################################
resgasup(temp==0) = 0;
resgasdown(temp==0) = 0;

###########################################################################################
# Saving data files:
###########################################################################################
shotname = strtok(datafile,"-"); #Taking tha name from datafile.

#Saving resistivities:
#Making the time from 0 to 1 vector:

t_un = linspace(0,1,length(resgasup(resgasup!=0)))'; %Linear values between initial melting and final ionization times

t_uniform = zeros(length(vmeas),1); %Vector of same length that other voltage measurements

t_uniform(resgasup!=0) = t_un; %Creating a time vector of the adequate length.

res = [tvec./us t_uniform resgasup resgasdown Rres];
redond = [4 3 4 4 4];
name = horzcat(shotname,"_res.dat"); #Adding the right sufix.
output = fopen(name,"w"); #Opening the file.
fdisp(output,"time(µs)	time_uniform	res.up(Ohm/m)	res.down(Ohm/m)	resistance(Ohms)"); #First line.
display_rounded_matrix(res, redond, output); 
fclose(output);

#Saving voltages:
voltages = [tvec./us temp vmeas vanode vres];
redond = [4 3 4 4 4];
name = horzcat(shotname,"_voltages.dat"); #Adding the right sufix.
output = fopen(name,"w"); #Opening the file.
fdisp(output,"time(µs)	Temper.(K)	Vmeas(V)	Vanode(V)	Vresis(V)"); #First line.
display_rounded_matrix(voltages, redond, output); 
fclose(output);

#Saving voltage and resistance:
combo = [tvec./us vres Rres];
redond = [4 4 4];
name = horzcat(shotname,"_combo.dat"); #Adding the right sufix.
output = fopen(name,"w"); #Opening the file.
fdisp(output,"time(µs)	Vresis(V)	Resis.(Ohm)"); #First line.
display_rounded_matrix(combo, redond, output); 
fclose(output);

#Saving parameters:
name = horzcat(shotname,"_param.txt"); #Adding the right suffix.
output = fopen(name, "w"); #Opening the file.
fdisp(output, [ "L_anode(Hr) " num2str(L2) "\n" ]); %Inductance of the anode
fdisp(output, [ "R_anode(Ohm) " num2str(R2) "\n"]); %Resistance of the anode
fdisp(output, "L and R errors:"); %Erros of the parameters.
fdisp(output, anode_error');
fdisp(output, [ "\n" "dark pause time interval(µs) " num2str(delta_t) "\n"]); %time interval of dark pause
fdisp(output, "Temperature parameters"); %Linear apporx. for temperture parameters (NEXT LINE)
fdisp(output, [alpha bet]); %As only numbers are given, they area automatically converted to strings and saved properly.
fdisp(output, ["\n" "moles " num2str(moles)]); %Mass in moles.
fclose(output);

toc;


#That`s that`s all folks!!!
