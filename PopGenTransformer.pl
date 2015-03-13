#!/usr/bin/perl
use warnings;
use strict;
print "
************************************************************************
************************************************************************
*****This PopGenetConverter is used to convert your genotyping data*****
*****   into format for multiple programs in Population Genetics   *****
************************************************************************
************************************************************************
\n";
###Created by Fu-Hao Lu
###Feel free to modify this script unless for commercial usage###
###Any problems or questions, please contact lufuhao@gmail.com
###Homepage: http://code.google.com/p/popgentransformer

our ($file_name, $mk_num, $ind_num, $project_name, @linearray, %pop_list, @pop_name_list, $pop_num, @mk_name_list, @program_choice_list, $program_choice, $script_cycle_control);

###--------------input customized data start---------------###

print "\nPlease input the name of TEXT file you want to convert:"."\n";
$file_name=<>;
chomp $file_name;
print "\nPlease indicate the marker number:", "\n";
$mk_num=<>;
chomp $mk_num;

print "\nPlease indicate the individual nnumber:\n";
$ind_num=<>;
chomp $ind_num;

print "\nplease indicate the job name:\n";
$project_name=<>;
chomp $project_name;

print "\nThere are $mk_num markers, and $ind_num individuals\n and the project name was assigned as: $project_name\n\n";
###--------------input customized data end---------------###





##########  Reading from the input file ...start  ##########
my($i,$j);
open (INPUTFILE, "$file_name") || die "Fail to open input file";

$i=1;
$j=3;

while (<INPUTFILE>){
	chomp $_;
	$linearray[$i] = [split(/\t/, $_)];
	$i++;
}

close INPUTFILE;

my ($temp_01, $temp_02); 
for ($temp_01=2; $temp_01 < ($ind_num + 2); $temp_01++){
	$pop_list{$linearray[$temp_01][1]}++;
}
@pop_name_list=keys(%pop_list);
$pop_num=scalar(@pop_name_list);



for ($temp_02=3; $temp_02 < ($mk_num * 2 + 3); $temp_02+=2){
	push(@mk_name_list, $linearray[1][$temp_02]);
}
##########  Reading from the input file ...End  ##########



$script_cycle_control=1;
while ($script_cycle_control == 1){
print "Format supported\n
*****************************************************************\n";
print "
1-Popgene,		2-PowerMarker,	3-Structure, 
4-GenePop,	5- Fstat,			6-Microsatellite_toolkit, 
7-GenAlEx,		8-SPAGeDi,			9-Arlequin
10-Genetix,		11-Convert,		12-TASSEL,
13-(on request)\n";
print 
"*****************************************************************\n
Please indicate the target format you want to convert: (1-11)\n
For example: 3(Enter for single convet)\n
Please input your choice:\n";

$program_choice=<>;
chomp $program_choice;
#@program_choice_list=<>;
#foreach $program_choice (@program_choice_list){
#	chomp $program_choice;
	if ($program_choice eq "1"){
		& popgene();
	} 
	if ($program_choice eq "2"){
		& powermarker();
	}	
	if ($program_choice eq "3"){
		& structure();
	}
	if ($program_choice eq "4"){
		& GenePop();
	}
	if ($program_choice eq "5"){
		& Fstat();
	}
	if ($program_choice eq "6"){
		& Microsatellite_toolkit();
	}
	if ($program_choice eq "7"){
		& GenAlEx();
	}
	if ($program_choice eq "8"){
		& SPAGeDi();
	}
	if ($program_choice eq "9"){
		& Arlequin();
	}
	if ($program_choice eq "10"){
		& Genetix();
	}
	if ($program_choice eq "11"){
		& MConvert();
	}
	if ($program_choice eq "12"){
		& Tassel();
	}
	if ($program_choice > 12){
		print "\nyou may send the program name that you want to the E-mail: \n"."lufuhao\@gmail.com\n";
		print "I will implement the function you needed when I get spare time\n";
		print "Thanks a lot for your interest in this Small Program. See ya...\n";
	}
#}
print "Do you want to convert to another format?\n1 for YES\n0 for No\n";
$script_cycle_control=<>;
chomp $script_cycle_control;
}




###---------------- PopGene conversion start... --------------------###
sub popgene {
my (@ind_genotype_matrix, %popgeneMK, $temp_03, $temp_04, $temp_05, $temp_06, $temp_07, $temp_08, $temp_09, $temp_10, $temp_11, $temp_12, $temp_13, $temp_14, $temp_15, %popgene_genotype_get, @keys, @popgene_mk_genotypes, $popgene_genotype, $ind_genotype_matrix, $popgene_ind_genotype, @ert);
@ind_genotype_matrix=();
%popgeneMK=();
%popgene_genotype_get=();
@keys=();
@popgene_mk_genotypes=();
@ert=();

#test#print @{$linearray[1]}, "\n"; #test#

print "PopGene conversion start...\n";
open(POPGENE, ">popgene_input.txt") || die "Can not create popgene_input.dat file";
print POPGENE 
"/* $project_name */
Number of populations = $pop_num
Number of loci = $mk_num
Locus name:
@mk_name_list\n";

@ind_genotype_matrix=();
%popgeneMK=();
$temp_15=0;
%popgene_genotype_get=();
for($temp_04=3; $temp_04<(2*$mk_num + 3); $temp_04++){
	for($temp_05=2; $temp_05 < $ind_num + 2; $temp_05++){
		if ($linearray[$temp_05][$temp_04] ne ("?"||".")){
			$popgene_genotype_get{$linearray[$temp_05][$temp_04]}++;
		}
	}
	if(($temp_04>3) && !($temp_04%2)){
		@keys=keys(%popgene_genotype_get);
		@popgene_mk_genotypes=sort (@keys);
		$temp_05="A";
		foreach $popgene_genotype (@popgene_mk_genotypes){
			$popgeneMK{$popgene_genotype}=$temp_05++;
		}
	
		$temp_14=$temp_04-1;
		for ($temp_06=2; $temp_06 < ($ind_num + 2); $temp_06++){
			if ($linearray[$temp_06][$temp_14] ne ("?"||".")){
				$temp_12=$linearray[$temp_06][$temp_14];
				$temp_08=$popgeneMK{$temp_12};
			}
			else {
			$temp_08=".";
			}
			if ($linearray[$temp_06][$temp_04] ne ("?"||".")){
				$temp_13=$linearray[$temp_06][$temp_04];
				$temp_09=$popgeneMK{$temp_13};
			}
			else {
			$temp_09=".";
			}
			
			$ind_genotype_matrix[$temp_06][$temp_15]="$temp_08"."$temp_09";
				
		}
		%popgene_genotype_get=();
		%popgeneMK=();
		$temp_15++;
	}
}

for ($temp_11=0; $temp_11<$pop_num; $temp_11++){
	print POPGENE "ID=",$pop_name_list[$temp_11], "\n";
	for ($temp_10=2; $temp_10<=($ind_num+1); $temp_10++){
		if ($linearray[$temp_10][1] eq $pop_name_list[$temp_11]){
			@ert=@{$ind_genotype_matrix[$temp_10]};
			$popgene_ind_genotype=join("\t", @ert);
			print POPGENE $popgene_ind_genotype, "\n";
		}
		if ($temp_10 == ($ind_num+1)){
			print POPGENE "\n";
		}
	}
}

close POPGENE;
print "PopGene conversion finished\n";
}
###---------------- PopGene conversion finish... --------------------###


###-------------- PowerMarker conversion start...  ------------------###
sub powermarker{
my (@powermarker_matrix, $temp_16, $temp_17, $temp_18, $temp_19, $temp_20, $temp_21, @powermarker_ind_list, @powermarker_ind_tab);
@powermarker_matrix=();
@powermarker_ind_list=();
@powermarker_ind_tab=();

print "Powermarker conversion start...\n";
open (POWERMARKER, ">PowerMarker_input.txt") || die "Failed";

$powermarker_matrix[1][0]="Ind_Cat";
$powermarker_matrix[1][1]="Geo_Cat";
$powermarker_matrix[1][2]="Struct_Cat";

for ($temp_16=2; $temp_16<($ind_num+2); $temp_16++){
	for ($temp_17=0; $temp_17<3; $temp_17++){
		$powermarker_matrix[$temp_16][$temp_17]=$linearray[$temp_16][$temp_17];
	}
}


for ($temp_18=1; $temp_18<($ind_num+2); $temp_18++){
	$temp_20=3;
	for ($temp_19=3; $temp_19<($mk_num*2+3); $temp_19+=2){
		if ($temp_18==1){
			$powermarker_matrix[$temp_18][$temp_20]=$linearray[$temp_18][$temp_19];
		}
		else {
			$powermarker_matrix[$temp_18][$temp_20]=$linearray[$temp_18][$temp_19]."/".$linearray[$temp_18][$temp_19+1];
		}
	$temp_20++;
	}
}

for ($temp_21=1; $temp_21<($ind_num+2); $temp_21++){
	@powermarker_ind_list=@{$powermarker_matrix[$temp_21]};
	$powermarker_ind_tab[$temp_21]=join("\t", @powermarker_ind_list);
	print POWERMARKER  $powermarker_ind_tab[$temp_21], "\n";
}

close POWERMARKER;
print "Powermarker conversion finished...\n";
}
###-------------- PowerMarker conversion finish...  ------------------###

###--------------  Structure conversion start...   ----------------------###
sub structure {
my ($temp_22, $temp_23, $temp_24, $temp_25, $temp_26, @structure_matrix,  $structure_pop, %structure_pop_conversion);
@structure_matrix=();
%structure_pop_conversion=();
print "Structure conversion start...\n";
open (STRUCTURE, ">structure_input.str") || die "Structure conversion failed";
print STRUCTURE "@mk_name_list"."\n";

$temp_23=1;
foreach $structure_pop (@pop_name_list){
	$structure_pop_conversion{$structure_pop}=$temp_23++;
}

for ($temp_22=2; $temp_22<($ind_num+2); $temp_22++){
	$structure_matrix[$temp_22*2-2][0]=$structure_matrix[$temp_22*2-1][0]=$linearray[$temp_22][0];
	$structure_matrix[$temp_22*2-2][1]=$structure_matrix[$temp_22*2-1][1]=$structure_pop_conversion{$linearray[$temp_22][1]};
}

for ($temp_24=3; $temp_24<($mk_num*2+3); $temp_24++){
	for ($temp_25=2; $temp_25<($ind_num+2); $temp_25++){
		if ($temp_24%2){
			if ($linearray[$temp_25][$temp_24] ne ("?" || ".")){
				$structure_matrix[$temp_25*2-2][($temp_24+1)/2]=$linearray[$temp_25][$temp_24];
			}
			else{
				$structure_matrix[$temp_25*2-2][($temp_24+1)/2]="-9";
			}
		}
		else {
			if ($linearray[$temp_25][$temp_24] ne ("?" || ".")){
				$structure_matrix[$temp_25*2-1][$temp_24/2]=$linearray[$temp_25][$temp_24];
			}
			else{
				$structure_matrix[$temp_25*2-1][$temp_24/2]="-9";
			}
		}
	}
}

for ($temp_26=2; $temp_26<($ind_num*2+2); $temp_26++){
	print STRUCTURE "@{$structure_matrix[$temp_26]}"."\n";
}

close STRUCTURE;
print "The population (col-2) was coded as follows:\n", "%structure_pop_conversion", "\n";
print "Structure conversion finished...\n";
}
###--------------  Structure conversion finish...   ----------------------###


###--------------   GenePop conversion start...   ----------------------###
sub GenePop {
my ($genepop_mk_name, $genepop_pop_name, @genepop_matrix, $temp_33, $temp_34, $temp_35, $temp_36, $temp_37, $temp_38);
@genepop_matrix=();
print "GenePop conversion started...\n";
open (GENEPOP, ">GenePop_input.gen") || die "GenePop conversion failed";
print GENEPOP $project_name ."\n";
foreach $genepop_mk_name (@mk_name_list){
	print GENEPOP $genepop_mk_name ."\n";
}

for ($temp_33=2; $temp_33<($ind_num+2); $temp_33++){
	$genepop_matrix[$temp_33][0]=$linearray[$temp_33][0];
	$genepop_matrix[$temp_33][1]=",";
}

foreach $genepop_pop_name (@pop_name_list){
	print GENEPOP "POP\n";
	for ($temp_34=2; $temp_34<($ind_num+2); $temp_34++){
		for ($temp_35=3; $temp_35<($mk_num*2+2); $temp_35+=2) {
			if ($linearray[$temp_34][$temp_35] eq ("?"||".")){
				$temp_36="000";
			}
			if (length($linearray[$temp_34][$temp_35]) == 1) {
				$temp_36="00".$linearray[$temp_34][$temp_35];
			}
			if (length($linearray[$temp_34][$temp_35]) == 2) {
				$temp_36="0".$linearray[$temp_34][$temp_35];
			}
			if (length($linearray[$temp_34][$temp_35]) == 3) {
				$temp_36=$linearray[$temp_34][$temp_35];
			}
			if (length($linearray[$temp_34][$temp_35]) > 3 || length($linearray[$temp_34][$temp_35]) <1) {
				print "the data at $temp_34 line and $temp_35 col have error";
				last;
			}
			
			$temp_37=$temp_35+1;
			if ($linearray[$temp_34][$temp_37] eq ("?"||".")){
				$temp_38="000";
			}
			if (length($linearray[$temp_34][$temp_37]) == 1) {
				$temp_38="00".$linearray[$temp_34][$temp_37];
			}
			if (length($linearray[$temp_34][$temp_37]) == 2) {
				$temp_38="0".$linearray[$temp_34][$temp_37];
			}
			if (length($linearray[$temp_34][$temp_37]) == 3) {
				$temp_38=$linearray[$temp_34][$temp_37];
			}
			if (length($linearray[$temp_34][$temp_37]) > 3 || length($linearray[$temp_34][$temp_37]) <1) {
				print "the data at $temp_34 line and $temp_37 col have error";
				last;
			}
			$genepop_matrix[$temp_34][($temp_35+1)/2]=$temp_36.$temp_38;
			
		}
		print GENEPOP "@{$genepop_matrix[$temp_34]}\n";
	}
}

close GENEPOP;
print "GenePop conversion finished...\n";
}
###--------------   GenePop conversion finish...   ----------------------###



###-----------------   Fstat conversion start...   ------------------------###
sub Fstat {
print "\nYou may use GenePop format for Fstat input\n";
print "In Fstat menu\\Utilities\\File_conversion\\Genepop->Fstat\n";
}
###----------------   Fstat conversion finish...   ------------------------###



###---------   Microsatellite_toolkit conversion start...  ----------------###
sub Microsatellite_toolkit {
my ($temp_39, $temp_40, $temp_41, $temp_42, @mstoolkit, $mstoolkit_tab);
@mstoolkit=();

print "Microsatellite_toolkit conversion started...";
open (MSTOOLKIT, ">Microsatellite_toolkit_input.txt") || die "Microsatellite_toolkit conversion failed...";
@mstoolkit=();
for ($temp_39=1; $temp_39<($ind_num+2); $temp_39++){
	$mstoolkit[$temp_39][0]=$linearray[$temp_39][0];
}
$mstoolkit[1][0]="";

for ($temp_40=1; $temp_40<($ind_num+2); $temp_40++){
	for ($temp_41=3; $temp_41<($mk_num*2+3); $temp_41++){
		if ($linearray[$temp_40][$temp_41] eq ("?"||".")){
			$mstoolkit[$temp_40][$temp_41-2]=0;
		}
		else {
		$mstoolkit[$temp_40][$temp_41-2]=$linearray[$temp_40][$temp_41];
		}
	}
}

for ($temp_42=1; $temp_42<($ind_num+2); $temp_42++){
	$mstoolkit_tab=join("\t", @{$mstoolkit[$temp_42]});
	print MSTOOLKIT $mstoolkit_tab ."\n";
}
close MSTOOLKIT;
print "Microsatellite_toolkit conversion finished...\n\n"."You may open the text file and copy it to Excel"."\n\n";
}
###---------   Microsatellite_toolkit conversion finish...   --------------###



###---------------   GenAlEx conversion start...   ----------------------###
sub GenAlEx {
my (@genalex_header, @genalex_matrix, $genalex_matrix_tab, $genalex_line1, $temp_43, $temp_44, $temp_45, $temp_46, $temp_47, $temp_48, $temp_49);
@genalex_header=();
@genalex_matrix=();
print "GenAlEx conversion started...\n";
open (GENALEX, ">GenAlEx_input.txt") || die "GenAlEx conversion failed...";

$genalex_header[0][0]=$mk_num;
$genalex_header[0][1]=$ind_num;
$genalex_header[0][2]=$pop_num;
$genalex_header[1][0]=$project_name;
$genalex_header[1][1]="";
$genalex_header[1][2]="";

foreach ($temp_43=0; $temp_43<$pop_num; $temp_43++){
	$genalex_header[1][$temp_43 +3]=$pop_name_list[$temp_43];
	$genalex_header[0][$temp_43 +3]=$pop_list{$pop_name_list[$temp_43]};
}

$temp_44=join("\t", @{$genalex_header[0]});
$temp_45=join("\t", @{$genalex_header[1]});
print GENALEX $temp_44 ."\n";
print GENALEX $temp_45 ."\n";

for ($temp_46=2; $temp_46<($ind_num+2); $temp_46++){
	for ($temp_47=0; $temp_47<($mk_num *2 +3); $temp_47++){
		if ($temp_47<2){
			$genalex_matrix[$temp_46][$temp_47]=$linearray[$temp_46][$temp_47];
		}
		if ($temp_47>2) {
			if ($linearray[$temp_46][$temp_47] ne ("?" || ".")){
				$genalex_matrix[$temp_46][$temp_47-1]=$linearray[$temp_46][$temp_47];
			}
			else {
				$genalex_matrix[$temp_46][$temp_47-1]="0";
			}
		}
	}
}

$genalex_line1=join("\t\t", @mk_name_list);
print GENALEX $linearray[1][0] . "\t" . $linearray[1][1] . "\t" . $genalex_line1 ."\t\n";

foreach $temp_49 (@pop_name_list){
	for ($temp_48=2; $temp_48<($ind_num+2); $temp_48++){
		if ($genalex_matrix[$temp_48][1] eq $temp_49){
			$genalex_matrix_tab=join("\t", @{$genalex_matrix[$temp_48]});
			print GENALEX $genalex_matrix_tab ."\n";
		}
	}
}
close GENALEX;
print "GenAlEx conversion finished...\n";
}
###--------------   GenAlEx conversion finish...   ----------------------###



###--------------   SPAGeDi conversion start...   ----------------------###
sub SPAGeDi {
	print "\nYou may use GenePop format for SPAGeDi input\n\n";
}
###-------------   SPAGeDi conversion finish...   ----------------------###



###--------------   Arlequin conversion start...   ----------------------###
sub Arlequin {
print "\n\nSupported later\n\n";
print "Alternative way:   You may use Convert program\\n\n";
}
###-------------   Arlequin conversion finish...   ----------------------###



###---------------   Genetix conversion start...   ----------------------###
sub Genetix {
print "\n\nSupported later\n\n";
}
###--------------   Genetix conversion finish...   ----------------------###



###----------------   Convert conversion start...   ----------------------###
sub MConvert {
my($mconvert_line1, @mconvert_matrix, $temp_50, $temp_51, $temp_52, $temp_53);
@mconvert_matrix=();
print "\nConvert conversion started...\n";
open (MCONVERT, ">convert_input.txt") || die "Convert conversion failed...\n";
print MCONVERT $project_name ."\n";
print MCONVERT "npops = " . $pop_num . "\n";
print MCONVERT "nloci = " . $mk_num . "\n";

$mconvert_line1=join("\t\t", @mk_name_list);
print MCONVERT "\t" . $mconvert_line1 . "\t\n";

foreach $temp_50 (@pop_name_list){
	print MCONVERT "pop = " . $temp_50 . "\n";
	for ($temp_51=2; $temp_51<($ind_num+2); $temp_51++){
		if ($linearray[$temp_51][1] eq $temp_50){
			$mconvert_matrix[$temp_51][0] =$linearray[$temp_51][0];
			for ($temp_52=3; $temp_52<($mk_num*2+3); $temp_52++){
				if ($linearray[$temp_51][$temp_52] eq ("?"||".")){
					$mconvert_matrix[$temp_51][$temp_52-2]="?";
				}
				else {
					$mconvert_matrix[$temp_51][$temp_52-2]=$linearray[$temp_51][$temp_52];
				}
			}
		$temp_53=join("\t", @{$mconvert_matrix[$temp_51]});
		print MCONVERT $temp_53 . "\n";
		}

	}
}

close MCONVERT;
print "Convert conversion finished...\n";
}
###---------------   Convert conversion finish...   ----------------------###



###---------------   TASSEL conversion start...   ----------------------###
sub Tassel {
my ($temp_27, $temp_28, $temp_29, $temp_30, $temp_31, $temp_32, $tassel_mk_tab, @tassel_matrix, $tassel_matrix_tab);
@tassel_matrix=();

print "TASSEL conversion start...\n";
open (TASSEL, ">Tassel_input.txt") || die "TASSEL conversion failed";
print TASSEL "$ind_num"."\t"."$mk_num".":2\n";
$tassel_mk_tab=join("\t", @mk_name_list);
print TASSEL $tassel_mk_tab ."\n";

for ($temp_27=2; $temp_27<($ind_num+2); $temp_27++){
	$tassel_matrix[$temp_27][0]=$linearray[$temp_27][0];
}

for ($temp_28=2; $temp_28<($ind_num+2); $temp_28++){
	for ($temp_29=3; $temp_29<($mk_num*2+2); $temp_29+=2){
		if ($linearray[$temp_28][$temp_29] eq ("?" || ".")){
			$temp_30="?";
		}
		else{
			$temp_30=$linearray[$temp_28][$temp_29];
		}
		if ($linearray[$temp_28][$temp_29+1] eq ("?" || ".")){
			$temp_31="?";
		}
		else {
			$temp_31=$linearray[$temp_28][$temp_29+1];
		}
		$tassel_matrix[$temp_28][($temp_29-1)/2]="$temp_30".":"."$temp_31";
	}
}

for ($temp_32=2; $temp_32<($ind_num+2); $temp_32++){
	$tassel_matrix_tab=join("\t", @{$tassel_matrix[$temp_32]});
	print TASSEL "$tassel_matrix_tab"."\n";
}

close TASSEL;
print "TASSEL conversion finished...\n";
}
###--------------   NTSYS conversion finish...   ----------------------###
