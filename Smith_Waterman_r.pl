#! /usr/bin/perl -w
use strict;
use Term::ANSIColor;
#use simple qw(outDeal timeCosted rev);
use Getopt::Long;
my ($protein);
GetOptions(
	"protein"=>\$protein
);


if(@ARGV==0){
	print STDERR<<USE;
usage: perl $0 [-p] query_sequence ref_sequence
option: -p : add this ,if you align the protein  reads
USE
	exit;
}
my$s_query=shift;
my$s_ref=shift;
my%h_score;
my%hash_out;
my%hash_out_rev;

my$matix="./align_Matrix/phrap_default";
$matix="./align_Matrix/BLOSUM-62" if(defined( $protein ));
&matrice($matix,\%h_score);#get score matrix

my($score1,$p_1)=split(/\n/,getTraceBack(fillTable($s_query,$s_ref),$s_query,$s_ref,\%hash_out) );#deal the qurey and the ref
my($score2,$p_2)=split(/\n/,getTraceBack(fillTable(&rev($s_query),$s_ref),&rev($s_query),$s_ref,\%hash_out_rev) );#deal the reversed query and the ref
my@out_align;
if($score1>$score2){
	&trace(\%hash_out);
}
elsif($score1<$score2){
	&trace(\%hash_out_rev);
}
else{
	&trace(\%hash_out);
	&trace(\%hash_out_rev);
}
sub rev{
	my $query=@_;
	my $rev_query=reverse($query);
	$rev_query=~tr/[ATCG]/[TAGC]/;

	return $rev_query;
}
sub trace{
	my($p_h)=@_;
	foreach my$p_a(keys %$p_h){
		my$p=$$p_h{$p_a};
		my($ref,$line)=($$p[0],$$p[1]);
		while($$p[-1]!=0){
			$p=$$p[-1];
			$ref=$ref.$$p[0];
			$line=$line.$$p[1];
		}
		&out_color($ref,$line);
	}
}
sub out_color{
	my($ref,$line)=@_;#line eq query
	my@ref_a=split("",uc($ref));
	$|=1;
	my@a=split("",uc($line));
	my$out;
	my$i=0;
	my($line1,$line2,$line3)=("","","");;
	while($i<scalar@a or $i<scalar @ref_a){
		if($i%100==0 and $i!=0){
			$out=$out."$line1\n$line2\n$line3\n";
			($line1,$line2,$line3)=("","","");
		}
		if(!defined($a[$i]) or !defined($ref_a[$i]) or uc($a[$i]) ne uc($ref_a[$i]) ){
			$line1=$line1.colored($a[$i],"bold red") if(defined($a[$i]));
			$line3=$line3.colored($ref_a[$i],"blue") if(defined($ref_a[$i]));
			$line2=$line2." ";
		}
		else{
			$line1.=$a[$i];
			$line2.="|";
			$line3=$line3.colored($ref_a[$i],"blue");
		}
		$i++;
	}
	$out=$out."$line1\n$line2\n$line3\n";
	print $out;
}
sub matrice{
	my($in,$p_hash)=@_;
	open IN,"$in"||die"$!";
	my$line=<IN>;#title
	$line=<IN>;
	chomp $line;
	my@a_base=split(/\s+/,$line);
	shift @a_base;
	while($line=<IN>){
		chomp $line;
		my@a=split(/\s+/,$line);
		my$base=shift@a;
		for(my$i=0;$i<scalar@a;$i++){
			$$p_hash{$base}{$a_base[$i]}=$a[$i];
		}
	}
	close IN;
}
sub max{
	my$i=$_[0];
	for(my$j=1;$j<@_;$j++){
		$i=$_[$j] if($_[$j]>$i);
	}
	return $i;
}
sub fillTable{
	my($s_query1,$s_ref1)=@_;
	my($s_query,$s_ref)=(uc($s_query1),uc($s_ref1));
	my@a_query=split("",$s_query);
	my@a_ref=split("",$s_ref);
	my@a_table;
	my($i,$j);
	
	for($i=0;$i<=length($s_query);$i++){
		$a_table[$i][0]=0;	
	}		
	for($j=0;$j<=length($s_ref);$j++){	
		$a_table[0][$j]=0;
	}
	for($i=1;$i<=length($s_query);$i++){
		for($j=1;$j<=length($s_ref);$j++){
			my($query,$ref)=($a_query[$i-1],$a_ref[$j-1]);
			unless(defined($h_score{$ref})){
				#print STDERR "error: unrecognized word :$ref been translated to space\n";
				$ref="*";
			}
			unless(defined($h_score{$query})){
				#print STDERR "error: unrecognized word :$query been translated to space\n";
				$query="*";
			}
			$a_table[$i][$j]=&max(0,&max( $a_table[$i-1][$j]+$h_score{$query}{"*"},$a_table[$i][$j-1]+$h_score{$ref}{"*"},$a_table[$i-1][$j-1]+$h_score{$ref}{$query} ));
		}
	}
=pod
	print STDERR " \t \t",join("\t",@a_ref),"\n";
	for($i=0;$i<=length($s_query);$i++){
		if($i==0){
			print STDERR " ";
		}
		else{
			print STDERR $a_query[$i-1];
		}
		for($j=0;$j<=length($s_ref);$j++){
			print STDERR "\t",$a_table[$i][$j];
		}
		print STDERR "\n";
	}
=cut
	return \@a_table;
}
sub getTraceBack{
	my($p_table,$s_query,$s_ref,$p_hash)=@_;
	my@align;#all the best alignment to return
	my@a_query1=split("",$s_query);
	my@a_ref1=split("",$s_ref);
	my@a_query=split("",uc($s_query));
	my@a_ref=split("",uc($s_ref));
	my($i,$j)=(0,0);
	my$i_max_row=$i;
	my$i_max_col=$j;
	my$max_score=$$p_table[$i_max_row][$i_max_col];
	$|=1;
	for($i=0;$i<=length($s_query);$i++){
		for($j=0;$j<=length($s_ref);$j++){
			if($max_score<$$p_table[$i][$j]){
			#	@align=([$a_ref1[$j-1],$a_query1[$i-1],0,0,0,$i,$j,0]);#
				@align=(["","",0,0,0,$i,$j,0]);#
				$max_score=$$p_table[$i][$j];
			}
			elsif($max_score==$$p_table[$i][$j] && $max_score!=0){
				#$align[@align]=[$a_ref1[$j-1],$a_query1[$i-1],0,0,0,$i,$j,0];
				$align[@align]=["","",0,0,0,$i,$j,0];
			}
		}
	}
#	print STDERR "the maxScore is: $max_score\n";
	foreach my$p_align(@align){
		my$p_next=$p_align;
		my$s_align_query="";
		my$s_align_ref="";
		($i,$j)=($$p_align[5],$$p_align[6]);
		for(my$m=$j;$m<scalar @a_ref;$m++){
			$s_align_ref=$s_align_ref.$a_ref1[$m];
		}
		for(my$m=$i;$m<scalar @a_query;$m++){
			$s_align_query=$s_align_query.$a_query1[$m];
		}
		$$p_align[0]=$$p_align[0].$s_align_ref;
		$$p_align[1]=$$p_align[1].$s_align_query;
		&TraceBack($i,$j,$p_align,$p_table,$s_query,$s_ref,$p_hash);
	}
	return $max_score."\n".\@align;
}
sub TraceBack{
	my($i,$j,$p_next,$p_table,$s_query,$s_ref,$p_hash)=@_;
	my@align;#all the best alignment to return
	my@a_query1=split("",$s_query);
	my@a_ref1=split("",$s_ref);
	my@a_query=split("",uc($s_query));
	my@a_ref=split("",uc($s_ref));
	if($$p_table[$i][$j]!=0){
		my($query,$ref)=($a_query[$i-1],$a_ref[$j-1]);
		$ref="*" unless(defined($h_score{$ref}));
		$query="*" unless(defined($h_score{$query}));
		if($$p_table[$i-1][$j]+$h_score{$query}{"*"} == $$p_table[$i][$j]){
			$$p_next[2]=["-",$a_query[$i-1],0,0,0,$p_next];
			$i--;
			&TraceBack($i,$j,$$p_next[2],$p_table,$s_query,$s_ref,$p_hash);
		}
		if($$p_table[$i][$j-1]+$h_score{$ref}{"*"} == $$p_table[$i][$j] ){
			$$p_next[3]=[$a_ref[$j-1],"-",0,0,0,$p_next];
			$j--;
			&TraceBack($i,$j,$$p_next[3],$p_table,$s_query,$s_ref,$p_hash);
		}
		if($$p_table[$i-1][$j-1]+$h_score{$ref}{$query} == $$p_table[$i][$j] ){
			$$p_next[4]=[$a_ref[$j-1],$a_query[$i-1],0,0,0,$p_next];
			$i--;
			$j--;
			&TraceBack($i,$j,$$p_next[4],$p_table,$s_query,$s_ref,$p_hash);
		}
	}
	else{
		my($s_align_ref,$s_align_query)=("","");
		for(my$m=$j-1;$m>=0;$m--){
			$s_align_ref=$a_ref1[$m].$s_align_ref;
		}
		for(my$m=$i-1;$m>=0;$m--){
			$s_align_query=$a_query1[$m].$s_align_query;
		}
		if($i>$j){
			$s_align_ref=" "x($i-$j).$s_align_ref;
		}
		elsif($i<$j){
			$s_align_query=" "x($j-$i).$s_align_query;
		}
		$$p_next[2]=[$s_align_ref,$s_align_query,0,0,0,$p_next];
		$$p_hash{$$p_next[2]}=$$p_next[2];
	}
}
