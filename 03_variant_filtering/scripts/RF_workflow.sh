#!/bin/bash

wrkdir=/PHShome/hb875/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/clean_pipeline/
cd $wrkdir

##see earler versions for creating required files##

##Generate initial datasets##
cat /data/talkowski/Samples/SFARI/deep_sv/asc_540/sv-pipeline/rdtest/split_beds/Phase1*>/PHShome/hb875/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/Phase1.bed

##Remove regions with poor sequencing just for training of data and get size to add to metrics##
##remove X && Y###
module load bedtools2/2.25.0
egrep -hv "^X|^Y" /PHShome/hb875/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/Phase1.bed| sort -k1,1 -k2,2n| coverageBed -a - -b /data/talkowski/rlc47/src/GRch37.segdups_gaps_abParts_heterochrom.lumpy.exclude.bed -sorted |awk '{print  $4 "\t" $3-$2 "\t" $NF}'|cat <(echo -e "name" '\t' "size" '\t'"poor_region_cov") ->filter_region_wSize.bed

##create overall metrics file with all metrics and variant types## 
cat /data/talkowski/Samples/SFARI/deep_sv/asc_540/sv-pipeline/03_genotyping/metrics/Phase1.*.metrics| \
    ##remove duplicate lines##	\
    awk '!seen[$0]++'| \
    ##create metric file remove variants where all samples are called to have a CNV## \
    fgrep  -v "All_samples_called_CNV_no_analysis"| \
    ####fix BAF p-value to be on log scale to match others### \
    awk '{if ($34!="NA" && $34!="BAF_KS_pval") {$34=-log($34)/log(10)} print}' | \
    ##replace inf from -log10 p-value with max -log10(P) of 300## \
    sed 's/inf/300/g' | \
    ##add poor region coverage and size as metrics at the end## \
    sort -k1,1 |join -a 1 -j 1 - filter_region_wSize.bed | \
    ##add size and coverage NA for none CNV sv types## \
    awk '{if (NF==41) print $0,"NA","NA";else print}'| \
    ##remove chr X and Y \
    awk -F"_" '{if ($3!="X" && $3!="Y" && $4!="X" && $4!="Y" ) print}' | \
    ##get rid of straggler header line## \
    tr ' ' '\t'|sort -k1,1|tail -n +2> $wrkdir/Phase1.all.metrics


##remove BAF failures (variants which BAF can not be assessed due to # of snps or ROH)##
cat $wrkdir/Phase1.all.metrics|fgrep DEL|awk '{if ($30 =="NA" || $31 =="NA") print $1 }'>$wrkdir/BAF/Phase1.BAF.noassess

cat $wrkdir/Phase1.all.metrics|fgrep DUP|awk '{if ($33 =="NA" || $34 =="NA") print $1}'>>$wrkdir/BAF/Phase1.BAF.noassess


###Start with BAF##
##Build Training set##
##del restrict to greater than 5 kb for training and in clean regions##
mkdir $wrkdir/BAF

fgrep DEL $wrkdir/Phase1.all.metrics|fgrep -wvf $wrkdir/BAF/Phase1.BAF.noassess|awk '{if ($NF<0.3 && $(NF-1)>5000 && $26<0.15) print $1,"Fail",$30,$31; else if ($NF<0.3 && $(NF-1)>5000 && $26>0.4) print $1,"Pass",$30,$31 }'|cat <(awk '{print "name","Status",$30,$31}' $wrkdir/Phase1.all.metrics|head -n 1) - |tr ' ' '\t'>$wrkdir/BAF/BAF.del.metrics

fgrep DEL $wrkdir/Phase1.all.metrics|fgrep -wvf $wrkdir/BAF/Phase1.BAF.noassess|cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/BAF/BAF.del.all.metrics

fgrep DUP $wrkdir/Phase1.all.metrics|fgrep -wvf $wrkdir/BAF/Phase1.BAF.noassess|awk '{if ($NF<0.3 && $(NF-1)>5000 && $26<0.15) print $1,"Fail",$33,$34; else if ($NF<0.3 && $(NF-1)>5000 && $26>0.4) print $1,"Pass",$33,$34 }'|cat <(awk '{print "name","Status",$33,$34}' $wrkdir/Phase1.all.metrics|head -n 1) - |tr ' ' '\t'>$wrkdir/BAF/BAF.dup.metrics

fgrep DUP $wrkdir/Phase1.all.metrics|fgrep -wvf $wrkdir/BAF/Phase1.BAF.noassess|cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/BAF/BAF.dup.all.metrics

##Run model##
##deletions##
Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/BAF/BAF.del.metrics $wrkdir/BAF/BAF.del.all.metrics 1343124 $wrkdir/BAF/BAF.del

##duplications##
Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/BAF/BAF.dup.metrics $wrkdir/BAF/BAF.dup.all.metrics 1343124 $wrkdir/BAF/BAF.dup

##PE##
##train PE on BAF results, remove any discordant##
mkdir $wrkdir/PE

##Determine training variants that are concordant with initial depth cutoff and new BAF results## 
##Remove depth only variants which have no PE metrics##
cat $wrkdir/BAF/BAF.del.pred $wrkdir/BAF/BAF.dup.pred|egrep -v "name|depth"|sort -k1,1|join -j 1 <(cat $wrkdir/BAF/BAF.del.metrics $wrkdir/BAF/BAF.dup.metrics|cut -f 1-2|sort -k1,1) -|awk '{if ($2=="Fail" && $3<0.5) print $1 "\t" $2; else if ($2=="Pass" && $3>=0.5) print $1 "\t" $2 }'>$wrkdir/PE/PE.assessment.txt

##create a training set with PE metrics##
awk '{print $1}' $wrkdir/PE/PE.assessment.txt|fgrep -wf - $wrkdir/Phase1.all.metrics|cut -f1,6-8|sort -k1,1|join -j 1  $wrkdir/PE/PE.assessment.txt -|cat <(awk '{print $1,"Status",$6,$7,$8}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/PE/PE.metrics

##print all variants to be assessed by PE metrics (exclude depth only and wham INV and BND)##
awk '{if ($1!~"depth" && !($1~"wham" && ($2~"INV" || $2~"BND"))) print}' $wrkdir/Phase1.all.metrics>$wrkdir/PE/PE.all.metrics

##Generate PE RF and make predictions##
Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/PE/PE.metrics $wrkdir/PE/PE.all.metrics 1343124 $wrkdir/PE/PE

##SR##
mkdir $wrkdir/SR

##require to pass both PE and BAF to be considered valid and fail both to be invalid###
cat $wrkdir/PE/PE.pred $wrkdir/BAF/BAF.del.pred $wrkdir/BAF/BAF.dup.pred|awk '{if ($2>=0.5) print $1 "\t" "Pass";else print $1 "\t" "Fail"}'|sort|uniq -c|awk '{if ($1==2) print $2 "\t" $3}'|sort -k1,1|join -j 1 - <(awk '{print $1,$11,$14,$17}' $wrkdir/Phase1.all.metrics|sort -k1,1)|cat <(awk '{print $1,"Status",$11,$14,$17}' $wrkdir/Phase1.all.metrics|head -n 1) - |tr ' ' '\t'>SR.metrics

##print all variants to be assessed by SR metrics (exclude depth only and wham INV and BND), NOTE: redundant with PE##

ln -s $wrkdir/PE/PE.all.metrics $wrkdir/SR/SR.all.metrics

##Run model##
Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/SR/SR.metrics $wrkdir/SR/SR.all.metrics 1343124 $wrkdir/SR/SR


##RD test##
##gt1kb & excluding depth based calls##
mkdir $wrkdir/RDgt1kb

##require PE or SR passing and BAF for training set, restrict data set to >1kb##
cat $wrkdir/PE/PE.pred $wrkdir/SR/SR.pred|sort -nrk2,2|awk '!seen[$1]++'|cat - $wrkdir/BAF/BAF.del.pred $wrkdir/BAF/BAF.dup.pred|awk '{if ($2>=0.5) print $1 "\t" "Pass";else print $1 "\t" "Fail"}'|sort|uniq -c|awk '{if ($1==2) print $2 "\t" $3}'|sort -k1,1|join -j 1 - <(awk '{if ($(NF-1)>1000) print}'  $wrkdir/Phase1.all.metrics|awk '{print $1,$26,$27,$28}'|sort -k1,1 )|cat <(awk '{print $1,"Status",$26,$27,$28}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/RDgt1kb/rd.gt1kb.nodepth.metrics

##print all variants that will be tested excluding depth based calls##
awk '{if ($(NF-1)>1000) print }' $wrkdir/Phase1.all.metrics|fgrep -v depth|egrep -w "DEL|DUP"|cat <(awk '{print}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/RDgt1kb/rd.all.gt1kb.nodepth.metrics

Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/RDgt1kb/rd.gt1kb.nodepth.metrics $wrkdir/RDgt1kb/rd.all.gt1kb.nodepth.metrics 1343124 $wrkdir/RDgt1kb/rd.gt1kb.nodepth 

##lt1kb##
mkdir $wrkdir/RDlt1kb

##require SR support and restrict to variants greater than 100 bp and less than 1000##

awk '{if ($2<0.5) print $1 "\t" "Fail";else if ($2>=0.5 ) print $1 "\t" "Pass"}' $wrkdir/SR/SR.pred|fgrep -v name|sort -k1,1|join -j 1 - <(awk '{if ($(NF-1)<=1000 && $(NF-1)>100 ) print $1,$26,$27,$28}' $wrkdir/Phase1.all.metrics|sort -k1,1)|cat <(awk '{print $1,"Status",$26,$27,$28}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/RDlt1kb/rd.lt1kb.nodepth.metrics

egrep -w "DEL|DUP" $wrkdir/Phase1.all.metrics|fgrep -v depth|awk '{if ($(NF-1)<=1000) print }'| cat <( head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/RDlt1kb/rd.all.lt1kb.nodepth

Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/RDlt1kb/rd.lt1kb.nodepth.metrics $wrkdir/RDlt1kb/rd.all.lt1kb.nodepth.metrics 1343124 $wrkdir/RDlt1kb/rd.lt1kb.nodepth 

##depth only##
mkdir $wrkdir/rd.depth.del

##DEL##
## training set from BAF variants that are made by a depth caller, variants must be in highly mappable region and greater than 5kb ##
cat $wrkdir/BAF/BAF.del.pred|fgrep depth|awk '{if ($2<0.5) print $1 "\t" "Fail";else if ($2>0.5) print $1 "\t" "Pass"}'|sort -k1,1|join -j 1 - <(sort -k1,1 $wrkdir/Phase1.all.metrics|awk '{if ($(NF-1)>=5000 && $NF<0.3) print $1,$26,$27,$28 }')|cat <(awk '{print $1,"Status",$26,$27,$28}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/rd.depth.del/rd.gt5000.depth.del.clean.metrics

## training set from BAF variants that are made by a depth caller, variants must be in poorly mapping region and greater than 5kb ##
cat $wrkdir/BAF/BAF.del.pred|fgrep depth|awk '{if ($2<0.5) print $1 "\t" "Fail";else if ($2>0.5) print $1 "\t" "Pass"}'|sort -k1,1|join -j 1 - <(sort -k1,1 $wrkdir/Phase1.all.metrics|awk '{if ($(NF-1)>=5000 && $NF>=0.3) print $1,$26,$27,$28 }')|cat <(awk '{print $1,"Status",$26,$27,$28}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/rd.depth.del/rd.gt5000.depth.del.poor.metrics

##Cleanly mapping depth variants gt 5kb##

egrep -w "DEL" $wrkdir/Phase1.all.metrics|fgrep depth|awk '{if ($(NF-1)>5000 && $NF<0.3) print }'|cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/rd.depth.del/rd.all.gt5000.depth.del.clean.metrics

##Poorly mapping depth variants gt 5kb##

egrep -w "DEL" $wrkdir/Phase1.all.metrics|fgrep depth|awk '{if ($(NF-1)>5000 && $NF>=0.3) print }'|cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/rd.depth.del/rd.all.gt5000.depth.del.poor.metrics

##Clean Random Forest##

Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/rd.depth.del/rd.gt5000.depth.del.clean.metrics $wrkdir/rd.depth.del/rd.all.gt5000.depth.del.clean.metrics 1343124 $wrkdir/rd.depth.del/rd.gt5000.depth.del.clean 

##Poor Random Forest##

Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/rd.depth.del/rd.gt5000.depth.del.poor.metrics $wrkdir/rd.depth.del/rd.all.gt5000.depth.del.poor.metrics 1343124 $wrkdir/rd.depth.del/rd.gt5000.depth.del.poor

##DUP##
mkdir $wrkdir/rd.depth.dup

## Hybrid approach: training set failing variants from BAF variants that are made by a depth caller, passing variants are gt5kb with PE/SR & BAF support  , variants must be in highly mappable region ##

fgrep depth $wrkdir/BAF/BAF.dup.pred|awk '{if ($2<0.5) print $1 "\t" "Fail"}'|sort -k1,1|join -j 1 - <(sort -k1,1 $wrkdir/Phase1.all.metrics|awk '{if ($(NF-1)>5000 && $NF<0.3) print $1,$26,$27,$28}')|tr ' ' '\t'|cat <(awk '{print $1,"Status",$26,$27,$28}' $wrkdir/Phase1.all.metrics|head -n 1) - <(fgrep Pass $wrkdir/RDgt1kb/rd.gt1kb.nodepth.metrics|awk '{print $1}'|fgrep -wf - $wrkdir/Phase1.all.metrics|fgrep DUP|awk '{if ($(NF-1)>5000 && $NF<0.3) print $1,"Pass",$26,$27,$28}')|tr ' ' '\t'>$wrkdir/rd.depth.dup/rd.gt5000.depth.dup.clean.metrics

## Hybrid approach: training set failing variants from BAF variants that are made by a depth caller, passing variants are gt5kb with PE/SR & BAF support  , variants must be in poorly mapping regions ##

fgrep depth $wrkdir/BAF/BAF.dup.pred|awk '{if ($2<0.5) print $1 "\t" "Fail"}'|sort -k1,1|join -j 1 - <(sort -k1,1 $wrkdir/Phase1.all.metrics|awk '{if ($(NF-1)>5000 && $NF>=0.3) print $1,$26,$27,$28}')|tr ' ' '\t'|cat <(awk '{print $1,"Status",$26,$27,$28}' $wrkdir/Phase1.all.metrics|head -n 1) - <(fgrep Pass $wrkdir/RDgt1kb/rd.gt1kb.nodepth.metrics|awk '{print $1}'|fgrep -wf - $wrkdir/Phase1.all.metrics|fgrep DUP|awk '{if ($(NF-1)>5000 && $NF>=0.3) print $1,"Pass",$26,$27,$28}')|tr ' ' '\t'>$wrkdir/rd.depth.dup/rd.gt5000.depth.dup.poor.metrics

##Cleanly mapping depth variants gt 5kb##

egrep -w "DUP" $wrkdir/Phase1.all.metrics|fgrep depth|awk '{if ($(NF-1)>5000 && $NF<0.3) print }'|cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/rd.depth.dup/rd.all.gt5000.depth.dup.clean.metrics

##Poorly mapping depth variants gt 5kb##

egrep -w "DUP" $wrkdir/Phase1.all.metrics|fgrep depth|awk '{if ($(NF-1)>5000 && $NF>=0.3) print }'|cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/rd.depth.dup/rd.all.gt5000.depth.dup.poor.metrics

##Clean Random Forest##

Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/rd.depth.dup/rd.gt5000.depth.dup.clean.metrics $wrkdir/rd.depth.dup/rd.all.gt5000.depth.dup.clean.metrics 1343124 $wrkdir/rd.depth.dup/rd.gt5000.depth.dup.clean 

##Poor Random Forest##

Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/rd.depth.dup/rd.gt5000.depth.dup.poor.metrics $wrkdir/rd.depth.dup/rd.all.gt5000.depth.dup.poor.metrics 1343124 $wrkdir/rd.depth.dup/rd.gt5000.depth.dup.poor 

##Put all prediction classes together###

cat $wrkdir/BAF/BAF.del.pred $wrkdir/BAF/BAF.dup.pred|awk '{print $1 "\t" $NF}'|sort -k 1,1 | \
	##join BAF and PE## \
	join -e NA -a 2 -a 1 -j 1 - <(sort -k1,1 $wrkdir/PE/PE.pred|awk '{print $1 "\t" $NF}') -o 1.1,2.1,1.2,2.2| \
	awk '{if ($1=="NA") print $2,$3,$4;else print $1,$3,$4 }'|sort -k1,1|\
	##Join BAF PE with SR## \
	join -e NA -a 2 -a 1 -j 1 -  <(sort -k1,1 $wrkdir/SR/SR.pred|awk '{print $1 "\t" $NF}') -o 1.1,2.1,1.2,1.3,2.2 | \
	awk '{if ($1=="NA") print $2,$3,$4,$5;else print $1,$3,$4,$5 }'|sort -k1,1| \
	##Add RD## \
	join -e NA -a 2 -a 1 -j 1 -  <(cat $wrkdir/RD*/rd.lt1kb.nodepth.pred \
	$wrkdir/RD*/rd.gt1kb.nodepth.pred $wrkdir/rd*/rd.gt5000.depth.*.pred| \
	awk '{print $1 "\t" $NF}'|sort -k1,1 ) -o 1.1,2.1,1.2,1.3,1.4,2.2| \
	awk '{if ($1=="NA") print $2,$3,$4,$5,$6;else print $1,$3,$4,$5,$6 }'|fgrep -v name | \
	##Add header## \
	awk '{if (NR==1) print "CNVID","BAF","PE","SR","RD" "\n" $0;else print $0}'| \
	tr ' ' '\t'>$wrkdir/metric.table

##combined p-value##

mkdir combined_prob

##Depth only clean & poor##
## if probability fails (p<0.5) use RD value, else add BAF minus 0.5 * 1-probRD to the probRD, which will either increase or decrease the probability ##

cat $wrkdir/rd.depth*/rd.gt5000.depth.*.clean.pred|awk '{print $1}'|fgrep -wf - $wrkdir/metric.table |awk '{if ($NF<.5 || $2=="NA") print $1 "\t" $NF ;else if ($2>=.5) print $1 "\t" $NF+(($2-.5)*(1-$NF));else if ($2<.5) print $1 "\t" $NF+(($2-0.5)*($NF-0.5)) }'>$wrkdir/combined_prob/Metric_table_depthonly_clean_combined.p

cat $wrkdir/rd.depth*/rd.gt5000.depth.*.poor.pred|awk '{print $1}'|fgrep -wf - $wrkdir/metric.table |awk '{if ($NF<.5 || $2=="NA") print $1 "\t" $NF ;else if ($2>=.5) print $1 "\t" $NF+(($2-.5)*(1-$NF));else if ($2<.5) print $1 "\t" $NF+(($2-0.5)*($NF-0.5)) }'>$wrkdir/combined_prob/Metric_table_depthonly_poor_combined.p

##PE/SR gt 1kb CNV##
##Pull out max PE/SR score add bonus for the other if > 0.5 and then average with RD + BAF bonus if probBAF >0.5##
awk '{print $1 }' $wrkdir/RDgt1kb/rd.gt1kb.nodepth.pred |fgrep -wf - $wrkdir/metric.table|awk '{ \
    if ($3>=$4 && $NF>=.5 && $4<.5 && ($2<.5 || $2=="NA")) print $1 "\t" ($3+$NF)/2 ; \
     else if ($4>$3 && $NF>=.5 && $3<.5 && ($2<.5 || $2=="NA")) print $1 "\t" ($4+$NF)/2 ; \
      else if ($3>=$4 && $NF>=.5 && $4>=.5 && ($2<.5 || $2=="NA")) print $1 "\t" (($3+($4-.5)*(1-$3))+$NF)/2 ; \
       else if ($4>$3 && $NF>=.5 && $3>=.5 && ($2<.5 || $2=="NA")) print $1 "\t" (($4+($3-.5)*(1-$4))+$NF)/2 ; \
        else if ($3>=$4 && $NF>=.5 && $4<.5 && $2>=.5) print $1 "\t" ($3+ ($NF+($2-.5)*(1-$NF)))/2 ; \
          else if ($4>$3 && $NF>=.5 && $3<.5 && $2>=.5) print $1 "\t" ($4+ ($NF+($2-.5)*(1-$NF)))/2 ; \
            else if ($3>=$4 && $NF>=.5 && $4<.5 && $2>=.5) print $1 "\t" ($3+ ($NF+($2-.5)*(1-$NF)))/2 ; \
             else if ($4>$3 && $NF>=.5 && $3<.5 && $2>=.5) print $1 "\t" ($4+ ($NF+($2-.5)*(1-$NF)))/2 ; \
              else if ($3>=$4 && $NF>=.5 && $4>=.5 && $2>=.5) print $1 "\t" (($3+($4-.5)*(1-$3))+ ($NF+($2-.5)*(1-$NF)))/2 ; \
               else if ($4>$3 && $NF>=.5 && $3>=.5 && $2>=.5) print $1 "\t" (($4+($3-.5)*(1-$4))+ ($NF+($2-.5)*(1-$NF)))/2 ; \
                }'>$wrkdir/combined_prob/Metric_table_PESR.rdgt1k_combined.p

##PE/SR support only no Depth##
##Pull out anything without read depth support (prob<0.5) and check if PE/SR pass, add bonus ##
cat $wrkdir/RDgt1kb/rd.gt1kb.nodepth.pred|awk '{if ($NF<0.5) print $1}'|fgrep -wf - $wrkdir/metric.table|awk '{if ($3>=$4 && $4<.5) print $1 "\t" $3 ; else if ($4>$3 && $3<.5) print $1 "\t" $4; else if ($3>=$4 && $4>=.5) print $1 "\t" ($3+($4-.5)*(1-$3));else if ($4>$3 && $3>=.5) print $1 "\t" ($4+($3-.5)*(1-$4)); }'> $wrkdir/combined_prob/Metric_table_PESR.only.gt1k_combined.p

#PE/SR lt 1kb CNV###
##Pull out max PE/SR score add bonus for the other if > 0.5 and then add an RD bonus if probRD >0.5##

awk '{print $1 }' $wrkdir/RDlt1kb/rd.lt1kb.nodepth.pred |fgrep -wf - $wrkdir/metric.table|awk '{ \
    if ($3>=$4 && $NF<.5 && $4<.5 ) print $1 "\t" $3 ; \
     else if ($4>$3 && $NF<.5 && $3<.5 ) print $1 "\t" $4 ; \
      else if ($3>=$4 && $NF<.5 && $4>=.5 ) print $1 "\t" ($3+($4-.5)*(1-$3)) ; \
       else if ($4>$3 && $NF<.5 && $3>=.5 ) print $1 "\t" ($4+($3-.5)*(1-$4)) ; \
        else if ($3>=$4 && $NF>=.5 && $4<.5 ) print $1 "\t" ($3+($NF-.5)*(1-$3)) ; \
         else if ($4>$3 && $NF>=.5 && $3<.5 ) print $1 "\t" ($4+($NF-.5)*(1-$4)) ; \
          else if ($3>=$4 && $NF>=.5 && $4>=.5 ) print $1 "\t" ($3+($NF-.5)*(1-$3)+($4-.5)*(1-$3)) ; \
           else if ($4>$3 && $NF>=.5 && $3>=.5 ) print $1 "\t" ($4+($NF-.5)*(1-$4)+($3-.5)*(1-$4)) ; \
                }'>$wrkdir/combined_prob/Metric_table_PESR.lt1k_combined.p

##BCA##
##Take max PE/SR and then add bonus##
egrep -wv "DEL|DUP|name" $wrkdir/Phase1.all.metrics|awk '{print $1}'|fgrep -wf - $wrkdir/metric.table|cut -f1,3,4|awk '{if ($2>=$3 && $3<.5) print $1 "\t" $2 ; else if ($3>$2 && $2<.5) print $1 "\t" $3; else if ($2>=$3 && $3>=.5) print $1 "\t" ($2+($3-.5)*(1-$2));else if ($3>$2 && $2>=.5) print $1 "\t" ($3+($2-.5)*(1-$3)); }' |cat <(echo -e "ID" '\t' "Probability") ->  $wrkdir/combined_prob/Metric_table_BCA_combined.p


cat Metric_table_PESR.lt1k_combined.p Metric_table_PESR.rdgt1k_combined.p Metric_table_depthonly_clean_combined.p Metric_table_depthonly_poor_combined.p|awk '{if ($NF>=.5) print}'|cat <(echo -e "ID" '\t' "Probability") ->combined_CNV.passing.p


cat Metric_table_PESR.only.gt1k_combined.p  $wrkdir/combined_prob/Metric_table_BCA_combined.p|awk '{if ($NF>=.5) print}'|cat <(echo -e "ID" '\t' "Probability") - >combined_BCA.passing.p
