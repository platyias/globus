b=$1
#Make folders
mkdir /data/results/$b
mkdir /data/results/$b/profile
mkdir /data/results/$b/blast

#Blast Swissprot
/data/usearch/usearch -ublast /data/genomes/$b -db /data/database/sprot_enzyme_reviewed_ec5.fasta.udb -evalue 5e-2 -userout /data/results/$b/$b.VsSprot2_5e-2_out6_seq.blastp -userfields query+target+ql+tl+id+mism+gaps+qlo+qhi+tlo+thi+evalue+bits+qrow+trow

#Blast Profile
for a in $(cat /data/database/71_Genomes.txt); do /data/usearch/usearch -ublast /data/genomes/$b -db /data/database/71Genomes/$a -evalue 5e-2 -userout /data/results/$b/profile/$b.TO.$a.blastp -userfields query+target+evalue; done

#Blast Ort Forward
for a in $(cat /data/database/SwissCluster_Genomes.txt); do /data/usearch/usearch -ublast /data/genomes/$b -db /data/database/ClusterGenomes/$a.pep -evalue 5e-2 -maxhits 5 -userout /data/results/$b/blast/$b.vs.$a.blastp -userfields query+target+ql+tl+id+mism+gaps+qlo+qhi+tlo+thi+evalue+bits; done

#Blast Ort Reverse
for a in $(cat /data/database/SwissCluster_Genomes.txt); do /data/usearch/usearch -ublast /data/database/ClusterGenomes/$a.pep -db /data/genomes/$b -evalue 5e-2 -maxhits 5 -userout /data/results/$b/blast/$a.vs.$b.blastp -userfields query+target+ql+tl+id+mism+gaps+qlo+qhi+tlo+thi+evalue+bits; done

#Parse SwissprotBlast
python /data/scripts/01-Parse_ublast_outfmt6_2.py $b 
#Make Gene file
python /data/scripts/02-homology_gene_ec_link.py $b

#Parse ProfileBlast
for a in $(cat /data/database/71_Genomes.txt); do python /data/scripts/03-Parse_71Genome_blastp_Phylo_corr.py $b $a; done
python /data/scripts/04-Align_Phylo_Profile.py $b
python /data/scripts/05-Calc_Phylo_corr_gene_in_network.py $b
python /data/scripts/06-Context_Corr2Zscore_2.py $b pc

#Get BBHs
for a in $(cat /data/database/SwissCluster_Genomes.txt); do python /data/scripts/07-Get_BBHs.py $b $a; done
cat /data/results/$b/blast/*.bbhs > /data/results/$b/$b.blast.bbh.out

# Make Orthology
python /data/scripts/08-Get_Orthology-updated.py $b

#Calculate clustering
perl /data/scripts/09-KeepGeneBBHs.pl $b
perl /data/scripts/10-Make_Dummy_Gene_Order.pl $b
java -classpath /data/scripts -Xmx2000M GeneClusteringSplitTL $b /data/results/$b/$b.dummy_order /data/database/genomes.txt /data/database/Gene_order113/ /data/results/$b/$b.blast.bbh.filtered.out > /data/results/$b/$b.gc
python /data/scripts/11-gc_context_corr2zscore_2.py $b

#Make metafile
perl /data/scripts/MakeMetafile.pl $b;
#Run GLOBUS
/data/CGibbs/CGibbs 1 /data/results/$b/$b.metafile.txt 2

#ParseResults
perl /data/scripts/12-All_Results+Identity+Sorted.pl $b > /data/results/$b/$b.All_P_sorted+Iden.txt
