b=$1
#Make folders
mkdir /data/results/$b
mkdir /data/results/$b/profile
mkdir /data/results/$b/blast

#Make diamond database
/data/diamond/diamond makedb --in /data/genomes/$b -d /data/genomes/$b

#Blast Swissprot
/data/diamond/diamond blastp -q /data/genomes/$b -d /data/database/sprot_enzyme_reviewed_ec5.fasta --evalue 5e-2 -o /data/results/$b/$b.VsSprot2_5e-2_out6_seq.blastp --outfmt 6 qseqid sseqid qlen slen pident mismatch gaps qstart qend sstart send evalue bitscore qseq sseq

#Blast Profile
for a in $(cat /data/database/71_Genomes.txt); do /data/diamond/diamond blastp -q /data/genomes/$b -d /data/database/71Genomes/$a --evalue 5e-2 -o /data/results/$b/profile/$b.TO.$a.blastp --outfmt 6 qseqid sseqid evalue; done

#Blast Ort Forward
for a in $(cat /data/database/SwissCluster_Genomes.txt); do /data/diamond/diamond blastp -q /data/genomes/$b -d /data/database/ClusterGenomes/$a.pep --evalue 5e-2 --max-target-seqs 5 -o /data/results/$b/blast/$b.vs.$a.blastp --outfmt 6 qseqid sseqid qlen slen pident mismatch gaps qstart qend sstart send evalue bitscore; done

#Blast Ort Reverse
for a in $(cat /data/database/SwissCluster_Genomes.txt); do /data/diamond/diamond blastp -q /data/database/ClusterGenomes/$a.pep -d /data/genomes/$b --evalue 5e-2 --max-target-seqs 5 -o /data/results/$b/blast/$a.vs.$b.blastp --outfmt 6 qseqid sseqid qlen slen pident mismatch gaps qstart qend sstart send evalue bitscore; done

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
