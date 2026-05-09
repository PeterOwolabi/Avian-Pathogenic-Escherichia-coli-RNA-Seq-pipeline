while read SRR
do
echo "Processing $SRR..."

prefetch $SRR
fasterq-dump $SRR --split-files -O raw_data/

fastqc raw_data/${SRR}_1.fastq raw_data/${SRR}_2.fastq -o qc_reports/

fastp -i raw_data/${SRR}_1.fastq -I raw_data/${SRR}_2.fastq -o trimmed/${SRR}_1_trim.fastq -O trimmed/${SRR}_2_trim.fastq -h qc_reports/${SRR}_fastp.html -j qc_reports/${SRR}_fastp.json

fastqc trimmed/${SRR}_1_trim.fastq trimmed/${SRR}_2_trim.fastq -o qc_reports/
echo "$SRR done."
done < samples.txt
