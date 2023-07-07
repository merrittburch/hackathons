# How do you unzip a fastq.gz file?
gunzip sample1_reverse.fastq.gz

# 2. How do you create a fastq.gz file and maintain the original uncompressed file?
touch fastq.gz

# 3. Make a directory using your name in Indigidata project directory on OSC
mkdir /fs/scratch/PAS2510/Indigidata_2023/merritt

# 4. How do you list all files in a directory, including hidden files?
ls -alh

# 5. How do you know when a file in a directory is a hidden file?
# If it has a '.' in front of it

# 6. How would you get the number of sequences in a fastq file?
wc -l /fs/ess/PAS2510/Indigidata_2023/raw-fastqc/MAY1717_A3KO_Repseq.fa

# 7. How would you read and print the Sample_location.csv file in terminal?
vim Sample_location.csv
head Sample_location.csv

# 8. How would you remove the misc. directory?
rm -r /misc

# 9. How would you move a file from ~/Raw-data to another directory?
mv /Raw-data /path/to/another/directory

# 10. What is the path for sample1_reverse.fastq.gz?
/fs/ess/PAS2510/Raw-data/

# 11. How would you create a new TAR archive containing the files in ~/Raw-data?
tar /Raw-data

# 12. How would you create a new file without using a text editor?
touch name_of_file.txt

# 13. What does the -n flag for the cat command do?
# number all output lines