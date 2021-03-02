#simply throw all single end alignment output (both chimeric and non-chimeric together).
# That should allow normalisation based on the

#set working directory
import os 
import sys

PATH="/home/wms_stephen/Adam/STAR_output" #set working directory
os.chdir(PATH)
#===========================================================================
#read both output sam file. Also create a output file to be written
list_of_arguments = sys.argv
file =list_of_args[1]
input1= file+"_chimera_Aligned.out.sam"
Paired_aligned_read = open (input1, 'r')

input3=file+"_chimera_Single_R1_Aligned.out.sam"
R1_Aligned_read = open (input3, 'r')
input4=file+"_chimera_Single_R2_flipped_Aligned.out.sam"
R2_Aligned_read = open (input4, 'r')
input5=file+"_chimera_Single_R1_Chimeric.out.sam"
R1_Chimeric_read = open (input5, 'r')
input6=file+"_chimera_Single_R2_flipped_Chimeric.out.sam"
R2_Chimeric_read = open (input6, 'r')
sam_ouput= file+"_chimera_Chimeric_Aligned_single_merged_.out.sam"
parsed_paired_aligned_read=[]
parsed_R1_Aligned_read=[]
parsed_R2_Aligned_read=[]
parsed_R1_Chimeric_read =[]
parsed_R2_Chimeric_read =[]
#define a function that read through both sam files. Return two list of reads
def read_sam(sam1,out1):
	for rows in sam1:
		out1.append(rows)
	sam1.close()     
	return out1[4:]
parsed_paired_aligned_read= read_sam(Paired_aligned_read ,parsed_paired_aligned_read)
parsed_R1_Aligned_read=read_sam(R1_Aligned_read ,parsed_R1_Aligned_read)
parsed_R2_Aligned_read=read_sam(R2_Aligned_read ,parsed_R2_Aligned_read)
parsed_R1_Chimeric_read=read_sam(R1_Chimeric_read ,parsed_R1_Chimeric_read)
parsed_R2_Chimeric_read=read_sam(R2_Chimeric_read ,parsed_R2_Chimeric_read)
print("parsing sam file completed")
#===========================================================================

       #copy content from the original chimeric_output file into another

temp_file=file+"_chimera_Single_R1_Chimeric.out.sam" #original chiermic paired output 
with open(sam_ouput, 'w') as f:
	temp=open(temp_file,"r")
	for line in temp:
		f.write(line)

with open(sam_ouput, 'a') as f:
	for item in parsed_R2_Chimeric_read:
		f.write(item)
	for item in parsed_R1_Aligned_read:
		f.write(item) 
	for item in parsed_R2_Aligned_read:
		f.write(item)   
	  


