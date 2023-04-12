#here is some of biological functions and run it using getopt function 
from Bio import pairwise2
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import AlignIO
import sys
import getopt

def GC():
    Seq=sys.argv[2]
    Seq=Seq.upper()
    counter = 0
    for i in Seq:
        if i=='G' or i=='C':
            counter+=1

    print("GC = ", counter/len(Seq)*100)



#function 2
def transcribe(Seq=''):
    Seq=sys.argv[2]
    Seq = Seq.upper()
    compseq=''
    for i in Seq:
        if i=='G':
            compseq+='C'
        elif i=='C':
            compseq+='G'
        elif i=='A':
            compseq+='U'
        elif i=='T':
            compseq+='A'
    print(compseq)



#function 3
def reverse_complement():
    Seq=sys.argv[2]
    Seq = Seq.upper()
    compseq = ''
    for i in Seq:
        if i == 'G':
            compseq += 'C'
        elif i == 'C':
            compseq += 'G'
        elif i == 'A':
            compseq += 'T'
        elif i == 'T':
            compseq += 'A'
    print(compseq[::-1])




#function 4
def calc_nbases():
    Seq=sys.argv[2]
    Seq = Seq.upper()
    counter=0
    for i in Seq:
        if i=='N':
            counter+=1
    print( counter)

#function 5
def is_valid():
    Seq=sys.argv[2]
    type=sys.argv[3]
    Seq = Seq.upper()
    flag=False
    protein=['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','W','X','Y','Z']
    dna=['A','T','G','C','N']
    rna=['A','U','G','C','N']
    is_protein=True
    is_dna = True
    is_rna = True
    for i in Seq:
        if i not in protein:
            is_protein=False
        if i not in dna:
            is_dna = False
        if i not in rna:
            is_rna = False
    if type=='protein' and is_protein:
        flag=True
    elif type=='dna' and is_dna:
        flag=True
    elif type=='rna' and is_rna:
        flag=True
    print(flag)


#function 6
def filter_nbases():
    Seq=sys.argv[2]
    Seq = Seq.upper()
    filteredseq=''
    for i in Seq:
        if i!='N':
            filteredseq+=i

    print(filteredseq)
#function 7
def seq_alignment(seq1,seq2,out=""):
    align = " "

    alignments = pairwise2.align.globalxx(seq1,seq2)
    for alignment in alignments:
        align=format_alignment(*alignment)

    if len(out)==0:
        print(align)
    else:
        print(align, file=open(out, 'w'))




#function 8
def seq_alignment_files (seq1,seq2,out=""):
    align = " "
    # seq1=SeqIO.parse(fileseq1,"fasta")
    # seq2 = SeqIO.parse(fileseq2, "fasta")
    alignments = pairwise2.align.globalxx(seq1, seq2)
    for alignment in alignments:
        align = format_alignment(*alignment)
    if len(out)==0:
        print(align)

    else:
        print(align, file=open(out, 'w'))


def online_alignment(seq,output=""):

    result_handle = NCBIWWW.qblast("blastn","nt", seq)
    blast_result=NCBIXML.read(result_handle)
    if len(output)==0:
        for align in blast_result.alignments:
            for hsp in align.hsps:
                print("======Alignment=====")
                print('\n')
                print("sequence ",align.title)
                print("length ",align.length)
                print("e-value ",hsp.expect)
                print(hsp.query)
                print(hsp.match)
                print(hsp.sbjct)


    else:

            file = open(output, "w")
            for align in blast_result.alignments:
                for hsp in align.hsps:
                    file.write("======Alignment=====")
                    file.write('\n')
                    file.write("sequence ")
                    file.write(align.title)
                    file.write('\n')
                    file.write("length ")
                    file.write(str(align.length))
                    file.write('\n')
                    file.write("e-value ")
                    file.write(str(hsp.expect))
                    file.write('\n')
                    file.write(hsp.query)
                    file.write(hsp.match)
                    file.write(hsp.sbjct)



def convert_to_fasta():
    file = sys.argv[2]

    with open(file) as input_handle, open("ls-orchid.fasta", "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "genbank")
        SeqIO.write(sequences, output_handle, "fasta")


def merge_fasta(output="",*files):
    #seq2=SeqIO.read(file2,"fasta")
    if len(output)>1:
        for i in files:
            with open(i) as input_handle,open(output,"a") as output_handle:
                seq=SeqIO.parse(input_handle,'fasta')
                SeqIO.write(seq,output_handle,'fasta')
    else:

        for i in files:
            for seq in SeqIO.parse(i,'fasta'):
                print(seq)


if __name__=='__main__':

    fun=["GC",'transcribe','reverse_complement',"calc_nbases","filter_nbases"]
    if sys.argv[1] in fun:
        if len(sys.argv[1:])<=1:
            print("Please, enter your sequence ")
        else:
            argv = sys.argv[2]
            globals()[sys.argv[1]]()

    elif sys.argv[1]=="is_valid":
        if len(sys.argv[1:]) == 1:
            print("please enter your sequence ")
        elif len(sys.argv[1:])==2:
            print("please enter the type of your sequence")
        else:
            is_valid()

    elif sys.argv[1] == "seq_alignment":
        if len(sys.argv[1:])<2:
            print("Please, enter two sequence ")
        else:
            argv = sys.argv[4:]
            seq1 = sys.argv[2]
            seq2 = sys.argv[3]

            try:
                opts, args = getopt.gnu_getopt(argv, "o:")
                if len(opts)==0:
                    out=""
                    seq_alignment(seq1,seq2,out)
                else:
                    seq_alignment(seq1,seq2,opts[0][1])

            except:
                print("ERROR, put correct arguments! ")

    elif sys.argv[1]=="seq_alignment_files":
        if len(sys.argv[1:]) < 3:
            print("Please,enter at least two files ")
        else:
            argv = sys.argv[4:]
            fileseq1 = sys.argv[2]
            fileseq2 = sys.argv[3]

            try:
                opts, args = getopt.gnu_getopt(argv, "o:")
                if len(opts)==0:
                    out=""
                    fileseq1 = SeqIO.read(fileseq1, "fasta")
                    fileseq2 = SeqIO.read(fileseq2, "fasta")
                    seq_alignment_files(fileseq1,fileseq2,out)
                else:
                    fileseq1 = SeqIO.read(fileseq1, "fasta")
                    fileseq2 = SeqIO.read(fileseq2, "fasta")
                    seq_alignment_files(fileseq1,fileseq2,opts[0][1])

            except:
                print("ERROR, put correct arguments! ")

    elif sys.argv[1] == "online_alignment":

        if len(sys.argv[1:])<2:
            print("Please, enter your sequence ")
        else:
            argv = sys.argv[3:]
            seq = sys.argv[2]
            try:
                opts, args = getopt.gnu_getopt(argv, "o:")
                if len(opts)==0:
                    out=""
                    online_alignment(seq,out)
                else:
                    online_alignment(seq,opts[0][1])

            except:
                print("ERROR, put correct arguments! ")

    elif sys.argv[1]=='merge_fasta':

        if len(sys.argv[1:]) < 3:
            print("Please,enter at least two files ")
        else:
            argv=sys.argv[2:]
    
            try:
                opts, args = getopt.gnu_getopt(argv, "o:")
                if len(opts)==0:
                    aa=''
                    merge_fasta(aa,*args)
                else:
    
                    merge_fasta(opts[0][1],*args)

            except:
                print("ERROR, put correct arguments! ")

    elif sys.argv[1] == "convert_to_fasta":
            if len(sys.argv[1:]) == 1:
                print("please enter yout genebank file ")
            else:
                convert_to_fasta()

    elif sys.argv[1] not in fun or sys.argv[1]!="online_alignment" or sys.argv[1]!="merge_fasta" or sys.argv[1]!="seq_alignment_files":
        print("please enter valid function")
    else:
        globals()[sys.argv[1]]()
        
