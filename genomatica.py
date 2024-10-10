from Bio import SeqIO
from Bio.Seq import Seq

# MANIPULATING A FASTA INPUT
def fastaInput():
    file = input("\nInsert the full absolute path to your FASTA (.fasta) file:\n")
    if file:
        for sequence in SeqIO.parse(file, "fasta"):
            print("\n-----------------------------------------------------")
            print("#####################################################")
            print("-----------------------------------------------------\n")
            print(f"-> SEQUENCE:\n{sequence.id}\n")
            translate(sequence.seq)
    id = None

# MANIPULATING A STRING INPUT
def stringInput():
    sequence = Seq(input("\nInsert the dNTP sequence you wish to translate:\n").upper())
    if sequence:
        translate(sequence)

# TRANSLATION ALGORITHM
def translate(sequence):
    max = 0
    
    ### READING THE 3 FRAMES ON THE FORWARD STRAND 
    # STARTING BY THE FIRST NUCLEOTIDE
    frame = 0
    i = 0
    # FOR EACH OF THE 3 FRAMES
    for frame in range(3):
        # STARTING BY THE NT POSITION CORRESPONDANT TO THE FRAME BEING TRANSLATED
        i = frame
        print(f"-> READING FRAME {frame+1} (Forward Strand):",  end=" ")
        # WHILE THE NUCLEOTIDE POSITION IS NOT THE LAST 2 OF THE SEQUENCE
        while i <= len(sequence) - 3:
            # A CODON STARTS IN THE i POSITION OF THE SEQUENCE AND ALSO HAS THE NEXT TWO NTs
            codon = sequence[i:i+3]
            found = 0
            # IF THIS IS A START CODON
            if codon == "ATG":
                # MARKS THAT A START CODON WAS FOUND ON THIS FRAME
                found = 1
                # COMPLETES THE SEQUENCE THAT STARTS HERE AS A MULTIPLE OF THREE SO IT DOESN'T RAISE AND EXCEPTION
                if len(sequence[i:]) % 3 == 1:
                    sequence += "NN"
                elif len(sequence[i:]) % 3 == 2:
                    sequence += "N"
                # TRANSLATES IT ACCORDING TO THE CHOSEN GENETIC CODE AND ONLY UNTIL IT FINDS THE FIRST STOP CODON IN FRAME
                translation = sequence[i:].translate(table=chosenTable, to_stop=True)
                # PRINTS IT MARKING A STOP CODON WITH AN ASTERISK
                print(translation + '*')
                # IF THE AA LENGTH OF THIS TRANSLATION IS THE GREATEST FOR THIS SEQUENCE, IT SUGGESTS THIS IS AS THE BEST PREDICTED SYNTHESIZED PROTEIN
                if len(translation) >= max:
                    max = len(translation)
                    suggestedTranslation = translation + '*'
                    suggestedFrame = frame+1
                break
            # BUT IF THIS IS NOT A START CODON, IT GOES TO THE NEXT CODON 
            else:
                i += 3
        # IF NO START CODON WAS FOUND IN THIS FRAME, IT PRINTS "NO ORF FOUND"
        if found == 0:
            print("No ORF found")

    ### READING THE 3 FRAMES ON THE REVERSE COMPLEMENT
    sequence = sequence.reverse_complement()

    # STARTING BY THE FIRST NUCLEOTIDE
    frame = 0
    i = 0
    # FOR EACH OF THE 3 FRAMES
    for frame in range(3):
        # STARTING BY THE NT POSITION CORRESPONDANT TO THE FRAME BEING TRANSLATED
        i = frame
        print(f"-> READING FRAME {frame+4} (Reverse Complement):",  end=" ")
        # WHILE THE NUCLEOTIDE POSITION IS NOT THE LAST 2 OF THE SEQUENCE
        while i <= len(sequence) - 3:
            # A CODON STARTS IN THE i POSITION OF THE SEQUENCE AND ALSO HAS THE NEXT TWO NTs
            codon = sequence[i:i+3]
            found = 0
            # IF THIS IS A START CODON
            if codon == "ATG":
                # MARKS THAT A START CODON WAS FOUND ON THIS FRAME
                found = 1
                # COMPLETES THE SEQUENCE THAT STARTS HERE AS A MULTIPLE OF THREE SO IT DOESN'T RAISE AND EXCEPTION
                if len(sequence[i:]) % 3 == 1:
                    sequence += "NN"
                elif len(sequence[i:]) % 3 == 2:
                    sequence += "N"
                # TRANSLATES IT ACCORDING TO THE CHOSEN GENETIC CODE AND ONLY UNTIL IT FINDS THE FIRST STOP CODON IN FRAME
                translation = sequence[i:].translate(table=chosenTable, to_stop=True)
                # PRINTS IT MARKING A STOP CODON WITH AN ASTERISK
                print(translation + '*')
                # IF THE AA LENGTH OF THIS TRANSLATION IS THE GREATEST FOR THIS SEQUENCE, IT SUGGESTS THIS IS AS THE BEST PREDICTED SYNTHESIZED PROTEIN
                if len(translation) >= max:
                    max = len(translation)
                    suggestedTranslation = translation + '*'
                    suggestedFrame = frame+4
                break
            # BUT IF THIS IS NOT A START CODON, IT GOES TO THE NEXT CODON 
            else:
                i += 3
        # IF NO START CODON WAS FOUND IN THIS FRAME, IT PRINTS "NO ORF FOUND"
        if found == 0:
            print("No ORF found")

    # PRINTS OUT THE BEST TRANSLATION FOUND FOR THIS SEQUENCE  
    print(f"""\n-> Most fitting translation (table {chosenTable}, {max} AAs, reading frame {suggestedFrame}):\n{suggestedTranslation}""")

# INITIALIZATION
print(f"""
                   IAGO MOCELIN 2024

                    *** ### ### ***
               *##                  ##*
           *##                          ##*
        *##                                ##*
      *##                                    ##*
    *##                                        ##*
   *##  *######*       *######*        *######*  ##*
  *##  # |   |  #     # |   |  #      # |   |  #   ##*
 *##  #  |   |   #   #  |   |   #    #  |   |   #   ##*
 *###    |   |    # #   |   |    #  #   |   |    # ###*
 *## #             #              #              #  ##*
 *####   |   |   #  #   |   |    #  #   |   |    # ###*
 *##  #  |   |  #    #  |   |   #    #  |   |   #   ##*
  *##  # |   | #      # |   |  #      # |   |  #   ##*
   *##  *######*       *######*        *######*  ##*
    *##                                        ##*
      *##                                    ##*
        *#                                 ##*
           *##                          ##*
               *##                  ##*
                    *** ### ### ***\n""")

# USER CHOOSES BY TYPING THE NUMBER CORRESPONDANT TO THE GENETIC CODE THEY'LL USE
# THE CHOSEN NUMBER IS STORED IN THE chosenTable variable
entry = 0
while entry not in [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]:
    entry = int(input(f"""Choose which of the following Genetic Codes you wish to use for your translation:
     1 - The Standard Code
     2 - The Vertebrate Mitochondrial Code
     3 - The Yeast Mitochondrial Code
     4 - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
     5 - The Invertebrate Mitochondrial Code
     6 - The Ciliate, Dasycladacean and Hexamita Nuclear Code
     9 - The Echinoderm and Flatworm Mitochondrial Code
    10 - The Euplotid Nuclear Code
    11 - The Bacterial, Archaeal and Plant Plastid Code
    12 - The Alternative Yeast Nuclear Code
    13 - The Ascidian Mitochondrial Code
    14 - The Alternative Flatworm Mitochondrial Code
    15 - Blepharisma Nuclear Code
    16 - Chlorophycean Mitochondrial Code
    21 - Trematode Mitochondrial Code
    22 - Scenedesmus obliquus Mitochondrial Code
    23 - Thraustochytrium Mitochondrial Code
    24 - Rhabdopleuridae Mitochondrial Code
    25 - Candidate Division SR1 and Gracilibacteria Code
    26 - Pachysolen tannophilus Nuclear Code
    27 - Karyorelict Nuclear Code
    28 - Condylostoma Nuclear Code
    29 - Mesodinium Nuclear Code
    30 - Peritrich Nuclear Code
    31 - Blastocrithidia Nuclear Code
    32 - Balanophoraceae Plastid Code
    33 - Cephalodiscidae Mitochondrial UAA-Tyr Code\n"""))
chosenTable = entry

# USER CHOOSES IF THEY'LL TRANSLATE ONE OR MULTIPLE SEQUENCES FROM A .fasta FILE AND THEN fastaInput() IS CALLED,
# OR IF THEY'LL MANUALLY INPUT A SEQUENCE AND stringInput() IS CALLED
entry = ''
while entry.lower() not in ['f', 's']:
    entry = input(f"""\nChoose the type of DNA input you want to use for translation by typing the correspondent letter on your keyboard.
        F/f - FASTA File (allows for the translation of multiple sequences)
        S/s - Typing a single NT Sequence as a String\n""")
if entry.lower() == 'f':
    fastaInput()
elif entry.lower() == 's':
    stringInput()