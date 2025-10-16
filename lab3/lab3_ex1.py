"""
The melting temperature is the temperature at which one halF of a particular DNA duplex
will dissociate and become a single stand of DNA. 
Primar length and sequence are of critical impmortance in the design
of the parameters of a succesful amplification (pcr).
The melting temperature of a nucleic acid duplex increases both with it`s length
and with increasing G+C content.
A simple formula used for the calculation of Tm (melting temperature) is:
Tm = 4*(G+C) + 2(A+T)

The actual Tm is influenced by the concentration of Mg2+, K+, and cosolvents.
An alternative formula is:
Tm = 81.5 + 16.6(log10[Na+]) + 0.41*(%G+C) - (600/length of duplex),
where [Na+] is the ion concentration of the solution and can take a value of  0.001.
"""

"""
Implement an application that calculates the melting temperature 
of a DNA sequence using one of these formulas or both.

Input = a string of DNA

Output = temperature in celsius


Note: Upload a single file called Project_L3.zip that contains:

a) ReadMe.txt file - Here you write your name if you worked alone on the project. If you worked in a team, put the authors in the order of their involvement in the project (one below the other).

b) Screenshot.jpg file - One or more screenshots of your GUI (graphical user interface) or your output console. Several screenshot files can be uploaded (Screenshot_1.jpg, Screenshot_2.jpg, Screenshot_3.jpg ...). 

c) L3.zip file - a compressed file that contains your project - the compiled project (if you are working in a programming language) and the source code of the project).
"""

import math

def calculate_tm_simple(dna):
    """
    Calculate melting temperature using the simple formula: 
    Tm = 4*(G+C) + 2*(A+T)
    """
    dna = dna.upper()
    a = dna.count('A')
    t = dna.count('T')
    g = dna.count('G')
    c = dna.count('C')
    return 4 * (g + c) + 2 * (a + t)

def calculate_tm_advanced(dna, na_conc=0.001):
    """
    Calculate melting temperature using the advanced formula:
    Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) - 600/length
    """
    dna = dna.upper()
    length = len(dna)
    if length == 0:
        return 0
    g = dna.count('G')
    c = dna.count('C')
    gc_percent = (g + c) / length * 100
    return 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - 600 / length

#this one was added as an alternative version of the advanced formula (after talking with the professor)
def calculate_tm_advanced2(dna, na_conc=0.001):
    """
    Calculate melting temperature using the advanced formula:
    Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) - 600/length
    """
    dna = dna.upper()
    length = len(dna)
    if length == 0:
        return 0
    g = dna.count('G')
    c = dna.count('C')
    gc_percent = (g + c) / length * 100
    return 81.5 + 0.41 * gc_percent - 675 / length 


def main():
    dna_sequence = input("Enter a DNA sequence: ").strip()
    if not dna_sequence:
        print("No DNA sequence provided.")
        return
    
    tm_simple = calculate_tm_simple(dna_sequence)
    tm_advanced = calculate_tm_advanced(dna_sequence)
    
    print(f"Simple formula Tm: {tm_simple:.2f} °C")
    print(f"Advanced formula Tm: {tm_advanced:.2f} °C")
    print(f"Second version of an Advanced formula for Tm: {calculate_tm_advanced2(dna_sequence):.2f} °C")

if __name__ == "__main__":
    main()

