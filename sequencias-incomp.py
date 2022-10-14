#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 13:58:58 2020

@author: miguelrocha
"""

def transcricao (dna):
    return dna.upper().replace("T","U")


# funcao traduz um codao (usando if's + expressoes regulares)
def traduzCodao (cod):
    tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A", 
      "TGT":"C", "TGC":"C",
      "GAT":"D", "GAC":"D",
      "GAA":"E", "GAG":"E",
      "TTT":"F", "TTC":"F",
      "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
      "CAT":"H", "CAC":"H",
      "ATA":"I", "ATT":"I", "ATC":"I",
      "AAA":"K", "AAG":"K",
      "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
      "ATG":"M", "AAT":"N", "AAC":"N",
      "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
      "CAA":"Q", "CAG":"Q",
      "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
      "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
      "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
      "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
      "TGG":"W",
      "TAT":"Y", "TAC":"Y",
      "TAA":"_", "TAG":"_", "TGA":"_"}
    if cod in tc:
        aa = tc[cod]
    else: aa = ""
    return aa

    
def valida (seqDNA):
    # completar

def frequencias(seqDNA):
    # completar

### PROGRAMA ####

seq_dna = input("Introduza sequencia:").upper()
if valida(seq_dna):
    print ("Transcricao: " + transcricao(seq_dna))
else:
    print("Sequencia nao valida")
