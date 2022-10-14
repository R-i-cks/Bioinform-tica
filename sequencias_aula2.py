
def transcricao (dna):
    """ Funcao que calcula a transcricao 
    de uma sequencia de DNA (argumento dna)
    Retorna uma string com a sequencia de RNA."""
    assert validaDNA(dna), "Sequencia invalida"
    return dna.upper().replace("T","U")

# funcao que calcula o complemento inverso de uma sequencia de DNA
def compinverso(dna):
    comp = ""
    for c in dna.upper():
        if c == 'A':
            comp = "T" + comp
        elif c == "T": 
            comp = "A" + comp
        elif c == "G": 
            comp = "C" + comp
        elif c== "C": 
            comp = "G" + comp
    return comp

# funcao que determina se sequencia de DNA e valida
def validaDNA (seqDNA):
    seqM = seqDNA.upper()
    for c in seqM:
        if c not in "ACGT": 
            return False
    
    return True


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


# funcao que traduz uma sequencia de DNA, a partir de uma posicao
def traduzSeq (seqDNA, iniPos= 0):
    seqM = seqDNA.upper()
    seqAA = ""
    for pos in range(iniPos,len(seqM)-2,3):
        codao = seqM[pos:pos+3]
        seqAA += traduzCodao(codao)
    return seqAA

# funcao que determina as 6 open reading frames associadas a uma sequencia de DNA 
def orfs (seqDNA):
    res = []
    res.append(traduzSeq(seqDNA,0))
    res.append(traduzSeq(seqDNA,1))
    res.append(traduzSeq(seqDNA,2))
    compinv = compinverso(seqDNA)
    res.append(traduzSeq(compinv,0))
    res.append(traduzSeq(compinv,1))
    res.append(traduzSeq(compinv,2))    
    return res

def todasProteinas (seqAA):
    seqAA = seqAA.upper() 
    protsAtuais = []
    proteinas = []
    for aa in seqAA:
        if aa == "_":
            if protsAtuais:
                for p in protsAtuais:
                    proteinas.append(p)
                protsAtuais = []
        else:
            if aa == "M":
                protsAtuais.append("")
            for i in range(len(protsAtuais)):
                protsAtuais[i] += aa
    return proteinas

def todasProteinasORFs (seqDNA):
    porfs = orfs(seqDNA) # lista de listas
    res = []
    for orf in porfs: # a cada lista 
        prots = todasProteinas(orf) # vamos traduzir cada lista em proteínas
        for p in prots: res.append(p) # Vamos tirar todas as proteínas na forma de uma lista
    return res  # É uma lista



def procuraPadrao (seq, pad):
    res = -1
    pos = 0
    while pos <= (len(seq) - len(pad)) and res < 0:
        if seq[pos:(pos+len(pad))] == pad:
            res = pos
        else: pos += 1
    return res
            
def procuraPadraoTodas (seq, pad):
    res = []
    for pos in range(len(seq) - len(pad) + 1):
        if seq[pos:(pos+len(pad))] == pad:
            res.append(pos)
    return res



# funcao para testar
def test():
    seq_dna = input("Introduza sequencia:").upper()

    if validaDNA(seq_dna):
        print ("Sequencia valida")
        print ("Transcricao: " + transcricao(seq_dna))
        print ("Comp. inverso:" + compinverso(seq_dna))
        print ("Traducao: " + traduzSeq(seq_dna))
        print ("ORFs:")
        for orf in orfs(seq_dna): print (orf)
        print ("Proteinas nas ORFs:")
        print (todasProteinasORFs(seq_dna))
    else:
        print ("Sequencia invalida")

if __name__ == '__main__': 
    test()

# --------------------------------------------------------------------------------------------------------------------- #
# Ficha 2 

# Exercício 1

def base_count(dna):
    bases = {'A':0 ,'C':0,'G':0,'T':0}
    for elem in dna.upper():
        if elem in bases:
            bases[elem] +=1
    return bases


# Exercício 2

def e_palindromo(palavra):
    i = 0
    inv =''
    while i < len(palavra):
        inv = palavra[i] + inv
        i += 1
    return palavra == inv

# Exercício 3 

def aminoacidos(dna):
    codoes = traduzSeq(dna,0)
    dic = {}
    print(codoes)
    for elem in codoes:
        if elem in dic:
            dic[elem] += 1
        else:
            dic[elem] = 1
    return dic
    

# Exercício 4

def stop_check(dna):
    codoes = traduzSeq(dna,0)
    return(codoes.find('_'))

# Exercício 5

def all_stop(dna):
    codoes = traduzSeq(dna,0)
    pos = []
    i = 0
    while i< len(codoes):
        if codoes[i] == '_':
            pos.append(i)
        i +=1
    return pos
   
#Exercício 6

def a_maior_p(amino): 
    j=0
    i=0
    enc = False
    while j< len(amino) and not enc:
        if amino[j] == '_':
            enc = j
        j+=1
    fim =[]
    while i < len(amino):
        if amino[i] =='_':
            fim.append(i)
        i+=1
    return(amino[j:max(fim)],enc,fim)


#Exercício 7

def ret_prt(dna):
    proteinas =[]
    pR = todasProteinasORFs(dna)  # lista com todas as proteínas possíveis
    return(sorted(lista_prot, key = len))

