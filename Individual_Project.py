#ANGELIQUE VELLA_ PROJET INDIVIDUEL

#CHARGER LE FICHIER
def load_fasta(nomfich):
    """fonction qui lit un fichier contenant une sequence au format fasta possiblement sur plusieurs lignes
    input : nom du fichier fasta avec l'extension du fichier
    output : string contenant la sequence en majuscule"""
    listseq=[] #initialisation des variables
    header = ''
    seq= ''
    f=open(nomfich,"r")
    ligne=f.readline() 
    while ligne != '': #pour parcourir tout le fichier
        if ligne[0] != '>': #la ligne ne commence pas par un chevron donc c'est la sequence
           ligne=ligne.strip() #enleve les sauts de ligne \n
           listseq.append(ligne)
        else :
             header=ligne[1:].strip()#enleve le chevron et le saut à la ligne
        ligne=f.readline()
    f.close()
    seq=''.join(listseq) #créer une chaine de caractère en collant les éléments en joignant les 'vide'
    return seq.upper()


#INVERSE COMPLEMENTAIRE
def inverse_complementaire(seq):
    """fonction qui associe la base nucleique complementaire à chaque base lue dans une liste
       puis retourne les indices de la liste pour la lire dans le sens 5' vers 3'
    input : string contenant la sequence en majuscule 
    output : string contenant la reverse transcrit de la sequence en majuscule"""
    complementaire = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}#dictionnaire qui associe chaque base à sa complémentaire
    seq_complementaire= [complementaire[i]for i in seq] #parcours les bases dans la séquence et les remplaces par les complémentaires
    seq_inverse_complementaire=''.join(seq_complementaire[::-1]) #inverse et reforme ma séquence
    return seq_inverse_complementaire


#RECHERCHER LES POSITIONS DE CODON
def rechREsimple(seq, codon):
    """Fonction qui renvoie les positions d'un codon dans une séquence en utilisant une expression régulière
    Input : seq est une chaîne de caractères contenant la séquence d'ADN en majuscules
            codon est une chaîne de caractères contenant le codon à rechercher, également en majuscules
    Output : une liste d'entiers contenant les indices des occurrences du codon dans la séquence"""
    import re
    correspondances = re.finditer(codon, seq)
    positions = [match.start() for match in correspondances] # Extraire l'indice du début du codon à partir de l'itérateur
    return positions


#Question 3 : RECHERCHE SELON CODON
def codon_in_seq(seq, codon):
    """Fonction qui détermine le nombre de codons ATG dans les six phases de lecture et récupère les positions dans des listes
    Input : seq est une chaîne de caractères contenant la séquence d'ADN en majuscules
            codon est une chaîne de caractères contenant le codon à rechercher, également en majuscules
    Output : un dictionnaire qui associe chaque phase de lecture à une liste de positions du codon"""
    import re 
    sortie = {}
    debut =[]
    debut_inverse = []

    réponse = re.finditer(codon, seq)#recherche de tous les codons dans la séquence
    for match in réponse:     
        debut.append(match.start()) #liste du début de tous les codons trouvés

    for i in range(0, 3): # Parcourir les trois décalages possibles dans la séquence directe
        phase = "F" + str(i + 1)# Créer le nom de la phase de lecture
        positions = []
        for nbr in debut :#Parcourir les postions trouvées dans la séquence
            if nbr %3==i:#Determine le cadre de lecture auquel la position appartient selon le modulo
                positions.append(nbr)
        sortie[phase] = positions 
    
    seq_inverse = inverse_complementaire(seq) # Obtenir la séquence inverse complémentaire en appelant la fonction inverse_complementaire
    réponse_inverse = re.finditer(codon, seq_inverse)#recherche de tous les codons dans la séquence inverse complementaire
    for match in réponse_inverse:     
        debut_inverse.append(match.start()) #liste du début de tous les codons trouvés dans la séquence inverse complementaire

    for i in range(0, 3): # Parcourir les trois décalages possibles dans la séquence inverse
        phase = "R" + str(i + 1)# Créer le nom de la phase de lecture
        positions_inverse = []
        for nbr in debut_inverse :#Parcourir les postions trouvées dans la séquence inverse
            if nbr %3==i:#Determine le cadre de lecture auquel la position appartient selon le modulo
                positions_inverse.append(nbr)
        sortie[phase] = positions_inverse       

    for phase, positions in sortie.items():
        print("%s : %s correspondances trouvées."%(phase,len(positions)))

    return sortie


#FUSIONNER LES LISTES
def fusion (dico1 , dico2, dico3):
    """Fonction qui fusionne et trie dans l'ordre croissant 3 dictionnaires qui associent
        les 6 phases de lecture à une liste de positions du codon.
    Input :  3 dictionnaires qui associent les 6 phases de lecture à une liste de positions du codon
    Output : 1 dictionnaire qui associe les 6 phases de lecture à une liste de positions du codon."""
    STOP ={}
    phases = ["F1","F2","F3","R1", "R2", "R3"] #liste des clés à fusionner
    for phase in phases :
        STOP[phase] = dico1[phase] + dico2[phase] + dico3[phase] #Fusionne les dictionnaires
        STOP[phase]= sorted (STOP[phase]) #Trie le dictionnaire par ordre croissant
    return STOP


#TROUVER UN ORF A PARTIR DE LA POSITION DES CODONS
def ORF_in_seq(START,STOP):
    """Fonction qui recherche les ORF de plus de 150 nucléotides dans une séquence donnée,
    Input : 2 dictionnaires qui associent les 6 phases de lecture à une liste de positions des codons START et STOP
    Output : un dictionnaire qui associent l'ORF à ses coordonnées dans la séquence. """
    ORF ={}
    phases = ["F1","F2","F3","R1", "R2", "R3"]
    for phase in phases : #parcourir les phases dans les dictionnaires 
        frame_start = START[phase]
        frame_stop = STOP[phase]
        k=0
        name_liste = 'ORF_in_' + str(phase)
        liste = []
        for i in frame_start: #parcourir les positions 
            j= k
            while i< j: #Selectionne le codon start après le dernier codon stop
                i+=1
            while j < len(frame_stop) and frame_stop[j] < i : #chercher le premier codon STOP après le codon START
                j+=1
            if j < len(frame_stop): #si on a trouvé un codon STOP
                k = j
                a=frame_stop[j] - i
                if a > 152: #si la distance entre STOP et START est supérieure à 3*50 (AA)+ 2 codons du START
                    liste.append([i,frame_stop[j]])
        if liste == []: #Si aucun ORF est trouvé 
            liste = None #Attribue la valeur None à la liste des positions
        liste = [x for i, x in enumerate(liste) if x[1] not in [y[1] for y in liste[:i]]] #Supprime les doublons d'ORF qui ont les mêmes STOP
        ORF[name_liste] = liste
    return ORF


#RECHERCHE DE MOTIFS
def rechRE(seq, motif):
    """Cette fonction recherche une expression régulière dans une séquence ADN et renvoie les coordonnées des correspondances
    Input : sequence de nucléotides, un motif sous forme d'une expression régulière (module re)
    Output : un dictionnaire qui contient le nombre de correspondances trouvées, ainsi que les coordonnées et la séquence de chaque correspondance.
    La clé est "Motif n°i" où i est le numéro de la correspondance, et la valeur est une liste de listes de la forme [[début, fin], séquence].  """
    import re
    Motif_in_seq ={}

    pattern = re.compile(motif)
    matches = pattern.finditer(seq)#Recherche toutes les occurrences du motif dans la séquence

    for match in matches: #Parcours la liste de correspondances
        positions = [match.start(), match.end()] #Détermine la position de début et de fin du motif
        sequence = match.group() #Récupère la séquence du motif
        Motif_in_seq[tuple(positions)]=sequence #associe les coordonnées à la séquence de nucléotides de chaque correspondances dans un dictionnaire
    return Motif_in_seq


################################################ MAIN #######################################
import os, sys, re #Importation des modules
entree =input("Entrez le nom du fichier: ")
try:
    open("%s" %entree,"r")
except:
    print("Le fichier %s est introuvable, veuillez vérifier le nom du fichier et réessayer."%(entree))
    sys.exit()
else:
    sequence = load_fasta(entree)

    G = "\033[1m" #Affichage en gras, code d'échappement ANSI
    N = "\033[0m" #Réinitialisation du format d'impression
    print (G+"\nRecherche pour le codon START"+N)
    START = codon_in_seq(sequence,'ATG')#Recherche du codon START
    print (G + "\nRecherche pour le codon STOP dont la séquence est TAA" +N)
    TAA = codon_in_seq(sequence,'TAA')#Recherche du condon TAA
    print (G+ "\nRecherche pour le codon STOP dont la séquence est TAG" + N)
    TAG = codon_in_seq(sequence,'TAG')#Recherche du condon TAG
    print (G + "\nRecherche pour le codon STOP dont la séquence est TGA" + N)
    TGA = codon_in_seq(sequence,'TGA')#Recherche du condon TGA

    STOP = fusion(TAA,TAG,TGA)#cré un seul dictionnaire de codons STOP
    print (G+ "\nRecherche d'ORF dans la séquence:"+ N)
    ORF = ORF_in_seq(START,STOP)
    for phase, positions in ORF.items():
        print(G+ "%s (%s correspondances trouvées) :"%(phase, len(positions))+ N)
        print(f"{positions}")
    
    print (G+ "\nRecherche d'un motif dans la séquence"+ N)
    tata = "TATATA.{10,20}(A{3,5})" #Motif en expression régulière
    TATA = rechRE(sequence, tata)
    for positions, sequence in TATA.items():
        print("%s : %s"%(positions,sequence))
