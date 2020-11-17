# Data assembly: using Trinity
# All samples or a subset?

# Are sequences paired or unpaired reads? Paired reads. 
# Trimomatic: élimination des débuts/fins de séquence de basse qualité
# A faire tourner à part puis vérifier qualité des données avant de lancer Trinity

# Read1: left, Read2: right
# SS lib type: a rediscuter au moment où on fera le code pour Trinity
# seqType: fq
# CPU: 4

# A ne lancer que sur un échantillon (on récupèrera les données assemblées plus tard)

# Trinity --seqType fq --left /R1.fq --right /R2.fq --max_memory 16G --CPU 4