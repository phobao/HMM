# The amino acid sequence
seq = 'DELIFFLIF'

# Initialize matrices
S = [[0] * (len(seq) + 1) for i in range(2)]
T = [[0] * (len(seq) + 1) for i in range(2)]
prev = [[0] * (len(seq) + 1) for i in range(2)]

# Initialize start probabilities
S[0][0] = 0.788136664480353
T[1][0] = 0.211863335519647

#emission probabilities
emissionPt=dict(A=0.08368888538652917, C=0.0265453849531121 , D=0.009107386101385808, E=0.009535174323700593 ,  F=0.1072397527834378, G=0.07634894009839129, H=0.007722703171261638,I=0.1287529973319524,   K=0.009208704364565626 , L=0.15801146022132412   , M=0.03163381328169854 ,   N=0.020927850139031173,P=0.028841932251854687    ,  Q=0.013936889979623772   ,      R=0.005944004773215954, S=0.06796203942406197     ,T=0.05233651172477457 ,   V=0.0982449425300296 , W=0.022571457519503766   , Y=0.041439169640545316)

emissionPs=dict(A=0.053551866265592565, C=0.01211090465613141,D=0.05813960526076878,E=0.06335074414579084     ,  F=0.046500789841608006  ,G=0.05123681327660193 ,H=0.022690545505165746,  I=0.058635904202199454    , K=0.07148823105741936 , L=0.08893737554698801   ,   M=0.021250068089793794 ,  N=0.05938943125351797 , P=0.043289977787596155   ,  Q=0.037068083741367726 , R=0.04995672515327769 , S=0.095576886995152 ,  T=0.06138067944535567   , V=0.05255321595661621 , W=0.01366032574157351 ,  Y=0.039231826077483156)

#transition probabilities
transitionP=dict(S=dict(S=0.9872140687721683 , T=0.012785931227831702), T=dict(T=0.9526608562484516, S=0.04733914375154832))

# Fill in matrices
for i in range(1, len(seq) + 1):
    # Calculate emission probabilities
    eS = 1
    eT = 1
    for aa in seq[i-1]:
        eS *= emissionPs[aa]
        eT *= emissionPt[aa]
    
    # Calculate probabilities for S state
    if S[0][i-1] * transitionP['S']['S'] * eS >= T[1][i-1] * transitionP['T']['S'] * eS:
        S[0][i] = S[0][i-1] * transitionP['S']['S'] * eS
        prev[0][i] = 0
    else:
        S[0][i] = T[1][i-1] * transitionP['T']['S'] * eS
        prev[0][i] = 1
    
    # Calculate probabilities for T state
    if T[1][i-1] * transitionP['T']['T'] * eT >= S[0][i-1] * transitionP['S']['T'] * eT:
        T[1][i] = T[1][i-1] * transitionP['T']['T'] * eT
        prev[1][i] = 1
    else:
        T[1][i] = S[0][i-1] * transitionP['S']['T'] * eT
        prev[1][i] = 0

# Determine most likely path
if S[0][len(seq)] >= T[1][len(seq)]:
    state = 0
else:
    state = 1
path = [state]
for i in range(len(seq), 0, -1):
    state = prev[state][i]
    path.append(state)
path.reverse()

# Print results
print('Viterbi matrix:')
print('S: ' + str(S))
print('T: ' + str(T))
print('Most likely state sequence: ' + ''.join(['S' if x==0 else 'T' for x in path[1:]]))

