from matrizestencil import banded
from numpy import concatenate, insert

def PiTildeDivergent(Pi, Nu, s):
    """Returns the PiTilde matrix given the Pi matrix, Nu Nullspace
    and the 's' stencil"""
    missingRowsCount = s.shape[0] - Pi.shape[0]

    A = banded(s, missingRowsCount)
    
    PiT = concatenate((Pi, A, N)).T

    return PiT

def PiTildeGradient(Pi, Nu, s):
    """Returns the PiTilde matrix given the Pi matrix, Nu Nullspace
    and the 's' stencil"""
    missingRowsCount = s.shape[0] - Pi.shape[0]

    A = banded(insert(s,0,0), missingRowsCount)
    
    PiT = concatenate((Pi, A, Nu)).T

    return PiT


# if operador == 'grad':
#     PI=np.concatenate((r,banded(int(k/2),np.insert(s,0,0)),N)).T
# elif operador == 'div':
#     PI=np.concatenate((r,banded(int(k/2)+1,s),N)).T
# 
# 
# 
# c = np.zeros(l)
# c[0] = -1
# #print("i j s[j] x")
# for i in range(int(k/2)+1,l):
#     x = 0
#     for j in range(i - int(k/2)):
#         x = x + s[j]
#         #print(i,j,s[j],-x)
#     c[i] = -x
# 
# Lprint("c = " + rowvector(rationalize(c)))
# 
# 
# # In[23]:
# 
# 
# PI_art = np.array([
#     [-11/12, 1/24, 0, 0, -1],
#     [ 17/24, -9/8, 1/24, 0 ,5],
#     [3/8, 9/8, -9/8, 1/24, -10],
#     [-5/24, -1/24, 9/8, -9/8, 10],
#     [1/24, 0, -1/24, 9/8, -5],
#     [0, 0, 0, -1/24, 1]
# ])
# 
# np.all((PI-PI_art) < 1e100)
# 
# 
# 
# PI_s = Matrix(list(map(mrationalize, PI)))
# c_s =  Matrix(rationalize(c))
# P_s = linsolve(PI_s.row_join(c_s))
# P_s
# 
# # ## Montando a matriz final
# # 
# # O último passo no algoritimo e calcular a matriz que representa o operador. Para qualquer aplicação deste algoritimo, é melhor não computar a matriz, mas simplesmente a coleção dos valores que unicamente definem os coeficientes de aproximação nos pontos de contorno, e no interior do dominio discretizado. Ainda assim para que este documento seja completo, vamos montar o operador discrto diferencial na forma matricial com o tamanho mínimo para que os valores das bordas não se cruzem.
