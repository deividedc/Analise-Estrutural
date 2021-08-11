
import numpy as np
no=[]
material=[]
elementos=[]
arquivo = open('entrada_portico.txt', 'r')
nos=[]
carregamento=[]
num_nos=0
num_materiais=0
num_elementos = 0
i=0
for linha in arquivo:
    #leitura dos dados de controle
    if(i==2):
        controle=linha.split()
        num_nos=controle[0]
        num_nos=int(num_nos)
        num_materiais=int(controle[1])
        num_elementos=int(controle[2])

     #leitura das informações dos nós
    if(4<i<5+num_nos):
        info_nos=linha.split()
        no.append(info_nos)

     #leitura das informações dos materiais
    if(6+num_nos<i<7+num_nos+num_materiais):
         info_materiais=linha.split()
         material.append(info_materiais)

     #leitura das informações dos elementos
    if(7+num_nos+num_materiais<i<8+num_nos+num_materiais+num_elementos):
        info_elementos = linha.split()
        elementos.append(info_elementos)
    if(i>8+num_nos+num_materiais+num_elementos):
        info_carregamento = linha.split()
        carregamento.append(info_carregamento)
    i+=1
arquivo.close()
for i in range(num_nos):
    nos.append([int(no[i][3]),int(no[i][4]),int(no[i][5])])
    #definindo o numero de equações
eq = 0
for i in range(num_nos):
    for j in range(3):
        if(nos[i][j]==0):
            eq+=1
            nos[i][j]=eq
        else:
            nos[i][j]=0
        constante = eq
for i in range(num_nos):
    for j in range(3):
        if(nos[i][j]==0):
            constante+=1
            nos[i][j]=constante

#definindo o vetor LM
LM=[]
k=0
i=0
j=0
for i in range(num_elementos):
    cont=[]
    for j in range(3):
        k = int(elementos[i][1])-1
        cont.append(nos[k][j])
    for j in range(3):
        k = int(elementos[i][2])-1
        cont.append(nos[k][j])

    LM.append(cont)
#calculando o valor de EA e EJ para cada barra
EA=[]
EJ=[]
k=0
for i in range (num_elementos):
    k=int(elementos[i][3])-1
    A = float(material[k][1])
    J = float(material[k][2])
    E = int(material[k][3])
    EA.append(E*A)
    EJ.append(E*J)
#calculando o comprimento de cada barra e seus coficientes de rotação
L=[] #vetor q vai receber o comprimento de cada barra
Cx=[] #vetor que vai receber o coeficente de rotação em x de cada barra
Cy = [] #vetor que cai receber o coeficente de rotação em x de cada barra
for i in range(num_elementos):
    k1 = int(elementos[i][1])-1 #obtem o no inicial da barra
    k2 = int(elementos[i][2])-1 #obtem o no final da barra
    x1= int(no[k1][1]) #pega a coordenada X do no inicial
    x2= int(no[k2][1]) #pega a coordenada X do no final
    y1= int(no[k1][2]) #pega a coordenada y do no inicial
    y2= int(no[k2][2]) #pega a coordenada y do no final

    l = ((x2-x1)**2+(y2-y1)**2)**(1/2) #calcula o comprimento da barra i
    Cx.append((x2-x1)/l) #calcula o coeficiente de rotação em x da barra i e armazena no vetor
    Cy.append((y2-y1)/l) #calcula o coeficiente de rotação em y da barra i e armazena do vetor
    L.append(l) #armazena o comprimento da barra i no vetor

#definindo a matriz de rotação para cada barra
RT=[]
for i in range(num_elementos):
    RT.append([[Cx[i],Cy[i],0,0,0,0],[-Cy[i],Cx[i],0,0,0,0],[0,0,1,0,0,0],[0,0,0,Cx[i],Cy[i],0],[0,0,0,-Cy[i],Cx[i],0],[0,0,0,0,0,1]])
#definindo as matrizes de rigides no referencial local
KL=[]
for i in range(num_elementos):
    KL.append([[EA[i]/L[i],0,0,-EA[i]/L[i],0,0],[0,12*EJ[i]/(L[i]**3),6*EJ[i]/(L[i]**2),0,-12*EJ[i]/(L[i]**3),6*EJ[i]/(L[i]**2)],[0,6*EJ[i]/L[i]**2,4*EJ[i]/L[i],0,-6*EJ[i]/L[i]**2,2*EJ[i]/L[i]],[-EA[i]/L[i],0,0,EA[i]/L[i],0,0],[0,-12*EJ[i]/(L[i])**3,-6*EJ[i]/L[i]**2,0,12*EJ[i]/(L[i])**3,-6*EJ[i]/L[i]**2],[0,6*EJ[i]/L[i]**2,2*EJ[i]/L[i],0,-6*EJ[i]/L[i]**2,4*EJ[i]/L[i]]])

#definindo a matriz de rigidez no referencial global
KG=[]
for i in range(num_elementos):
    KG.append(np.dot(np.dot(np.transpose(RT[i]),KL[i]),RT[i]))


#criando a matriz K de valores zeros
K=[]
for i in range(eq):
    K.append([0]*eq)

#busca nas matrizes KG os valres correspondentes para cada posição da matriz K
for linha_eq in range(eq):
    for coluna_eq in range(eq):
        for i in range(num_elementos):
            for linha_lm in range(len(LM[i])):
                if((linha_eq+1)==LM[i][linha_lm]):
                    for coluna_lm in range(len(LM[i])):
                        if((coluna_eq+1)==LM[i][coluna_lm]):
                            K[linha_eq][coluna_eq]=K[linha_eq][coluna_eq]+KG[i][linha_lm][coluna_lm]
#criando a matriz K_RD de valores zeros
aux = constante-eq
K_RD=[]
for i in range(eq):
    K_RD.append([0]*aux)

#busca nas matrizes KG os valres correspondentes para cada posição da matriz K
for linha_eq in range(eq):
    for coluna_eq in range(aux):
        for i in range(num_elementos):
            for linha_lm in range(len(LM[i])):
                if((linha_eq+1)==LM[i][linha_lm]):
                    for coluna_lm in range(len(LM[i])):
                        if((coluna_eq+(eq+1))==LM[i][coluna_lm]):
                            K_RD[linha_eq][coluna_eq]=K_RD[linha_eq][coluna_eq]+KG[i][linha_lm][coluna_lm]
K_RD = np.transpose(K_RD)
#Criando o vetor f_N
f_N = [0]*constante
for i in range(num_nos):
    for j in range(len(nos[i])):
        aux = nos[i][j]-1
        forca = int(no[i][j+6])
        f_N[aux] = f_N[aux]+forca
#calculo de A_ML
A_ML=[]
q=0
V=0
M=0
barra=0
no1=0
no2=0
for i in range(len(carregamento)):
    q=int(carregamento[i][1])
    barra=int(carregamento[i][2])-1
    A_ml=[0]*6
    if(carregamento[i][0]=='1'): #carga distribuida retangular
        for j in range(len(LM[barra])):
            V=q*L[barra]/2
            M=q*L[barra]**2/12
            if(j==1):
                A_ml[1]=A_ml[1]-V
            if(j==2):
                A_ml[2]=A_ml[2]-M
            if(j==4):
                A_ml[4]=A_ml[4]-V
            if(j==5):
                A_ml[5]=A_ml[5]+M
                
    if(carregamento[i][0]=='2'): #carga distribuida triangular começando da extremidade A
        for j in range(len(LM[barra])):

            if(j==1):
                A_ml[1]=A_ml[1]-3*q*L[barra]/20
            if(j==2):
                A_ml[2]=A_ml[2]-q*L[barra]**2/30
            if(j==4):
                A_ml[4]=A_ml[4]-7*q*L[barra]/20
            if(j==5):
                A_ml[5]=A_ml[5]+q*L[barra]**2/20
            
    if(carregamento[i][0]=='3'): #carga distribuida triangular começando da extremidade B
        for j in range(len(LM[barra])):

            if(j==1):
                A_ml[1]=A_ml[1]-7*q*L[barra]/20
            if(j==2):
                A_ml[2]=A_ml[2]-q*L[barra]**2/20
            if(j==4):
                A_ml[4]=A_ml[4]-3*q*L[barra]/20
            if(j==5):
                A_ml[5]=A_ml[5]+q*L[barra]**2/30
    A_ML.append(A_ml)

#definicção do vetro{f}e
f_E = [0]*constante
for i in range(num_elementos):
    aux = np.dot((A_ML[i]),RT[i])
    aux = aux*-1
    for j in range(len(f_E)):
        for k in range(len(LM[i])):
            if(j==(LM[i][k]-1)):
                f_E[j]=f_E[j]+aux[k]

f=[]
for i in range(len(f_E)):
    f.append(f_E[i]+f_N[i])
    
f_RL=[]
for i in range(constante-eq):
    f_RL.append(-f[i+eq])


#calculando os deslocamentos
f_d = []
for i in range(eq):
    f_d.append(f[i])

d=np.linalg.solve(K,f_d)
di=[]
for i in range(num_elementos):
    d_aux=[0]*6
    for k in range(eq):
        for j in range(len(LM[i])):
            if((k+1)==(LM[i][j])):
                d_aux[j]=d[k]
    di.append(d_aux)
                
#Ações de extremidades na barra
A_M=[]
for i in range(num_elementos):
    A_M.append(np.dot(np.dot(KL[i],RT[i]),di[i]))
 
for i in range(num_elementos):
    A_M[i]=A_M[i]+A_ML[i]

f_R=f_RL+np.dot(K_RD,d)
#Saida dos dados
apoios=[]
apoio=0
n_apoios=0
no_apoio=[]
tipo_apoio = []
for i in range(len(no)):
    apoio=0
    if(no[i][3]=='1'):
        apoio+= 1
    if(no[i][4]=='1'):
        apoio+= 1
    if(no[i][5]=='1'):
        apoio+= 1
    if(apoio>0):
        n_apoios+=1
        no_apoio.append(int(no[i][0]))
        tipo_apoio.append(apoio)


for i in range(n_apoios):
    apoios.append([0]*3)

aux = []
for i in range(len(f_R)):
    aux.append(f_R[i])
for i in range(n_apoios):
    if(tipo_apoio[i]==3):
        apoios[i][0]=aux[0]
        apoios[i][1]=aux[1]
        apoios[i][2]=aux[2]
        del(aux[0])
        del(aux[0])
        del(aux[0])
    if(tipo_apoio[i]==2):
        apoios[i][0]=aux[0]
        apoios[i][1]=aux[1]
        del(aux[0])
        del(aux[0])
    if(tipo_apoio[i]==1):
        if(no[no_apoio[i]-1][3]=='1'):
            apoios[i][0]=aux[0]
            del(aux[0])
        if(no[no_apoio[i]-1][4]=='1'):
            apoios[i][1]=aux[0]
            del(aux[0])
arquivo = open('resultados.txt', 'w')
des=[]
for i in range(num_nos):
    des.append([0]*3)

for i in range(len(nos)):
    for k in range(3):
        for j in range(len(d)):
            if(nos[i][k]-1==j):
                des[i][k]=d[j]

arquivo.write("Esforços nas estremidades das barras:\n")
arquivo.write("------------------------------------------------------------------------------\n")
arquivo.write("Elemento Ni[N] Vi[N] Mi[Nm] Nf[N] Vf[N]Mf[Nm]\n")
for i in range(num_elementos):
    arquivo.write(" ")
    arquivo.write(str(i+1))
    arquivo.write(" ")

    for j in range(6):
        arquivo.write(str('{:e}'.format(A_M[i][j])))
        if(A_M[i][j]>0):
            arquivo.write(" ")
        else:
            arquivo.write(" ")
    arquivo.write("\n")

aux=[]
arquivo.write('\nVetor de deslocamentos nodais\n')
arquivo.write('---------------------------------------------------------------------------\n')
arquivo.write('nó dx[m] dy[m] dz[rad]\n')
for i in range(num_nos):
    arquivo.write('')
    arquivo.write(str(i+1))
    arquivo.write(' ')
    arquivo.write(str('{:e}'.format(des[i][0])))
    arquivo.write(' ')
    arquivo.write(str('{:e}'.format(des[i][1])))
    arquivo.write(' ')
    arquivo.write(str('{:e}'.format(des[i][2])))
    arquivo.write('\n')
    
arquivo.write("\nReações de apoio\n")
arquivo.write("------------------------------------------------------------------------------\n")
arquivo.write("Nó Rx[N] Ry[N] Rz[Nm] \n")
for i in range(n_apoios):
    arquivo.write('')
    arquivo.write(str(no_apoio[i]))
    arquivo.write(' ')
    arquivo.write(str('{:e}'.format(apoios[i][0])))
    arquivo.write(' ')
    arquivo.write(str('{:e}'.format(apoios[i][1])))
    arquivo.write(' ')
    arquivo.write(str('{:e}'.format(apoios[i][2])))
    arquivo.write('\n')
arquivo.close()