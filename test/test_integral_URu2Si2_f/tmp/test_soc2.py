import numpy as np
from module.mod_dump import dump_4dc, dump_2dc


soc_mat = np.zeros((10,10),dtype=np.complex128)
with open('soc.dat','r') as f:
    fileread = f.readlines()
    for cnt, line in enumerate(fileread):
        x,y, tmp1, tmp2 = line.strip().split()
        soc_mat[int(x)-1,int(y)-1] = complex(float(tmp1),float(tmp2))

umat = np.zeros((10,10),dtype=np.complex128)
umat_o = np.zeros((5,5),dtype=np.complex128)
su2 = np.identity(2,dtype=np.complex128)
umat_o[0,4] = 1.0
umat_o[1,2] = 1.0
umat_o[2,1] = 1.0
umat_o[3,3] = 1.0
umat_o[4,0] = 1.0
umat = np.kron(umat_o,su2)

soc_new = np.dot(np.transpose(np.conjugate(umat)),np.dot(soc_mat,umat))
umat2 = np.zeros((10,10),dtype=np.complex128)
umat2[0,0] = 1.0
umat2[2,1] = 1.0
umat2[4,2] = 1.0
umat2[6,3] = 1.0
umat2[8,4] = 1.0
umat2[1,5] = 1.0
umat2[3,6] = 1.0
umat2[5,7] = 1.0
umat2[7,8] = 1.0
umat2[9,9] = 1.0
soc_new = np.dot(np.transpose(np.conjugate(umat2)),np.dot(soc_new,umat2))


dump_2dc(10,10,soc_new,path='test_soc_zjwang.dat',prec=1.0e-6)

# transform in columns
op1 = np.array([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])
op2 = np.array([[1,0,0,0,0],[0,-1,0,0,0],[0,0,-1,0,0],[0,0,0,1,0 ],[0,0,0,0,1 ]])
op3 = np.array([[1,0,0,0,0],[0,0,-1,0,0],[0,1,0,0,0 ],[0,0,0,-1,0],[0,0,0,0,-1]])
op4 = np.array([[1,0,0,0,0],[0,0,1,0,0 ],[0,-1,0,0,0],[0,0,0,-1,0],[0,0,0,0,-1]])
op5 = np.array([[1,0,0,0,0],[0,1,0,0,0 ],[0,0,-1,0,0],[0,0,0,1,0 ],[0,0,0,0,-1]])
op6 = np.array([[1,0,0,0,0],[0,-1,0,0,0],[0,0,1,0,0 ],[0,0,0,1,0 ],[0,0,0,0,-1]])
op7 = np.array([[1,0,0,0,0],[0,0,-1,0,0],[0,-1,0,0,0],[0,0,0,-1,0],[0,0,0,0,1 ]])
op8 = np.array([[1,0,0,0,0],[0,0,1,0,0 ],[0,1,0,0,0 ],[0,0,0,-1,0],[0,0,0,0,1 ]])
op9  = op1
op10 = op2
op11 = op3
op12 = op4
op13 = op5
op14 = op6
op15 = op7
op16 = op8
#su1  = np.array([[1,0],[0,1]])
#su2  = np.array([[-1j,0],[0,1j]])
#su3  = np.array([[(1-1j)*np.sqrt(2)/2,0],[0,(1+1j)*np.sqrt(2)/2]])
#su4  = np.array([[(1+1j)*np.sqrt(2)/2,0],[0,(1-1j)*np.sqrt(2)/2]])
#su5  = np.array([[0,-1],[1,0]])
#su6  = np.array([[0,-1j],[-1j,0]])
#su7  = np.array([[0,-(1+1j)*np.sqrt(2)/2],[(1-1j)*np.sqrt(2)/2,0]])
#su8  = np.array([[0,-(1-1j)*np.sqrt(2)/2],[(1+1j)*np.sqrt(2)/2,0]])
su1  = np.array([[1,0],[0,1]])
su2  = np.array([[-1j,0],[0,1j]])
su3  = np.array([[(1-1j)*np.sqrt(2)/2,0],[0,(1+1j)*np.sqrt(2)/2]])
su4  = np.array([[(1+1j)*np.sqrt(2)/2,0],[0,(1-1j)*np.sqrt(2)/2]])
su5  = np.array([[0,-1],[1,0]])
su6  = np.array([[0,-1j],[-1j,0]])
su7  = np.array([[0,-(1+1j)*np.sqrt(2)/2],[(1-1j)*np.sqrt(2)/2,0]])
su8  = np.array([[0,-(1-1j)*np.sqrt(2)/2],[(1+1j)*np.sqrt(2)/2,0]])
su9  = su1 
su10 = su2
su11 = su3
su12 = su4
su13 = su5
su14 = su6
su15 = su7
su16 = su8
oop1 = np.kron(np.transpose(op1),np.transpose(su1))
oop2 = np.kron(np.transpose(op2),np.transpose(su2))
oop3 = np.kron(np.transpose(op3),np.transpose(su3))
oop4 = np.kron(np.transpose(op4),np.transpose(su4))
oop5 = np.kron(np.transpose(op5),np.transpose(su5))
oop6 = np.kron(np.transpose(op6),np.transpose(su6))
oop7 = np.kron(np.transpose(op7),np.transpose(su7))
oop8 = np.kron(np.transpose(op8),np.transpose(su8))
oop9 = np.kron(np.transpose(op9),np.transpose(su9))
oop10= np.kron(np.transpose(op10),np.transpose(su10))
oop11= np.kron(np.transpose(op11),np.transpose(su11))
oop12= np.kron(np.transpose(op12),np.transpose(su12))
oop13= np.kron(np.transpose(op13),np.transpose(su13))
oop14= np.kron(np.transpose(op14),np.transpose(su14))
oop15= np.kron(np.transpose(op15),np.transpose(su15))
oop16= np.kron(np.transpose(op16),np.transpose(su16))
su17 =-su1 
su18 =-su2 
su19 =-su3 
su20 =-su4 
su21 =-su5 
su22 =-su6 
su23 =-su7 
su24 =-su8 
su25 =-su9 
su26 =-su10
su27 =-su11
su28 =-su12
su29 =-su13
su30 =-su14
su31 =-su15
su32 =-su16
oop17=-oop1
oop18=-oop2
oop19=-oop3
oop20=-oop4 
oop21=-oop5 
oop22=-oop6 
oop23=-oop7 
oop24=-oop8 
oop25=-oop9 
oop26=-oop10
oop27=-oop11
oop28=-oop12
oop29=-oop13
oop30=-oop14
oop31=-oop15
oop32=-oop16

oop = [oop1,oop2,oop3,oop4,oop5,oop6,oop7,oop8,oop9,oop10,\
        oop11,oop12,oop13,oop14,oop15,oop16,oop17,oop18,oop19,oop20,\
        oop21,oop22,oop23,oop24,oop25,oop26,oop27,oop28,oop29,oop30,oop31,oop32]
ssu = [su1, su2, su3, su4, su5, su6, su7, su8, su9, su10,\
       su11,su12,su13,su14,su15,su16,su17,su18,su19,su20,\
       su21,su22,su23,su24,su25,su26,su27,su28,su29,su30,\
       su31,su32]

# no inversion
#oop = [oop1,oop2,oop3,oop4,oop5,oop6,oop7,oop8,\
#        oop17,oop18,oop19,oop20,\
#        oop21,oop22,oop23,oop24]

# no C4z
#oop = [oop1,oop2,oop5,oop6,oop9,oop10,\
#        oop13,oop14,oop17,oop18,\
#        oop21,oop22,oop25,oop26,oop29,oop30]
#oop = [oop3,oop5,oop9]
#oop = [oop7,oop8,oop15,oop16,oop23,oop24,oop31,oop32]
#oop = [oop1,oop2,oop3,oop4,oop5,oop6,oop9,oop10,\
#        oop11,oop12,oop13,oop14,oop17,oop18,oop19,oop20,\
#        oop21,oop22,oop25,oop26,oop27,oop28,oop29,oop30]
# test
##oop = [oop1,oop2,oop5,oop6,oop9,oop10,\
##        oop3,oop4,\
# oop7,oop8
##        oop11,oop12,\
##        oop13,oop14,oop17,oop18,\
# oop15,oop16,\
##        oop19,oop20,\
# oop23,oop24,\
##        oop27,oop28,\
# oop31,oop32,\
##        oop21,oop22,oop25,oop26,oop29,oop30]

#for i in range(int(len(ssu)/2)):
#    for j in range(len(ssu)):
#        if np.sum(np.abs(np.dot(ssu[i],ssu[j])-su7)) < 1.0e-6:
#            print('i = ',i,'j = ',j)
for i in range(len(oop)):
    for j in range(len(oop)):
        if np.sum(np.abs(np.dot(oop[i],oop[j])-oop7)) < 1.0e-6:
            print('i = ',i+1,'j = ',j+1)
print('su3:',su3)
print('su6:',su6)
print('su7:',su7)
print('np.dot(su3,su6)=',np.dot(su3,su6))
print('np.dot(su3,su6)-su7',np.dot(su3,su6)-su7)

soc_symm = np.zeros((10,10),dtype=np.complex128)
for i in range(len(oop)):
    soc_symm += np.dot(np.dot(np.conjugate(oop[i]),soc_mat),np.transpose(oop[i]))

soc_symm = soc_symm / float(len(oop))

error = np.sum(np.abs(soc_symm - soc_mat))
print('error = ',error)

print('test the rank:\n')
#for i in range():
# test c_4^4 = \bat{E}
#test_mat1 = np.dot(oop3,oop3)
#test_mat2 = np.dot(test_mat1,test_mat1)

# test c_4 * c_2y
#print('opp3:\n',oop3)
#print('opp5:\n',oop5)
#print('opp7:\n',oop7)
#print('c_4 * c_2y (1) =\n', np.dot(oop5,oop3)-oop7)
#print('op3\n', op3)
#print('su3\n', su3)
#print('op5\n', op5)
#print('su5\n', su5)
#print('op7\n', op7)
#print('su7\n', su7)
#print('c_4 * c_2y (2) =\n', np.kron(np.dot(op5,op3),np.dot(su5,su3)) - oop7 )


