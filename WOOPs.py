#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import linecache
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import h5py
import sys

#-------------------------------subroutine used-------------------------------------
def read_input():

    ''' a subroutine to get AO and MO number and number of kpoints and otehr information from input.woops'''

    dataset={}
    file = open('input.woops', "r")
    data = file.readlines()
    cal = "get_WOOP"
    cprec = 1e-4
    bprec = 1e-4
    pprec = 1e-2
    for line in data:
        key, value = line.split("=")
        dataset[key.strip()] = value.strip()
    if "cal" in dataset.keys():
        cal = str(dataset["cal"])
    number_MO = int(dataset["num_MO"])
    number_AO = int(dataset["num_AO"])
    number_kpts = int(dataset["num_kpts"])
    cell_dim = [float(dataset["cell_param"].split(" ")[0]),float(dataset["cell_param"].split(" ")[1]),float(dataset["cell_param"].split(" ")[2])]
    dim = str(dataset["cell_dim"])
    readmode = str(dataset["readmode"])
    if "cprec" in dataset.keys():
        cprec = float(dataset["cprec"])
    if "bprec" in dataset.keys():
        bprec = float(dataset["bprec"])
    if "pprec" in dataset.keys():
        pprec = float(dataset["pprec"])
    return number_MO, number_AO, number_kpts, cell_dim, dim, cal, readmode, cprec, bprec, pprec

def get_u_matrix(file_name,dimension_fix,dimension,num_kpoints):

    ''' a subroutine to get u matrix from seedname_u.mat, optional padding provided '''
    ''' dimension_fix -> total number of bloch bands(can be extended by padding zeros), dimension -> number of wannier functions, num_kpoints -> number of kpoints '''

    #readin data to data[]
    file = open(file_name,"r")
    content = [x.rstrip("\n") for x in file]
    data = [x.split()[:4] for x in content[3:]]
    #initialise array
    matrix = np.zeros((dimension_fix**2*num_kpoints,2)) #initialize an empty array with zeros.
    matrixs  = np.zeros((num_kpoints,dimension_fix,dimension,2))
    kpoints = np.zeros((num_kpoints,3))
    #labeling each dataset
    #######################################
    #the format of seedname_u.mat is
    # U_band1_wann1[rel] U_band1_wann1[img]
    # U_band2_wann1[rel] U_band1_wann1[img]
    # U_band3_wann1[rel] U_band1_wann1[img]
    # ...
    # U_bandn_wannm[rel] U_bandn_wannm[img]
    #######################################
    #reading Kpoints
    for k in range(num_kpoints):
        kpoints[k]=data[k*(dimension*dimension+2)]
    #differentiating img and rel part of data
    for k in range(num_kpoints):
        for i in range(k*(dimension**2+2)+1,k*(dimension**2+2)+(dimension**2+1)):
            for dat in range(2):
                matrix[i-(k*(dimension**2+2))-1][dat]=data[i][dat]
    #labeling bloch band and wannier band number
        for num_wann in range(dimension):
            for num_bands in range(dimension):
                for dat in range(2):
#                    matrixs[k][num_wann][num_bands][dat]= matrix[dimension*num_bands+num_wann][dat] # Manual tells us "row first", and its a LIE
                    matrixs[k][num_bands][num_wann][dat]= matrix[dimension*num_wann+num_bands][dat]
    return matrixs,kpoints

def get_AO(file_name,num_bands,num_wann,num_kpoints):

    ''' a subrouting to reconstruct AO(i+Jj) from get_u_matrix '''
    ''' num_bands -> number of bloch bands, num_wann -> number of wannier functions, num_kpoints -> number of kpoints'''

    #Using get_u_matrix to get wannier orbital and kpoints
    AO_raw, kpoints = get_u_matrix(file_name,num_bands,num_wann,num_kpoints)
    #initialize data array
    AO_rel = np.zeros((num_kpoints,num_wann,num_wann))
    AO_img = np.zeros((num_kpoints,num_wann,num_wann))
    AO = np.zeros((num_kpoints,num_bands,num_wann),dtype=complex)
    #Filtering useless infitestimally small data
    for k in range(num_kpoints):
        for i in range(num_wann):
            for j in range(num_wann):
                if abs(AO_raw[k][i][j][0]) > 1e-100:
                    AO_rel[k][i][j] = AO_raw[k][i][j][0]
                if abs(AO_raw[k][i][j][1]) > 1e-100:
                    AO_img[k][i][j] = AO_raw[k][i][j][1]
                #recombining data to its complex form
                AO[k][i][j] = AO_rel[k][i][j] + 1j*AO_img[k][i][j]
                #k -> kpoints, i -> bloch band number, j -> wannier_band number
    return AO,kpoints

def get_C_matrix(AO, MO, AO_wann, MO_wann, kpts, R_nb):

    ''' a subroutine to get coefficient constant from AO_U and MO_U by projecting <AO|MO> '''
    ''' AO MO -> complex matrixs, AO_wann MO_wann -> number of AO and MO, kpts -> kpoints, R_nb -> neighbour length to origin cell'''

    #initialize data array
    C_nt = np.zeros((MO_wann,AO_wann),dtype=complex)
    #Calculating overlap matrix
    for n in range(MO_wann):
        for t in range(AO_wann):
            C_nt_tmp = 0+1j*0
            #change len(kpts) to 1 to use only gamma
            for k in range(len(kpts)):
                C_nt_tmp_inner = 0+1j*0
                #Change MO_wann to AO_wann and change dimension_fix to AO_wann to use all bloch function
                for m in range(AO_wann):
                    C_nt_tmp_inner += complex(np.dot(MO[k][m][n],np.conj(AO[k][m][t])))
                # In PRB 91.195120, phase factor e^(-ik*R) is used to translate wannierfunction to each cell, R_nb is the lenght from home cell to other cell(in crystal coordinate)
                C_nt_tmp += np.dot(np.exp(1j*2*np.pi*np.dot(kpts[k],R_nb)),C_nt_tmp_inner)
            C_nt[n][t]=C_nt_tmp/len(kpts)
            #C_nt[n][t]: n -> MO_band, t -> AO_band
    return C_nt

def get_WOOP(AO, kpts, C_ij, AO_wann, MO_wann, use_nrpt, R_place):

    ''' a subroutine to get WOOP by <C_in|S_mn|C_im> '''
    ''' AO -> complex matrixs, AO_wann MO_wann -> number of AO and MO, kpts -> kpoints, C_ij -> Projection coefficient'''
    #initialize data array
    B_iml=np.zeros((MO_wann,AO_wann,AO_wann,len(use_nrpt),len(use_nrpt)),dtype=complex)
    #    S_iml = np.zeros((AO_wann,AO_wann),dtype=complex)
    #Calcualte WOOP
    for I in range(MO_wann): # Chose one MO
        for nrpt_1 in range(len(use_nrpt)): #use_nrpt[i] where i labels from  [0,0,0] to  [-1,0,0].
            for m in range(AO_wann):
                for nrpt_2 in range(len(use_nrpt)):
                    for l in range(AO_wann):
                        S_iml = 0+1j*0
                        for k in range(len(kpts)):
                            S_iml_tmp = 0+1j*0
                            for n in range(AO_wann): #n stands for bloch band number in the u matrix.
                                S_iml_tmp += np.dot(np.conj(AO[k][n][m]),AO[k][n][l])
                            S_iml += np.dot(np.exp(1j*2*np.pi*np.dot(kpts[k],(R_place[int(use_nrpt[nrpt_1])]-R_place[int(use_nrpt[nrpt_2])]))),S_iml_tmp)
                        S_iml /= len(kpts)
                        B_iml[I][m][l][nrpt_1][nrpt_2] = np.dot(np.dot(np.conj(C_ij[nrpt_1][I][m]),S_iml),C_ij[nrpt_2][I][l])
                        # B_iml in: I -> MO, m -> AO1, l -> AO2, nrpt_1 -> position of AO1, nrpt_2 -> position of AO2
    return B_iml

def get_r_matrix(file_name, num_wann):

    ''' a subroutine to get <ul|r|um> matrix from seedname_r.mat '''
    ''' num_wann -> number of MO'''

    for num_f in range(len(open(file_name,'r').readlines())):
        if linecache.getline(file_name,num_f) == "\n":
            start_line = num_f
            #print(num_f)
            break
    file = open(file_name, "r")
    content = [x.rstrip("\n") for x in file]
    nrpts = int(content[5])
    mat_j = np.zeros((nrpts,num_wann,num_wann,3),dtype=complex)
    R_place = np.zeros((nrpts,3))
    data_start = start_line-1+((num_wann*num_wann+2)*int(nrpts))-1
    for nrpt in range(nrpts):
        R_place[nrpt][0] = content[data_start+nrpt*(num_wann*num_wann+2)+2].split()[:6][0]
        R_place[nrpt][1] = content[data_start+nrpt*(num_wann*num_wann+2)+2].split()[:6][1]
        R_place[nrpt][2] = content[data_start+nrpt*(num_wann*num_wann+2)+2].split()[:6][2]
        for ao_home in range(num_wann):
            for ao_other in range(num_wann):
                mat_j[nrpt][ao_home][ao_other][0]=float(content[data_start+nrpt*(num_wann*num_wann+2)+3+ao_home+ao_other*num_wann].split()[:16][2])+1j*float(content[data_start+nrpt*(num_wann*num_wann+2)+3+ao_home+ao_other*num_wann].split()[:16][3])
                mat_j[nrpt][ao_home][ao_other][1]=float(content[data_start+nrpt*(num_wann*num_wann+2)+3+ao_home+ao_other*num_wann].split()[:16][4])+1j*float(content[data_start+nrpt*(num_wann*num_wann+2)+3+ao_home+ao_other*num_wann].split()[:16][5])
                mat_j[nrpt][ao_home][ao_other][2]=float(content[data_start+nrpt*(num_wann*num_wann+2)+3+ao_home+ao_other*num_wann].split()[:16][6])+1j*float(content[data_start+nrpt*(num_wann*num_wann+2)+3+ao_home+ao_other*num_wann].split()[:16][7])
    # mat_j is in mat_j[nrpt][l][m] = <ul_origin|r|um_nrpt>
    return mat_j, R_place, nrpts

def get_WOPP(rmat, kpts,C_ij, AO_wann, MO_wann, use_nrpt, R_place, cell_dim):

    ''' a subroutine to get WOPP by <C_in|r_mn|C_im> '''

    #initialize data array
    D_imn=np.zeros((3,MO_wann,AO_wann,AO_wann,len(use_nrpt),len(use_nrpt)),dtype=complex)
    rmat_tmp = np.zeros(3,dtype=complex)
    #Calculating WOOP
    for dir in range(3):
        for MO in range(MO_wann):
            for nrpt_1 in range(len(use_nrpt)):
                for m in range(AO_wann):
                    for nrpt_2 in range(len(use_nrpt)):
                        for n in range(AO_wann):
                            #if nrpt_2 is not home cell and nrpt_2 == nrpt_1 then both AO are in cell other than home, translation neede.
                            rmat_tmp_0 = rmat[R_place.tolist().index((-R_place[int(use_nrpt[nrpt_1])]+R_place[int(use_nrpt[nrpt_2])]).tolist())][m][n][:]#[dir] + R_place[int(use_nrpt[nrpt_1])][dir]*cell_dim[dir]
                            if nrpt_1 != 0 and nrpt_1 == nrpt_2 and m == n:# and np.real(C_ij[nrpt_1][MO][m]) > 0.1:
                                rmat_tmp[0] = rmat_tmp_0[0] + R_place[int(use_nrpt[nrpt_1])][0]*cell_dim[0]
                                rmat_tmp[1] = rmat_tmp_0[1] + R_place[int(use_nrpt[nrpt_1])][1]*cell_dim[1]
                                rmat_tmp[2] = rmat_tmp_0[2] + R_place[int(use_nrpt[nrpt_1])][2]*cell_dim[2]
                            else:
                                rmat_tmp[0] = rmat_tmp_0[0]
                                rmat_tmp[1] = rmat_tmp_0[1]
                                rmat_tmp[2] = rmat_tmp_0[2]
                            D_imn[dir][MO][m][n][nrpt_1][nrpt_2] = np.dot(np.dot(np.conj(C_ij[nrpt_1][MO][m]),rmat_tmp[dir]),C_ij[nrpt_2][MO][n])
    return D_imn


#===================================================================================#
#                                                                                   #
#                         Here comes the main program                               #
#                                                                                   #
#===================================================================================#
print(" _    _  _____  ___________\n"
      "| |  | ||  _  ||  _  | ___ \\\n"
      "| |  | || | | || | | | |_/ /__\n"
      "| |/\\| || | | || | | |  __/ __|\n"
      "\\  /\\  /\\ \\_/ /\\ \\_/ / |  \\__ \\\n"
      " \\/  \\/  \\___/  \\___/\\_|  |___/\n"
      "                 version 1.0.0\n")


# User input value start
number_MO, number_AO, number_kpts, cell_dim, dim, cal, readmo, cprec, bprec, pprec=read_input()
MO_filename="wannier90_u_MO.mat" #sys.argv[4]
AO_filename="wannier90_u_AO.mat" #sys.argv[5]
AO_r_filename="wannier90_tb.dat" #sys.argv[6]
# User input value end
cal_orb, cal_comple, cal_c, cal_charg, cal_woop, cal_wopp, readmode = False, False, False, False, False, False, False
if cal == "WOPP":
    cal_orb, cal_comple, cal_c, cal_charg, cal_woop, cal_wopp =  True, True, True, True, True, True
elif cal == "WOOP":
    cal_orb, cal_comple, cal_c, cal_charg, cal_woop = True, True, True, True, True
elif cal == "get_charge":
    cal_orb, cal_comple, cal_c, cal_charg = True, True, True, True
elif cal == "get_c_mat":
    cal_orb, cal_comple, cal_c = True, True, True
elif cal == "check_completeness":
    cal_orb, cal_comple = True, True
elif cal == "get_orbital":
    cal_orb = True
else:
    print("cal tag neede.")
    exit()
if readmo == "True":
    readmode = True
################################################################
################################################################
## Reading MO and AO and r_matrix
if cal_orb == True and readmode == False:
    sys.stdout.write('reading AO, MO, R_mat \t')
    sys.stdout.flush()
    MO, kpts = get_AO(MO_filename,number_AO,number_MO,number_kpts)
    AO, kpts = get_AO(AO_filename,number_AO,number_AO,number_kpts)
    rmat, R_place,nrpts=get_r_matrix(AO_r_filename,number_AO)
    sys.stdout.write('done \n')
    sys.stdout.flush()
    #SAVE
    with h5py.File('WAN_MAT.h5', 'w') as hf:
        hf.create_dataset("MO", data=MO)
        hf.create_dataset("AO", data=AO)
        hf.create_dataset("kpts", data=kpts)
        hf.create_dataset("rmat", data=rmat)
        hf.create_dataset("R_place", data=R_place)
        hf.create_dataset("nrpts", shape=(1,), data=nrpts)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
elif cal_orb == True and readmode == True:
    # Read save file
    f = h5py.File("WAN_MAT.h5", "r")
    MO = f['MO'][:]
    AO = f['AO'][:]
    kpts = f['kpts'][:]
    rmat = f['rmat'][:]
    R_place = f['R_place'][:]
    nrpts = f['nrpts'][0]
    f.close()
################################################################
################################################################

#get nearest neighbour nrpt number
if dim == "0D":
    nn=1
    use_nrpt=np.zeros(1)
    for nrpt in range(nrpts):
        if np.all(R_place[nrpt]==[0.,0.,0.]):
            use_nrpt[0]=nrpt

if dim == "1D":
    nn=3
    use_nrpt=np.zeros(3)
    for nrpt in range(nrpts):
        if np.all(R_place[nrpt]==[0.,0.,0.]):
            use_nrpt[0]=nrpt
        elif np.all(R_place[nrpt]==[1.,0.,0.]):
            use_nrpt[1]=nrpt
        elif np.all(R_place[nrpt]==[-1.,0.,0.]):
            use_nrpt[2]=nrpt

if dim == "2D":
    nn=5
    use_nrpt=np.zeros(nn)
    for nrpt in range(nrpts):
        if np.all(R_place[nrpt]==[0.,0.,0.]):
            use_nrpt[0]=nrpt
        elif np.all(R_place[nrpt]==[1.,0.,0.]):
            use_nrpt[1]=nrpt
        elif np.all(R_place[nrpt]==[0.,1.,0.]):
            use_nrpt[2]=nrpt
        elif np.all(R_place[nrpt]==[-1.,0.,0.]):
            use_nrpt[3]=nrpt
        elif np.all(R_place[nrpt]==[0.,-1.,0.]):
            use_nrpt[4]=nrpt
elif dim == "3D":
    nn=7
    use_nrpt=np.zeros(nn)
    for nrpt in range(nrpts):
        if np.all(R_place[nrpt]==[0.,0.,0.]):
            use_nrpt[0]=nrpt
        elif np.all(R_place[nrpt]==[1.,0.,0.]):
            use_nrpt[1]=nrpt
        elif np.all(R_place[nrpt]==[0.,1.,0.]):
            use_nrpt[2]=nrpt
        elif np.all(R_place[nrpt]==[0.,0.,1.]):
            use_nrpt[3]=nrpt
        elif np.all(R_place[nrpt]==[-1.,0.,0.]):
            use_nrpt[4]=nrpt
        elif np.all(R_place[nrpt]==[0.,-1.,0.]):
            use_nrpt[5]=nrpt
        elif np.all(R_place[nrpt]==[0.,0.,-1.]):
            use_nrpt[6]=nrpt

################################################################
################################################################
if cal_c == True and readmode == False:
    sys.stdout.write("Calculating C_mat \t")
    sys.stdout.flush()
    #Calculating C matrix
    C_nt=[]
    for i in range(len(use_nrpt)):
        C_nt.append(get_C_matrix(AO, MO, number_AO, number_MO, kpts, R_place[int(use_nrpt[i])]))  #GAM
    #C_nt now is  C_nt[nrpt][n][t] where nrpt -> cell position to origin, n -> MO, t -> AO
    sys.stdout.write("done \n")
    sys.stdout.flush()
    #SAVE
    with h5py.File('WAN_MAT.h5', 'r+') as hf:
        hf.create_dataset("C_nt", data=C_nt)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
elif cal_c == True and readmode == True:
    #READ
    cf = h5py.File("WAN_MAT.h5", "r")
    C_nt = cf['C_nt'][:]
    cf.close()
################################################################
################################################################

################################################################
################################################################
if cal_woop == True and readmode == False:
    sys.stdout.write("Calculating WOOP \t")
    sys.stdout.flush()
    #Calculating B_iml aka WOOP
    WOOP=[]
    WOOP=get_WOOP(AO, kpts,C_nt, number_AO, number_MO, use_nrpt, R_place)
    sys.stdout.write("done \n")
    sys.stdout.flush()
    # SAVE
    with h5py.File('WAN_MAT.h5', 'r+') as hf:
        hf.create_dataset("WOOP", data=WOOP)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
elif cal_woop == True and readmode == True:
    #READ
    opf = h5py.File("WAN_MAT.h5", "r")
    WOOP = opf['WOOP'][:]
    opf.close()
################################################################
################################################################

################################################################
################################################################
if cal_wopp == True and readmode == False:
    sys.stdout.write("Calculating WOPP \t")
    sys.stdout.flush()
    # Calculating D_iml aka WOPP
    WOPP = get_WOPP(rmat, kpts, C_nt, number_AO, number_MO, use_nrpt, R_place, cell_dim)
    sys.stdout.write("done \n")
    sys.stdout.flush()
    # SAVE
    with h5py.File('WAN_MAT.h5', 'r+') as hf:
        hf.create_dataset("WOPP", data=WOPP)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
elif cal_wopp == True and readmode == True:
    #READ
    opf = h5py.File("WAN_MAT.h5", "r")
    WOPP = opf['WOPP'][:]
    opf.close()
################################################################

#Printing procedure

################################################################

if cal_orb == True:
    with open('MO.txt', 'w') as f:
        print("*######################################################*",file=f)
        print("#                                                      #",file=f)
        print("#                                                      #",file=f)
        print("#                      Printing MO                     #",file=f)
        print("#                                                      #",file=f)
        print("#                                                      #",file=f)
        print("*######################################################*",file=f)

        #print only 0th wannierfunction at gamma point
        #for k in range(1):
        #    print("Kpoints= ",k,file=f)
        #    for n in range(number_MO):
        #        print("band_number= ",n,"MO_number= 1","MO=",MO[k][n][0],file=f)

        for k in range(number_kpts):
            print("*######################################################*",file=f)
            for m in range(number_MO):
                for n in range(number_AO):
                    #print("kpts_num=",k,"kpts=", kpts[k], "MO_number= ",m,"band_number= ",n,"MO=",MO[k][n][m],file=f)
                    print("kpts_num= {0:5d} | kpts= {1: 5f},{2: 5f},{3: 5f} | wann_number= {4:3d} | bloch_band_number= {5:3d} | U_mat = [{6.real: 5f} + {6.imag: 5f} i]".format(k,kpts[k][0],kpts[k][1],kpts[k][2],m,n,MO[k][n][m]),file=f)

    with open('AO.txt', 'w') as f:
        print("*######################################################*",file=f)
        print("#                                                      #",file=f)
        print("#                                                      #",file=f)
        print("#                      Printing AO                     #",file=f)
        print("#                                                      #",file=f)
        print("#                                                      #",file=f)
        print("*######################################################*",file=f)

        #print only 0th wannierfunction at gamma point
        #for k in range(1):
        #    print("Kpoints= ",k,file=f)
        #    for n in range(number_AO):
        #        print("band_number= ",n,"AO_number=9","MO=",AO[k][n][0],file=f)

        for k in range(number_kpts):
            print("*######################################################*",file=f)
            for m in range(number_AO):
                for n in range(number_AO):
                    #print("kpts_num=",k,"kpts=", kpts[k], "AO_number= ",m,"band_number= ",n,"AO=",AO[k][n][m],file=f)
                    print("kpts_num= {0:5d} | kpts= {1: 5f},{2: 5f},{3: 5f} | wann_number= {4:3d} | bloch_band_number= {5:3d} | U_mat = [{6.real: 5f} + {6.imag: 5f} i]".format(k,kpts[k][0],kpts[k][1],kpts[k][2],m,n,AO[k][n][m]),file=f)


    with open('r_mat.txt', 'w') as f:
        print("*######################################################*",file=f)
        print("#                                                      #",file=f)
        print("#                                                      #",file=f)
        print("#                  Printing r_mat                      #",file=f)
        print("#                                                      #",file=f)
        print("#                                                      #",file=f)
        print("*######################################################*",file=f)

        for nrptt in range(nrpts):
            for i in range(number_AO):
                for j in range(number_AO):
                    #print("R_place= ",R_place[nrptt],"AO_number_i= ",i,"AO_number_j= ",j,"rmat= ", rmat[nrptt][i][j],file=f)
                    print("R_place= {0: 5f},{1: 5f},{2: 5f} | AO_number_i= {3:3d} | AO_number_j= {4:3d} | r_mat= [{5.real: 5f} + {5.imag: 5f}],  [{6.real: 5f} + {6.imag: 5f}],  [{7.real: 5f} + {7.imag: 5f}] ".format(float(R_place[nrptt][0]),float(R_place[nrptt][1]),float(R_place[nrptt][2]),i,j,rmat[nrptt][i][j][0],rmat[nrptt][i][j][1],rmat[nrptt][i][j][2]),file=f)


if cal_comple == True:
    with open('completeness.txt', 'w') as f:
        print("*######################################################*",file=f)
        print("#                                                      #",file=f)
        print("#                 Checking Completeness                #",file=f)
        print("#                    Using <W_I|W_I>                   #",file=f)
        print("#                                                      #",file=f)
        print("*######################################################*",file=f)

        for n in range(number_MO):
            cpletns = 0+1j*0
            for k in range(number_kpts):
                for m in range(number_MO):
                    cpletns += np.dot(np.matrix(MO[k][m][n]).getH(),MO[k][m][n])
            cpletns /= number_kpts
            print("MO_number= [{0.real:5s}+{0.imag:5s}]".format(cpletns[0][0]),file=f)

################################################################
if cal_c == True:
    with open('C_mat.txt', 'w') as f:
        print("*######################################################*",file=f)
        print("#                                                      #",file=f)
        print("#                                                      #",file=f)
        print("#                    Printing C_nt                     #",file=f)
        print("#                                                      #",file=f)
        print("#                                                      #",file=f)
        print("*######################################################*",file=f)

        for Rl in range(nn):
            for j in range(number_MO):
                print("*##################       ",R_place[int(use_nrpt[Rl])],"       ###################*",file=f)
                for i in range(number_AO):
                    if  abs(np.real(C_nt[Rl][j][i])) > cprec:
                        print("MO= {0:5d} | AO= {1:5d} | C_nt= {2:10f}".format(j,i,np.real(C_nt[Rl][j][i])),file=f) #,j,"  AO=",i,"C_nt=",np.real(C_nt[Rl][j][i]))

if cal_woop == True:
    with open('WOOP.txt', 'w') as f:
        print("*######################################################*",file=f)
        print("#                                                      #",file=f)
        print("#                                                      #",file=f)
        print("#                   Printing B_iml                     #",file=f)
        print("#                     aka WOOP                         #",file=f)
        print("#                                                      #",file=f)
        print("*######################################################*",file=f)

        for I in range(number_MO):
            for nrpt_1 in range(len(use_nrpt)):
                for m in range(number_AO):
                    for nrpt_2 in range(len(use_nrpt)):
                        for l in range(number_AO):
                            if np.abs(np.real(WOOP[I][m][l][nrpt_1][nrpt_2])) > bprec:
                                #to ingest spin degeneracy here we muptiply 2 into the B_iml
                                print("MO= {0:2d} | AO_m= {1:2d} [{2:2.0f},{3:2.0f},{4:2.0f}] | AO_l= {5:2.0f} [{6:2.0f},{7:2.0f},{8:2.0f}] | B_iml= {9:10f}".format(I,m,R_place[int(use_nrpt[nrpt_1])][0],R_place[int(use_nrpt[nrpt_1])][1],R_place[int(use_nrpt[nrpt_1])][2],l,R_place[int(use_nrpt[nrpt_2])][0],R_place[int(use_nrpt[nrpt_2])][1],R_place[int(use_nrpt[nrpt_2])][2],np.real(WOOP[I][m][l][nrpt_1][nrpt_2])*2),file=f)

if cal_charg == True:
    with open('charge.txt', 'w') as f:
        print("*######################################################*",file=f)
        print("#                                                      #",file=f)
        print("#                                                      #",file=f)
        print("#                  Printing Charge                     #",file=f)
        print("#                                                      #",file=f)
        print("#                                                      #",file=f)
        print("*######################################################*",file=f)

        #initialize data array
        CHG = np.zeros((number_MO))
        #Calculating total charge
        for I in range(number_MO):
            for nrpt_1 in range(len(use_nrpt)):
                for m in range(number_AO):
                    for nrpt_2 in range(len(use_nrpt)):
                        for l in range(number_AO):
                            CHG[I] += np.real(WOOP[I][m][l][nrpt_1][nrpt_2])*2 #only WOOP in orig cell
            print("MO_number= {0:2d} | totalcharge= {1:10f}".format(I,CHG[I]),file=f)#,I, "totalcharge= ", CHG[I])

if cal_wopp == True:
    with open('WOPP.txt', 'w') as f:
        print("*######################################################*",file=f)
        print("#                                                      #",file=f)
        print("#                                                      #",file=f)
        print("#                  Printing D_iml                      #",file=f)
        print("#                    aka WOPP                          #",file=f)
        print("#                                                      #",file=f)
        print("*######################################################*",file=f)

        for dir in range(dim):
            total_moment = 0
            for i in range(number_MO):
                total_part = 0
                for nrpt_1 in range(len(use_nrpt)):
                    for m in range(number_AO):
                            for nrpt_2 in range(len(use_nrpt)):
                                for n in range(number_AO):
                                    total_part += np.real(WOPP[dir][i][m][n][nrpt_1][nrpt_2])
                                    total_moment += np.real(WOPP[dir][i][m][n][nrpt_1][nrpt_2])
                                    if np.abs(2*np.real(WOPP[dir][i][m][n][nrpt_1][nrpt_2])) > pprec:
                                        print("Direction= {0:2d} MO_number_i= {1:2d} AO_number_m= {2:2d}  [{3:2.0f},{4:2.0f},{5:2.0f}] AO_number_n= {6:2d} [{7:2.0f},{8:2.0f},{9:2.0f}] WOPP= {10:10f}".format(dir,i,m,R_place[int(use_nrpt[nrpt_1])][0],R_place[int(use_nrpt[nrpt_1])][1],R_place[int(use_nrpt[nrpt_1])][2],n,R_place[int(use_nrpt[nrpt_2])][0],R_place[int(use_nrpt[nrpt_2])][1],R_place[int(use_nrpt[nrpt_2])][2],-2.*np.real(WOPP[dir][i][m][n][nrpt_1][nrpt_2])),file=f)
                print("number of MO {0:2d} total_moment= {1:10f} ".format(i,total_part*-2.),file=f)
            print("Direction {0:2d} total_moment= {1:10f} ".format(dir,total_moment*-2.),file=f)

print("All calculations done, see you next time :)")
