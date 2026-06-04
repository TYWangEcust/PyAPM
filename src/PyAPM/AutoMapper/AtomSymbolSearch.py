#The purpose of this code is to find possible atom pairs involved in the reaction 
#and identify their atom IDs before and after the reaction for use in automap. 
#It does not support other reaction types, such as atom deletion or addition.
#Added AtomSymbolSearch.py, modified parts of PathSearch.py and AutoMapper.py
#Example: 
#python AutoMapper.py . newmap cleanedpre_reaction.data cleanedpost_reaction.data --save_name pre-molecule.data post-molecule.data --ba N C --ebt H H C C N O O S
#-ba means the element between the react bond, order does not matter
#-ebt means the element in cleanedpre_reaction.data and cleanedpost_reaction.data

import os
import sys
from glob import glob
from LammpsSearchFuncs import get_data, find_sections, get_neighbours
from LammpsTreatmentFuncs import clean_data,refine_data
from MapProcessor import map_processor


def get_neighbours_element(atomIDs,potencial_atomID,types):
    neighboursindexDicts = []
    atomIDs_keys = list(atomIDs.keys())
    potencial_stand_by_atomID = []
    neighboursindexDicks = {}
    for n in range(0,len(atomIDs)):
        for i in range(0,len(potencial_atomID)):
            if potencial_atomID[i][0] == atomIDs_keys[n]:
                potencial_stand_by_atomID.append(potencial_atomID[i][0])
                neighboursindexDicks = {key: atomIDs[key] for key in potencial_stand_by_atomID}  

    #Get all atoms ID probably react and their neighbors: neighboursindexDicks{key 1：atomIDs[2，18，19]}
    potencial_atomIDs_values_len=[]
    neighboursindexDicts_sb = []
    potencial_atomIDs_keys = list(neighboursindexDicks.keys())
    potencial_atomIDs_values =list(neighboursindexDicks.values())
    for i in range(0,len(potencial_atomIDs_keys)):
        for n in potencial_atomIDs_values[i]:
            for m in range(0,len(types)):
                if types[m][0] == n:
                    neighboursindexDicts_sb.append(types[m][1])
        potencial_atomIDs_values_len.append(len(potencial_atomIDs_values[i]))
        neighboursindexDicts = split_list(neighboursindexDicts_sb,potencial_atomIDs_values_len)
    return neighboursindexDicts, neighboursindexDicks

def split_list(lst, chunk_sizes):  

    start = 0  
    chunks = []  
    for size in chunk_sizes:  
        if start + size <= len(lst):  
            chunks.append(lst[start:start+size])  
            start += size  
        else:  
            chunks.append(lst[start:])  
            break  
    return chunks  

def  get_atomID(potencial_elementID,types):
    potencial_atomIDs = [] 
    for i in range(len(types)):
        for n in range(len(potencial_elementID)):
            if potencial_elementID[n] == types[i][1]:
                potencial_atomIDs.append(types[i])
                continue
    return potencial_atomIDs 


def remove_common_sublists(list1, list2):  

    list3 =list1[:]
    list4 =list2[:]   
    while i < len(list3):  
        for j in range(len(list4)):  
            if list3[i] == list4[j]:  
                del list3[i]  
                del list4[j]  
                break  
        else:  
            i += 1  
        
    return list3, list4  

def compare_atomID_difference(
    potencial_pre_atomID,
    potencial_post_atomID,
    pre_react_1,
    pre_react_2,
    preneighboursindex_Dict,
    postneighboursindex_Dict,
    preindexdict,
    postindexdict,
    potencial_1_moleculeID,
    potencial_2_moleculeID):
    #potencial_pre_atomID is the index and element of atoms that may react [1, 5]
    #pre_react_1,2 are the indices of confirmed elements before the reaction [5]
    #preneighboursindex_dict is the element types of neighbors of atoms that may react [4,1,1]
    #preindexdict is the atom IDs of an atom and its neighbors {1: [2, 18, 19]}
    #preatomIDdict is the atom IDs of atoms that may react [1, 2, 45]
    #potencial_1_moleculeID and potencial_2_moleculeID are the atom IDs within the two molecules in the system before the reaction [1...29]

    #find atomID of post-reaction atoms
    postindex_dict_key = list(postindexdict.keys())
    postindex_dict_value = list(postindexdict.values())
    post_react_atomID_1_sb = []
    post_react_atomID_1 = []
    post_react_atomID_2 = []
    indexs = []
    atomID_at_indics = []
    for i in range(len(postindex_dict_key)):
        for m in range(len(pre_react_1)):
            if potencial_post_atomID[i][1] == pre_react_1[m]:
                post_react_atomID_1_sb.append(potencial_post_atomID[i][0])
                indexs.append(i)   
    elements_at_indices = [postneighboursindex_Dict[index] for index in indexs]
    atomID_at_indics = [postindex_dict_value[index] for index in indexs]
    for p in range(len(elements_at_indices)):
        for n in range(len(elements_at_indices[p])):
            for q in range(len(pre_react_2)):
                if elements_at_indices[p][n] == pre_react_2[q]:
                    post_react_atomID_1.append(post_react_atomID_1_sb[p])
                    post_react_atomID_2.append(atomID_at_indics[p][n])

    #find atomID of pre-reaction atoms
    preindex_dict_key = list(preindexdict.keys())
    preindex_dict_value = list(preindexdict.values())
    pre_react_atomID_1_sb = []
    pre_react_atomID_2_sb = []
    pre_react_atomID_1 = []
    pre_react_atomID_2 = []  
    indexs = []
    atomID_at_indics = []
    for i in range(len(preindex_dict_key)):
        for m in range(len(pre_react_1)):
            if potencial_pre_atomID[i][1] == pre_react_1[m]:
                pre_react_atomID_1_sb.append(potencial_pre_atomID[i][0])
                indexs.append(i)
    
    elements_at_indices = [preneighboursindex_Dict[index] for index in indexs]
    atomID_at_indics = [preindex_dict_value[index] for index in indexs]
    pre_atom_divorce = []
    i = 0
    for p in range(len(elements_at_indices)):
        for n in range(len(elements_at_indices[p])):
            for q in range(len(pre_react_2)):
                if elements_at_indices[p][n] == pre_react_2[q]:
                    pre_atom_divorce.append([])
                    pre_react_atomID_1.append(pre_react_atomID_1_sb[p])
                    pre_react_atomID_2.append(atomID_at_indics[p][n])
                    pre_atom_divorce[i].append(pre_react_atomID_1_sb[p])
                    pre_atom_divorce[i].append(atomID_at_indics[p][n])
                    i = i+1
                    #此处的原子对是在反应前就有的配对，等会配对的时候不需要再配对
    print(pre_atom_divorce)#no use

    #注意此处pre和post开始有所不同
    for i in range(len(preindex_dict_key)):
        for m in range(len(pre_react_2)):
            if potencial_pre_atomID[i][1] == pre_react_2[m]:
                pre_react_atomID_2_sb.append(potencial_pre_atomID[i][0])
    pre_atomID_1_1 =[]
    pre_atomID_1_2 =[]
    pre_atomID_2_1 =[]
    pre_atomID_2_2 =[]

    for i in range(len(pre_react_atomID_1_sb)):
        for m in range(len(potencial_1_moleculeID)):
            if pre_react_atomID_1_sb[i] == potencial_1_moleculeID[m]:
                pre_atomID_1_1.append(pre_react_atomID_1_sb[i])
    for i in range(len(pre_react_atomID_1_sb)):
        for m in range(len(potencial_2_moleculeID)):
            if pre_react_atomID_1_sb[i] == potencial_2_moleculeID[m]:
                pre_atomID_1_2.append(pre_react_atomID_1_sb[i])
    for i in range(len(pre_react_atomID_2_sb)):
        for m in range(len(potencial_1_moleculeID)):
            if pre_react_atomID_2_sb[i] == potencial_1_moleculeID[m]:
                pre_atomID_2_1.append(pre_react_atomID_2_sb[i])
    for i in range(len(pre_react_atomID_2_sb)):
        for m in range(len(potencial_2_moleculeID)):
            if pre_react_atomID_2_sb[i] == potencial_2_moleculeID[m]:
                pre_atomID_2_2.append(pre_react_atomID_2_sb[i])

    #Output the atom pairs before and after the reaction
    pre_sets = [] 
    post_sets = [] 
    i = 0
    for m in range(len(pre_atomID_1_1)):#pre_atomID_1
        for n in range(len(pre_atomID_2_2)):#pre_atomID_2
            pre_sets.append([])
            pre_sets[i].append(pre_atomID_1_1[m])
            pre_sets[i].append(pre_atomID_2_2[n])
            i = i+1            

    for i in range((len(post_react_atomID_1))):
        post_sets.append([])
        post_sets[i].append(post_react_atomID_1[i])
        post_sets[i].append(post_react_atomID_2[i])

    return pre_sets,post_sets
 
def try_map(
    directory, 
    preDataFileName, 
    postDataFileName, 
    preMoleculeFileName, 
    postMoleculeFileName, 
    pre_sets,post_sets, 
    deleteAtoms, 
    elementsByType, 
    createAtoms, 
    debug=False):

    filename = 'automap.data'  
    file_path = os.path.join(directory, filename)  
    m=0
    for i in range(len(pre_sets)):
        for n in range(len(post_sets)):
            try:
                command = map_processor(directory, preDataFileName, postDataFileName, preMoleculeFileName, postMoleculeFileName, pre_sets[i], post_sets[n], deleteAtoms, elementsByType, createAtoms, debug=False)          
                if os.path.exists(file_path):                    
                    m = m+1
                    if m == 1:
                        with open('output.txt', 'w') as f:  
                            f.write("新的map文件被生成：")
                        with open('output.txt', 'a') as f:
                            f.write(f"{pre_sets[i]}\n")
                        with open('output.txt', 'a') as f:                              
                            f.write(f"{post_sets[n]}\n")                              
                        os.rename('automap.data', 'automap1') 
                        os.rename('pre-molecule.data', 'pre-molecule1.data')
                        os.rename('post-molecule.data', 'post-molecule1.data')                             
                    else:
                        print(pre_sets[i])
                        print(post_sets[n])
                        with open('output.txt', 'a') as f:  
                            f.write("存在多个map文件被生成，请手动判断map文件是否合理")
                        with open('output.txt', 'a') as f:
                            f.write(f"{pre_sets[i]}\n")
                        with open('output.txt', 'a') as f:
                            f.write("post_sets:")    
                        with open('output.txt', 'a') as f:                              
                            f.write(f"{post_sets[n]}\n") 
                        new_automap = 'automap{}'.format(m) 
                        new_pre = 'pre-molecule{}.data'.format(m) 
                        new_post  = 'post-molecule{}.data'.format(m)                                                                             
                        os.rename('automap.data', new_automap) 
                        os.rename('pre-molecule.data', new_pre)
                        os.rename('post-molecule.data', new_post)

            except Exception as e:    
                print(f"发生了一个未预期的异常：{e}")               
    return  m

def get_elements(neighbourIDs, elementDict):
    return [elementDict[atomID]for atomID in neighbourIDs]

def lammps_to_molecule_masses(fileName, element_1, element_2, validIDSet=None, renumberedAtomDict=None):

    tidiedLines = clean_data(fileName)

    sectionIndexList = find_sections(tidiedLines)

    masses = get_data('Masses', tidiedLines, sectionIndexList)
    masses = refine_data(masses, 0, validIDSet, renumberedAtomDict)
    noselecttypes = get_data('Atoms',tidiedLines, sectionIndexList)
    types = [[] for _ in range(len(noselecttypes))]  
    for i in range(len(noselecttypes)):
        types[i].append(noselecttypes[i][0])
        types[i].append(noselecttypes[i][2])

    atoms = [row[0] for row in types]
    bonds = get_data('Bonds',tidiedLines, sectionIndexList)

    potencial_elementID = []
    mass_str=[]
    pre_react_1=[]
    pre_react_2=[]

    for i in range(0,len(masses)):
        mass = masses[i][1] 
        if mass == '1.008':
            mass_str = 'H'
            if element_1 == mass_str:
                potencial_elementID.append(masses[i][0])
                pre_react_1.append(masses[i][0])

        elif mass == '12.011':    
            mass_str = 'C'
            if element_1 == mass_str:
                potencial_elementID.append(masses[i][0])
                pre_react_1.append(masses[i][0])

        elif mass == '14.007':
            mass_str = 'N'
            if element_1 == mass_str:
                potencial_elementID.append(masses[i][0])
                pre_react_1.append(masses[i][0])
 
        elif mass == '15.999':
            mass_str = 'O'
            if element_1 == mass_str:
                potencial_elementID.append(masses[i][0])
                pre_react_1.append(masses[i][0])
 
        elif mass == '32.06':
            mass_str = 'S'
            if element_1 == mass_str:
                potencial_elementID.append(masses[i][0])
                pre_react_1.append(masses[i][0])


    for i in range(0,len(masses)):
        mass = masses[i][1] 
        if mass == '1.008':
            mass_str = 'H'
            if element_1 == mass_str:
                potencial_elementID.append(masses[i][0])
                pre_react_2.append(masses[i][0])
              
        elif mass == '12.011':    
            mass_str = 'C'
            if element_2 == mass_str:
                potencial_elementID.append(masses[i][0])
                pre_react_2.append(masses[i][0])

        elif mass == '14.007':
            mass_str = 'N'
            if element_2 == mass_str:
                potencial_elementID.append(masses[i][0])
                pre_react_2.append(masses[i][0])
 
        elif mass == '15.999':
            mass_str = 'O'
            if element_2 == mass_str:
                potencial_elementID.append(masses[i][0])
                pre_react_2.append(masses[i][0])
 
        elif mass == '32.06':
            mass_str = 'S'
            if element_2 == mass_str:
                potencial_elementID.append(masses[i][0])
                pre_react_2.append(masses[i][0])
    return potencial_elementID,pre_react_1,pre_react_2,atoms,bonds,types

def get_moleculeID(fileName,validIDSet=None, renumberedAtomDict=None):
    
    tidiedLines = clean_data(fileName)

    sectionIndexList = find_sections(tidiedLines)
    potencial_1_moleculeID = []
    potencial_2_moleculeID = []
    atoms = get_data('Atoms', tidiedLines, sectionIndexList)
    atoms = refine_data(atoms, 0, validIDSet, renumberedAtomDict)
    for i in range(0,len(atoms)):
        molecule = atoms[i][1]
        if molecule == '1':
            potencial_1_moleculeID.append(atoms[i][0])
        else:
            potencial_2_moleculeID.append(atoms[i][0])

    return potencial_1_moleculeID,potencial_2_moleculeID

def ask_user_for_confirmation(): 
    #user_input = input("Do you want to confirm? (yes/no): ")  
    #user_input = user_input.lower().strip()  
    user_input = 'yes'
    if user_input == 'yes':  
        print("Continuing...")  

    elif user_input == 'no':  
        print("Exiting...")  
        sys.exit() 
    else:  
        print("Invalid input. Please enter 'yes' or 'no'.")

def delete_redundant_mapfile(correct_times):
    files_to_delete =[]
    for i in range(2,correct_times+1):
        files_to_delete.append('automap{}'.format(i)) 
    file_to_rename = 'automap1'  
    new_name = 'automap.data'  

    for file_name in files_to_delete:  
        try:  
            if os.path.exists(file_name):  
                os.remove(file_name)  
                print(f"Deleted {file_name}")  
            else:  
                print(f"{file_name} does not exist.")  
        except OSError as e:  
            print(f"Error: {e.strerror}")  
 
    try:  
        if os.path.exists(file_to_rename):  
            os.rename(file_to_rename, new_name)  
            print(f"Renamed {file_to_rename} to {new_name}")  
        else:  
            print(f"{file_to_rename} does not exist.")  
    except OSError as e:  
        print(f"Error: {e.strerror}")

def delete_redundant_prefile(correct_times):
    files_to_delete =[]
    for i in range(2,correct_times+1):
        files_to_delete.append('pre-molecule{}.data'.format(i)) 
    file_to_rename = 'pre-molecule1.data'  
    new_name = 'pre-molecule.data'  

    for file_name in files_to_delete:  
        try:  
            if os.path.exists(file_name):  
                os.remove(file_name)  
                print(f"Deleted {file_name}")  
            else:  
                print(f"{file_name} does not exist.")  
        except OSError as e:  
            print(f"Error: {e.strerror}")  
 
    try:  
        if os.path.exists(file_to_rename):  
            os.rename(file_to_rename, new_name)  
            print(f"Renamed {file_to_rename} to {new_name}")  
        else:  
            print(f"{file_to_rename} does not exist.")  
    except OSError as e:  
        print(f"Error: {e.strerror}")

def delete_redundant_postfile(correct_times):
    files_to_delete =[]
    for i in range(2,correct_times+1):
        files_to_delete.append('post-molecule{}.data'.format(i)) 
    file_to_rename = 'post-molecule1.data'  
    new_name = 'post-molecule.data'  

    for file_name in files_to_delete:  
        try:  
            if os.path.exists(file_name):  
                os.remove(file_name)  
                print(f"Deleted {file_name}")  
            else:  
                print(f"{file_name} does not exist.")  
        except OSError as e:  
            print(f"Error: {e.strerror}")  

    try:  
        if os.path.exists(file_to_rename):  
            os.rename(file_to_rename, new_name)  
            print(f"Renamed {file_to_rename} to {new_name}")  
        else:  
            print(f"{file_to_rename} does not exist.")  
    except OSError as e:  
        print(f"Error: {e.strerror}")

def check_map(correct_times):
    if correct_times == 0:
        print("没有map生成，请检查输入的交联键两方元素的正确性，如仍无法生成map请依据原著手动输入生成map【注意：如果map过程中设计删除副产物原子也请手动生成map】")
        sys.exit()
    elif correct_times == 1:
        print("有且只有一个map生成，无需担心")
        delete_redundant_mapfile(correct_times)
        delete_redundant_prefile(correct_times)
        delete_redundant_postfile(correct_times)        
    elif correct_times >= 1:
        print("有{}个map生成，请检查原始结构是否有相同的交联位点，例如4官能度的和2官能度的两种反应物有2x4=8, 8种map方法，如果生成的map数量正确则无需担心，回复yes。如果生成多余数量的map,请根据output.txt查看原子对手动生成map，回复no".format(correct_times))
        ask_user_for_confirmation()
        delete_redundant_mapfile(correct_times)
        delete_redundant_prefile(correct_times)
        delete_redundant_postfile(correct_times)  
    return


def new_map_processor(directory, preDataFileName, postDataFileName, preMoleculeFileName, postMoleculeFileName, preBondingAtomselement, postBondingAtomselement, deleteAtoms, elementsByType, createAtoms, debug=False):
    #获得所有可能反应元素的对应的序号【5，3，4】，所有原子的原子序号和原子类型列表potencial_pre_atomID【1，5】
    with open(postDataFileName, 'r') as f:
        lines = f.readlines()
    potencial_post_elementID,post_react_1_no_use,post_react_2_no_use,atomIDs_2,bonds_2,types_2 = lammps_to_molecule_masses(lines, preBondingAtomselement, postBondingAtomselement, validIDSet=None, renumberedAtomDict=None)
    potencial_post_atomID = get_atomID(potencial_post_elementID,types_2)

    #获得所有可能反应元素的对应的序号【5，3，4】，所有原子的原子序号和原子类型列表potencial_pre_atomID【1，5】
    with open(preDataFileName, 'r') as f:
        lines = f.readlines()
    potencial_pre_elementID,pre_react_1,pre_react_2,atomIDs_1,bonds_1,types_1 = lammps_to_molecule_masses(lines, preBondingAtomselement, postBondingAtomselement, validIDSet=None, renumberedAtomDict=None)
    potencial_pre_atomID = get_atomID(potencial_pre_elementID,types_1)

    #确定在反应前原子各自处于哪个分子之中
    potencial_1_moleculeID,potencial_2_moleculeID = get_moleculeID(lines,validIDSet=None, renumberedAtomDict=None)

    #获得所有原子的邻居的atomID{1:[2,18,19]}
    potencial_pre_neighbor_atomID = get_neighbours(atomIDs_1,bonds_1)
    potencial_post_neighbor_atomID = get_neighbours(atomIDs_2,bonds_2)
    #获得所有可能反应邻居的原子序号和邻居的原子种类preindexdict和preneighboursindexDict
    preneighboursindexDict,preindexdict = get_neighbours_element(potencial_pre_neighbor_atomID,potencial_pre_atomID,types_1)
    postneighboursindexDict,postindexdict = get_neighbours_element(potencial_post_neighbor_atomID,potencial_post_atomID,types_2)

    pre_sets,post_sets = compare_atomID_difference(potencial_pre_atomID,potencial_post_atomID,pre_react_1,pre_react_2,preneighboursindexDict,postneighboursindexDict,preindexdict,postindexdict,potencial_1_moleculeID,potencial_2_moleculeID)

    correct_times = try_map(directory, preDataFileName, postDataFileName, preMoleculeFileName, postMoleculeFileName, pre_sets,post_sets, deleteAtoms, elementsByType, createAtoms, debug=False)
     
    final_check = check_map(correct_times)
