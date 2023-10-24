import sys
import pybel
import yaml

f=open(sys.argv[2],'w')
f.write(sys.argv[3]+'\n')
for mol in pybel.readfile("sdf", sys.argv[1]):
    f.write(mol.title+','+mol.data["Score"]+'\n')

f.close()

