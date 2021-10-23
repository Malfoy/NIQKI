import sys
import matplotlib.pyplot as plt

# To use:
# python parse.py file_dashing file_nihm file_output_pdf.pdf

def build_dashing_dict(file_name):
    preambule=True
    dd = dict()
    f = open(file_name, "r")
    for x in f:
        if not preambule:
            li_values = sorted([y for y in x.strip().split('\t')[1:] if len(y)>0])
            for col in range(line+1,len(li_values)):
                dd[li_names[line],li_names[col]] = float(li_values[col])
            line += 1

        if x[:2] == "##":
            preambule = False
            li_names = sorted([y for y in x.strip().split('\t')[1:] if len(y)>0])
            line = 0

    return dd

def build_nihm_dict(file_name):
    dn = dict()
    f = open(file_name, "r")
    for x in f:
        if x[:2] == "##":
            li_names = sorted([y for y in x.strip().split('\t')[1:] if len(y)>0])
            line = 0
        li_values = sorted([y for y in x.strip().split('\t')[1:] if len(y)>0])
        for col in range(line+1,len(li_values)):
          dn[li_names[line],li_names[col]] = float(li_values[col])
        line += 1
    return dn

file_name_dashing = sys.argv[1]
file_name_nihm = sys.argv[2]
file_name_figoutput = sys.argv[3]

di_dashing = build_dashing_dict(file_name_dashing)
di_nihm = build_dashing_dict(file_name_nihm)


#print("---- Print dashing dictionary ----")
#i = 0
#for key, value in di_dashing.items() :
#    print (key, value)
#    i += 1
#print(i)
#i = 0
#print("---- Print nihm dictionary ----")
#for key, value in di_nihm.items() :
#    print (key, value)
#    i += 1
#print(i)

if di_dashing.keys() != di_nihm.keys():
    print("Not the same set of files")
    exit(0)

scores_dashing = []
scores_nihm = []

for (f1,f2) in di_dashing:
    scores_dashing.append(di_dashing[(f1,f2)])
    scores_nihm.append(di_nihm[(f1,f2)])
plt.scatter(scores_dashing, scores_nihm, marker='o')


plt.savefig(file_name_figoutput)





#
