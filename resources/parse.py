import sys
import matplotlib.pyplot as plt

# To use:
# python parse.py file_dashing file_nihm file_output_pdf.pdf

def build_dict(file_name):
    preambule=True
    d = dict()

    f = open(file_name, "r")
    for x in f:
        if not preambule:
            li_values = [y for y in x.strip().split('\t')[1:] if len(y)>0]
            for col in range(line+1,len(li_values)):
                d[li_names[line],li_names[col]] = float(li_values[col])
            line += 1

        if x[:2] == "##":
            preambule = False
            li_names = [y for y in x.strip().split('\t')[1:] if len(y)>0]
            line = 0

    return d

file_name_dashing = sys.argv[1]
file_name_nihm = sys.argv[2]
file_name_figoutput = sys.argv[3]

di_dashing = build_dict(file_name_dashing)
di_nihm = build_dict(file_name_nihm)

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
