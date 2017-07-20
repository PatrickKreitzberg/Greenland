linCMFile = open('./bathCPT.txt', 'r')
out = open('./b.py','w')
out.write('[')
for line in linCMFile:
    vals = line.split()
    s = '[' + str(vals[1]) + ',' + str(vals[2]) + ','+ str(vals[3]) + '],\n'
    out.write(s)
    # v0 = []
    # v0.append(int(vals[1]))
    # v0.append(int(vals[2]))
    # v0.append(int(vals[3]))

out.write(']')


linCMFile.close()
out.close()
